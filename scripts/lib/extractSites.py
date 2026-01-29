#!/usr/bin/env python2.6

import math
import os
import sys
from collections import defaultdict
from operator import attrgetter

import numpy as np
from Bio.PDB import PDBIO, Chain, Model, NeighborSearch, Residue, Structure
from Bio.PDB import Atom as BioAtom
from biopython import MyProtein
from lib.metalSite import MetalSite
from lib.nomenclature import metalList

# ===============================================================================
"""CLASS"""
# ===============================================================================


class SourcePDB:
    def __init__(self, pathPdb, maxDist=5.0):
        self.maxDist = maxDist

        self.code = pathPdb.split("/")[-1].replace(".pdb", "").lower()
        print(self.code)

        self.metals = []
        self.donors = {}

        self.sites = {}
        self.ligands_site = {}

        self.sitesSource = {}
        self.atoms4site = {}  # ? vuole essere un dizionario di liste key = site e value = lista atomi (oggetto)

    def findMetals(self, atoms_list):  # atoms_list added
        # 'Take out from pdb all the metal atoms found in list of metals'
        print("findMetals")
        self.metals = [x for x in atoms_list if x.aa in metalList]
        if not self.metals:
            raise Exception("Unable to find specified metal in the PDB file.")
            sys.exit()
        else:
            for met in self.metals:
                met.beta = 50.00
        # ? Test sostituzione findMetals() con Oggetto Atom created with Biopython

    def findDonors(self, atom_list):
        not_donors = ["C", "H", "D"]
        # distance = 3.0
        distance = 2.8
        print("findDonors")
        ligands_names = {}

        neighbor_search = NeighborSearch(atom_list)
        for metal in self.metals:
            self.donors[metal] = []
            metal_coords = np.array((metal.x, metal.y, metal.z), dtype="d")
            neighbor = neighbor_search.search(metal_coords, distance)
            # print(neighbor)
            # filter / remove metal and not donors from neighbor
            for atom in neighbor:
                if (atom.element not in not_donors) and (atom.element not in metalList):
                    atom.beta = 40.00
                    self.donors[metal].append(atom)
            if not self.donors[metal]:
                ligands_names[metal] = "-"
            else:
                temp = []
                for lig in self.donors[metal]:
                    temp.append(f"{lig.aa}_{lig.resid}_{lig.chain}")
                    ligands_names[metal] = list(set(temp))
            #! DEBUG =====================
            # for x in self.donors[metal]:
            #     print(vars(x))
            #! ===========================
        return ligands_names

    def findSites(self, ligands_names):
        "Cluster metals into sites"
        print("findSites")
        maxDist = 3.5
        unassigned = list(self.metals)

        site_numb = 0
        while unassigned:
            site_numb = site_numb + 1
            self.sites[site_numb] = [unassigned.pop(0)]
            new_member = True
            while new_member:
                new_member = False
                to_add = []
                for metal in unassigned:
                    for mt in self.sites[site_numb]:
                        dx = mt.x - metal.x
                        dy = mt.y - metal.y
                        dz = mt.z - metal.z
                        distance = math.sqrt(dx * dx + dy * dy + dz * dz)
                        ligands_in_common = set(ligands_names[metal]) & set(
                            ligands_names[mt]
                        )
                        if ligands_in_common == {"no_ligands"} or ligands_in_common == {
                            "-"
                        }:
                            ligands_in_common = []
                        if distance < maxDist or ligands_in_common:
                            to_add.append(metal)
                            break

                if to_add:
                    for metal in to_add:
                        self.sites[site_numb].append(metal)
                        unassigned.remove(metal)
                    new_member = True

    # Cloning module
    def Clone_atoms_4_site(self, atom_list):
        for site in self.sites:
            self.atoms4site[site] = [x.clone() for x in atom_list]

    # I expect as many keys as there are sites, and each site contains a list of atoms that is a copy of all the atoms, so that they can be modified in isolation.

    def findFirstSphere(self, userSites=None):
        print("findFirstSphere")

        if userSites:
            self.sites = userSites

        for site in self.sites:
            self.sitesSource[site] = []
            self.ligands_site[site] = []
            self.sitesSource[site].extend(self.sites[site])  # 50.00 - metals
            donors = []

            for metal in self.sites[site]:
                donors_temp = list(self.donors[metal])
                donors.extend(donors_temp)

            donors = list(set(donors))

            atomsFromSameDonorResidue_temp = []

            donor_keys = {(x.aa, x.resid, x.chain) for x in donors}
            atom_list = self.atoms4site[site]
            # ? DONORS has been created before cloning so the object's ID is different and we must exclude by attributes' comparison
            filtered_atom_list = [
                x
                for x in atom_list
                if ((x.aa, x.resid, x.chain) in donor_keys)
                and all(vars(x) != vars(y) for y in donors)
            ]
            #! DEBUG =======================
            # donors_debug_ = {
            #     (x.aa, x.resid, x.chain, x.element) for x in self.donors[self.sites[site][0]]
            # }
            # debug_check_ = [x for x in filtered_atom_list if ((x.aa, x.resid, x.chain, x.element) in donors_debug_)]
            # for x in debug_check_:
            #    print(vars(x))
            #! =============================
            for x in filtered_atom_list:
                x.beta = 30.00
                atomsFromSameDonorResidue_temp.append(x)
            # ---------------------------------------------------------
            atomsFromSameDonorResidues = list(set(atomsFromSameDonorResidue_temp))

            self.ligands_site[site].extend(donors)  # 40.00 - donors
            self.ligands_site[site].extend(
                atomsFromSameDonorResidues
            )  # 30.00 - ligands 1st sphere

            self.sitesSource[site].extend(donors)
            self.sitesSource[site].extend(atomsFromSameDonorResidues)

    def findSecondSphere(self, sites=None):
        print("findSecondSphere")

        for site in sites:
            allAtomsOfSecondSphere_temp = []
            secondSphereDonors_temp = []
            maxi_temp = []
            excluder = []
            atom_list = self.atoms4site[site]
            for atom in self.ligands_site[site]:
                excluder.append(atom)
                # ? ATTENTION: the donors atom inside excluder and inside ligands_site has been collecter before cloning so the id is different from those found throught neighbor search, so they won't be excluded by the filter because of the diffent ID. Again we must execute by attributes' comparison.
                neighbour_search = NeighborSearch(atom_list)
                atom_coords = np.array([atom.x, atom.y, atom.z], dtype="d")
                temp = neighbour_search.search(atom_coords, self.maxDist)
                maxi_temp.extend(temp)

            maxi_temp = list(set(maxi_temp))  # remove duplets
            filter_temp = [
                x
                for x in maxi_temp
                if (x not in excluder)
                and all(vars(x) != vars(y) for y in self.metals)
                and all(
                    vars(x) != vars(y)
                    for metal in self.sites[site]
                    for y in self.donors[metal]
                )
            ]
            for atm in filter_temp:
                if atm.type_ != "H" and atm.type_ != "D":
                    atm.beta = 20.00
                    secondSphereDonors_temp.append(atm)

            secondSphereDonors = list(set(secondSphereDonors_temp))
            self.sitesSource[site].extend(secondSphereDonors)

            # searching for atoms further than 5A but belonging to res within 5A.
            filter_dict = {(x.aa, x.resid, x.chain) for x in self.sitesSource[site]}
            atoms_in_sec = [
                x
                for x in atom_list
                if ((x.aa, x.resid, x.chain) in filter_dict)
                and (x not in self.sitesSource[site])
                and (x.type_ != "HETATM")
                and all(
                    vars(x) != vars(y)
                    for metal in self.sites[site]
                    for y in self.donors[metal]
                )
            ]
            excluder = [
                x
                for x in atoms_in_sec
                if ((x in secondSphereDonors) or (x in self.ligands_site[site]))
            ]
            #! DEBUG ===========================
            # donors_debug_ = {
            #     (x.aa, x.resid, x.chain, x.element)
            #     for metal in self.sites[site] for x in self.donors[metal]
            # }
            # debug_check_ = [x for x in atoms_in_sec if ((x.aa, x.resid, x.chain, x.element) in donors_debug_)]
            # for x in debug_check_:
            #     print(x)
            #! =================================
            atoms_in_sec = [x for x in atoms_in_sec if x not in excluder]
            for atom in atoms_in_sec:
                atom.beta = 10.00
                allAtomsOfSecondSphere_temp.append(atom)

            allAtomsOfSecondSphere = list(set(allAtomsOfSecondSphere_temp))
            self.sitesSource[site].extend(allAtomsOfSecondSphere)

    def findSitesInPDB(self, atoms_list):
        try:
            self.findMetals(atoms_list=atoms_list)
            ligands_names = self.findDonors(atom_list=atoms_list)
            self.findSites(ligands_names)
            self.Clone_atoms_4_site(atom_list=atoms_list)
        except Exception:
            raise Exception("Unable to extract sites")

    def getSites(self):
        self.findFirstSphere()
        self.findSecondSphere(sites=self.sites)

        metalSites = []

        # TEST SOSTITUZIONE MODULO PROTEIN P3D
        # for site in self.sites:
        #     try:
        #         path_to_site_file = self.dumpSiteToFile(site)
        #         inSite = Protein(path_to_site_file)
        #     except:
        #         raise Exception('Unable to open %s' % self.code)

        #     try:
        #         metalSite = MetalSite(inSite)
        #     except Exception:
        #         os.remove(path_to_site_file)
        #         continue
        #     else:
        #         metalSites.append(metalSite)
        #         os.remove(path_to_site_file)

        # return metalSites
        for site in self.sites:
            atom_list = self.sitesSource.get(site, [])
            try:
                path_to_site_file = self.dumpSiteToFile(site, atoms=atom_list)
            except Exception:
                raise Exception(f"Unable to open {self.code}")

            try:
                metalsite = MetalSite(
                    rawSite=site, atoms=atom_list, filename=path_to_site_file
                )
            except Exception:
                os.remove(path_to_site_file)
                continue
            else:
                metalSites.append(metalsite)
                os.remove(path_to_site_file)

        return metalSites

        # TODO: Continuare con la parte di lettura del sito.pdb e return del metalsite

    def dumpSiteToFile(self, site, atoms):
        print("Site:", site)
        atoms = self.sitesSource.get(site, [])
        print("Atoms in sitesSource:", len(atoms))
        #! DEBUG ==============================
        # donors_debug_ = {
        #     (x.aa, x.resid, x.chain, x.element, x.x, x.y, x.z) for metal in self.sites[site] for x in self.donors[metal]
        # }
        # debug_check_ = [x for x in atoms if ((x.aa, x.resid, x.chain, x.element, x.x, x.y, x.z) in donors_debug_)]
        # for x in debug_check_:
        #     print(vars(x))
        #! ====================================
        metal_sites_dir = os.path.dirname(os.path.abspath(__file__))
        output_name = f"{self.code}_{site}.site.pdb"
        mSitePath = os.path.join(metal_sites_dir, output_name)

        # 1. Create empty Structure
        structure = Structure.Structure(f"{self.code}_{site}")
        model = Model.Model(0)
        structure.add(model)

        # 2. Group for chain and residue
        chains = defaultdict(lambda: defaultdict(list))

        for atom in atoms:
            chains[atom.chain][(atom.aa, atom.resid)].append(atom)

        # 3. Hierarchy construction
        serial = 1

        for chain_id in sorted(chains.keys()):
            chain = Chain.Chain(chain_id)
            model.add(chain)

            for (resname, resid), atom_list in chains[chain_id].items():
                if resname in metalList:
                    residue = Residue.Residue((f"H_{resname}", resid, " "), resname, "")
                else:
                    residue = Residue.Residue(
                        (" ", resid, " "),
                        resname,
                        "",  # residue id PDB
                    )
                chain.add(residue)

                for atom in sorted(atom_list, key=attrgetter("idx")):
                    bio_atom = BioAtom.Atom(
                        name=atom.atype,
                        coord=np.array([atom.x, atom.y, atom.z], dtype=float),
                        bfactor=float(atom.beta),
                        occupancy=1.0,
                        altloc=" ",
                        fullname=f"{atom.atype:>4}",
                        serial_number=serial,
                        element=atom.element,
                    )
                    serial += 1
                    residue.add(bio_atom)

        # 4. Write PDB
        io = PDBIO()
        io.set_structure(structure)
        io.save(mSitePath)

        return mSitePath


# ===============================================================================
"""METHODS"""
# ===============================================================================


def getSitesFromPdbFile(pathPdb, metalID=None, maxDist=5.0):
    # ? create atom object
    protein = MyProtein.MyProtein(pathPdb)
    atoms = protein.atoms_raw
    sourcePDB = SourcePDB(pathPdb, maxDist)
    try:
        sourcePDB.findSitesInPDB(atoms_list=atoms)
        if metalID:
            metalSites = sourcePDB.getSiteWithMetal(metalID)
        else:
            metalSites = sourcePDB.getSites()
    except Exception:
        sys.exit()

    print(metalSites)
    return metalSites
