#!/usr/bin/env python2.6

import os
import sys

sys.path.append(os.path.join(os.path.expanduser("~"), "development"))

import glob
import shutil
from collections import defaultdict
from operator import attrgetter

import numpy as np
from Bio.PDB import PDBIO, Chain, Model, Residue, Structure
from Bio.PDB import Atom as BioAtom
from lib.nomenclature import metalList
from lib.sequence import writeSequenceToFile


class ScoreReport:
    "Collection of data for a score reports"

    def __init__(self, site1, site2):
        self.getNames(site1, site2)

        self.rmsd = None

        self.p1 = None
        self.p2 = None
        self.p3 = None

        self.totalScore = None
        self.percentage = None

        self.metals1 = ""
        self.metals2 = ""
        self.ligands1 = ""
        self.ligands2 = ""

        self.metals1 = ", ".join(
            f"{atom.aa}_{atom.resid}_{atom.chain}" for atom in site1.metals
        )
        self.metals2 = ", ".join(
            f"{atom.aa}_{atom.resid}_{atom.chain}" for atom in site2.metals
        )
        self.ligands1 = ", ".join(
            f"{atom.aa}_{atom.resid}_{atom.chain}" for atom in site1.donorAtoms
        )
        self.ligands2 = ", ".join(
            f"{atom.aa}_{atom.resid}_{atom.chain}" for atom in site2.donorAtoms
        )

        self.initial_score = None
        self.min_initial_score = None
        self.max_initial_score = None
        self.threshold_percentage = None

    def getNames(self, site1, site2):
        self.name1 = site1.name
        self.name2 = site2.name

    def setScores(self, rmsd, p1, p2, p3, totalScore, percentage):
        self.rmsd = rmsd

        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

        self.totalScore = totalScore
        self.percentage = percentage


def storeAlignment(site1, site2, outputPath, writeSites=True):
    if not os.path.exists(outputPath):
        try:
            os.mkdir(outputPath)
        except OSError:
            print(
                "Alignment is not complete. Check if you have permissions to create directories."
            )
            sys.exit()

    rootDir = os.path.join(outputPath, f"{site1.name}_vs_{site2.name}")
    try:
        os.mkdir(rootDir)
    except Exception:
        # print '\nOutput directory "%s" already exists.\nDo you want to overwrite it? (y/n)\n\nWarning! Deleting a directory will also delete all contents within the directory.\n' % (sitesDirRoot)
        overwriteDir = True

        if overwriteDir:  # The directory is to be overwritten
            shutil.rmtree(rootDir)
            os.mkdir(rootDir)
        if not overwriteDir:
            sys.exit()

    if writeSites:
        sitesDir = os.path.join(rootDir, "sites")
        os.mkdir(sitesDir)

        writeSiteToFile(sitesDir, f"{site1.name}_query.site.pdb", site1)
        writeSiteToFile(sitesDir, f"{site2.name}_target.site.pdb", site2)
        visualizeAlignment(sitesDir)

    return rootDir


def reportAlignment(site1, site2, scoreReport, outputPath, maxDist, scoreOnly=False):
    if scoreOnly:
        scorePath = os.path.join(outputPath, f"score_{site1.name}_vs_{site2.name}.txt")
    else:
        scorePath = os.path.join(outputPath, "score.txt")
    writeScoresToFile(scoreReport, scorePath)

    if not scoreOnly:
        sequencePath = os.path.join(outputPath, "sequence.txt")
        writeSequenceToFile(site1, site2, sequencePath, maxDist)


def writeScoresToFile(scoreReport, path):
    out_file = open(path, "w")

    name1 = f"{scoreReport.name1}_{scoreReport.metals1}"
    fst_line = "{:<10}{:<25}\n".format("Name 1:", name1)
    out_file.write(fst_line)
    name2 = f"{scoreReport.name2}_{scoreReport.metals2}"
    snd_line = "{:<10}{:<25}\n".format("Name 2:", name2)
    out_file.write(snd_line)

    out_file.write("\n")

    rmsd_line = "{:<27}{:<25}\n".format("RMSD over aligned region:", scoreReport.rmsd)
    out_file.write(rmsd_line)

    out_file.write("\n")

    p1_line = "{:<27}{:<25}\n".format("Relative coverage term:", scoreReport.p1)
    out_file.write(p1_line)
    p2_line = "{:<27}{:<25}\n".format("Chemical similarity term:", scoreReport.p2)
    out_file.write(p2_line)
    p3_line = "{:<27}{:<25}\n".format("Continuity term:", scoreReport.p3)
    out_file.write(p3_line)

    out_file.write("\n")

    total_line = "{:<27}{:<25}\n".format("Total score:", scoreReport.totalScore)
    out_file.write(total_line)

    out_file.close()


def visualizeAlignment(path):
    fileName = "visualisation.pml"
    filePath = os.path.join(path, fileName)

    out = open(filePath, "w")

    outList = []
    for fname in glob.iglob(os.path.join(path, "*.pdb")):
        nm = f"load {fname.split('/')[-1]}\n"
        outList.append(nm)

    outList.append("deselect\n")
    outList.append("hide everything\nshow cartoon\n")
    outList.append(f"sele me,resn {metalList[0]}")
    for i in range(1, len(metalList)):
        outList.append(f" or resn {metalList[i]}")

    outList.append(
        "\nshow spheres,me\nset sphere_scale, 0.65\ncolor red, me\ndeselect\n"
    )
    out.writelines(outList)

    out.close()


def writeSiteToFile(path, name, site):
    "Writes an ASCII text representation of an MetalSite object to a file"

    # Defines output name
    trgPath = os.path.join(path, name)

    # 1. Create empy structure
    structure = Structure.Structure(f"{name}")
    model = Model.Model(0)
    structure.add(model)
    # 2. Group for chain and residue
    chains = defaultdict(lambda: defaultdict(list))
    whole_atoms = []
    whole_atoms.extend(site.metals)
    whole_atoms.extend(site.donorAtoms)
    whole_atoms.extend(site.approxSphereAtoms)
    for atom in whole_atoms:
        chains[atom.chain][(atom.aa, atom.resid)].append(atom)
    # 3. Hierarchy collection
    serial = 1
    for chain_id in sorted(chains.keys()):
        chain = Chain.Chain(chain_id)
        model.add(chain)

        for (resname, resid), atom_list in chains[chain_id].items():
            if resname in metalList:
                residue = Residue.Residue((f"H {resname}", resid, " "), resname, "")
            else:
                residue = Residue.Residue((" ", resid, " "), resname, "")
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
    # 4. write PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(trgPath)

    return trgPath
