"""
Leonardo Galeotti 11-12-2025
Creazione della classe MyProtein che sostituisce la obsoleta Protiein di p3d per l'estrapolazione di informazioni strutturali e non essenziali all'allineamento dei siti.
"""

import os

from Bio.PDB import PDBParser

# import sys

metalsTempFactor = 50.00
donorsTempFactor = 40.00


class Atom:
    def __init__(
        self, aa, resid, chain, x, y, z, atype, type_, element, element_type, beta, idx
    ):
        self.aa = aa
        self.resid = resid
        self.chain = chain

        self.x = x
        self.y = y
        self.z = z

        self.atype = atype
        self.type_ = type_
        self.element = element
        self.element_type = element
        self.beta = beta
        self.idx = idx

    def clone(self):
        return Atom(
            aa=self.aa,
            resid=self.resid,
            chain=self.chain,
            x=self.x,
            y=self.y,
            z=self.z,
            atype=self.atype,
            type_=self.type_,
            element=self.element,
            element_type=self.element_type,
            beta=self.beta,
            idx=self.idx,
        )


class MyProtein:
    def __init__(self, pdbfile):
        if not os.path.exists(pdbfile):
            raise FileNotFoundError("File does not exist!")

        parser = PDBParser(QUIET=True)
        self.filename = pdbfile

        # Carica struttura biopython
        self.structure = parser.get_structure("protein", pdbfile)

        # Lista di nostri atomi compatibili con MetalS2
        self.atoms_raw = []

        # Itera atomi BioPython → Atom personalizzato
        for i, atom in enumerate(self.structure.get_atoms()):
            aa = atom.get_parent().get_resname()  # residuo
            resid = atom.get_parent().id[1]  # numero residuo
            chain = atom.get_parent().get_parent().id  # es: "A"
            x, y, z = atom.get_coord()  # coordinate
            atype = atom.get_name()  # nome atomo

            # Se è un HETATM: resname + flag residuo
            type_ = "HETATM" if "H_" in atom.get_parent().id[0].strip() else "ATOM"

            # Filter / remove waters
            if aa == "HOH":
                continue

            element = atom.element
            element_type = element
            beta = atom.get_bfactor()

            # Crea Atom compatibile con MetalS2
            self.atoms_raw.append(
                Atom(
                    aa,
                    resid,
                    chain,
                    x,
                    y,
                    z,
                    atype,
                    type_,
                    element,
                    element_type,
                    beta,
                    idx=i,
                )
            )


class InputError(Exception):
    pass


if __name__ == "__main__":
    print("Test")
    pdb_file_path_test = "/home/leonardo/Scrivania/metals2_2.0/data/1EN7.pdb"
    protein = MyProtein(pdbfile=pdb_file_path_test)
    atoms = protein.atoms_raw
    print(vars(atoms[1]))
    # Estraiamo le coordinati di uno ione metallico
    metals = ["ZN", "CA"]
    ions = [x.type_ for x in atoms if x.aa in metals]
    print(ions)
    # Test clone
    # pippo = [x.clone() for x in atoms]
    # for x in pippo:
    #     x.beta = 2
    # vars_pippo = [vars(x) for x in pippo[:10]]
    # print(vars_pippo)
    # vars_native = [vars(x) for x in atoms[:10]]
    # print(vars_native)
