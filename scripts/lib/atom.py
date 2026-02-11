#!/usr/bin/env python

from lib.constants import *
from lib.nomenclature import *


def atomIsFromProtein(atom):
    is_protein = False
    res_name = atom.aa.upper().strip()
    if res_name in aminoAcids:
        is_protein = True

    return is_protein


def atomIsFromBackbone(atom):
    roleName = atom.atype.upper().strip()
    if (roleName == roleName_proteinBackbone) or (roleName == roleName_nucleicBackbone):
        return True
    else:
        return False


def atomIsFromSideChain(atom):
    roleName = atom.atype.upper().strip()
    resName = atom.aa.upper().strip()
    if (
        (roleName == roleName_proteinSideChain)
        or ((resName in C_U_T) and (roleName == roleName_nucleicPurines))
        or ((resName in A_G) and (roleName == roleName_nucleicPyrimidines))
    ):
        return True
    else:
        return False
    
def atomIsDonor(atom):
    if atom.beta == donorsTempFactor:
        return True
    else:
        return False

def atomIsFromLigand(atom):
    if atom.beta == ligandsTempFactor:
        return True
    else:
        return False


def atomIsNeighbour(atom):
    if (atom.beta == neighborDonorsTempFactor) or (
        atom.beta == neighborLigandsTempFactor
    ):
        return True
    else:
        return False


def bothAtomsAreFromSameMolecule(atom1, atom2):
    if atomIsFromProtein(atom1) and atomIsFromProtein(atom2):
        return True
    else:
        return False


def getOneLetterResudueCode(atom):
    res = atom.aa.upper().strip()

    if res in standardAminoAcids:
        code = aminoacidCode.get(res)

    elif res in nonStandardAminoAcids.keys():
        res_sub = nonStandardAminoAcids[res]
        code = aminoacidCode.get(res_sub)

    elif res in standardNucleicAcids:
        if len(res) > 1:
            code = res[-1]
        else:
            code = res

    elif res in nonStandardNucleicAcids.keys():
        res_sub = nonStandardNucleicAcids[res]
        if len(res_sub) > 1:
            code = res_sub[-1]
        else:
            code = res_sub

    else:
        code = "X"

    return code
