#!/usr/bin/env python

import copy

from lib.atom import (
    atomIsFromBackbone,
    atomIsFromLigand,
    atomIsFromProtein,
    atomIsFromSideChain,
    atomIsNeighbour,
)
from lib.constants import *


class MetalSite:
    """
        Single metal site object shifted to the origin of coordinates

    usage MetalSite(metals, donors, sphere)

        metals:   list of metal atoms
        donors:   list of donor atoms (bounded to a metal directly)
        sphere:   list of atoms populating approximate sphere around metals
    """

    def __init__(self, atoms=None, filename=None, rawSite=None):
        if rawSite != None:
            # self.nameExt = rawSite.filename.split('/')[-1]
            # self.name = rawSite.filename.split('/')[-1].split('.')[0]
            self.nameExt = filename.split("/")[-1]
            self.name = filename.split("/")[-1].split(".")[0]

            metals = [
                atom for atom in atoms if (atom.beta == metalsTempFactor)
            ]  # da constants
            donors = [atom for atom in atoms if (atom.beta == donorsTempFactor)]
            sphere = [
                atom
                for atom in atoms
                if ((atom.beta != metalsTempFactor) and (atom.beta != donorsTempFactor))
            ]

            if not metals:
                print("\nNo metals are found or file in a bad format.")
                raise Exception
            if not donors:
                print(
                    f"\n {self.name} has no ligands or file in a bad format. Site will not be processed."
                )
                raise Exception
            else:
                self.gcX, self.gcY, self.gcZ = self.getGeometricalCentre(metals)
                self.gcX_old, self.gcY_old, self.gcZ_old = self.gcX, self.gcY, self.gcZ

                self.metals = self.shiftCoordinates(metals)
                self.donorAtoms = self.shiftCoordinates(donors)
                self.approxSphereAtoms = self.shiftCoordinates(sphere)

                # dictionaries
                self.getAtomsForMatching()

                # lists
                self.getBackboneAndSidechain()

        else:
            self.name = "none"
            self.metals = []
            self.donorAtoms = []
            self.approxSphereAtoms = []

    def getGeometricalCentre(self, metals):
        "Calculates geometric centre of metals and returns coordinates of the geometric centre"

        gcX = 0.0
        gcY = 0.0
        gcZ = 0.0

        nAtoms = len(metals)

        gcX = float(sum([metal.x for metal in metals])) / nAtoms
        gcY = float(sum([metal.y for metal in metals])) / nAtoms
        gcZ = float(sum([metal.z for metal in metals])) / nAtoms

        return gcX, gcY, gcZ

    def shiftCoordinates(self, atoms):
        "Shifts all atoms by coordinates of geometric center"

        for atom in atoms:
            atom.x = atom.x - self.gcX
            atom.y = atom.y - self.gcY
            atom.z = atom.z - self.gcZ

        return atoms

    def getAtomsForMatching(self):
        self.backboneDonors = {}
        self.backboneNeighbours = {}

        self.sideChainDonors = {}
        self.sideChainNeighbours = {}

        for atom in self.approxSphereAtoms:
            if atomIsFromBackbone(atom):
                if atomIsFromLigand(atom):
                    if atomIsFromProtein(atom):
                        handle = "{0}{1}{2}".format(
                            atom.chain.upper().strip(), "_", atom.resid
                        )
                        if handle not in self.backboneDonors.keys():
                            self.backboneDonors[handle] = atom
                elif atomIsNeighbour(atom):
                    if atomIsFromProtein(atom):
                        handle = "{0}{1}{2}".format(
                            atom.chain.upper().strip(), "_", atom.resid
                        )
                        if handle not in self.backboneNeighbours.keys():
                            self.backboneNeighbours[handle] = atom

            elif atomIsFromSideChain(atom):
                if atomIsFromLigand(atom):
                    if atomIsFromProtein(atom):
                        handle = "{0}{1}{2}".format(
                            atom.chain.upper().strip(), "_", atom.resid
                        )
                        if (
                            handle in self.backboneDonors.keys()
                            and handle not in self.sideChainDonors.keys()
                        ):
                            self.sideChainDonors[handle] = atom
                elif atomIsNeighbour(atom):
                    if atomIsFromProtein(atom):
                        handle = "{0}{1}{2}".format(
                            atom.chain.upper().strip(), "_", atom.resid
                        )
                        if (
                            handle in self.backboneNeighbours.keys()
                            and handle not in self.sideChainNeighbours.keys()
                        ):
                            self.sideChainNeighbours[handle] = atom

    def getBackboneAndSidechain(self):
        # self.backbone = self.backboneDonors.items() + self.backboneNeighbours.items() python2
        # self.sideChain = self.sideChainDonors.items() + self.sideChainNeighbours.items() python2
        self.backbone = [*self.backboneDonors.items(), *self.backboneNeighbours.items()]
        self.sideChain = [
            *self.sideChainDonors.items(),
            *self.sideChainNeighbours.items(),
        ]


# ===============================================================================
"""METHODS"""
# ===============================================================================


def multipleCoordination(site1, site2):
    if len(site1.donorAtoms) > 1 and len(site2.donorAtoms) > 1:
        return True
    else:
        return False


def singleCoordination(site1, site2):
    if len(site1.donorAtoms) == 1 and len(site2.donorAtoms) == 1:
        return True
    else:
        return False


def getCoordinatesOfAtomsInList(atmList):
    "Extracts coordinates of atoms"

    coordinates = [[atom.x, atom.y, atom.z] for atom in atmList]

    return coordinates


def copySite(site, tempSite):
    """
    Copies a given site to a tempSite

    Returns copied site
    """
    name = site.name

    metals = []
    for atom in site.metals:
        tempAtom = copy.copy(atom)
        tempAtom.x = atom.x
        tempAtom.y = atom.y
        tempAtom.z = atom.z
        metals.append(tempAtom)

    donorAtoms = []
    for atom in site.donorAtoms:
        tempAtom = copy.copy(atom)
        tempAtom.x = atom.x
        tempAtom.y = atom.y
        tempAtom.z = atom.z
        donorAtoms.append(tempAtom)

    approxSphereAtoms = []
    for atom in site.approxSphereAtoms:
        tempAtom = copy.copy(atom)
        tempAtom.x = atom.x
        tempAtom.y = atom.y
        tempAtom.z = atom.z
        approxSphereAtoms.append(tempAtom)

    tempSite.name = name

    tempSite.metals = metals
    tempSite.donorAtoms = donorAtoms
    tempSite.approxSphereAtoms = approxSphereAtoms

    tempSite.getAtomsForMatching()
    tempSite.getBackboneAndSidechain()

    return tempSite
