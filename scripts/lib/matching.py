#!/usr/bin/env python

import math

import numpy
from lib.atom import bothAtomsAreFromSameMolecule  # atomIsFromBackbone,
from scipy.spatial import cKDTree


def matchAtoms(site1, site2, maxDist):
    """
    For an each CA/C1 atom from a site1 finds the closest neighbor in a kd-tree generated based on a site2

    Returns a list of the closest neighbors and the biggest possible number of pairs
    """

    backbone1 = site1.backboneDonors.copy()
    sidechain1 = site1.sideChainDonors.copy()
    backbone2 = site2.backboneDonors.copy()
    sidechain2 = site2.sideChainDonors.copy()

    (
        donors_pairs,
        res_backbone1,
        res_sidechain1,
        res_backbone2,
        res_sidechain2,
    ) = getNearestNeighborsInPairs(backbone1, sidechain1, backbone2, sidechain2, 15.0)

    backbone1 = site1.backboneNeighbours.copy()
    if res_backbone1:
        backbone1.update(res_backbone1)
    sidechain1 = site1.sideChainNeighbours.copy()
    if res_sidechain1:
        sidechain1.update(res_sidechain1)
    backbone2 = site2.backboneNeighbours.copy()
    if res_backbone2:
        backbone2.update(res_backbone2)
    sidechain2 = site2.sideChainNeighbours.copy()
    if res_sidechain2:
        sidechain2.update(res_sidechain2)

    (
        neighbours_pairs,
        res_backbone1,
        res_sidechain1,
        res_backbone2,
        res_sidechain2,
    ) = getNearestNeighborsInPairs(
        backbone1, sidechain1, backbone2, sidechain2, maxDist
    )

    pairs = donors_pairs + neighbours_pairs

    return pairs


def getNearestNeighborsInPairs(backbone1, sidechain1, backbone2, sidechain2, maxDist):
    (
        backbone1Coords,
        sidechain1Coords,
        backbone2Coords,
        sidechain2Coords,
    ) = getCoordinatesForNNS(backbone1, sidechain1, backbone2, sidechain2)

    keys = findNearestNeighbors(backbone1Coords, backbone2Coords, maxDist)
    pairs = getAtomsByKeys(keys, backbone1, sidechain1, backbone2, sidechain2, maxDist)

    res_backbone1 = backbone1.copy()
    res_sidechain1 = sidechain1.copy()
    res_backbone2 = backbone2.copy()
    res_sidechain2 = sidechain2.copy()
    for pair_keys in keys:
        key1 = pair_keys[0]
        key2 = pair_keys[1]
        res_backbone1.pop(key1, None)
        res_sidechain1.pop(key1, None)
        res_backbone2.pop(key2, None)
        res_sidechain2.pop(key2, None)

    return pairs, res_backbone1, res_sidechain1, res_backbone2, res_sidechain2


def getAtomsByKeys(keys, source1bb, source1sc, source2bb, source2sc, maxDist):
    atomPairs = []

    for atomsKeys in keys:
        atom1_bb = source1bb[atomsKeys[0]]
        atom2_bb = source2bb[atomsKeys[1]]

        # check if residues belong to the same molecul type to avoid matching amino acids with nucleic acids
        if bothAtomsAreFromSameMolecule(atom1_bb, atom2_bb):
            atomPairs.append([atom1_bb, atom2_bb])

            # check if side chain atoms have a pair
            if (
                source1sc.get(atomsKeys[0]) is not None
                and source2sc.get(atomsKeys[1]) is not None
            ):
                atom1_sc = source1sc[atomsKeys[0]]
                atom2_sc = source2sc[atomsKeys[1]]
                dx = atom1_sc.x - atom2_sc.x
                dy = atom1_sc.x - atom2_sc.y
                dz = atom1_sc.z - atom2_sc.z
                distance = math.sqrt(dx * dx + dy * dy + dz * dz)
                if distance < maxDist:
                    atomPairs.append([atom1_sc, atom2_sc])

    return atomPairs


def findNearestNeighbors(dict1, dict2, maxDist):
    """
    Takes two dictionaries as an input, represents second dictionary as a kd-tree, and identifies the nearest neighbors points in target dictionary to the query.

    Returns:
    --------
    keyPairs: a list of paired keys from two dictionaries
    """

    keys1 = list(dict1.keys())
    set1 = numpy.array(list(dict1.values()), dtype=float)
    keys2 = list(dict2.keys())
    set2 = numpy.array(list(dict2.values()), dtype=float)

    keyPairs = []
    setPairs = []

    while (len(set1) > 0) and (len(set2) > 0):
        set2tree = cKDTree(set2)
        dist, ind2 = set2tree.query(set1, k=1, eps=0, p=2, distance_upper_bound=maxDist)

        ind1 = dist.argmin()

        if dist[ind1] > maxDist:
            break

        keyPairs.append([keys1[ind1], keys2[ind2[ind1]]])
        setPairs.append([set1[ind1], set2[ind2[ind1]]])

        keys1.pop(ind1)
        set1 = numpy.delete(set1, ind1, axis=0)
        keys2.pop(ind2[ind1])
        set2 = numpy.delete(set2, ind2[ind1], axis=0)

    return keyPairs


def getCoordinatesForNNS(backbone1, sidechain1, backbone2, sidechain2):
    # Find pairs among residues
    backbone1coords = {}
    sidechain1coords = {}
    backbone2coords = {}
    sidechain2coords = {}

    # Derive coordinates of backbone and side-chain atoms

    # Site 1 -> backbone
    for key in backbone1:
        elem = backbone1[key]
        backbone1coords[key] = [elem.x, elem.y, elem.z]

    # Site 1 -> side chain
    for key in sidechain1:
        elem = sidechain1[key]
        sidechain1coords[key] = [elem.x, elem.y, elem.z]

    # Site 2 -> backbone
    for key in backbone2:
        elem = backbone2[key]
        backbone2coords[key] = [elem.x, elem.y, elem.z]

    # Site 2 -> side chain
    for key in sidechain2:
        elem = sidechain2[key]
        sidechain2coords[key] = [elem.x, elem.y, elem.z]

    return backbone1coords, sidechain1coords, backbone2coords, sidechain2coords
