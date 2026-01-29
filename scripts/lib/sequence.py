#!/usr/bin/env python

import copy

from lib.atom import (
    atomIsFromBackbone,
    atomIsFromLigand,
    atomIsFromProtein,
    getOneLetterResudueCode,
)
from lib.matching import matchAtoms


def getPenalties(atomsList):
    P = 0
    N = len(atomsList)

    atom1_1 = atomsList[0][0]
    atom2_1 = atomsList[0][1]

    id1 = atom1_1.resid
    id2 = atom2_1.resid

    n = 1

    for i in range(1, len(atomsList)):
        atom1_i = atomsList[i][0]
        atom2_i = atomsList[i][1]

        if ((atom1_i.resid != id1 + 1) and (atom1_i.resid != id1 - 1)) or (
            (atom2_i.resid != id2 + 1) and (atom2_i.resid != id2 - 1)
        ):
            P = P + 1.0 / n
            n = 1

        else:
            n += 1

        id1 = atom1_i.resid
        id2 = atom2_i.resid

    P = P + 1.0 / n

    return P, N


def penalizeSequence(pairsOfAtoms):
    reducedPairsOfAtoms = eliminateSideChains(pairsOfAtoms)
    atomsList = quicksortPairsByChainThenID(reducedPairsOfAtoms, 0)

    atomsListProtein = []
    atomsListDNAorRNA = []
    for atomPair in atomsList:
        if atomIsFromProtein(atomPair[0]):
            atomsListProtein.append(atomPair)

    P = 0
    N = 0
    if atomsListProtein:
        P_prot, N_prot = getPenalties(atomsListProtein)
        P += P_prot
        N += N_prot
    if atomsListDNAorRNA:
        P_nucl, N_nucl = getPenalties(atomsListDNAorRNA)
        P += P_nucl
        N += N_nucl

    return P, N


def writeSequenceToFile(site1, site2, outputPath, maxDist):
    pairsOfAtoms = matchAtoms(site1, site2, maxDist)
    reducedPairsOfAtoms = eliminateSideChains(pairsOfAtoms)

    atomsList = quicksortPairsByChainThenID(reducedPairsOfAtoms, 0)

    atomsListProtein = []
    atomsListDNAorRNA = []
    for atomPair in atomsList:
        if atomIsFromProtein(atomPair[0]):
            atomsListProtein.append(atomPair)

    sequence_list = []
    sequence_list_prot = []
    sequence_list_nucl = []

    outFile = open(outputPath, "w")
    #
    querySiteName = "{:<7}{:<4}".format("Query:", site1.name)
    sbjctSiteName = "{:<7}{:<4}".format("Sbjct:", site2.name)
    #
    outFile.writelines(querySiteName)
    outFile.writelines("\n")
    outFile.writelines(sbjctSiteName)
    outFile.writelines("\n")

    if atomsListProtein:
        outFile.writelines("\n")
        outFile.writelines("Protein alignment:\n")

        sequence_list = writeLinesOfSequenceAlignment(
            atomsListProtein, outFile, sequence_list
        )
        sequence_list_prot = copy.deepcopy(sequence_list)
        sequence_list = []

    if atomsListDNAorRNA:
        outFile.writelines("\n")
        outFile.writelines("DNA/RNA alignment:\n")

        sequence_list = writeLinesOfSequenceAlignment(
            atomsListDNAorRNA, outFile, sequence_list
        )
        sequence_list_nucl = copy.deepcopy(sequence_list)
    outFile.close()

    return sequence_list_prot, sequence_list_nucl


def getLinesOfSequenceAlignment(atomsList):
    #
    atom1_1 = atomsList[0][0]
    chain1 = atom1_1.chain
    id1 = atom1_1.resid
    aa1 = getOneLetterResudueCode(atom1_1)
    #
    atom2_1 = atomsList[0][1]
    chain2 = atom2_1.chain
    id2 = atom2_1.resid
    aa2 = getOneLetterResudueCode(atom2_1)
    #

    #
    chainList1 = []
    chainList1.append(f"{chain1:<4}")
    idList1 = []
    idList1.append(f"{id1:<4}")
    aaList1 = []
    if atomIsFromLigand(atom1_1):
        aaList1.append("{:<1}{:<3}".format(aa1, "*"))
    else:
        aaList1.append(f"{aa1:<4}")
    #
    pipeList = []
    #
    chainList2 = []
    chainList2.append(f"{chain2:<4}")
    idList2 = []
    idList2.append(f"{id2:<4}")
    aaList2 = []
    if atomIsFromLigand(atom2_1):
        aaList2.append("{:<1}{:<3}".format(aa2, "*"))
    else:
        aaList2.append(f"{aa2:<4}")
    #

    if aaList1[0] == aaList2[0]:
        pipeList.append("|  ")
    else:
        pipeList.append("   ")

    k = 0
    for i in range(1, len(atomsList)):
        atom1_i = atomsList[i][0]
        atom2_i = atomsList[i][1]

        if (
            getOneLetterResudueCode(atom1_i) is not None
            and getOneLetterResudueCode(atom2_i) is not None
        ):
            flag = False
            k = k + 1
            # The first atom from a pair
            if atom1_i.chain != chain1:
                chainList1.append(" | ")
                idList1.append(" | ")
                aaList1.append(" | ")
                #
                pipeList.append("   ")
                #
                chainList2.append(" | ")
                idList2.append(" | ")
                aaList2.append(" | ")
                chainList1.append(f"{atom1_i.chain:<4}")
                chain1 = atom1_i.chain
                flag = True
                k = k + 1
            else:
                chainList1.append(" " * 4)

            if (atom1_i.resid != id1 + 1) and (atom1_i.resid != id1 - 1):
                if flag is False:
                    temp = chainList1[k]
                    chainList1[k] = "   "
                    chainList1.append(temp)
                    idList1.append(" | ")
                    aaList1.append(" | ")
                    #
                    pipeList.append("   ")
                    #
                    chainList2.append("   ")
                    idList2.append(" | ")
                    aaList2.append(" | ")
                    k = k + 1
                    flag = True

                idList1[k - 2] = f"{atomsList[i - 1][0].resid:<4}"
                idList1.append(f"{atom1_i.resid:<4}")
                id1 = atom1_i.resid
            else:
                idList1.append(f"{atom1_i.resid:<4}")
                id1 = atom1_i.resid

            if atomIsFromLigand(atom1_i):
                aaList1.append(
                    "{:<1}{:<3}".format(getOneLetterResudueCode(atom1_i), "*")
                )
            else:
                aaList1.append(f"{getOneLetterResudueCode(atom1_i):<4}")

            # The second atom from a pair
            if atom2_i.chain != chain2:
                # If chain is changed in second sequence, pipe is inserted in both sequenses and first sequence is shifted
                if flag is False:
                    temp = chainList1[k]
                    chainList1[k] = " | "
                    chainList1.append(temp)
                    temp = idList1[k]
                    idList1[k] = " | "
                    idList1.append(temp)
                    temp = aaList1[k]
                    aaList1[k] = " | "
                    aaList1.append(temp)
                    chainList2.append(" | ")
                    idList2.append(" | ")
                    aaList2.append(" | ")
                    pipeList.append("   ")
                    k = k + 1
                    flag = True
                else:
                    chainList1[k - 1] = " | "
                    chainList2[k - 1] = " | "

                chainList2.append(f"{atom2_i.chain:<4}")
                chain2 = atom2_i.chain
            else:
                chainList2.append(" " * 4)

            aaTemp = aaList1[k]
            if (atom2_i.resid != id2 + 1) and (atom2_i.resid != id2 - 1):
                if flag is False:
                    temp = chainList1[k]
                    chainList1[k] = "   "
                    chainList1.append(temp)
                    temp = idList1[k]
                    idList1[k] = " | "
                    idList1.append(temp)
                    temp = aaList1[k]
                    aaList1[k] = " | "
                    aaList1.append(temp)
                    temp = chainList2[k]
                    chainList2[k] = "   "
                    chainList2.append(temp)
                    idList2.append(" | ")
                    aaList2.append(" | ")
                    pipeList.append("   ")
                    k = k + 1
                    flag = True

                idList2[k - 2] = f"{atomsList[i - 1][1].resid:<4}"
                idList2.append(f"{atom2_i.resid:<4}")
                id2 = atom2_i.resid
            else:
                idList2.append(f"{atom2_i.resid:<4}")
                id2 = atom2_i.resid

            if atomIsFromLigand(atom2_i):
                aaList2.append(
                    "{:<1}{:<3}".format(getOneLetterResudueCode(atom2_i), "*")
                )
            else:
                aaList2.append(f"{getOneLetterResudueCode(atom2_i):<4}")

            if aaTemp == aaList2[-1]:
                pipeList.append(" |  ")
            else:
                pipeList.append("    ")

    return chainList1, idList1, aaList1, chainList2, idList2, aaList2, pipeList


def writeLinesOfSequenceAlignment(atomsList, outFile, sequence_list):
    (
        chainList1,
        idList1,
        aaList1,
        chainList2,
        idList2,
        aaList2,
        pipeList,
    ) = getLinesOfSequenceAlignment(atomsList)
    #

    sequence_list.append(
        (chainList1, idList1, aaList1, pipeList, aaList2, idList2, chainList2)
    )

    outFile.writelines("\n")

    chainList1.insert(0, "Chain:  | ")
    outFile.writelines(chainList1)
    outFile.writelines("\n")
    idList1.insert(0, "ResID:  | ")
    outFile.writelines(idList1)
    outFile.writelines("\n")
    aaList1.insert(0, "Query:  | ")
    outFile.writelines(aaList1)
    outFile.writelines("\n")

    pipeList.insert(0, "          ")
    outFile.writelines(pipeList)
    outFile.writelines("\n")

    aaList2.insert(0, "Target: | ")
    outFile.writelines(aaList2)
    outFile.writelines("\n")
    idList2.insert(0, "ResID:  | ")
    outFile.writelines(idList2)
    outFile.writelines("\n")
    chainList2.insert(0, "Chain:  | ")
    outFile.writelines(chainList2)
    outFile.writelines("\n")

    return sequence_list


def eliminateSideChains(listOfAtomPairs):
    "Takes a list of all pairs and remains carbon alpha of protein residues and C1 of DNA/RNA residues"

    reducedListOfAtomPairs = [
        pairOfAtoms
        for pairOfAtoms in listOfAtomPairs
        if atomIsFromBackbone(pairOfAtoms[0])
    ]

    return reducedListOfAtomPairs


# LIST OF PAIRS SORTING


def quicksortPairsByChainThenID(pairList, pairID):
    chainIDsortedPairs = quicksortPairsByChain(pairList, pairID)

    blocksOfChains = []
    chainBlock = []

    prev_item = chainIDsortedPairs[0]
    for next_item in chainIDsortedPairs:
        if prev_item[pairID].chain == next_item[pairID].chain:
            chainBlock.append(next_item)
            prev_item = next_item
        else:
            sequence = quicksortPairsByResID(chainBlock, pairID)
            blocksOfChains.append(sequence)
            chainBlock = []
            chainBlock.append(next_item)
            prev_item = next_item

    sequence = quicksortPairsByResID(chainBlock, pairID)
    blocksOfChains.append(sequence)

    pairsSortedByChainThenID = []
    for block in blocksOfChains:
        pairsSortedByChainThenID.extend(block)

    return pairsSortedByChainThenID


def quicksortPairsByResID(pairList, pairID):
    "Takes a list of atom pairs and sorts all atoms by residue ID of atoms with a given ordinal number"

    if len(pairList) <= 1:
        return pairList

    pivot = pairList.pop()
    before = [pair for pair in pairList if pair[pairID].resid <= pivot[pairID].resid]
    after = [pair for pair in pairList if pair[pairID].resid > pivot[pairID].resid]

    return (
        quicksortPairsByResID(before, pairID)
        + [pivot]
        + quicksortPairsByResID(after, pairID)
    )


def quicksortPairsByChain(pairList, pairID):
    if len(pairList) <= 1:
        return pairList

    pivot = pairList.pop()
    before = [pair for pair in pairList if pair[pairID].chain <= pivot[pairID].chain]
    after = [pair for pair in pairList if pair[pairID].chain > pivot[pairID].chain]

    return (
        quicksortPairsByChain(before, pairID)
        + [pivot]
        + quicksortPairsByChain(after, pairID)
    )
