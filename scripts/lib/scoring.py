#!/usr/bin/env python2.6

import numpy
from lib.atom import atomIsFromBackbone, atomIsFromProtein
from lib.matching import matchAtoms
from lib.nomenclature import nonStandardAminoAcids, scoringMatrixBlosum62
from lib.sequence import penalizeSequence


class Scoring:
    def __init__(self, site1, site2, maxDist, mode):
        self.maxDist = maxDist
        self.maxS = self.getMaxSequenceScore(site1, site2, mode)
        self.N = self.getMaxNumberOfPairs(site1, site2, mode)

        self.lengthSite1 = len(site1.backbone) + len(site1.sideChain)
        self.lengthSite2 = len(site2.backbone) + len(site2.sideChain)

        self.poses_selected = None

        self.rmsd = None

        self.p1 = None
        self.p2 = None
        self.p3 = None

        self.totalScore = None
        self.percentage = None

    def getScores(self):
        return [self.rmsd, self.p1, self.p2, self.p3, self.totalScore, self.percentage]

    def getTotalScore(self):
        return self.totalScore

    def scoreSites(self, site1, site2):
        "Finds correspondence of atoms between the two sites and calculates score."

        self.pairs = matchAtoms(site1, site2, self.maxDist)

        if not self.pairs:
            self.totalScore = float("inf")
        else:
            w1 = numpy.float32(1.500)
            w2 = numpy.float32(1.000)
            w3 = numpy.float32(2.350)
            w4 = numpy.float32(0.001)

            sum_d = 0
            S = 0
            count = 0
            total_count = 0

            n = len(self.pairs)

            for pair in self.pairs:
                # RMSD
                d = numpy.square(
                    numpy.linalg.norm(
                        numpy.subtract(
                            [pair[0].x, pair[0].y, pair[0].z],
                            [pair[1].x, pair[1].y, pair[1].z],
                        )
                    )
                )
                sum_d = sum_d + d

                # Sequence similarity term
                atom1 = pair[0]
                atom2 = pair[1]

                if atomIsFromBackbone(atom1) and atomIsFromBackbone(atom2):
                    total_count += 1

                    if atomIsFromProtein(atom1) and atomIsFromProtein(atom2):
                        aa1 = atom1.aa.upper().strip()
                        aa2 = atom2.aa.upper().strip()
                        aa_score = scoreAminoAcids(aa1, aa2)
                        S = S + aa_score
                        if aa1 == aa2:
                            count += 1

            self.rmsd = numpy.sqrt(sum_d / n)

            # Relative coverage term
            sizePart = numpy.log(numpy.float32(self.N) / numpy.float32(n))
            self.p1 = w1 * sizePart

            # Chemical term
            chemPart = 1 - numpy.float32(S) / numpy.float32(self.maxS)
            self.p2 = w2 * chemPart

            # Continuity term
            P, N_p = penalizeSequence(self.pairs)
            self.p3 = w3 * numpy.float32(P) / numpy.float32(N_p)

            # rmsd
            self.p4 = w4 * self.rmsd

            # Total
            self.totalScore = self.p1 + self.p2 + self.p3 + self.p4

            # Percentage
            self.percentage = (float(count) / total_count) * 100

    def getMaxNumberOfPairs(self, site1, site2, mode):
        if mode == "pairwise":
            N = len(site1.backbone) + len(site1.sideChain)
            if N > len(site2.backbone) + len(site2.sideChain):
                N = len(site2.backbone) + len(site2.sideChain)
        elif mode == "database":
            N = len(site1.backbone) + len(site1.sideChain)

        return N

    def getMaxSequenceScore(self, site1, site2, mode):
        S_max = 0
        if mode == "pairwise":
            backbone = site1.backbone
            if len(site1.backbone) >= len(site2.backbone):
                backbone = site2.backbone

        elif mode == "database":
            backbone = site1.backbone

        for item in backbone:
            atom = item[1]
            residueCode = atom.aa.upper().strip()

            if atomIsFromProtein(atom):
                aa_score = scoreAminoAcids(residueCode, residueCode)
                S_max = S_max + aa_score

        return S_max


# ===============================================================================
"""METHODS"""
# ===============================================================================


def scoreAminoAcids(aa1, aa2):
    if (aa1 in scoringMatrixBlosum62.keys()) and (aa2 in scoringMatrixBlosum62.keys()):
        aa_score = (scoringMatrixBlosum62.get(aa1)).get(aa2)

    elif (aa1 in nonStandardAminoAcids.keys()) and (
        aa2 in scoringMatrixBlosum62.keys()
    ):
        aa1_sub = nonStandardAminoAcids[aa1]
        aa_score = (scoringMatrixBlosum62.get(aa1_sub)).get(aa2)

    elif (aa2 in nonStandardAminoAcids.keys()) and (
        aa1 in scoringMatrixBlosum62.keys()
    ):
        aa2_sub = nonStandardAminoAcids[aa2]
        aa_score = (scoringMatrixBlosum62.get(aa1)).get(aa2_sub)

    elif (aa1 in nonStandardAminoAcids.keys()) and (
        aa2 in nonStandardAminoAcids.keys()
    ):
        aa1_sub = nonStandardAminoAcids[aa1]
        aa2_sub = nonStandardAminoAcids[aa2]
        aa_score = (scoringMatrixBlosum62.get(aa1_sub)).get(aa2_sub)
    else:
        aa_score = 0

    return aa_score


def scoreIsValid(score):
    # Score equal to zero may exist for two equal sites
    if score is not None and score != float("inf"):
        return True
    else:
        return False
