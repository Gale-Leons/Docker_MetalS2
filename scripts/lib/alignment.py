from operator import attrgetter

import numpy
from lib.fitting import Fitting
from lib.matching import matchAtoms
from lib.metalSite import MetalSite, copySite, getCoordinatesOfAtomsInList
from lib.reporting import ScoreReport
from lib.scoring import Scoring, scoreIsValid
from lib.transformation import (
    Rotations,
    applyRotationMatrix,
    applyRotationMatrix2,
    applyTransformationToSite,
)


def alignSites(site1, site2, maxDist=2.0, mode="pairwise", pose_threshold=40):
    maxDist = float(maxDist)

    scoreReport = ScoreReport(site1, site2)
    scoresObj = Scoring(site1, site2, maxDist, mode)

    rotations = Rotations(site1, site2)
    poses = rotations.getPoses()

    initial = scorePoses(site1, site2, poses, rotations.coordination, scoresObj)

    if not initial:
        raise Exception(
            "I could not align any of residues. Make sure you try to align sites of the same type of molecule or try to change the RMSD threshold.\n"
        )

    else:
        min_score_pose = min(initial, key=attrgetter("score"))
        scoreReport.min_initial_score = min_score_pose.score
        max_score_pose = max(initial, key=attrgetter("score"))
        scoreReport.max_initial_score = max_score_pose.score

        refineBeforeFitting(initial, pose_threshold)
        scoreReport.poses_selected = len(initial)
        fitInitialAlignments(site1, site2, initial, scoresObj, maxDist)

        min_init_score_pose = min(initial, key=attrgetter("score"))
        sitePositionInit = initial.index(min_init_score_pose)
        bestInitScore = min_init_score_pose.score

        min_opt_score_pose = min(initial, key=attrgetter("optScore"))
        sitePositionOpt = initial.index(min_opt_score_pose)
        bestOptScore = min_opt_score_pose.optScore

        if bestOptScore < bestInitScore:
            alignment = initial[sitePositionOpt]
            rotSite = applyRotationMatrix(site2, alignment.rotationMatrix)
            alignedSite = applyTransformationToSite(
                rotSite, alignment.optRotMatrix, alignment.optTranslVector
            )

        else:
            alignment = initial[sitePositionInit]
            alignedSite = applyRotationMatrix(site2, alignment.rotationMatrix)

        scoresObj.scoreSites(site1, alignedSite)
        scores = scoresObj.getScores()
        scoreReport.setScores(*scores)

        return alignedSite, scoreReport


def refineBeforeFitting(poses, percentage):
    from operator import attrgetter

    max_score_pose = max(poses, key=attrgetter("score"))
    min_score_pose = min(poses, key=attrgetter("score"))

    threshold = (
        min_score_pose.score
        + (max_score_pose.score - min_score_pose.score) * percentage / 100
    )

    poses[:] = [pose for pose in poses if pose.score <= threshold]

    return poses


def optimizeRMSD(site1, site, maxDist):
    pairs = matchAtoms(site1, site, maxDist)

    site1_mathed = []
    site2_mathed = []
    for pair in pairs:
        site1_mathed.append(pair[0])
        site2_mathed.append(pair[1])

    set1 = getCoordinatesOfAtomsInList(site1_mathed)
    set2 = getCoordinatesOfAtomsInList(site2_mathed)

    me = [0.0, 0.0, 0.0]
    set1.append(me)
    set2.append(me)

    x = numpy.array(set1, "f")
    y = numpy.array(set2, "f")

    sup = Fitting()
    sup.set(x, y)
    sup.run()
    rot, tran = sup.get_rotran()

    return rot, tran


def metalsAreClose(site1, site2, nuclearity):
    gcX_1, gcY_1, gcZ_1 = site1.getGeometricalCentre(site1.metals)
    gcX_2, gcY_2, gcZ_2 = site2.getGeometricalCentre(site2.metals)

    distance = numpy.square(
        numpy.linalg.norm(numpy.subtract([gcX_1, gcY_1, gcZ_1], [gcX_2, gcY_2, gcZ_2]))
    )

    if nuclearity == "mono":
        threshold = 1.75
    else:
        threshold = 10.0

    if distance < threshold:
        return True
    else:
        return False


def fitInitialAlignments(site1, site2, initial_alignments, scoresObj, maxDist):
    tempSite = MetalSite()

    if len(site1.metals) == 1 and len(site2.metals) == 1:
        nuclearity = "mono"
    else:
        nuclearity = "poly"

    for alignment in initial_alignments:
        tempSite = copySite(site2, tempSite)
        site = applyRotationMatrix2(site2, tempSite, alignment.rotationMatrix)

        rot, tran = optimizeRMSD(site1, site, maxDist)
        fitted_alignment = applyTransformationToSite(site, rot, tran)

        if not metalsAreClose(site1, fitted_alignment, nuclearity):
            alignment.optScore = alignment.score
            alignment.optRotMatrix = None
            alignment.optTranslVector = None

            continue

        else:
            scoresObj.scoreSites(site1, fitted_alignment)
            optimal_score = scoresObj.getTotalScore()
            current_score = alignment.score

            if current_score < optimal_score:
                alignment.optScore = current_score
                alignment.optRotMatrix = None
                alignment.optTranslVector = None

            else:
                alignment.optScore = optimal_score
                alignment.optRotMatrix = rot
                alignment.optTranslVector = tran

                while True:
                    current_score = optimal_score

                    rot, tran = optimizeRMSD(site1, fitted_alignment, maxDist)
                    fitted_alignment = applyTransformationToSite(
                        fitted_alignment, rot, tran
                    )
                    scoresObj.scoreSites(site1, fitted_alignment)
                    optimal_score = scoresObj.getTotalScore()

                    if (optimal_score == current_score) and metalsAreClose(
                        site1, fitted_alignment, nuclearity
                    ):
                        alignment.optScore = optimal_score
                        alignment.optRotMatrix = numpy.dot(alignment.optRotMatrix, rot)
                        alignment.optTranslVector = (
                            numpy.dot(alignment.optTranslVector, rot) + tran
                        )
                        break

                    elif (optimal_score < current_score) and metalsAreClose(
                        site1, fitted_alignment, nuclearity
                    ):
                        alignment.optScore = optimal_score
                        alignment.optRotMatrix = numpy.dot(alignment.optRotMatrix, rot)
                        alignment.optTranslVector = (
                            numpy.dot(alignment.optTranslVector, rot) + tran
                        )

                    else:
                        break


def scorePoses(site1, site2, poses, coordination, scores):
    if coordination == "multiple":
        for pose in poses:
            rotated_site = applyRotationMatrix(site2, pose.rotationMatrix)
            scores.scoreSites(site1, rotated_site)
            score = scores.getTotalScore()

            pose.score = score

    poses[:] = [pose for pose in poses if scoreIsValid(pose.score)]

    return poses
