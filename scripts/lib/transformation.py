import copy

import numpy
from lib.metalSite import MetalSite, getCoordinatesOfAtomsInList, multipleCoordination
from lib.poses import PoseWithPlanarPatterns

# ===============================================================================
"""CONSTANTS"""
# ===============================================================================

# Epsilon for testing whether a number == close to zero
_EPS = numpy.finfo(float).eps * 4.0

# ===============================================================================
"""CLASS"""
# ===============================================================================


class Rotations:
    """
    Rotations are defined as all rigid body movements minimizing a value of rmsd between each pair of local patterns
    Rotations are represented by rotation matrices derived from superposition of pars of local patterns
    """

    def __init__(self, site1, site2):
        self.site1 = site1
        self.site2 = site2

        if multipleCoordination(self.site1, self.site2):
            self.coordination = "multiple"
        else:
            self.coordination = "single"

    def getPoses(self):
        self.getLocalPatterns()

        poses = []

        for query_pattern in self.patterns1:
            query_set = getCoordinatesOfAtomsInList(query_pattern)

            for target_pattern in self.patterns2:
                if self.coordination == "multiple":
                    pose = PoseWithPlanarPatterns()
                else:
                    # pose = PoseWithLinearPatterns()
                    raise Exception

                target_set = getCoordinatesOfAtomsInList(target_pattern)
                transformationMatrix = findBestTransformation(query_set, target_set)

                pose.rotationMatrix = extractRotationMatrix(transformationMatrix)
                pose.pattern1 = query_set
                pose.pattern2 = target_set

                poses.append(pose)

        return poses

    def getLocalPatterns(self):
        if multipleCoordination(self.site1, self.site2):
            # Both sites have more than one donor atom
            triangles1 = self.findPlanarLocalPatterns(self.site1)
            triangles2 = self.findPlanarLocalPatterns(self.site2)

            permutations = self.generatePermutations(triangles2)

            self.patterns1 = triangles1
            self.patterns2 = triangles2 + permutations

        else:
            # At least one site is a single donor site
            self.patterns1 = self.findLinearLocalPatterns(self.site1)
            self.patterns2 = self.findLinearLocalPatterns(self.site2)

    def generatePermutations(self, triangles):
        """
        Explicitly sets the permutations for triangles
        Takes a triangle: a list of atoms - > [atomA', atomB']
        Returns the pair of atoms reversed -> [atomB', atomA']
        """

        permutations = []
        for triangle in triangles:
            permutations.append([triangle[1], triangle[0]])

        return permutations

    def findPlanarLocalPatterns(self, site):
        """
        Composes a metal atom and two donor atoms in all possible triangles originated from the metal
        site: object of MetalSite class
        Returns a list of atom pairs -> [[donor_atom1,donor_atom2],[donor_atom1, donor_atom2],...]
        """

        localPatterns = []

        for i in range(0, len(site.donorAtoms) - 1):
            for j in range(i + 1, len(site.donorAtoms)):
                localPatterns.append([site.donorAtoms[i], site.donorAtoms[j]])

        return localPatterns


# ===============================================================================
"""METHODS"""
# ===============================================================================


def extractRotationMatrix(M):
    """
    Extracts a rotation matrix from a tranformation matrix

    M: tranformation matrix [4x4]

    Returns a rotation matrix [3x3]
    """

    rotM = [numpy.take(M[i], [0, 1, 2]) for i in range(0, 3)]

    return rotM


def applyRotationMatrix(site, rotM):
    """
    Rotates a site using a rotation matrix M

    site: MetalSite object
    rotM: rotation matrix [3x3]

    Returns a site with changed coordinates (creates a new site object)
    """

    rotSite = MetalSite()

    name = site.name

    metals = []
    for atom in site.metals:
        v = [atom.x, atom.y, atom.z]
        v = numpy.array(v, dtype=numpy.float64, copy=True)
        rot_v = numpy.dot(rotM, v)
        tempAtom = copy.copy(atom)
        tempAtom.x = rot_v[0]
        tempAtom.y = rot_v[1]
        tempAtom.z = rot_v[2]
        metals.append(tempAtom)

    donorAtoms = []
    for atom in site.donorAtoms:
        v = [atom.x, atom.y, atom.z]
        v = numpy.array(v, dtype=numpy.float64, copy=True)
        rot_v = numpy.dot(rotM, v)
        tempAtom = copy.copy(atom)
        tempAtom.x = rot_v[0]
        tempAtom.y = rot_v[1]
        tempAtom.z = rot_v[2]
        donorAtoms.append(tempAtom)

    approxSphereAtoms = []
    for atom in site.approxSphereAtoms:
        v = [atom.x, atom.y, atom.z]
        v = numpy.array(v, dtype=numpy.float64, copy=True)
        rot_v = numpy.dot(rotM, v)
        tempAtom = copy.copy(atom)
        tempAtom.x = rot_v[0]
        tempAtom.y = rot_v[1]
        tempAtom.z = rot_v[2]
        approxSphereAtoms.append(tempAtom)

    rotSite.name = name
    rotSite.metals = metals
    rotSite.donorAtoms = donorAtoms
    rotSite.approxSphereAtoms = approxSphereAtoms

    rotSite.getAtomsForMatching()
    rotSite.getBackboneAndSidechain()

    return rotSite


def applyRotationMatrix2(site, tempSite, rotM):
    """
    Rotates a site using a rotation matrix M

    site: MetalSite object
    rotM: rotation matrix [3x3]

    Returns a site with changed coordinates (changes tempSite)
    """

    v = numpy.array([0, 0, 0], dtype=numpy.float64, copy=True)
    i = 0
    for atom in site.metals:
        v[0] = atom.x
        v[1] = atom.y
        v[2] = atom.z

        rot_v = numpy.dot(rotM, v)
        tempSite.metals[i].x = rot_v[0]
        tempSite.metals[i].y = rot_v[1]
        tempSite.metals[i].z = rot_v[2]
        i = i + 1

    i = 0
    for atom in site.donorAtoms:
        v[0] = atom.x
        v[1] = atom.y
        v[2] = atom.z
        rot_v = numpy.dot(rotM, v)
        tempSite.donorAtoms[i].x = rot_v[0]
        tempSite.donorAtoms[i].y = rot_v[1]
        tempSite.donorAtoms[i].z = rot_v[2]
        i = i + 1

    i = 0
    for atom in site.approxSphereAtoms:
        v[0] = atom.x
        v[1] = atom.y
        v[2] = atom.z
        rot_v = numpy.dot(rotM, v)
        tempSite.approxSphereAtoms[i].x = rot_v[0]
        tempSite.approxSphereAtoms[i].y = rot_v[1]
        tempSite.approxSphereAtoms[i].z = rot_v[2]
        i = i + 1

    return tempSite


def findBestTransformation(v0, v1):
    """
    Takes as input a pair of triangles and finds superpositions between the two, optimizing a value of RMSD
    v0: vector of coordinates
    v1: vector of coordinates

    Returns a transformation matrix [4x4]
    """

    v0 = numpy.array(v0, dtype=numpy.float64, copy=True)
    v1 = numpy.array(v1, dtype=numpy.float64, copy=True)

    # [A'-A:[x'-x, y'-y, z'-z], B'-B:[x'-x, y'-y, z'-z]]
    delta_m = v1 - v0
    # [A'+A:[x'+x, y'+y, z'+z], B'+B:[x'+x, y'+y, z'+z]]
    delta_p = v1 + v0

    # [(A'-A)^2:[(x'-x)^2, (y'-y)^2, (z'-z)^2], (B'-B)^2:[(x'-x)^2, (y'-y)^2, (z'-z)^2]]
    sq_delta_m = delta_m * delta_m
    # [(A'+A)^2:[(x'+x)^2, (y'+y)^2, (z'+z)^2], (B'+B)^2:[(x'+x)^2, (y'+y)^2, (z'+z)^2]]
    sq_delta_p = delta_p * delta_p

    xx_m, yy_m, zz_m = numpy.sum(sq_delta_m, axis=0)
    xx_p, yy_p, zz_p = numpy.sum(sq_delta_p, axis=0)

    xm_ym, ym_zm, zm_xm = numpy.sum(delta_m * numpy.roll(delta_m, -1, axis=1), axis=0)
    xp_yp, yp_zp, zp_xp = numpy.sum(delta_p * numpy.roll(delta_p, -1, axis=1), axis=0)

    xm_yp, ym_zp, zm_xp = numpy.sum(delta_m * numpy.roll(delta_p, -1, axis=1), axis=0)
    xm_zp, ym_xp, zm_yp = numpy.sum(delta_m * numpy.roll(delta_p, -2, axis=1), axis=0)

    # Computes symmetric matrix N
    N = [
        [xx_m + yy_m + zz_m, zm_yp - ym_zp, xm_zp - zm_xp, ym_xp - xm_yp],
        [zm_yp - ym_zp, yy_p + zz_p + xx_m, xm_ym - xp_yp, zm_xm - zp_xp],
        [xm_zp - zm_xp, xm_ym - xp_yp, xx_p + zz_p + yy_m, ym_zm - yp_zp],
        [ym_xp - xm_yp, zm_xm - zp_xp, ym_zm - yp_zp, xx_p + yy_p + zz_m],
    ]

    # Finds eignvalues and eignvectors
    w, V = numpy.linalg.eigh(N)

    # Gives an eigenvector corresponding to the smallest positive eigenvalue
    q = V[:, numpy.argmin(w)]

    # Unit quaternion
    q /= vectorNorm(q)

    # Homogeneous transformation matrix
    M = matrixOutOfQuaternion(q)

    return M


def vectorNorm(data, axis=None, out=None):
    """
    Returns length, i.e. eucledian norm, of ndarray along axis.
    """

    data = numpy.array(data, dtype=numpy.float64, copy=True)
    if out is None:
        if data.ndim == 1:
            return numpy.sqrt(numpy.dot(data, data))
        data *= data
        out = numpy.atleast_1d(numpy.sum(data, axis=axis))
        numpy.sqrt(out, out)
        return out
    else:
        data *= data
        numpy.sum(data, axis=axis, out=out)
        numpy.sqrt(out, out)


def matrixOutOfQuaternion(quaternion):
    """
    Returns a homogeneous transformation matrix from quaternion.

    quaternion: quaternion vector [1x4]

    Returns a transformation matrix [4x4]
    """

    q = numpy.array(quaternion, dtype=numpy.float64, copy=True)
    n = numpy.dot(q, q)
    if n < _EPS:
        return numpy.identity(4)
    q *= numpy.sqrt(2.0 / n)
    q = numpy.outer(q, q)

    return numpy.array(
        [
            [1.0 - q[2, 2] - q[3, 3], q[1, 2] - q[3, 0], q[1, 3] + q[2, 0], 0.0],
            [q[1, 2] + q[3, 0], 1.0 - q[1, 1] - q[3, 3], q[2, 3] - q[1, 0], 0.0],
            [q[1, 3] - q[2, 0], q[2, 3] + q[1, 0], 1.0 - q[1, 1] - q[2, 2], 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )


def applyTransformationToSite(site, rotM, tranM):
    for atom in site.metals:
        v = [atom.x, atom.y, atom.z]
        v = numpy.array(v, dtype=numpy.float64, copy=True)
        transf_v = numpy.dot(v, rotM) + tranM

        atom.x = transf_v[0]
        atom.y = transf_v[1]
        atom.z = transf_v[2]

    for atom in site.donorAtoms:
        v = [atom.x, atom.y, atom.z]
        v = numpy.array(v, dtype=numpy.float64, copy=True)
        transf_v = numpy.dot(v, rotM) + tranM

        atom.x = transf_v[0]
        atom.y = transf_v[1]
        atom.z = transf_v[2]

    for atom in site.approxSphereAtoms:
        v = [atom.x, atom.y, atom.z]
        v = numpy.array(v, dtype=numpy.float64, copy=True)
        transf_v = numpy.dot(v, rotM) + tranM

        atom.x = transf_v[0]
        atom.y = transf_v[1]
        atom.z = transf_v[2]

    site.getAtomsForMatching()
    site.getBackboneAndSidechain()

    return site
