from numpy import dot, transpose
from numpy.linalg import det, svd


class Fitting(object):
    def __init__(self):
        self._clear()

    def _clear(self):
        self.reference_coords = None
        self.coords = None
        self.transformed_coords = None
        self.rot = None
        self.tran = None
        self.rms = None
        self.init_rms = None

    def set(self, reference_coords, coords):
        """
        Set the coordinates to be superimposed.
        coords will be put on top of reference_coords.

        o reference_coords: an NxDIM array
        o coords: an NxDIM array

        DIM is the dimension of the points, N is the number
        of points to be superimposed.
        """
        # clear everything from previous runs
        self._clear()
        # store cordinates
        self.reference_coords = reference_coords
        self.coords = coords
        n = reference_coords.shape
        m = coords.shape
        if n != m or not (n[1] == m[1] == 3):
            raise Exception("Coordinate number/dimension mismatch.")
        self.n = n[0]

    def run(self):
        "Superimpose the coordinate sets."
        if self.coords is None or self.reference_coords is None:
            raise Exception("No coordinates set.")
        coords = self.coords
        reference_coords = self.reference_coords
        # center on centroid
        av1 = sum(coords) / self.n
        av2 = sum(reference_coords) / self.n
        coords = coords - av1
        reference_coords = reference_coords - av2
        # correlation matrix
        a = dot(transpose(coords), reference_coords)
        u, d, vt = svd(a)
        self.rot = transpose(dot(transpose(vt), transpose(u)))
        # check if we have found a reflection
        if det(self.rot) < 0:
            vt[2] = -vt[2]
            self.rot = transpose(dot(transpose(vt), transpose(u)))
        self.tran = av2 - dot(av1, self.rot)

    def get_rotran(self):
        "Right multiplying rotation matrix and translation."
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        return self.rot, self.tran
