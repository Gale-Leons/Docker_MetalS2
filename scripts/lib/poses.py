# ===============================================================================
"""CLASS"""
# ===============================================================================


class Pose(object):
    def __init__(self):
        self.rotationMatrix = None
        self.pattern1 = None
        self.pattern2 = None
        self.score = None
        self.site = None


class PoseWithLinearPatterns(Pose):
    pass


class PoseWithPlanarPatterns(Pose):
    pass
