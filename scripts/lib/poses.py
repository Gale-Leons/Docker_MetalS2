# ===============================================================================
"""CLASS"""
# ===============================================================================


class Pose:
    def __init__(self):
        self.rotationMatrix = None
        self.pattern1 = None
        self.pattern2 = None
        self.score = None
        self.site = None


class PoseWithLinearPatterns(Pose):
    
    def getAxis(self):
        axis = self.pattern1[0]
        return axis


class PoseWithPlanarPatterns(Pose):
    pass
