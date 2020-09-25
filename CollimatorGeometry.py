#! usr/bin/env python
import math

"""
Class for describing the geometry of the spot-size collimator.
You can set the values of the parameters yourself to make a custom geometry,
or use the methods useStraightHoleGeometry() or useTippedHoleGeometry() for 
pre-buit geometries.
"""


class CollimatorGeometry:

    # class variables
    # stopping target position in Mu2e coords given in Mu2e-doc-26586
    inchToCm = 2.54
    stoppingTargetCenterPosition = [0.0, 0.0, -1*(4064.0-(627.1+547.1)/2)]
    stoppingTargetLength = 80.0
    stoppingTargetRadius = 7.5

    def __init__(self):

        # units are CENTIMETERS!
        self.collimatorThickness = None
        self.leftHoleCenterFrontCoords = None
        self.leftHoleRadius = None
        self.rightHoleCenterFrontCoords = None
        self.rightHoleRadius = None
        self.leftHoleCenterBackCoords = None
        self.rightHoleCenterBackCoords = None
        self.name = None

    @staticmethod
    def vectorToStoppingTargetCenter(position):
        dx = CollimatorGeometry.stoppingTargetCenterPosition[0] - position[0]
        dy = CollimatorGeometry.stoppingTargetCenterPosition[1] - position[1]
        dz = CollimatorGeometry.stoppingTargetCenterPosition[2] - position[2]
        return [dx, dy, dz]

    def angleHoles(self, angleX=0, angleY=0):
        # ANGLES ARE IN RADIANS!
        # This might get a little weird for large angles, so keep them small!
        # takes whatever geometry you already have and angles the holes
        # Just the holes, not the whole block
        # Leaves the center back coordinates alone, changes the front ones
        xShift = self.collimatorThickness * math.tan(angleX)
        yShift = self.collimatorThickness * math.tan(angleY)
        self.leftHoleCenterFrontCoords = [self.leftHoleCenterFrontCoords[0]+xShift, self.leftHoleCenterFrontCoords[1]+yShift, self.leftHoleCenterFrontCoords[2]]
        self.rightHoleCenterFrontCoords = [self.rightHoleCenterFrontCoords[0]+xShift, self.rightHoleCenterFrontCoords[1]+yShift, self.rightHoleCenterFrontCoords[2]]

    def setAngleOffset(self, angleX=0, angleY=0):
        # ANGLES ARE IN RADIANS!
        # This might get a little weird for large angles, so keep them small!
        # keep back coords the same
        # assume that without this angle offset everything would be parallel
        # to the z axis
        xShift = self.collimatorThickness * math.tan(angleX)
        yShift = self.collimatorThickness * math.tan(angleY)
        self.leftHoleCenterFrontCoords = [self.leftHoleCenterBackCoords[0]+xShift, self.leftHoleCenterBackCoords[1]+yShift, self.leftHoleCenterBackCoords[2]]
        self.rightHoleCenterFrontCoords = [self.rightHoleCenterBackCoords[0]+xShift, self.rightHoleCenterBackCoords[1]+yShift, self.rightHoleCenterBackCoords[2]]



    def useStraightHoleGeometry(self):
        self.collimatorThickness = 9.144
        self.leftHoleCenterBackCoords = [-4.06, 0.0, 0.0]
        self.rightHoleCenterBackCoords = [4.06, 0.0, 0.0]

        self.leftHoleCenterFrontCoords = [-4.06, 0, -1*self.collimatorThickness]
        self.leftHoleRadius = 1.128/2.0 #0.5642
        self.rightHoleCenterFrontCoords = [4.06, 0, -1*self.collimatorThickness]
        self.rightHoleRadius = 3.568/2.0
        self.name = "StraightHole"

    def useThickerBlockGeometry(self):
        self.collimatorThickness = 6.0*self.inchToCm
        self.leftHoleCenterBackCoords = [-4.06, 0.0, 0.0]
        self.rightHoleCenterBackCoords = [4.06, 0.0, 0.0]

        self.leftHoleCenterFrontCoords = [-4.06, 0, -1*self.collimatorThickness]
        self.leftHoleRadius = 1.128/2.0 #0.5642
        self.rightHoleCenterFrontCoords = [4.06, 0, -1*self.collimatorThickness]
        self.rightHoleRadius = 3.568/2.0
        self.name = "SixInchTungsten"

    # angle in radians
    def useTippedHoleGeometry(self, angle=None):
        self.collimatorThickness = 9.144
        self.leftHoleCenterBackCoords = [-4.06, 0.0, 0.0]
        self.rightHoleCenterBackCoords = [4.06, 0.0, 0.0]
        self.leftHoleRadius = 1.128/2.0 #0.5642
        self.rightHoleRadius = 3.568/2.0

        if angle is None:
            _leftVec = self.vectorToStoppingTargetCenter(self.leftHoleCenterBackCoords)
            # collimator thickness is our z component
            _leftFrac = abs(self.collimatorThickness / _leftVec[2])
            _leftHoleCenterFrontX = self.leftHoleCenterBackCoords[0] + (_leftFrac * _leftVec[0])
            _leftHoleCenterFrontY = self.leftHoleCenterBackCoords[1] + (_leftFrac * _leftVec[1])
            _leftHoleCenterFrontZ = self.leftHoleCenterBackCoords[2] + (_leftFrac * _leftVec[2])
            self.leftHoleCenterFrontCoords = [_leftHoleCenterFrontX, _leftHoleCenterFrontY, _leftHoleCenterFrontZ]

            _rightVec = self.vectorToStoppingTargetCenter(self.rightHoleCenterBackCoords)
            # collimator thickness is our z component
            _rightFrac = abs(self.collimatorThickness / _rightVec[2])
            _rightHoleCenterFrontX = self.rightHoleCenterBackCoords[0] + (_rightFrac * _rightVec[0])
            _rightHoleCenterFrontY = self.rightHoleCenterBackCoords[1] + (_rightFrac * _rightVec[1])
            _rightHoleCenterFrontZ = self.rightHoleCenterBackCoords[2] + (_rightFrac * _rightVec[2])
            self.rightHoleCenterFrontCoords = [_rightHoleCenterFrontX, _rightHoleCenterFrontY, _rightHoleCenterFrontZ]
            self.name = "TippedHole"
        else:
            # angle is in x only
            xShift = self.collimatorThickness * math.tan(angle)
            # shift each hole closer to axis
            self.leftHoleCenterFrontCoords = [self.leftHoleCenterBackCoords[0]+xShift, self.leftHoleCenterBackCoords[1], self.leftHoleCenterBackCoords[2]]
            self.rightHoleCenterFrontCoords = [self.rightHoleCenterBackCoords[0]-xShift, self.rightHoleCenterBackCoords[1], self.rightHoleCenterBackCoords[2]]
            self.name = "Tipped{0}".format(angle)
    