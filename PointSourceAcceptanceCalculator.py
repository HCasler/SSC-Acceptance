#! usr/bin/env python
from CollimatorGeometry import CollimatorGeometry
import circleOverlap
import math

class PointSourceAcceptanceCalculator:

    def __init__(self, collimatorGeometry):
        self.geometry = collimatorGeometry

    def getNonNormalizedAcceptance(self, sourceCoords, leftOrRight):
        if leftOrRight is not "left" and leftOrRight is not "right":
            raise RuntimeError("parameter leftOrRight must be either \"left\" or \"right\"")
        
        xfc = None
        yfc = None
        zfc = None
        xbc = None
        ybc = None
        zbc = None
        rc = None
        if leftOrRight == "left":
            xfc = self.geometry.leftHoleCenterFrontCoords[0]
            yfc = self.geometry.leftHoleCenterFrontCoords[1]
            zfc = self.geometry.leftHoleCenterFrontCoords[2]
            xbc = self.geometry.leftHoleCenterBackCoords[0]
            ybc = self.geometry.leftHoleCenterBackCoords[1]
            zbc = self.geometry.leftHoleCenterBackCoords[2]
            rc = self.geometry.leftHoleRadius
        elif leftOrRight == "right":
            xfc = self.geometry.rightHoleCenterFrontCoords[0]
            yfc = self.geometry.rightHoleCenterFrontCoords[1]
            zfc = self.geometry.rightHoleCenterFrontCoords[2]
            xbc = self.geometry.rightHoleCenterBackCoords[0]
            ybc = self.geometry.rightHoleCenterBackCoords[1]
            zbc = self.geometry.rightHoleCenterBackCoords[2]
            rc = self.geometry.rightHoleRadius
        # we now have the coordinates for the collimator hole ends.
        # get the center and radius of the circle of particles that pass though
        # the front of the collimator (the "projected" curcle)
        xpc = xfc + self.geometry.collimatorThickness*(xfc - sourceCoords[0])/(zfc - sourceCoords[2])
        ypc = yfc + self.geometry.collimatorThickness*(yfc - sourceCoords[1])/(zfc - sourceCoords[2])
        zpc = zbc
        rp = rc*(1 + self.geometry.collimatorThickness/(zfc - sourceCoords[2]))
        # OK. Of all the particles that pass into the front of the collimator, how
        # many make it out the back? To find that, we take the overlap of the 
        # projected circle above and the circle representing the back and of 
        # the collimator.
        centerP = [xpc, ypc, zpc]
        centerB = [xbc, ybc, zbc]

        acceptArea = circleOverlap.circleOverlap(centerP, rp, centerB, rc)
        return acceptArea

    def getFullSphericalArea(self, sourceCoords, testPosition):
        rFromSourceSqr = (testPosition[0] - sourceCoords[0])**2 + (testPosition[1] - sourceCoords[1])**2 + (testPosition[2] - sourceCoords[2])**2
        fullArea = 4 * math.pi * rFromSourceSqr
        return fullArea

    def getGeometricAcceptance(self, sourceCoords, leftOrRight):
        acceptArea = self.getNonNormalizedAcceptance(sourceCoords, leftOrRight)

        testPosition = [0,0,0]
        if leftOrRight == "left":
            testPosition = self.geometry.leftHoleCenterBackCoords
        elif leftOrRight == "right":
            testPosition = self.geometry.rightHoleCenterBackCoords

        # finally: the fraction of all particles from the source
        fullArea = self.getFullSphericalArea(sourceCoords, testPosition)

        return acceptArea/fullArea


if __name__ == "__main__":
    from matplotlib import pyplot
    geom = CollimatorGeometry()
    geom.useStraightHoleGeometry()
    calc = PointSourceAcceptanceCalculator(geom)
    print "Quick validation.."
    # a little validation
    # put the point source directly in front of one of the collimator holes
    zDist = 200.0
    sourcePos = [geom.leftHoleCenterBackCoords[0], geom.leftHoleCenterBackCoords[1], geom.leftHoleCenterBackCoords[2]-zDist]
    # we should get pi*leftHoleRadius^2 / ( 4 pi zDist^2 ) 
    #               = leftHoleRadius^2 / ( 4 zDist^2 ) as our acceptance
    expectedAnswer = geom.leftHoleRadius**2 / (4*zDist**2)
    actualAnswer = calc.getGeometricAcceptance(sourcePos, "left")
    diff = abs(expectedAnswer - actualAnswer)
    if diff > 0.0:
        print "FAIL: placing point source directly in front of left collimator"\
        " hole should have resulted in acceptance {0}. Instead got {1}".format(\
            expectedAnswer, actualAnswer)
    else:
        print "PASS: point source directly in front of left collimator hole"

    # put it in front of a collimator hole, less than  radius from the center
    xShift = geom.leftHoleRadius*0.9
    sourcePos = [geom.leftHoleCenterBackCoords[0]+xShift, geom.leftHoleCenterBackCoords[1], geom.leftHoleCenterBackCoords[2]-zDist]
    expectedAnswer = math.pi*geom.leftHoleRadius**2
    actualAnswer = calc.getNonNormalizedAcceptance(sourcePos, "left")
    diff = abs(expectedAnswer - actualAnswer)
    if diff > 0.0:
        print "FAIL: placing point source in front of left collimator within"\
        " 1 radius of center should have resulted in non-normalized "\
        "acceptance area {0}. Instead got {1}".format(expectedAnswer, actualAnswer)
    else:
        print "PASS: point source in front of left collimator hole within 1"\
        " radius of center"

    # check the x-position of the source which gets the acceptance to 0
    zZero = -200.0
    # got the following by setting xpc + rp < xbc - rb, substituting the expression
    # for the projected circle's coordinates, and solving for the source x-position
    xZero = geom.leftHoleCenterFrontCoords[0] + geom.leftHoleRadius - (geom.leftHoleCenterFrontCoords[2] - zZero)*(geom.leftHoleCenterBackCoords[0] - geom.leftHoleCenterFrontCoords[0] - 2*geom.leftHoleRadius)/geom.collimatorThickness
    yZero = 0.0
    sourcePos = [xZero, yZero, zZero]
    expectedAnswer = 0.0
    actualAnswer = calc.getGeometricAcceptance(sourcePos, "left")
    diff = abs(expectedAnswer - actualAnswer)
    if diff > 0.0:
        print "FAIL: placing point source at special position ({0}, {1}, {2})"\
        " should have resulted in acceptance {3}. Instead got {4}"\
        .format(xZero, yZero, zZero, expectedAnswer, actualAnswer)
    else:
        print "PASS: point source located for zero acceptance"

    print "Validation done"


    # distance of source from collimator
    meters = 4
    print "Source distance: ", meters, "meters"
    # start it on the left, move it to the right, plot
    
    xSize = 8.0
    sourcePos = [-1.0*xSize, 0.0, -100.0*meters]
    stepSize = 0.1
    lefts = []
    rights = []
    poses = []
    while sourcePos[0] < xSize:
        lefts.append(calc.getGeometricAcceptance(sourcePos, "left"))
        rights.append(calc.getGeometricAcceptance(sourcePos, "right"))
        poses.append(sourcePos[0])
        sourcePos[0] += stepSize

    veryMaxLeft = max(lefts)
    veryMaxRight = max(rights)
    print "\nMax value on the left: ", veryMaxLeft
    print "\nMax value on the right: ", veryMaxRight

    pyplot.plot(poses, lefts)
    pyplot.title("left collimator, point source at {0} meters fraction acceptance".format(meters))
    pyplot.xlabel("point source x from center (cm)")
    pyplot.ylabel("fraction acceptance")
    pyplot.show()
    pyplot.plot(poses, rights)
    pyplot.title("right collimator, point source and {0} meters fraction acceptance".format(meters))
    pyplot.xlabel("point source x from center (cm)")
    pyplot.ylabel("fraction acceptance")
    pyplot.show()


    # OK let's be serious now
    # start 1/8 inch to the left, move till 1/8 inch to the right
    radius = 0.125*geom.inchToCm/2.0
    sourcePos = [-1*radius, 0.0, -1*100.0*meters]
    numSteps = 10
    stepSize = 2*radius/numSteps
    poses = []
    lefts = []
    rights = []
    for i in range(0, numSteps+1):
        poses.append(sourcePos[0])
        lefts.append(calc.getGeometricAcceptance(sourcePos, "left"))
        rights.append(calc.getGeometricAcceptance(sourcePos, "right"))
        sourcePos[0] += stepSize
    pyplot.plot(poses, lefts)
    pyplot.title("left collimator, point source ({0} m) fraction acceptance".format(meters))
    pyplot.xlabel("point source x from center (cm)")
    pyplot.ylabel("fraction acceptance")
    pyplot.show()

    pyplot.plot(poses, rights)
    pyplot.title("right collimator, point source ({0} m) fraction acceptance".format(meters))
    pyplot.xlabel("point source x from center (cm)")
    pyplot.ylabel("fraction acceptance")
    pyplot.show()

    # save and show some actually-useful numbers
    centerLeft = 0.0
    centerright = 0.0
    minLeft = min(lefts)
    minRight = min(rights)
    maxLeft = max(lefts)
    maxRight = max(rights)

    for i in range(0, len(poses)):
        if poses[i] == 0.0:
            centerLeft = lefts[i]
            centerRight = rights[i]

    
    print "\n\nWhen the source is in the very center, the left collimator gets ", centerLeft # to make sure I got it
    print "Range in the left in this x-range: (", minLeft, ", ", maxLeft, ")"
    print "When the source is in the very center, the right collimator gets ", centerRight
    print "Range in the right in this x-range: (", minRight, ", ", maxRight, ")"
    print "Range of fractions in left collimator:", (maxLeft - minLeft)
    print "Range of fractions in right collimator:", (maxRight - minRight)
    leftPercent = 100.0*(maxLeft - centerLeft)/centerLeft
    rightPercent = 100.0*(maxRight - centerRight)/centerRight
    print "Moving the source across the extended surface causes the acceptance in the left  collimator to vary by +-", leftPercent, "%"
    print "Moving the source across the extended surface causes the acceptance in the right collimator to vary by +-", rightPercent, "%"

