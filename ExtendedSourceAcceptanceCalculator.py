#! usr/bin/env python
#import calcGeomAcceptancePointSource as ptSrc
from PointSourceAcceptanceCalculator import PointSourceAcceptanceCalculator
import math
from scipy import integrate
from CollimatorGeometry import CollimatorGeometry

# make some plots for point and extended sources

class ExtendedSourceAcceptanceCalculator:
    inchToCm = 2.54

    # using default straight-hole geometry, or you can pass one in
    def __init__(self, collimatorGeometry=None):
        self.sourceRadius = 0.5*0.125*self.inchToCm # diameter 1/8"
        self.sourceHalfDepth = 0.09*0.5*self.inchToCm
        if collimatorGeometry is None:
            self.geometry = CollimatorGeometry()
            self.geometry.useStraightHoleGeometry()
        else:
            self.geometry = collimatorGeometry
        self.ptSrc = PointSourceAcceptanceCalculator(self.geometry)


    def numeratorLeft(self, z, y, x):
        return self.ptSrc.getNonNormalizedAcceptance([x, y, z], "left")
    
    def numeratorRight(self, z, y, x):
        return self.ptSrc.getNonNormalizedAcceptance([x, y, z], "right")
    
    def denominatorLeft(self, z, y, x):
        return self.ptSrc.getFullSphericalArea([x, y, z], self.geometry.leftHoleCenterBackCoords)
    
    def denominatorRight(self, z, y, x):
        return self.ptSrc.getFullSphericalArea([x, y, z], self.geometry.rightHoleCenterBackCoords)
    @staticmethod
    def errPropDivision(num, numErr, denom, denomErr):
        resErr1 = numErr/denom
        resErr2 = -0.5*num*denomErr/ denom**2
        totalErr = math.sqrt(resErr1**2 + resErr2**2)
        return totalErr


    def _geomAccept_pt(self, sourceCenterPos, leftOrRight):
        #print "doing point source"
        accept = self.ptSrc.getGeometricAcceptance(sourceCenterPos, leftOrRight)
        err = 0
        return (accept, err)

    def _geomAccept_zLine(self, sourceCenterPos, leftOrRight):
        #print "doing line source"
        x0 = sourceCenterPos[0]
        y0 = sourceCenterPos[1]
        z0 =sourceCenterPos[2]
        zLowVal = z0 - self.sourceHalfDepth
        zHighVal = z0 + self.sourceHalfDepth
        num = None
        numErr = None
        denom = None
        denomErr = None
        if leftOrRight == "left":
            num, numErr = integrate.quad(self.numeratorLeft, zLowVal, zHighVal, args=(y0, x0))
            denom, denomErr = integrate.quad(self.denominatorLeft, zLowVal, zHighVal, args=(y0, x0))
        if leftOrRight == "right":
            num, numErr = integrate.quad(self.numeratorRight, zLowVal, zHighVal, args=(y0, x0))
            denom, denomErr = integrate.quad(self.denominatorRight, zLowVal, zHighVal, args=(y0, x0))
        accept = num/denom
        err = self.errPropDivision(num, numErr, denom, denomErr)
        return (accept, err)

    def _geomAccept_xyDisc(self, sourceCenterPos, leftOrRight):
        #print "doing disc source"
        accept = None 
        err = None
        x0 = sourceCenterPos[0]
        y0 = sourceCenterPos[1]
        z0 =sourceCenterPos[2]
        xLow = x0 - self.sourceRadius
        xHigh = x0 + self.sourceRadius
        def yLow(x):
            return y0 - math.sqrt(self.sourceRadius**2 - (x - x0)**2)
        def yHigh(x):
            return y0 + math.sqrt(self.sourceRadius**2 - (x - x0)**2)

        def numLeft2D(y, x):
            return self.ptSrc.getNonNormalizedAcceptance([x, y, z0], "left")
        def numRight2D(y, x):
            return self.ptSrc.getNonNormalizedAcceptance([x, y, z0], "right")
        def denomLeft2D(y, x):
            return self.ptSrc.getFullSphericalArea([x, y, z0], self.geometry.leftHoleCenterBackCoords)
        def denomRight2D(y, x):
            return self.ptSrc.getFullSphericalArea([x, y, z0], self.geometry.rightHoleCenterBackCoords)

        if leftOrRight is "left":
            leftTop = integrate.dblquad(numLeft2D, xLow, xHigh, yLow, yHigh)
            leftBottom = integrate.dblquad(denomLeft2D, xLow, xHigh, yLow, yHigh)
            accept = leftTop[0] / leftBottom[0]
            err = self.errPropDivision(leftTop[0], leftTop[1], leftBottom[0], leftBottom[1])
        elif leftOrRight is "right":
            rightTop = integrate.dblquad(numRight2D, xLow, xHigh, yLow, yHigh)
            rightBottom = integrate.dblquad(denomRight2D, xLow, xHigh, yLow, yHigh)
            accept = rightTop[0] / rightBottom[0]
            err = self.errPropDivision(rightTop[0], rightTop[1], rightBottom[0], rightBottom[1])
        else:
            raise RuntimeError("second argument to getGeometricAcceptance() must be either \"left\" or \"right\"")
        return (accept, err)


    def _geomAccept_xyzCyl(self, sourceCenterPos, leftOrRight):
        #print "doing cylinder source"
        accept = None 
        err = None 
        x0 = sourceCenterPos[0]
        y0 = sourceCenterPos[1]
        z0 =sourceCenterPos[2]

        xLow = x0 - self.sourceRadius
        xHigh = x0 + self.sourceRadius

        def yLow(x):
            return y0 - math.sqrt(self.sourceRadius**2 - (x - x0)**2)
        
        def yHigh(x):
            return y0 + math.sqrt(self.sourceRadius**2 - (x - x0)**2)
        
        def zLow(x, y):
            return z0 - self.sourceHalfDepth
        
        def zHigh(x, y):
            return z0 + self.sourceHalfDepth

        if leftOrRight is "left":
            leftTop = integrate.tplquad(self.numeratorLeft, xLow, xHigh, yLow, yHigh, zLow, zHigh)
            leftBottom = integrate.tplquad(self.denominatorLeft, xLow, xHigh, yLow, yHigh, zLow, zHigh)
            accept = leftTop[0] / leftBottom[0]
            err = self.errPropDivision(leftTop[0], leftTop[1], leftBottom[0], leftBottom[1])
        elif leftOrRight is "right":
            rightTop = integrate.tplquad(self.numeratorRight, xLow, xHigh, yLow, yHigh, zLow, zHigh)
            rightBottom = integrate.tplquad(self.denominatorRight, xLow, xHigh, yLow, yHigh, zLow, zHigh)
            accept = rightTop[0] / rightBottom[0]
            err = self.errPropDivision(rightTop[0], rightTop[1], rightBottom[0], rightBottom[1])
        else:
            raise RuntimeError("second argument to getGeometricAcceptance() must be either \"left\" or \"right\"")
        return (accept, err)


    def getGeometricAcceptance(self, sourceCenterPos, leftOrRight):
        if self.sourceRadius == 0 and self.sourceHalfDepth == 0:
            # this is just a point source then
            return self._geomAccept_pt(sourceCenterPos, leftOrRight)
        elif self.sourceRadius == 0 and self.sourceHalfDepth != 0:
            # this is a line in the z direction
            return self._geomAccept_zLine(sourceCenterPos, leftOrRight)
        elif self.sourceHalfDepth == 0 and self.sourceRadius != 0:
            # this is a flat disc in the x-y plane
            return self._geomAccept_xyDisc(sourceCenterPos, leftOrRight)
        else:
            # this is a 3d cylinder
            return self._geomAccept_xyzCyl(sourceCenterPos, leftOrRight)

        

if __name__ == "__main__":
    from matplotlib import pyplot
    from scipy import special
    import numpy as np
    import sys

    print "Validation checks for ExtendedSourceAcceptanceCalculator"
    calc = ExtendedSourceAcceptanceCalculator()
    geom = calc.geometry

    print "Validation for point sources"

    # OK, step 1: put the point source directly in front of a collimator hole
    zDists = [9.144+1.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0]

    calc.sourceRadius = 0
    calc.sourceHalfDepth = 0 # point source
    collRadius = geom.leftHoleRadius
    collPosBack = geom.leftHoleCenterBackCoords
    gottenAnswers = []
    gottenErrs = []
    expectedAnswers = []
    diffsPercent = []
    for zDist in zDists:
        expectedAnswer = geom.leftHoleRadius**2 / (4*zDist**2)
        calcSourcePos = [collPosBack[0], collPosBack[1], collPosBack[2]-zDist]
        calced, err = calc.getGeometricAcceptance(calcSourcePos, "left")
        expectedAnswers.append(expectedAnswer)
        gottenAnswers.append(calced)
        gottenErrs.append(err)
        diff = 100.0*(calced - expectedAnswer)/expectedAnswer
        diffsPercent.append(diff)

    pyplot.plot(zDists, diffsPercent)
    pyplot.suptitle("Varying source distance, point source centered on small collimator")
    pyplot.title("Diff between analytic calculation and this integrator")
    pyplot.xlabel("point source distance (cm)")
    pyplot.ylabel("difference in acceptance (%)")
    pyplot.show()


    # Step 2: place the point source somewhere that it *should* have an
    # acceptance of 0. See that it does.
    calc.sourceRadius = 0
    calc.sourceHalfDepth = 0 # point source
    collRadius = geom.leftHoleRadius
    collPosBack = geom.leftHoleCenterBackCoords
    gottenAnswers = []
    gottenErrs = []
    expectedAnswers = []
    diffs = []
    for zDist in zDists:
        zZero = collPosBack[2]-zDist
        # got the following by setting xpc + rp < xbc - rb, substituting the expression
        # for the projected circle's coordinates, and solving for the source x-position
        xZero = geom.leftHoleCenterFrontCoords[0] + collRadius - (geom.leftHoleCenterFrontCoords[2] - zZero)*(geom.leftHoleCenterBackCoords[0] - geom.leftHoleCenterFrontCoords[0] - 2*collRadius)/geom.collimatorThickness
        yZero = 0.0
        sourcePos = [xZero, yZero, zZero]
        calced, err = calc.getGeometricAcceptance(sourcePos, "left")
        expected = 0
        diff = calced - expected
        gottenAnswers.append(calced)
        gottenErrs.append(err)
        expectedAnswers.append(expected)
        diffs.append(diff)

    pyplot.plot(zDists, diffs)
    pyplot.suptitle("Varying source distance, point source, located at x position where acceptance should be zero")
    pyplot.title("Diff between analytic calculation and this integrator")
    pyplot.xlabel("point source distance (cm)")
    pyplot.ylabel("difference in acceptance (raw)")
    pyplot.show()


    print "Validation for 2D disc sources"

    def validationFuncCoaxialDisc(r_source, r_collimator, z_dist):
        # from L. Ruby, J.B. Rechen,
        # A simpler approach to the geometrical efficiency of a parallel-disk source and detector system,
        # Nuclear Instruments and Methods,
        # Volume 58, Issue 2,
        # 1968,
        # Pages 345-346,
        # ISSN 0029-554X,
        # https://doi.org/10.1016/0029-554X(68)90491-6.
        prefactor = r_collimator / r_source
        def integrand(k):
            fac1 = math.exp(-1*k*z_dist) / k
            fac2 = special.jn(1, k*r_source)
            fac3 = special.jn(1, k*r_collimator)
            return fac1*fac2*fac3
        integral, err = integrate.quad(integrand, 0, np.inf)
        return prefactor * integral

    

    discRadii_cm = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0]
    discRadiiZoomed_cm = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.55625, 0.5625, 0.56875, 0.575, 0.58125, 0.5875, 0.59375, 0.6]
    discDist_cm = [9.144+1.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0]
    discDistZoomed_cm = [100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0]
    discThickness_cm = [0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5]

    
    # OK, step 1: at a fixed distance (1 m), get the acceptance for different radii
    collRadius = geom.leftHoleRadius
    collPosBack = geom.leftHoleCenterBackCoords
    overallZ = 100.0
    calcSourcePos = [collPosBack[0], collPosBack[1], collPosBack[2]-overallZ]
    radiusAcceptRubyRechen = []
    radiusAcceptCalc = []
    diffsPercent = []
    for radius in discRadii_cm:
        analytic = validationFuncCoaxialDisc(radius, collRadius, overallZ)
        calc.sourceRadius = radius
        calc.sourceHalfDepth = 0
        calced, err = calc.getGeometricAcceptance(calcSourcePos, "left")
        radiusAcceptRubyRechen.append(analytic)
        radiusAcceptCalc.append(calced)
        diff = 100.0*(calced - analytic)/analytic
        diffsPercent.append(diff)

    verticalLineInchX = (1.0/16.0)*calc.inchToCm
    yRange = max(diffsPercent) - min(diffsPercent) 
    verticalLineYMax = max(diffsPercent) + 0.1*yRange
    verticalLineYMin = min(diffsPercent) - 0.1*yRange

    verticalLineInchX2 = collRadius

    pyplot.plot(discRadii_cm, diffsPercent)
    pyplot.plot([verticalLineInchX, verticalLineInchX], [verticalLineYMin, verticalLineYMax])
    pyplot.plot([verticalLineInchX2, verticalLineInchX2], [verticalLineYMin, verticalLineYMax])
    pyplot.suptitle("Varying source radius, disc source, 1 m distance")
    pyplot.title("Diff between Ruby+Rechen and this integrator")
    pyplot.xlabel("disc source radius (cm)")
    pyplot.ylabel("difference in acceptance (%)")
    pyplot.show()

    # OK, step 2: at a fixed distance (1 m), get the acceptance for different radii
    # BUT ZOOM IN
    collRadius = geom.leftHoleRadius
    collPosBack = geom.leftHoleCenterBackCoords
    overallZ = 100.0
    calcSourcePos = [collPosBack[0], collPosBack[1], collPosBack[2]-overallZ]
    radiusAcceptRubyRechen = []
    radiusAcceptCalc = []
    diffsPercent = []
    for radius in discRadiiZoomed_cm:
        analytic = validationFuncCoaxialDisc(radius, collRadius, overallZ)
        calc.sourceRadius = radius
        calc.sourceHalfDepth = 0
        calced, err = calc.getGeometricAcceptance(calcSourcePos, "left")
        radiusAcceptRubyRechen.append(analytic)
        radiusAcceptCalc.append(calced)
        diff = 100.0*(calced - analytic)/analytic
        diffsPercent.append(diff)

    verticalLineInchX = (1.0/16.0)*geom.inchToCm
    yRange = max(diffsPercent) - min(diffsPercent) 
    verticalLineYMax = max(diffsPercent) + 0.1*yRange
    verticalLineYMin = min(diffsPercent) - 0.1*yRange
    verticalLineInchX2 = collRadius

    pyplot.plot(discRadiiZoomed_cm, diffsPercent)
    pyplot.plot([verticalLineInchX, verticalLineInchX], [verticalLineYMin, verticalLineYMax])
    pyplot.plot([verticalLineInchX2, verticalLineInchX2], [verticalLineYMin, verticalLineYMax])
    pyplot.suptitle("Varying source radius, disc source, 1 m distance -- zoomed in")
    pyplot.title("Diff between Ruby+Rechen and this integrator")
    pyplot.xlabel("disc source radius (cm)")
    pyplot.ylabel("difference in acceptance (%)")
    pyplot.show()

    # OK, step 3: using our real source radius, but a flat disc: see how we 
    # compare to Ruby-Rechen as a function of distance from the collimator
    collRadius = geom.leftHoleRadius
    collPosBack = geom.leftHoleCenterBackCoords
    distAcceptRubyRechen = []
    distAcceptCalc = []
    diffsPercent = []
    calc.sourceRadius = 0.5*0.125*geom.inchToCm # diameter 1/8"
    for zDist in discDist_cm:
        analytic = validationFuncCoaxialDisc(calc.sourceRadius, collRadius, zDist)
        sourcePos = [collPosBack[0], collPosBack[1], collPosBack[2]-zDist]
        calced, err = calc.getGeometricAcceptance(sourcePos, "left")
        distAcceptRubyRechen.append(analytic)
        distAcceptCalc.append(calced)
        diff = 100.0*(calced - analytic)/analytic
        diffsPercent.append(diff)

    pyplot.plot(discDist_cm, diffsPercent)
    pyplot.suptitle("Varying source distance, disc source, 1/16 inch radius")
    pyplot.title("Diff between Ruby+Rechen and this integrator")
    pyplot.xlabel("disc source distance (cm)")
    pyplot.ylabel("difference in acceptance (%)")
    pyplot.show()

    # pyplot.plot(discDist_cm, distAcceptCalc)
    # pyplot.suptitle("Varying source distance, disc source, 1/16 inch radius")
    # pyplot.title("Acceptance using this integrator")
    # pyplot.xlabel("disc source distance (cm)")
    # pyplot.ylabel("acceptance")
    # pyplot.show()

    # pyplot.plot(discDist_cm, distAcceptRubyRechen)
    # pyplot.suptitle("Varying source distance, disc source, 1/16 inch radius")
    # pyplot.title("Acceptance using this Ruby-Rechen")
    # pyplot.xlabel("disc source distance (cm)")
    # pyplot.ylabel("acceptance")
    # pyplot.show()


    # OK, step 4: same as step 3, but focus on distances >= 1 m
    collRadius = geom.leftHoleRadius
    collPosBack = geom.leftHoleCenterBackCoords
    distAcceptRubyRechen = []
    distAcceptCalc = []
    diffsPercent = []
    calc.sourceRadius = 0.5*0.125*geom.inchToCm # diameter 1/8"
    for zDist in discDistZoomed_cm:
        analytic = validationFuncCoaxialDisc(calc.sourceRadius, collRadius, zDist)
        sourcePos = [collPosBack[0], collPosBack[1], collPosBack[2]-zDist]
        calced, err = calc.getGeometricAcceptance(sourcePos, "left")
        distAcceptRubyRechen.append(analytic)
        distAcceptCalc.append(calced)
        diff = 100.0*(calced - analytic)/analytic
        diffsPercent.append(diff)

    pyplot.plot(discDistZoomed_cm, diffsPercent)
    pyplot.suptitle("Varying source distance, disc source, 1/16 inch radius -- zoomed in")
    pyplot.title("Diff between Ruby+Rechen and this integrator")
    pyplot.xlabel("disc source distance (cm)")
    pyplot.ylabel("difference in acceptance (%)")
    pyplot.show()

    # OK step 5: does the result for a cylinder approach the result for a disc
    # as the cylinder gets thinner and thinner?
    # Pick a single distance -- 1 m
    # Pick a single radius -- our real source radius (1/16 inch)
    collRadius = geom.leftHoleRadius
    collPosBack =geom.leftHoleCenterBackCoords
    overallZ = 100.0
    calcSourcePos = [collPosBack[0], collPosBack[1], collPosBack[2]-overallZ]
    calc.sourceRadius = 0.5*0.125*geom.inchToCm
    radiusAcceptRubyRechen = []
    radiusAcceptCalc = []
    diffsPercent = []
    for thickness in discThickness_cm:
        analytic = validationFuncCoaxialDisc(radius, collRadius, overallZ)
        calc.sourceHalfDepth = thickness/2
        calced, err = calc.getGeometricAcceptance(calcSourcePos, "left")
        radiusAcceptRubyRechen.append(analytic)
        radiusAcceptCalc.append(calced)
        diff = 100.0*(calced - analytic)/analytic
        diffsPercent.append(diff)

    verticalLineInchX = (0.09)*geom.inchToCm
    yRange = max(diffsPercent) - min(diffsPercent) 
    verticalLineYMax = max(diffsPercent) + 0.1*yRange
    verticalLineYMin = min(diffsPercent) - 0.1*yRange

    pyplot.plot(discThickness_cm, diffsPercent)
    pyplot.plot([verticalLineInchX, verticalLineInchX], [verticalLineYMin, verticalLineYMax])
    pyplot.suptitle("Varying source thickness, cylinder source, 1 m distance")
    pyplot.title("Diff between Ruby+Rechen and this integrator")
    pyplot.xlabel("disc source thickness (cm)")
    pyplot.ylabel("difference in acceptance (%)")
    pyplot.show()
