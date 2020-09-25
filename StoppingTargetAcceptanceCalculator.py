#! usr/bin/env python
#import calcGeomAcceptancePointSource as ptSrc
from PointSourceAcceptanceCalculator import PointSourceAcceptanceCalculator
from CollimatorGeometry import CollimatorGeometry
import math
from scipy import integrate


class StoppingTargetAcceptanceCalculator:

    def __init__(self, collimatorGeometry=None):
        # default values from 26586 (stopping target CRR, June 2019)
        self.numFoils = 37
        self.foilSpacing = 2.22 # cm! Everything else I've been doing is in cm
        self.r_in = 2.1
        self.r_out = 7.5
        self.foilThick = 0.01
        # TODO: default position!
        if collimatorGeometry is None:
            self.geometry = CollimatorGeometry()
            self.geometry.useStraightHoleGeometry()
        else:
            self.geometry = collimatorGeometry
        self.ptSrc = PointSourceAcceptanceCalculator(self.geometry)


    #@staticmethod
    def numeratorLeft(self, z, y, x):
        return self.ptSrc.getNonNormalizedAcceptance([x, y, z], "left")
    #@staticmethod
    def numeratorRight(self, z, y, x):
        return self.ptSrc.getNonNormalizedAcceptance([x, y, z], "right")
    #@staticmethod
    def denominatorLeft(self, z, y, x):
        return self.ptSrc.getFullSphericalArea([x, y, z], self.geometry.leftHoleCenterBackCoords)
    #@staticmethod
    def denominatorRight(self, z, y, x):
        return self.ptSrc.getFullSphericalArea([x, y, z], self.geometry.rightHoleCenterBackCoords)
    @staticmethod
    def errPropDivision(num, numErr, denom, denomErr):
        resErr1 = numErr/denom
        resErr2 = -0.5*num*denomErr/ denom**2
        totalErr = math.sqrt(resErr1**2 + resErr2**2)
        return totalErr
    @staticmethod
    def errPropSum(errsList):
        sqrSum = 0.0
        for err in errsList:
            sqrSum += err**2
        return math.sqrt(sqrSum)

    def singleFoilNumerator(self, centerPos, leftOrRight):
        accept = None 
        err = None 
        x0 = centerPos[0]
        y0 = centerPos[1]
        z0 = centerPos[2]

        xLowOut = x0 - self.r_out
        xHighOut = x0 + self.r_out
        xLowIn = x0 - self.r_in
        xHighIn = x0 + self.r_in

        def yLowOut(x):
            return y0 - math.sqrt(self.r_out**2 - (x - x0)**2)
        
        def yHighOut(x):
            return y0 + math.sqrt(self.r_out**2 - (x - x0)**2)

        def yLowIn(x):
            return y0 - math.sqrt(self.r_in**2 - (x - x0)**2)
        
        def yHighIn(x):
            return y0 + math.sqrt(self.r_in**2 - (x - x0)**2)
        
        def zLow(x, y):
            return z0 - self.foilThick/2
        
        def zHigh(x, y):
            return z0 + self.foilThick/2

        # do an integral over the outer radius, subtract the interal over the
        # inner radius
        if leftOrRight is "left":
            unSubbed = integrate.tplquad(self.numeratorLeft, xLowOut, xHighOut, yLowOut, yHighOut, zLow, zHigh)
            center = integrate.tplquad(self.numeratorLeft, xLowIn, xHighIn, yLowIn, yHighIn, zLow, zHigh)
            accept = unSubbed[0] - center[0]
            err = self.errPropSum([unSubbed[1], center[1]])
        elif leftOrRight is "right":
            unSubbed = integrate.tplquad(self.numeratorRight, xLowOut, xHighOut, yLowOut, yHighOut, zLow, zHigh)
            center = integrate.tplquad(self.numeratorRight, xLowIn, xHighIn, yLowIn, yHighIn, zLow, zHigh)
            accept = unSubbed[0] - center[0]
            err = self.errPropSum([unSubbed[1], center[1]])

        return (accept, err)

    def singleFoilDenominator(self, centerPos, leftOrRight):
        accept = None 
        err = None 
        x0 = centerPos[0]
        y0 = centerPos[1]
        z0 = centerPos[2]

        xLowOut = x0 - self.r_out
        xHighOut = x0 + self.r_out
        xLowIn = x0 - self.r_in
        xHighIn = x0 + self.r_in

        def yLowOut(x):
            return y0 - math.sqrt(self.r_out**2 - (x - x0)**2)
        
        def yHighOut(x):
            return y0 + math.sqrt(self.r_out**2 - (x - x0)**2)

        def yLowIn(x):
            return y0 - math.sqrt(self.r_in**2 - (x - x0)**2)
        
        def yHighIn(x):
            return y0 + math.sqrt(self.r_in**2 - (x - x0)**2)
        
        def zLow(x, y):
            return z0 - self.foilThick/2
        
        def zHigh(x, y):
            return z0 + self.foilThick/2

        # do an integral over the outer radius, subtract the interal over the
        # inner radius
        if leftOrRight is "left":
            unSubbed = integrate.tplquad(self.denominatorLeft, xLowOut, xHighOut, yLowOut, yHighOut, zLow, zHigh)
            center = integrate.tplquad(self.denominatorLeft, xLowIn, xHighIn, yLowIn, yHighIn, zLow, zHigh)
            accept = unSubbed[0] - center[0]
            err = self.errPropSum([unSubbed[1], center[1]])
        elif leftOrRight is "right":
            unSubbed = integrate.tplquad(self.denominatorRight, xLowOut, xHighOut, yLowOut, yHighOut, zLow, zHigh)
            center = integrate.tplquad(self.denominatorRight, xLowIn, xHighIn, yLowIn, yHighIn, zLow, zHigh)
            accept = unSubbed[0] - center[0]
            err = self.errPropSum([unSubbed[1], center[1]])

        return (accept, err)


    def getGeometricAcceptance(self, sourceCenterPos, leftOrRight):
        if leftOrRight is not "left" and leftOrRight is not "right":
            raise RuntimeError("argument leftOrRight must be either \"left\" or \"right\"")
        zLow = sourceCenterPos[2] - (self.numFoils-1)*self.foilSpacing/2.0
        # currently assuming stopping target is aligned with "my" z axis
        positions = []
        for i in range(0, self.numFoils):
            z = zLow + i*self.foilSpacing
            newPosition = [sourceCenterPos[0], sourceCenterPos[1], z]
            positions.append(newPosition)

        numSum = 0.0
        numErrs = []
        denomSum = 0.0
        denomErrs = []
        for position in positions:
            accept1, err1 = self.singleFoilNumerator(position, leftOrRight)
            accept2, err2 = self.singleFoilDenominator(position, leftOrRight)
            numSum += accept1
            denomSum += accept2
            numErrs.append(err1)
            denomErrs.append(err2)
        numErr = self.errPropSum(numErrs)
        denomErr = self.errPropSum(denomErrs)

        accept = numSum / denomSum
        err = self.errPropDivision(numSum, numErr, denomSum, denomErr)
        return accept, err
    

if __name__ == "__main__":

    # give it a little test, just to be sure we don't error out anywhere
    calc = StoppingTargetAcceptanceCalculator()
    #calc.numFoils = 3
    print "num foils: ", calc.numFoils

    centerPos = [0.0, 0.0, -2400.0]
    centerPos = CollimatorGeometry.stoppingTargetCenterPosition
    print "getting left..."
    leftResult = calc.getGeometricAcceptance(centerPos, "left")
    print "left hole result: ", leftResult
    print "getting right..."
    rightResult = calc.getGeometricAcceptance(centerPos, "right")
    print "right result: ", rightResult

