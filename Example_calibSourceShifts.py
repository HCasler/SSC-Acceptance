#! usr/bin/env python
from CollimatorGeometry import CollimatorGeometry
from ExtendedSourceAcceptanceCalculator import ExtendedSourceAcceptanceCalculator
from matplotlib import pyplot


print "Running example script..."
chosenZ = -200.0
calibSourceCenterPosition = [0.0, 0.0, chosenZ]

shiftValuesSmall = []
minShift = -0.2 #cm, so -2 mm
maxShift = 0.2
numSteps = 20
stepSize = (maxShift - minShift)/numSteps
for i in range(0, numSteps+1):
    val = minShift + (i*stepSize)
    shiftValuesSmall.append(val)

leftAcceptances = []

geom = CollimatorGeometry()
geom.useThickerBlockGeometry()

calculator = ExtendedSourceAcceptanceCalculator(geom)

count = 1
for xShift in shiftValuesSmall:
    print "calculating acceptance {0}/{1}".format(count, len(shiftValuesSmall))
    position = [calibSourceCenterPosition[0]+xShift, calibSourceCenterPosition[1], calibSourceCenterPosition[2]]
    accept = calculator.getGeometricAcceptance(position, "left")[0]
    leftAcceptances.append(accept)
    count += 1

pyplot.plot(shiftValuesSmall, leftAcceptances)
pyplot.title("Geometric acceptance in left collimator vs source X shift")
pyplot.xlabel("source shift (cm)")
pyplot.ylabel("acceptance (frac)")
pyplot.show()

