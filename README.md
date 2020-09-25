# Spot-size collimator geometric acceptance

This code calculates the fraction of particles emitted by a source that make it through the STM spot-size collimator. It assumes no other hindrances exist between the source and the collimator, or at the very least, that these hindrances do not block any particles that would have made it to the collimator. So, if using the stopping target as a source, the assumption is that the field-of-view collimator is wide enough, and well-aligned enough, to not change the acceptance in the spot-size collimator.

This is the code used to calculate the results in Mu2e-doc-35098. That document also gives a basic overall desription of how the calculation works. In short, it calculates the geometric acceptance for a point source, and then integrates that over the volume of the full radiation source.

## Note on the geometry
The coordinate system used in these scripts palces the origin at the *downstream end* of the spot-size collimator. The *z* axis is the same as the *z* axis of the stopping target in the usual Mu2e geometry. The positive *z* direction is from the stopping target, toward the STM. This is why the *z* positions of the calibration source, FOV collimator, etc are all negative.

Length units are centimeters.

The radiation source is generally assumed to be a cylinder, whose central axis is parallel to the *z* axis of the coordinate system. The source is also assumed to emit uniformly and isotropically.

## Installation (such as it is)
Just copy the scripts to the directory where you will use them. 

## Dependencies
The code relies on scipy.integrate. The example script also uses matplotlib for plotting.

## Usage

An example script is given in **Example_calibSourceShifts.py**. The script works as follows:

    from CollimatorGeometry import CollimatorGeometry
    from ExtendedSourceAcceptanceCalculator import ExtendedSourceAcceptanceCalculator

CollimatorGeometry is a class that lets you describe the positions and sizes of the collimator holes, as well as the thickness of the collimator itself. ExtendedSourceAcceptanceCalculator does the actual acceptance calculation.

    geom = CollimatorGeometry()
    geom.useThickerBlockGeometry()

Create a CollimatorGeometry instance, and use the method useThickerBlockGeometry() to set the required internal values. If you open up **CollimatorGeometry.py**, you can see which values are set, and which geometries can be created "automatically" via methods like useThickerBlockGeometry(). You can also set these value individually by hand, to make a custom geometry.

    calculator = ExtendedSourceAcceptanceCalculator(geom)

Create an ExtendedSourceAcceptanceCalculator instance, passing the collimator geometry as an argument to the constructor. If you do not include a constructor argument, it defaults to the old, 9-cm-thick geometry, where one collimator hole is 1 cm^2 and the other if 10 cm^2. 

You can optionally set the source geometry as well. By default, ExtendedSourceAcceptanceCalculator assumes you are working with the calibration source, which is a cylinder 1/8" in diameter and 0.09" in length. To change this, use the following:

    calculator = ExtendedSourceAcceptanceCalculator(geom)
    calculator.sourceRadius = 0.5*someDiameterInCm
    calculator.sourceHalfDepth = 0.5*someLengthInCm

The source is assumed to be a cylinder, whose central axis is parallel to the *z* axis of the coordinate system. You can get a point source by setting the radius and half-depth to both be zero. Set just the half-depth to zero for a flat disc, or the radius to zero for a line source.

    accept = calculator.getGeometricAcceptance(position, "left")[0]

The method getGeometricAcceptance() takes two arguments. The first is the position of the radiation source, in the form of a three-element list [*x*, *y*, *z*]. The second is either the string "left" or "right". These strings refer to which collimator aperture you want. The "left" aperture means the one on the negative-*x* side of the origin -- that is "left" = "HPGe". The "right" one refers to the positive-*x* side, or the LaBr. The method returns a pair of values:

    (acceptance, error)

This "error" is the error reported by the integrator. It's strictly about the computation of the integral, and **should *not* be taken to refer to physics-related errors coming from count rate, photopeak acceptance, etc.**

## Stopping target as a source
The class StoppingTargetAcceptanceCalculator lets you calculate the geometric acceptance using the stopping target as a source. It operates the same way as ExtendedSourceAcceptanceCalculator, but instead of treating the source as a single cylinder, it treats it as a number of annular foils. The default source geometry is the one given for the stopping target in Mu2e-doc-26586. The method getGeometricAcceptance() works the same as well -- you pass in a position and "left" or "right", and get back an acceptance and an error.

## Other stuff
CollimatorGeometry contains a static method called vectorToStoppingTargetCenter(). Pass it a position and it gives you a vector from that position to the center of the stopping target.