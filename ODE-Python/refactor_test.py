import numpy as np


import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad

# Standard procedure: define a path, then define a thickness profile using the
#                     length of that path, then pass both of those objects into
#                     the overall spring object that contains deformation
#                     methods

testPath = PATHDEF.Minimal_Polynomial_Definition(n=1, fullParamLength=4,
                                       radii = np.array([1,(1+2.4)/2,2.4]),
                                       betaAngles = np.array([0,45,90])*deg2rad,
                                       XYFactors = np.array([0.5]))
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                                        IcPts = np.array([.004/2, .0004, .004/2]),
                                        IcParamLens = np.array([0.5]))
# fuck this
testPath.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.Maraging300Steel)

#just try to print it out
# testSpring.plot_spring(showBoolFalse)

res, SF, divergeFlag, i = testSpring.deform_by_torque(4550, testSpring.deform_ODE)
print(testSpring.solnerr)
print(SF)
testSpring.plot_deform()