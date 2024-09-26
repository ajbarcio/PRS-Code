import numpy as np
from scipy import optimize as opt

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad
import matplotlib.pyplot as plt

IR = 1
OR = 2

testPath = PATHDEF.Minimal_Polynomial_Definition4(n=2, fullParamLength=4,
                                       radii = np.array([IR,(IR+OR)/2*1.15,OR]),
                                       ffradii = np.array([1.14, 2.113]),
                                       alphaAngles = np.array([45,45])*deg2rad,
                                       betaAngles = np.array([0,87.5,175])*deg2rad,
                                       XYFactors = np.array([0.5]))
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath, 
                                        IcPts=np.array([0.0025,0.0025]),
                                        IcParamLens=np.array([]))
# fuck this
testPath.get_crscRef(testCrsc)
testSpring = Spring(testCrsc, materials.Maraging300Steel)

testTorque = 3000

startingIndex = 1

err, res = testSpring.optimization_integration(testSpring.full_ODE_oneSided,
                                               np.array([10,10,testTorque]),
                                               np.array([0.1,0.12]),testTorque,
                                               startingIndex)
print(err)
plt.figure()
plt.plot(testSpring.res[1,:],testSpring.res[2,:])
plt.figure()
plt.plot(testSpring.ximesh[startingIndex:],testSpring.res[0,:])
plt.figure()
plt.plot(testSpring.ximesh[startingIndex:],testSpring.res[3,:])
plt.plot(testSpring.ximesh[startingIndex:],testSpring.res[4,:])
plt.show()
