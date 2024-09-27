import numpy as np
from scipy import optimize as opt

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad
import matplotlib.pyplot as plt

from utils import numerical_fixed_mesh_diff

IR = 1
OR = 2


testPath = PATHDEF.Minimal_Polynomial_Definition4(n=2, fullParamLength=4,
                                       radii = np.array([IR,(IR+OR)/2*1.15,OR]),
                                       ffradii = np.array([1.14, 2.113]),
                                       alphaAngles = np.array([45,45])*deg2rad,
                                       betaAngles = np.array([0,87.5,175])*deg2rad,
                                       XYFactors = np.array([0.5]))
testPath2 = PATHDEF.Linear_Rn_Spiral()
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath2,
                                        IcPts=np.array([0.0025,0.0025]),
                                        IcParamLens=np.array([]))
# fuck this
testPath.get_crscRef(testCrsc)
testPath2.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.Maraging300Steel, resolution=1000)

testTorque = 3000
startingIndex = 0

# res, solnSF, divergeFlag, i = testSpring.opt_deform_by_torque(testTorque, testSpring.full_ODE,SF=np.array([0,0,0]), plotBool=True)
# testSpring.detailed_deform_regression(testTorque, testSpring.full_ODE, 5, 1)

# print(i, divergeFlag, solnSF)
# print(testSpring.solnerrVec)
plt.show()
err, res = testSpring.optimization_integration(testSpring.full_ODE,
                                               np.array([30,50,testTorque]),
                                               np.array([0.05,0.1]),testTorque,
                                               startingIndex)

rn = testSpring.path.get_rn(testSpring.ximesh[startingIndex:])
la = res[3,:]
lb = res[4,:]
stressInner = testSpring.E**res[0,:]*rn*(la/(rn-la))
stressOuter = testSpring.E**res[0,:]*rn*(lb/(rn+lb))

print(err)
plt.figure("stress success")
plt.plot(testSpring.ximesh[startingIndex:],stressInner)
plt.plot(testSpring.ximesh[startingIndex:],stressOuter)
plt.figure("deformed neutral surface")
plt.plot(testSpring.path.get_neutralSurface(200)[:,0],testSpring.path.get_neutralSurface(200)[:,1],)
plt.plot(testSpring.res[1,:],res[2,:])
plt.figure("gamma vs arc length")
plt.plot(testSpring.ximesh[startingIndex:],res[0,:])
plt.figure("la, lb")
plt.plot(testSpring.ximesh[startingIndex:],-la)
plt.plot(testSpring.ximesh[startingIndex:],lb)
plt.figure("la, lb derivs")
plt.plot(testSpring.ximesh[startingIndex:],numerical_fixed_mesh_diff(res[3,:],testSpring.ximesh[startingIndex:]))
plt.plot(testSpring.ximesh[startingIndex:],numerical_fixed_mesh_diff(res[4,:],testSpring.ximesh[startingIndex:]))
plt.show()
