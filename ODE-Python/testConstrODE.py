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
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                                        IcPts=np.array([0.0025,0.0025]),
                                        IcParamLens=np.array([]))
# fuck this
testPath.get_crscRef(testCrsc)
testPath2.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.TestMaterial, resolution=100)

testTorque = 3000
testSF     = np.array([50,-150,testTorque])
print(testSpring.designStress)
h0 = np.sqrt(6*testTorque/(testSpring.t*testSpring.designStress))
print(h0)
la0 = h0/2
print(la0)
confirmStress = testTorque*la0/(testSpring.t*h0**3/12)
print(confirmStress)


err, res = testSpring.forward_integration_constr(testSpring.constr_deform_ODE,
                                                testSF,
                                                testTorque,
                                                la0,la0)

print(testSpring.designStress)

plt.figure()
plt.plot(testSpring.xiData, [x * -1 for x in testSpring.laData])
plt.plot(testSpring.xiData, testSpring.lbData)

plt.figure()
plt.plot(testSpring.xiData, testSpring.stressData)

# plt.show()
#res: gamma, x, y, la, lb
# rn = testSpring.path.get_rn(testSpring.ximesh)
# la = res[3,:]
# lb = res[4,:]
# dgdxi = numerical_fixed_mesh_diff(res[0,:], testSpring.ximesh)
# dxids = testSpring.path.get_dxi_n(testSpring.ximesh)
# dgds = dgdxi*dxids
# stressInner = abs(testSpring.E*dgds*rn*(la/(rn-la)))
# stressOuter = abs(testSpring.E*dgds*rn*(lb/(rn+lb)))

# cmap = plt.get_cmap('rainbow')

# print(err)
# plt.figure("stress success")
# plt.plot(testSpring.ximesh[startingIndex:],stressInner, label = "stress la")
# plt.plot(testSpring.ximesh[startingIndex:],stressOuter, label = "stress lb")
# plt.legend()
plt.figure("deformed neutral surface")
plt.plot(testSpring.path.get_neutralSurface(200)[:,0],testSpring.path.get_neutralSurface(200)[:,1], label = "original surface")
plt.plot(testSpring.res[1,:],res[2,:], label = "deformed surface")
plt.legend()
# plt.figure("gamma vs arc length")
# plt.plot(testSpring.ximesh[startingIndex:],res[0,:])
# plt.figure("la, lb")
# plt.plot(testSpring.ximesh[startingIndex:],-la, label = "la")
# plt.plot(testSpring.ximesh[startingIndex:],lb, label = "lb")
# plt.legend()
# plt.figure("la, lb derivs")
# plt.plot(testSpring.ximesh[startingIndex:],numerical_fixed_mesh_diff(res[3,:],testSpring.ximesh[startingIndex:])*dxids, label = "dlads")
# plt.plot(testSpring.ximesh[startingIndex:],numerical_fixed_mesh_diff(res[4,:],testSpring.ximesh[startingIndex:])*dxids, label = "dlbds")
# plt.legend()
plt.show()
