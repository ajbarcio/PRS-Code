import numpy as np
import numpy.linalg as lin
from scipy import optimize as opt

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad
import matplotlib.pyplot as plt

from utils import numerical_fixed_mesh_diff

IR = 1
OR = 2

testTorque = 2500

testSF     = np.array([0, 0,testTorque])

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


preDesignStress = materials.TestMaterial.yieldStress*0.9
# print(testSpring.designStress)
h0 = np.sqrt(6*testTorque/(testCrsc.t*preDesignStress))
# print(h0)
la0 = h0/2
# print(la0)
confirmStress = testTorque*la0/(testCrsc.t*h0**3/12)
# print(confirmStress)

print("outside H0", h0)
print("outside design stress", preDesignStress)

constCrsc = CRSCDEF.Constant_Ic(pathDef=testPath,t=0.375,h0=h0)
# fuck this
testPath.get_crscRef(testCrsc)
testPath2.get_crscRef(testCrsc)

testSpring = Spring(constCrsc, materials.TestMaterial, torqueCapacity=testTorque, resolution=5000)

res, solnSFstraight, divergeFlag, i = testSpring.deform_by_torque(testTorque, testSpring.deform_ODE, SF=testSF)

err, res = testSpring.forward_integration(testSpring.constr_deform_ODE,
                                                solnSFstraight,
                                                testTorque)

print(lin.norm(err))

# res, solnSF, divergeFlag, i = testSpring.deform_by_torque(testTorque, testSpring.constr_deform_ODE, SF=solnSFstraight)

# print(testSpring.solnerr)
print(solnSFstraight)


# fxmesh = np.linspace(-3000,3000,100)
# fymesh = np.linspace(-3000,3000,100)
# errs = np.empty([len(fxmesh),len(fymesh)])
# errs2 = np.empty([len(fxmesh),len(fymesh)])
# for indexa, i in enumerate(fxmesh):
#     for indexb, j in enumerate(fymesh):
#         err, res = testSpring.forward_integration(testSpring.deform_ODE,
#                                                   np.array([i,j,testTorque]),testTorque)
#         err2, res2 = testSpring.forward_integration(testSpring.constr_deform_ODE,
#                                                   np.array([i,j,testTorque]),testTorque)
#         errs[indexa,indexb] = lin.norm(err)
#         errs2[indexa,indexb] = lin.norm(err2)
#         print(indexa, indexb)
# I, J = np.meshgrid(fxmesh, fymesh)
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# ax.plot_surface(I,J,errs)
# ax.plot_surface(I,J,errs2)
# ax.plot_surface(I,J,np.zeros_like(J))

# print(np.min(errs))
# print(np.min(errs2))

# plt.show()

# print(testSpring.designStress)

# res, solnSF, divergeFlag, i = testSpring.deform_by_torque(testTorque, testSpring.constr_deform_ODE, testSF)

plt.figure()
plt.plot(testSpring.xiData, [x for x in testSpring.laData], label="la")
plt.plot(testSpring.xiData, [x for x in testSpring.lbData], label="lb")
plt.plot(testSpring.ximesh, testSpring.crsc.get_lalb(testSpring.ximesh).T, "--")

plt.figure()
plt.plot(testSpring.xiData, testSpring.stressData)

plt.figure()
plt.plot(testSpring.xiData, testSpring.IcData)
plt.plot(testSpring.ximesh, testSpring.crsc.get_Ic(testSpring.ximesh))

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
plt.figure("gamma vs arc length")
plt.plot(testSpring.ximesh,res[0,:])
# plt.figure("la, lb")
# plt.plot(testSpring.ximesh[startingIndex:],-la, label = "la")
# plt.plot(testSpring.ximesh[startingIndex:],lb, label = "lb")
# plt.legend()
# plt.figure("la, lb derivs")
# plt.plot(testSpring.ximesh[startingIndex:],numerical_fixed_mesh_diff(res[3,:],testSpring.ximesh[startingIndex:])*dxids, label = "dlads")
# plt.plot(testSpring.ximesh[startingIndex:],numerical_fixed_mesh_diff(res[4,:],testSpring.ximesh[startingIndex:])*dxids, label = "dlbds")
# plt.legend()
plt.show()
