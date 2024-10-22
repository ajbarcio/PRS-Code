import numpy as np
import numpy.linalg as lin
from scipy import optimize as opt

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad
import matplotlib.pyplot as plt

from utils import numerical_fixed_mesh_diff

# TO TEST OUR WEIGHTED ODE WE:

# Define some paths to test
IR = 1
OR = 2

# some polynomial-based spring with local straightness at the ends
testPath = PATHDEF.Minimal_Polynomial_Definition4(n=2, fullParamLength=4,
                                       radii = np.array([IR,(IR+OR)/2*1.15,OR]),
                                       ffradii = np.array([1.14, 2.113]),
                                       alphaAngles = np.array([45,45])*deg2rad,
                                       betaAngles = np.array([0,87.5,175])*deg2rad,
                                       XYFactors = np.array([0.5]))
# some spiral that linearly increases r_n
testPath2 = PATHDEF.Linear_Rn_Spiral(start=5, end=5)
# design a cross section to test: this one uses a constant I_c at all points
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                                        IcPts=np.array([0.0025,0.0025]),
                                        IcParamLens=np.array([]))

# Determine loading conditions
testTorque = 2500

testSF     = np.array([150, -50, testTorque])

testResl   = 1000

# Determine initial geometry based on torque
preDesignStress = materials.TestMaterial.yieldStress*0.9
# print(testSpring.designStress)
h0 = np.sqrt(6*testTorque/(testCrsc.t*preDesignStress))
# print(h0)
la0 = h0/2

# One more cross section for testing: this one uses whatever I_c is needed at 
# the start to acheive target stress, and uses it for the whole length of the beam
constCrsc = CRSCDEF.Constant_Ic(pathDef=testPath2,t=0.375,h0=h0)
# fuck this
testPath.get_crscRef(testCrsc)
testPath2.get_crscRef(testCrsc)

print("outside H0", h0)
print("outside design stress", preDesignStress)

# create the test spring
testSpring = Spring(constCrsc, materials.TestMaterial, torqueCapacity=testTorque, resolution=testResl)

# startingIndex = 50

err, res = testSpring.forward_integration(testSpring.weighted_ODE,
                                                testSF,
                                                testTorque)

print(lin.norm(err))

gamma = res[0,:]
x     = res[1,:]
y     = res[2,:]
la    = res[3,:]
lb    = res[4,:]
lalb  = np.vstack((np.atleast_2d(la),np.atleast_2d(lb)))

dgds = numerical_fixed_mesh_diff(gamma,testSpring.ximesh)

print(lalb)

# print(la, lb)

plt.figure("lb-la")
plt.plot(testSpring.ximesh, lb-la)

plt.figure("eqn 1 (rn) error")
rn = testSpring.path.get_rn(testSpring.ximesh)
# plt.plot(testSpring.ximesh, rn)
# plt.plot(testSpring.ximesh, (la+lb)/(np.log((rn+lb)/(rn-la))))
plt.plot(testSpring.ximesh, rn-(la+lb)/(np.log((rn+lb)/(rn-la))))

plt.figure("eqn 2 (lb^2-la^2) error")
Fx = testSF[0]
Fy = testSF[1]
sigma = preDesignStress
fsigd = np.maximum(la/(rn-la),lb/(rn+lb))
FStar = Fy*(x-testSpring.x0)-Fx*(y-testSpring.y0)
diffTarget = 2*FStar/(testSpring.t*(sigma/fsigd-testSpring.dgds0))
plt.plot(testSpring.ximesh, (lb**2-la**2)-diffTarget)

plt.figure("weights")
plt.plot(testSpring.ximesh, testSpring.lambda1(rn))
plt.plot(testSpring.ximesh, testSpring.lambda2(rn))

plt.figure("la, lb")
# plt.plot(testSpring.xiData, [x for x in testSpring.laData], label="la")
# plt.plot(testSpring.xiData, [x for x in testSpring.lbData], label="lb")
plt.plot(testSpring.ximesh, lalb.T, label="ODE SOLN")
# plt.plot(testSpring.ximesh, testSpring.crsc.get_lalb(testSpring.ximesh).T, "--", label="CONSTANT IC SOLN")

plt.figure("stress is correct")
stressFactor = np.maximum(la/(rn-la),lb/(rn+lb))
print(dgds[0], testSpring.dgds0)
stress = testSpring.E*rn*dgds*stressFactor
# plt.plot(testSpring.ximesh, stress)
plt.plot(testSpring.xiData, testSpring.stressData)

# plt.figure("beam geometry")
# arrayla = np.array(testSpring.laData)
# arraylb = np.array(testSpring.lbData)
# e = (arraylb-arrayla)/2
# plt.plot(testSpring.xiData, -arrayla-e, label="a")
# plt.plot(testSpring.xiData, arraylb-e, label="b")
# plt.plot(testSpring.xiData, -e, label="rn")

# plt.figure("thickness")
# h=arrayla+arraylb
# plt.plot(testSpring.xiData, h)


# plt.figure("stress success")
# plt.plot(testSpring.xiData, testSpring.stressData)

# plt.figure("Ic")
# plt.plot(testSpring.xiData, testSpring.IcData)
# plt.plot(testSpring.ximesh, testSpring.crsc.get_Ic(testSpring.ximesh))

# plt.figure("deformed neutral surface")
# plt.plot(testSpring.path.get_neutralSurface(200)[:,0],testSpring.path.get_neutralSurface(200)[:,1], label = "original surface")
# plt.plot(testSpring.res[1,:],res[2,:], label = "deformed surface")
# plt.legend()
# plt.figure("gamma vs arc length")
# plt.plot(testSpring.ximesh,res[0,:])

# plt.show()
