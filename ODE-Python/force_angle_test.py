import numpy as np
from scipy.integrate import solve_ivp as variable_mesh_solve
import numpy.linalg as lin
import matplotlib.pyplot as plt
import time
from stiffness_library import *
from StatProfiler import SSProfile

then = time.monotonic()
# Give the spring some parameters that hopefully do not affect anything about the math

## this is a spring we know doesn't work

deg2rad = np.pi/180
n = 2
fullArcLength = 5

E = 27.5*10**6
self.t = .375

globalInnerRadiusLimit = 0.75
globalOuterRadiusLimit = 6/2

R0 = 2.2/2
R3 = 5.9/2

R1 = (R0+R3)/2+.26
R2 = (R0+R3)/2-.25

funnyAngle = 150

betaB = funnyAngle/3*deg2rad*.5
betaC = 2*funnyAngle/3*deg2rad
betaD = funnyAngle*deg2rad
beta0 = betaD

x0 = R0
y0 = 0

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
cIs  = np.array([.008, .001, .008])

ctrlcIs      = np.array([0,fullArcLength*.35,fullArcLength])
ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

# maxTorque = 4554.5938
maxTorque = 0
maxDBeta  = 0.087266463

dragVector0 = [R0, R1, R2, R3, \
               betaB, betaC, beta0, \
               cIs[0], cIs[1], cIs[2], \
               ctrlcIs[1], \
               ctrlLengths[1], ctrlLengths[2],
               fullArcLength]

geometryDef, smesh = drag_vector_spring(dragVector0)

funnyResults = np.zeros([4,10])

# for i in range(10):
#     SF = np.array([0,0,i+1])
#     res,SF,divergeFlag = spring_deform_eval(SF, geometryDef, fullArcLength, deform_ODE)
#     dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(coord(fullArcLength, geometryDef[1]),coord(fullArcLength, geometryDef[0])))
#     torque = SF[2]
#     forceAngle = np.arctan2(SF[1],SF[0])/deg2rad
#     Fx=SF[0]
#     Fy =SF[1]
#     plt.figure("Force Vectors")
#     plt.arrow(i,0,SF[0]/lin.norm([SF[0],SF[1]]),SF[1]/lin.norm([SF[0],SF[1]]))
#     print(i,divergeFlag)
#     funnyResults[:,i]=[torque,forceAngle,Fx,Fy]
# plt.axis('equal')
# print(funnyResults)
# print(np.transpose(funnyResults))
# plt.figure("results")
# plt.plot(funnyResults[0,:],np.transpose(funnyResults[1:4,:]))
# plt.legend(("angle", "Fx", "Fy"))
# print("beta0:", np.arctan2(coord(fullArcLength, geometryDef[1]),coord(fullArcLength, geometryDef[0])))

res, SF, divergeFlag = deform_spring_by_torque2(4554.5938,geometryDef,fullArcLength,deform_ODE)

xorg = coord(smesh, geometryDef[0])
yorg = coord(smesh, geometryDef[1])

# plt.show()

# this is a spring that should work

deg2rad = np.pi/180
n = 2
fullArcLength = 5

E = 27.5*10**6
self.t = .375

globalInnerRadiusLimit = 0.75
globalOuterRadiusLimit = 6/2

R0 = 2.2/2
R3 = 5.9/2

R1 = (R0+R3)/2+.26
R2 = (R0+R3)/2-.25

funnyAngle = 150

betaB = funnyAngle/3*deg2rad*.5
betaC = 2*funnyAngle/3*deg2rad
betaD = funnyAngle*deg2rad
beta0 = betaD

x0 = R0
y0 = 0

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
cIs  = np.array([.005, .0001, .002])*.4

ctrlcIs      = np.array([0,fullArcLength*.6,fullArcLength])
ctrlLengths = np.array([0,fullArcLength*0.4,fullArcLength*0.667,fullArcLength])

# maxTorque = 4554.5938
maxTorque = 0
maxDBeta  = 0.087266463

dragVector0 = [R0, R1, R2, R3, \
               betaB, betaC, beta0, \
               cIs[0], cIs[1], cIs[2], \
               ctrlcIs[1], \
               ctrlLengths[1], ctrlLengths[2],
               fullArcLength]

geometryDef, smesh = drag_vector_spring(dragVector0)

res, SF, divergeFlag = deform_spring_by_torque2(4554.5938,geometryDef,fullArcLength,deform_ODE)

xorg = coord(smesh, geometryDef[0])
yorg = coord(smesh, geometryDef[1])
