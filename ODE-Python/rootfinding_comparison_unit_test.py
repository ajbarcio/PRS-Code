from mpl_toolkits import mplot3d
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from materials import *
from openpyxl import *
from scipy.integrate import solve_ivp as ODE45
from scipy import optimize as op
import time
from stiffness_library import *
from copy import deepcopy as dc
import winsound as ws

deg2rad = np.pi/180
n = 2
finiteDifferenceLength = 0.01
finiteDifferenceAngle  = 5*deg2rad
finiteDifferenceForce  = 0.5
finiteDifferenceTorque = 1
finiteDifferenceCI     = 0.000001

straights = []

fullArcLength = 5.2
globalRes = 200 
globalLen = globalRes+1
globalStep = fullArcLength/globalRes
globalMaxIndex = globalLen-1

E = 27.5*10**6
outPlaneThickness = .375

# profile design variables:
#   R0: inner radius
#   R3: outer radius
#   R1: radius ctrl. pt. 1
#   R2: radius ctrl. pt. 2
#   beta0/betaD: angle of final point
#   betaB:       angle ctrl. pt. 1
#   betaC:       angle ctrl. pt. 2
#   hs:          spline with 3 controlled in-plane thicknesses
#   ctrlHs:      arc-length positions of thicknesses: only ctrlHs[1] is a parameter
#   ctrlLengths: arc-length positions of ctrl points: only ctrlLengths[1] and ctrlLengths[2] are parameters

# Total # of design variables: 13 :((((((((((

globalInnerRadiusLimit = 0.73
globalOuterRadiusLimit = 5.9/2

R0 = 2.1/2
R3 = 5.7/2*.9

R1 = (R0+R3)/2+.27
R2 = (R0+R3)/2-.26

betaB = 152/3*deg2rad*.5
betaC = 2*153/3*deg2rad
betaD = 158*deg2rad
beta0 = betaD

x0 = R0
y0 = 0

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
cIs  = np.array([.008, .004, .008])

ctrlcIs      = np.array([0,fullArcLength*.5,fullArcLength])
ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

maxTorque = 4554.5938
maxDBeta  = 0.087266463

dragVector0 = [R0, R1, R2, R3, \
               betaB, betaC, beta0, \
               cIs[0], cIs[1], cIs[2], \
               ctrlcIs[1], \
               ctrlLengths[1], ctrlLengths[2],
               fullArcLength]
print("initial initial guess", dragVector0)

geometryDef = form_spring(pts, cIs, ctrlLengths, ctrlcIs)
smesh = np.linspace(0,fullArcLength, globalLen)

print("start AB method")
start = time.time()
a = np.empty(len(smesh))
b = np.empty(len(smesh))
hAB = np.empty(len(smesh))
for i in range(len(smesh)):
    AB = a_b_rootfinding(smesh[i], geometryDef[0], geometryDef[1], geometryDef[2], False)
    # print(AB)
    a[i] = AB[0]
    b[i] = AB[1]
    if np.isinf(a[i]) or np.isinf(b[i]):
        hAB[i] = np.cbrt(12*cI_s(smesh[i],geometryDef[2])/outPlaneThickness)
    else:
        hAB[i] = b[i]-a[i]
end = time.time()
print("end AB method")
abtime = end-start
print("AB time:",end-start)

print("start lAB method")
start = time.time()
la = np.empty(len(smesh))
lb = np.empty(len(smesh))
hLALB = np.empty(len(smesh))
lABPrev = [0, 0]
for i in range(len(smesh)):
    lAB = l_a_l_b_rootfinding(smesh[i], lABPrev, geometryDef[0], geometryDef[1], geometryDef[2], False)
    # print(AB)
    la[i] = lAB[0]
    lb[i] = lAB[1]
    hLALB[i] = lb[i]+la[i]
    lABPrev = lAB
end = time.time()
print("end lAB method")
labtime = end-start
print("lAB time:",end-start)
print("ratio:", labtime/abtime)
plt.plot(smesh, hAB, label='old way')
plt.plot(smesh, hLALB, label='new way')
plt.legend()
print(lin.norm(hLALB-hAB))
plt.show()