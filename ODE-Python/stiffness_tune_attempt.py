from mpl_toolkits import mplot3d
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp as ODE45
from scipy import optimize as op
import time
from stiffness_library import *
from copy import deepcopy as dc
import winsound as ws # comment out these lines if on linux

deg2rad = np.pi/180
n = 2
finiteDifferenceLength = 0.01
finiteDifferenceAngle  = 5*deg2rad
ffForce  = 0.5
finiteDifferenceTorque = 1
finiteDifferenceCI     = 0.000001

straights = []

fullArcLength = 5.2
globalRes = 200 
globalLen = globalRes+1
globalStep = fullArcLength/globalRes
globalMaxIndex = globalLen-1

E = 27.5*10**6
self.t = .375

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

globalInnerRadiusLimit = 0.75
globalOuterRadiusLimit = 6/2

R0 = 2.2/2
R3 = 5.9/2

R1 = (R0+R3)/2*1.3
R2 = (R0+R3)/2*1.1

funnyAngle = 160

betaB = funnyAngle/3*deg2rad*.7
betaC = 2*funnyAngle/3*deg2rad*1.2
betaD = funnyAngle*deg2rad
beta0 = betaD


x0 = R0
y0 = 0

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
cIs  = np.array([.01, .002, .01])

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

# geometryDef = form_spring(pts, cIs, ctrlLengths, ctrlcIs)
# res = deform_spring_by_torque(maxTorque/2, geometryDef)
# print(res[0,-1])
# res = deform_spring_by_torque(maxTorque, geometryDef)
# print(res[0,-1])
# smesh = np.linspace(0,fullArcLength, globalLen)
# for value in res[0]:
#     print(value)
# print(geometryDef[2])
# assert(False)
# ws.Beep(831,200)

discludeVector = [ 0, 0, 0, 0,       # Gains for Radii
                   0, 0, 0,          # Gains for Angles 
                   1, 1, 1,          # Gains for Control cIs
                   0,                # Gain for cI minimum point
                   0, 0,             # Gains for length control points
                   0      ]      # Gain for overall length]

print("FIRST ATTEMPT; DRAG ONLY IC")
print("target stiffness:", maxTorque/maxDBeta)
stiffness, res, dragVector, dragVector0 = tune_stiffness(maxTorque/maxDBeta, maxDBeta, dragVector0, discludeVector)

# stiffness, maxStress, res, dragVector, dragVector0 = stress_stiffness_tuning(maxTorque/maxDBeta, maxDBeta, 247760, dragVector0, discludeVector)
geometryDef, smesh = drag_vector_spring(dragVector)
lAB0 = np.cbrt(12*PPoly_Eval(0, geometryDef[2])/self.t)/2

print("overall relative change", dragVector-dragVector0)
print("target stiffness:", maxTorque/maxDBeta)
print("acheived stiffness:", stiffness)
print("original guess", dragVector0)
print("final guess", dragVector)

ws.Beep(831,200)

assert(not violates_bounds(dragVector))

xorg = coord(smesh, geometryDef[0])
yorg = coord(smesh, geometryDef[1])

R0 = dragVector[0]
R1 = dragVector[1]
R2 = dragVector[2]
R3 = dragVector[3]
fullArcLength = dragVector[-1]
# print("geometryDef", geometryDef)

plt.figure("geometry results")
plt.plot(xorg, yorg)
plt.plot(res[1,:],res[2,:])
theta = np.linspace(0, 2*np.pi, 100)
outerCircleX = R3*np.cos(theta)
outerCircleY = R3*np.sin(theta)
innerCircleX = R0*np.cos(theta)
innerCircleY = R0*np.sin(theta)
plt.plot(outerCircleX,outerCircleY)
plt.plot(innerCircleX,innerCircleY)
plt.axis('equal')
rn = np.empty(len(smesh))
Ic = np.empty(len(smesh))
start = time.time()
for i in range(len(smesh)):
    rn[i] = r_n(smesh[i], geometryDef[0], geometryDef[1])
    Ic[i] = PPoly_Eval(smesh[i], geometryDef[2])
plt.figure(0)
plt.plot(np.transpose(res), label='results')
plt.plot(xorg)
plt.plot(yorg)
plt.legend()
la = np.empty(len(smesh))
lb = np.empty(len(smesh))
h = np.empty(len(smesh))
start = time.time()
lABPrev = [0, 0]
for i in range(len(smesh)):
    # print(i)
    lAB = l_a_l_b_rootfinding(smesh[i], lABPrev, geometryDef[0], geometryDef[1], geometryDef[2], False)
    # print(AB)
    la[i] = lAB[0]
    lb[i] = lAB[1]
    h[i] = lb[i]+la[i]
    lABPrev = lAB
end=time.time()
ecc = Ic/(self.t*h*rn)
print("rootfinding time,", end-start)
# print(smesh)
# print(rn)

xb = -lb*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
xa = -la*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
yb = lb*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
ya = la*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

xrc = ecc*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
yrc = ecc*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

# print(xrc, yrc)

plt.figure("geometry results")
plt.plot(xorg+xb,yorg+yb)
plt.plot(xorg-xa,yorg-ya)
plt.plot(-(xorg+xb),-(yorg+yb))
plt.plot(-(xorg-xa),-(yorg-ya))
plt.plot(xorg+xrc, yorg+yrc, label="AB rootfinding centroidal axis")
plt.legend()

plt.figure(99)
plt.plot(smesh, PPoly_Eval(smesh, geometryDef[2])*1000)
plt.plot(smesh, h)
plt.plot(smesh,la)
plt.plot(smesh,lb)
plt.plot(smesh,rn)

# res, SF = deform_spring_by_torque(maxTorque, geometryDef)

plt.show()