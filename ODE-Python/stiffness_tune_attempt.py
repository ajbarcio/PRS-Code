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

globalInnerRadiusLimit = 0.75
globalOuterRadiusLimit = 6/2

R0 = 2.2/2
R3 = 5.9/2*.9

R1 = (R0+R3)/2+.26
R2 = (R0+R3)/2-.25

betaB = 150/3*deg2rad*.5
betaC = 2*150/3*deg2rad
betaD = 160*deg2rad
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

discludeVector = [False,False,False,False,False,False,False, \
                  True,True,True, \
                  False,False,False,False]

print("FIRST ATTEMPT; DRAG ONLY IC")
stiffness, res, dragVector, dragVector0 = tune_stiffness(maxTorque/maxDBeta, maxDBeta, dragVector0, discludeVector)
geometryDef, smesh = drag_vector_spring(dragVector)

print("overall relative change", dragVector/dragVector0)
print("target stiffness:", maxTorque/maxDBeta)
print("acheived stiffness:", stiffness)
print("original guess", dragVector0)
print("final guess", dragVector)

# dragVector0 = dc(dragVector)
# assert(np.all(dragVector0==dragVector))
# print("STARTING OVER, USING BEST PREVIOUS TO START")
# stiffness, res, dragVector, dragVector0 = tune_stiffness(maxTorque/maxDBeta, maxDBeta, dragVector0, discludeVector)
# geometryDef, smesh = drag_vector_spring(dragVector)
# # ws.Beep(831,200)
# print("overall relative change", dragVector/dragVector0)
# print("target stiffness:", maxTorque/maxDBeta)
# print("acheived stiffness:", stiffness)
# print("original guess", dragVector0)
# print("final guess", dragVector)

# discludeVector = [True,True,True,True, \
#                   False,False,False,False,False,False,False,False,False,False]
# dragVector0 = dc(dragVector)
# print("SECOND ATTEMPT; DRAG ONLY Radii")
# stiffness, res, dragVector, dragVector0 = tune_stiffness(maxTorque/maxDBeta, maxDBeta, dragVector0, discludeVector)
# geometryDef, smesh = drag_vector_spring(dragVector)

# print("overall relative change", dragVector/dragVector0)
# print("target stiffness:", maxTorque/maxDBeta)
# print("acheived stiffness:", stiffness)
# print("original guess", dragVector0)
# print("final guess", dragVector)

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
# rc = np.empty(len(smesh))
rn = np.empty(len(smesh))
Ic = np.empty(len(smesh))
start = time.time()
for i in range(len(smesh)):
    # rc[i] = rc_rootfinding(smesh[i], geometryDef[0], geometryDef[1], geometryDef[2], True)
    rn[i] = r_n(smesh[i], geometryDef[0], geometryDef[1])
    Ic[i] = cI_s(smesh[i], geometryDef[2])
# e = rc-rn
# end = time.time()
# print("RC Rootfinding", end-start)
# print(e)
# rn = r_n(smesh, geometryDef[0], geometryDef[1])
# index_max = max(range(len(rn[1:-2])), key=rn[1:-2].__getitem__)
# plt.arrow(0,0,xorg[index_max],yorg[index_max])

# for i in range(len(e)):
#     if np.isnan(e[i]):
#         e[i] = 0

# rcx = xorg-e*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
# rcy = yorg+e*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
# # print("sin",np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1])))
# # print("cos",np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1])))
# plt.plot(rcx,rcy, label='rc_rootfinding centroidal axis')

plt.figure(0)
plt.plot(np.transpose(res), label='results')
plt.plot(xorg)
plt.plot(yorg)
plt.legend()
# for i in range(len(smesh)):
#     print(xorg[i], yorg[i], res[1,i], res[2,i])
la = np.empty(len(smesh))
lb = np.empty(len(smesh))
h = np.empty(len(smesh))
start = time.time()
lABPrev = [0, 0]
for i in range(len(smesh)):
    lAB = l_a_l_b_rootfinding(smesh[i], lABPrev, geometryDef[0], geometryDef[1], geometryDef[2], False)
    # print(AB)
    la[i] = lAB[0]
    lb[i] = lAB[1]
    h[i] = lb[i]+la[i]
end=time.time()
ecc = Ic/(outPlaneThickness*h*rn)
print("rootfinding time,", end-start)
# print(ecc)

# nb = h/2+ecc
# na = h/2-ecc

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
plt.plot(smesh, cI_s(smesh, geometryDef[2])*1000)
plt.plot(smesh, h)
plt.plot(smesh,la)
plt.plot(smesh,lb)
plt.plot(smesh,rn)

# for i in range(len(smesh)):
#     print(coord(smesh[i],geometryDef[0]), coord(smesh[i],geometryDef[1]))
res, SF = deform_spring_by_torque(maxTorque, geometryDef)

plt.show()