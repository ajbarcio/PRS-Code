import numpy as np
from scipy.integrate import solve_ivp as variable_mesh_solve
import numpy.linalg as lin
import matplotlib.pyplot as plt
import time
from stiffness_library import *
from StatProfiler import SSProfile

then = time.monotonic()

deg2rad = np.pi/180
n = 2
fullArcLength = 5.2

E = 27.5*10**6
outPlaneThickness = .375

# These won't change ^^^

# globalInnerRadiusLimit = 0.75
# globalOuterRadiusLimit = 6/2

# R0 = 2.2/2
# R3 = 5.9/2

# R1 = (R0+R3)/2
# R2 = (R0+R3)/2

# funnyAngle = 165

# betaB = funnyAngle/3*deg2rad
# betaC = 2*funnyAngle/3*deg2rad
# betaD = funnyAngle*deg2rad
# beta0 = betaD

# x0 = R0
# y0 = 0

# pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
# cIs  = np.array([.001, .0001, .001])

# ctrlcIs      = np.array([0,fullArcLength*.5,fullArcLength])
# ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

# # maxTorque = 4554.5938
# maxTorque = 0
# maxDBeta  = 0.087266463

# dragVector0 = [R0, R1, R2, R3, \
#                betaB, betaC, beta0, \
#                cIs[0], cIs[1], cIs[2], \
#                ctrlcIs[1], \
#                ctrlLengths[1], ctrlLengths[2],
#                fullArcLength]

# geometryDef, smesh = drag_vector_spring(dragVector0)

# dxdxi = d_coord_d_s(smesh, geometryDef[0])
# dydxi = d_coord_d_s(smesh, geometryDef[1])
# dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
# dxds  = dxdxi*dxids
# dyds  = dydxi*dxids

# res, SF, divergeFlag = deform_spring_by_torque(4554.5938,geometryDef,fullArcLength,deform_ODE)
# print(SF)

# xorg = coord(smesh, geometryDef[0])
# yorg = coord(smesh, geometryDef[1])

# beta0 = np.arctan2(yorg[-1],xorg[-1])
# alpha0 = alpha_xy(fullArcLength, geometryDef[0], geometryDef[1])
# print("some stuff,",beta0+alpha0)
# print(beta0)
# print(alpha0)
# print(np.arctan2(d_coord_d_s(fullArcLength,geometryDef[1]),d_coord_d_s(fullArcLength,geometryDef[0])))

# print("fullArcLength Outside:", fullArcLength)
# dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
# print("graphed dbeta:",dBeta*180/np.pi)
# print("graphed gamma final:",res[0,-1]*180/np.pi)
# print(res[0,-1])

# R0 = dragVector0[0]
# R1 = dragVector0[1]
# R2 = dragVector0[2]
# R3 = dragVector0[3]
# fullArcLength = dragVector0[-1]

# plt.figure("geometry results")
# plt.plot(xorg, yorg)
# plt.plot(res[1,:],res[2,:])
# theta = np.linspace(0, 2*np.pi, 100)
# outerCircleX = R3*np.cos(theta)
# outerCircleY = R3*np.sin(theta)
# innerCircleX = R0*np.cos(theta)
# innerCircleY = R0*np.sin(theta)
# plt.plot(outerCircleX,outerCircleY)
# plt.plot(innerCircleX,innerCircleY)
# plt.axis('equal')

# rn = r_n(smesh,geometryDef[0],geometryDef[1])
# Ic = cI_s(smesh,geometryDef[2])

# la = np.empty(len(smesh))
# lb = np.empty(len(smesh))
# h = np.empty(len(smesh))
# start = time.time()
# lABPrev = [0, 0]
# for i in range(len(smesh)):
#     lAB = l_a_l_b_rootfinding(smesh[i], lABPrev, geometryDef[0], geometryDef[1], geometryDef[2], False)

#     la[i] = lAB[0]
#     lb[i] = lAB[1]
#     h[i] = lb[i]+la[i]
#     lABPrev = lAB
# end=time.time()
# ecc = Ic/(outPlaneThickness*h*rn)
# print("rootfinding time,", end-start)

# a = rn-la
# b = rn+lb

# dgdxi = np.zeros(len(res[0,:]))
# dgds = np.zeros(len(res[0,:]))
# step = smesh[1]-smesh[0]
# for i in range(len(smesh)):
#     if i==0:
#         dgdxi[i] = (res[0,i+1]-res[0,i])/step
#         dgds[i] = dgdxi[i]*dxids[i]
#     if i==len(dgdxi)-1:
#         dgdxi[i] = (res[0,i]-res[0,i-1])/step
#         dgds[i] = dgdxi[i]*dxids[i]
#     else:
#         dgdxi[i] = (res[0,i+1]-res[0,i-1])/(2*step)
#         dgds[i] = dgdxi[i]*dxids[i]

# allowableStress = 270000
# stressInner = abs(E*(1-rn/a)*rn*dgds) # O SHIT THIS IS right??
# stressOuter = abs(E*(1-rn/b)*rn*dgds)
# maxInnerStress = np.nanmax(stressInner)
# maxOuterStress = np.nanmax(stressOuter)

# maxInnerStress = np.nanmax(stressInner)
# maxOuterStress = np.nanmax(stressOuter)

# nomalizedInnerStress = stressInner/allowableStress
# nomalizedOuterStress = stressOuter/allowableStress


# maxStress = np.nanmax([stressInner, stressOuter])
# print(maxStress)

# xb = -lb*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
# xa = -la*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
# yb = lb*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
# ya = la*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

# xrc = ecc*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
# yrc = ecc*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

# xbdef = -lb*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1])+res[0,:])
# xadef = -la*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1])+res[0,:])
# ybdef = lb*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1])+res[0,:])
# yadef = la*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1])+res[0,:])

# figure = plt.figure("geometry results")

# colorline(res[1,:]+xbdef,res[2,:]+ybdef,nomalizedInnerStress,cmap=plt.get_cmap('rainbow'))
# colorline(res[1,:]-xadef,res[2,:]-yadef,nomalizedOuterStress,cmap=plt.get_cmap('rainbow'))

# colorline(-(res[1,:]+xbdef),-(res[2,:]+ybdef),nomalizedInnerStress,cmap=plt.get_cmap('rainbow'))
# colorline(-(res[1,:]-xadef),-(res[2,:]-yadef),nomalizedOuterStress,cmap=plt.get_cmap('rainbow'))

# colorline(np.ones(101)*4,np.linspace(-1,1,101),np.linspace(0,1,101),cmap=plt.get_cmap('rainbow'),linewidth=10)

# plt.legend()

# now = time.monotonic()
# print("Total time, spring 1:", now-then)

#### SPRING 2

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

dxdxi = d_coord_d_s(smesh, geometryDef[0])
dydxi = d_coord_d_s(smesh, geometryDef[1])
dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
dxds  = dxdxi*dxids
dyds  = dydxi*dxids


res, SF, divergeFlag = deform_spring_by_torque(4554.5938,geometryDef,fullArcLength,deform_ODE)
print(SF)

xorg = coord(smesh, geometryDef[0])
yorg = coord(smesh, geometryDef[1])

beta0 = np.arctan2(yorg[-1],xorg[-1])
alpha0 = alpha_xy(fullArcLength, geometryDef[0], geometryDef[1])
print("some stuff,",beta0+alpha0)
print(beta0)
print(alpha0)
print(np.arctan2(d_coord_d_s(fullArcLength,geometryDef[1]),d_coord_d_s(fullArcLength,geometryDef[0])))

print("fullArcLength Outside:", fullArcLength)
dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
print("graphed dbeta:",dBeta*180/np.pi)
print("graphed gamma final:",res[0,-1]*180/np.pi)
print(res[0,-1])

R0 = dragVector0[0]
R1 = dragVector0[1]
R2 = dragVector0[2]
R3 = dragVector0[3]
fullArcLength = dragVector0[-1]

plt.figure("geometry results 2")
# plt.plot(xorg, yorg)
plt.plot(res[1,:],res[2,:])
theta = np.linspace(0, 2*np.pi, 100)
outerCircleX = R3*np.cos(theta)
outerCircleY = R3*np.sin(theta)
innerCircleX = R0*np.cos(theta)
innerCircleY = R0*np.sin(theta)
plt.plot(outerCircleX,outerCircleY,'k')
plt.plot(innerCircleX,innerCircleY,'k')
plt.axis('equal')

rn = r_n(smesh,geometryDef[0],geometryDef[1])
Ic = cI_s(smesh,geometryDef[2])

plt.figure("funny")
plt.plot(smesh, Ic)

la = np.empty(len(smesh))
lb = np.empty(len(smesh))
h = np.empty(len(smesh))
start = time.time()
lABPrev = [0, 0]
for i in range(len(smesh)):
    lAB = l_a_l_b_rootfinding(smesh[i], lABPrev, geometryDef[0], geometryDef[1], geometryDef[2], False)

    la[i] = lAB[0]
    lb[i] = lAB[1]
    h[i] = lb[i]+la[i]
    lABPrev = lAB
end=time.time()
ecc = Ic/(outPlaneThickness*h*rn)
print("rootfinding time,", end-start)

a = rn-la
b = rn+lb

dgdxi = np.zeros(len(res[0,:]))
dgds = np.zeros(len(res[0,:]))
step = smesh[1]-smesh[0]
for i in range(len(smesh)):
    if i==0:
        dgdxi[i] = (res[0,i+1]-res[0,i])/step
        dgds[i] = dgdxi[i]*dxids[i]
    if i==len(dgdxi)-1:
        dgdxi[i] = (res[0,i]-res[0,i-1])/step
        dgds[i] = dgdxi[i]*dxids[i]
    else:
        dgdxi[i] = (res[0,i+1]-res[0,i-1])/(2*step)
        dgds[i] = dgdxi[i]*dxids[i]

allowableStress = 270000
stressInner = abs(E*(1-rn/a)*rn*dgds) # O SHIT THIS IS right??
stressOuter = abs(E*(1-rn/b)*rn*dgds)
maxInnerStress = np.nanmax(stressInner)
maxOuterStress = np.nanmax(stressOuter)

maxInnerStress = np.nanmax(stressInner)
maxOuterStress = np.nanmax(stressOuter)

normalizedInnerStress = stressInner/allowableStress
normalizedOuterStress = stressOuter/allowableStress

maxStressList = np.empty(len(normalizedInnerStress))
for i in range(len(normalizedInnerStress)):
    if normalizedInnerStress[i] > normalizedOuterStress[i]:
        maxStressList[i] = normalizedInnerStress[i]
    else:
        maxStressList[i] = normalizedOuterStress[i]
maxStress = np.nanmax([stressInner, stressOuter])
print(maxStress)

xb = -lb*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
xa = -la*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
yb = lb*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
ya = la*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

xrc = ecc*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
yrc = ecc*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

xbdef = -lb*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1])+res[0,:])
xadef = -la*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1])+res[0,:])
ybdef = lb*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1])+res[0,:])
yadef = la*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1])+res[0,:])



figure = plt.figure("geometry results 2")

colorline(res[1,:],res[2,:],maxStressList,cmap=plt.get_cmap('rainbow'))
colorline(-res[1,:],-res[2,:],maxStressList,cmap=plt.get_cmap('rainbow'))

colorline(res[1,:]+xbdef,res[2,:]+ybdef,normalizedInnerStress,cmap=plt.get_cmap('rainbow'))
colorline(res[1,:]-xadef,res[2,:]-yadef,normalizedOuterStress,cmap=plt.get_cmap('rainbow'))

colorline(-(res[1,:]+xbdef),-(res[2,:]+ybdef),normalizedInnerStress,cmap=plt.get_cmap('rainbow'))
colorline(-(res[1,:]-xadef),-(res[2,:]-yadef),normalizedOuterStress,cmap=plt.get_cmap('rainbow'))

colorline(np.ones(101)*4,np.linspace(-1,1,101),np.linspace(0,1,101),cmap=plt.get_cmap('rainbow'),linewidth=10)

# ax = plt.gca()
# ax.set_facecolor('#cad3e3')
# plt.legend()

now = time.monotonic()
print("Total time, spring 1:", now-then)

print()
plt.show()