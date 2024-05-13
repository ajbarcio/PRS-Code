import numpy as np
from scipy.integrate import solve_ivp as variable_mesh_solve
import numpy.linalg as lin
import matplotlib.pyplot as plt
import time
from stiffness_library import *
from StatProfiler import SSProfile
from copy import deepcopy as dc
from drawnow import drawnow, figure

def reform_drag_vector(length):
    deg2rad = np.pi/180
    n = 2
    fullArcLength = length

    E = 27.5*10**6
    outPlaneThickness = .375

    globalInnerRadiusLimit = 0.75
    globalOuterRadiusLimit = 6/2

    R0 = 2.2/2
    R3 = 5.9/2

    R1 = (R0+R3)/2
    R2 = (R0+R3)/2

    fullAngleDegrees = 160
    betaD = fullAngleDegrees*deg2rad
    beta0 = betaD

    betaB = fullAngleDegrees/3*deg2rad*.5
    betaC = 2*fullAngleDegrees/3*deg2rad


    x0 = R0
    y0 = 0

    pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
    cIs  = np.array([.008, .0001, .008])

    ctrlcIs      = np.array([0,fullArcLength*.5,fullArcLength])
    ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

    # maxTorque = 4554.5938
    maxTorque = 5000
    maxDBeta  = 0.087266463

    dragVector0 = [R0, R1, R2, R3, \
                betaB, betaC, beta0, \
                cIs[0], cIs[1], cIs[2], \
                ctrlcIs[1], \
                ctrlLengths[1], ctrlLengths[2],
                fullArcLength]
    
    return dragVector0

def xi_s_unification(dragVector0):
    err = 1
    dragVector = dc(dragVector0)
    length = dragVector[-1]
    while err>10**-6:
        # get error at point
        geometryDef, smesh = drag_vector_spring(dragVector)
        dxdxi = d_coord_d_s(smesh, geometryDef[0])
        dydxi = d_coord_d_s(smesh, geometryDef[1])
        dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
        ref = np.ones(len(dxids))
        err = lin.norm(dxids-ref)
        print(err)
        # get derivative
        dragVectorF = dc(dragVector)
        print(dragVectorF)
        lengthF = length+1
        dragVectorF = reform_drag_vector(lengthF)
        print(dragVectorF)
        geometryDef, smesh = drag_vector_spring(dragVector)
        dxdxiF = d_coord_d_s(smesh, geometryDef[0])
        dydxiF = d_coord_d_s(smesh, geometryDef[1])
        dxidsF = 1/np.sqrt(dxdxiF**2+dydxiF**2)
        errF = lin.norm(dxidsF-ref)
        print(errF)
        Jac = (errF-err)/0.01
        print(Jac)
        length = length-err/Jac
        dragVector = reform_drag_vector(length)
        print("got through")
    return dragVector

deg2rad = np.pi/180
n = 2
fullArcLength = 3

E = 27.5*10**6
outPlaneThickness = .375

globalInnerRadiusLimit = 0.75
globalOuterRadiusLimit = 6/2

R0 = 2.2/2
R3 = 5.9/2

R1 = (R0+R3)/2
R2 = (R0+R3)/2

fullAngleDegrees = 160
betaD = fullAngleDegrees*deg2rad
beta0 = betaD

betaB = fullAngleDegrees/3*deg2rad*.5
betaC = 2*fullAngleDegrees/3*deg2rad

x0 = R0
y0 = 0

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
cIs  = np.array([.008, .0001, .008])

ctrlcIs      = np.array([0,fullArcLength*.5,fullArcLength])
ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

# maxTorque = 4554.5938
maxTorque = 5000
maxDBeta  = 0.087266463

dragVector0 = [R0, R1, R2, R3, \
               betaB, betaC, beta0, \
               cIs[0], cIs[1], cIs[2], \
               ctrlcIs[1], \
               ctrlLengths[1], ctrlLengths[2],
               fullArcLength]

# dragVector0 = [1.10000000e+00, 2.13750000e+00, 1.62750000e+00, 2.65500000e+00,
#  4.36332313e-01, 1.74532925e+00, 2.79252680e+00, 
#  4.74085857e-03, 7.40837797e-04, 4.74084339e-03, 
#  2.60000000e+00,
#  1.73160000e+00, 3.46840000e+00, 
#  3]
# fullArcLength = dragVector0[-1]
# dragVector0 = [1.10000000e+00, 2.13750000e+00, 1.62750000e+00, 2.65500000e+00,
#  4.36332313e-01, 1.74532925e+00, 2.79252680e+00, 
#  2*4.74085857e-03, 7.40837797e-04, 5*4.74085857e-03, 
#  2.60000000e+00,
#  1.73160000e+00, 3.46840000e+00, 
#  5.20000000e+00]

# dragVector0 = xi_s_unification(dragVector0)
# print(dragVector0[-1])
# fullArcLength = dragVector0[-1]

geometryDef, smesh = drag_vector_spring(dragVector0)

start = 0.01
end = 6
lengthMesh = np.linspace(start,end,1001)

errVec = np.ones(len(lengthMesh))

step = lengthMesh[1]-lengthMesh[0]
fig = plt.figure("dxids")
def draw_the_damn_thing():
    plt.plot(smesh,dxids,color=c)

for i in lengthMesh:
    dragVectorCycle = reform_drag_vector(i)
    geometryDef, smesh = drag_vector_spring(dragVectorCycle)

    dxdxi = d_coord_d_s(smesh, geometryDef[0])
    dydxi = d_coord_d_s(smesh, geometryDef[1])
    dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
    cmap=plt.get_cmap('rainbow')
    c=cmap(i/(end-start))
    # plt.plot(smesh,dxids,color=c)
    drawnow(draw_the_damn_thing)

    ref = np.ones(len(dxids))
    err = lin.norm(dxids-ref)
    print(err)
    errVec[int((i-1)/step)] = err
plt.figure("arc length parameter tuning")
plt.plot(lengthMesh, errVec)




dxdxi = d_coord_d_s(smesh, geometryDef[0])
dydxi = d_coord_d_s(smesh, geometryDef[1])
dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
dxds  = dxdxi*dxids
dyds  = dydxi*dxids

ref = np.ones(len(dxids))
err = lin.norm(dxids-ref)
print("xi/s err",err)

plt.figure("derivatives")
plt.plot(smesh, dxdxi, label = "dxdxi")
plt.plot(smesh, dydxi, label = "dydxi")
plt.plot(smesh, dxds, label = "dxds")
plt.plot(smesh, dyds, label = "dyds")
plt.plot(smesh,dxids, label="dxids")
plt.legend()
plt.figure("angles")
plt.plot(smesh, np.arctan2(dydxi, dxdxi), label="angle of xi")
plt.plot(smesh, np.arctan2(dyds, dxds), label="angle of s")
plt.legend()

rn = np.empty(len(smesh))
Ic = np.empty(len(smesh))
start = time.time()
for i in range(len(smesh)):
    rn[i] = r_n(smesh[i], geometryDef[0], geometryDef[1])
    Ic[i] = cI_s(smesh[i], geometryDef[2])
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
ecc = Ic/(outPlaneThickness*h*rn)
print("rootfinding time,", end-start)

cICoeffs = geometryDef[2]
xCoeffs = geometryDef[0]
yCoeffs=geometryDef[1]

plt.figure("polynomials")
plt.plot(smesh, coord(smesh,xCoeffs),label="x(xi)")
plt.plot(smesh, coord(smesh,yCoeffs),label="y(xi)")
# plt.plot(smesh, cI_s(smesh, cICoeffs),label="cI(xi)")
plt.plot(smesh, d_coord_d_s(smesh,xCoeffs),label="dxdxi")
plt.plot(smesh, d_coord_d_s(smesh,yCoeffs),label="dydxi")

plt.plot(smesh, d2_coord_d_s2(smesh,xCoeffs),label="d2xdxi2")
plt.plot(smesh, d2_coord_d_s2(smesh,yCoeffs),label="d2ydxi2")
plt.plot(smesh, rn, label="rn")
plt.legend()

print(d2_coord_d_s2(fullArcLength,xCoeffs), d2_coord_d_s2(fullArcLength,yCoeffs))

plt.figure("geometry resultants")
plt.plot(smesh, h,label="thicknesses")
plt.plot(smesh, ecc,label="eccentricities")
plt.legend()


xorg = coord(smesh, geometryDef[0])
yorg = coord(smesh, geometryDef[1])

xb = -lb*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
xa = -la*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
yb = lb*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
ya = la*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

xrc = ecc*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
yrc = ecc*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

plt.figure("geometry results")
# plt.xlim(-3.5,3.5)
plt.plot(xorg, yorg)
theta = np.linspace(0, 2*np.pi, 100)
outerCircleX = R3*np.cos(theta)
outerCircleY = R3*np.sin(theta)
innerCircleX = R0*np.cos(theta)
innerCircleY = R0*np.sin(theta)
plt.plot(outerCircleX,outerCircleY)
plt.plot(innerCircleX,innerCircleY)
plt.axis('equal')
# plt.xlim(-3.5,3.5)
plt.figure("geometry results")
plt.plot(xorg+xb,yorg+yb)
plt.plot(xorg-xa,yorg-ya)
plt.plot(-(xorg+xb),-(yorg+yb))
plt.plot(-(xorg-xa),-(yorg-ya))
plt.plot(xorg+xrc, yorg+yrc, label="AB rootfinding centroidal axis")
plt.legend()
print(xCoeffs)
plt.show()