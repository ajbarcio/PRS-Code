from mpl_toolkits import mplot3d
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp as ODE45
from scipy import optimize as op
import os
import random
import math
import time
from copy import deepcopy as dc
from StatProfiler import SSProfile
import matplotlib.collections as mcoll
import matplotlib.path as mpath
from stiffness_library import *

deg2rad = np.pi/180

maxTorque = 5000
maxDBeta  = 0.087266463
maxStress = 270000

lABPrevOuter = [0,0]
# maxTorque = 4554.5938

dragVector0 = [R0, R1, R2, R3, \
               betaB, betaC, beta0, \
               cIs[0], cIs[1], cIs[2], \
               ctrlcIs[1], \
               ctrlLengths[1], ctrlLengths[2],
               fullArcLength]

geometryDef , smesh = drag_vector_spring(dragVector0)
cI_dumper = np.empty(len(smesh))

def l_a_l_b_rootfinding_with_direct_params(lABPrev, rn, cI, printBool):
    def func(x, rn, cI):
        f1 = (x[0]+x[1])/(np.log((rn+x[1])/(rn-x[0])))-rn
        # f2 = (x[0]+x[1])*(outPlaneThickness*(x[1]-(x[1]+x[0])/2)*rn)-cI
        f2 = outPlaneThickness*rn*(x[0]+x[1])*(x[1]/2-x[0]/2)-cI
        return np.array([f1, f2])
    def jac(x, rn, cI):
        return np.array([[1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn-x[0])), \
                          1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn+x[1]))], \
                         [-rn*outPlaneThickness*x[0], rn*outPlaneThickness*x[1]]])
    # print(rn)
    if not np.isinf(rn):
        if lABPrev[0]==lABPrev[1]:
            x0 = [lABPrev[0], lABPrev[1]+0.001]
        else:
            x0 = lABPrev
        err = 1
    else:
        err = 0
        l = np.cbrt(12*cI/outPlaneThickness)/2
        x0 = [l, l]
        # print("entered")
    # if (x0[0]>=abs(rn)):
    #     x0[0] = abs(rn)-10e-5
    x = x0
    # print("x0:",x0)
    iii = 0
    then = time.monotonic()
    print("gonna start rootfinding                                                           ", end="\r")
    while err > 10**-6 and iii <500 and lin.norm(x)<1:
        # print("entered while")
        xprev = x
        # if (x[0]>=abs(rn)):
        #     x[0] = abs(rn)-10e-5
        x = x - np.transpose(lin.inv(jac(x, rn, cI)).dot(func(x, rn, cI)))
        # print(x)
        err = lin.norm(x-xprev)
        if not all(np.isfinite(x)):
            print(x, xprev, rn, iii)
            x = xprev
            err = 0
        assert(all(np.isfinite(x)))
        # err = lin.norm(func(x, rn, cI))
        iii+=1
    print("done rootfinding                                                                           ", end="\r")
    if(lin.norm(x)>lin.norm(lABPrev)*100):
        l = np.cbrt(12*cI/outPlaneThickness)/2
        x = [l,l]
    if(printBool):
        print(x0)
        print(x)
        print(iii)
        print(rn)
        print(err)
        if iii > 2999:
            print("--------------------DID NOT CONVERGE-------------------------")

    # print("DONE WITH S = ",s)
    return x    # here x is [la, lb]

def ODE_with_rootfinding_for_stress(s, p, *args):
    # print("ODE called")
    Fx          = args[0][0][0][0]
    Fy          = args[0][0][0][1]
    geometryDef = args[0][0][0][2]
    dgds0       = args[0][0][0][3]
    # print("pass check:", dgds0)

    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]

    # Use pretty large guess for cI at first
    cI_initial = 0.005
    cI = cI_initial

    finiteDifferenceCI=0.0005

    stressErr = 1
    iiii = 0
    print("gonna start rootfinding for stress EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE", end="\r")
    while abs(stressErr) > 10**-6 and iiii<50:
        Mdim = E*cI        
        dxds = d_coord_d_s(s, xCoeffs)
        dyds = d_coord_d_s(s, yCoeffs)
        # evaluate diffEQ
        LHS, stressErr = diffEq(s, p, xCoeffs, yCoeffs, cI, dgds0, Fy, Fx, Mdim, dyds, dxds)
        # approsimate the derivative
        LHSG, stressErrF = diffEq(s, p, xCoeffs, yCoeffs, cI+finiteDifferenceCI, dgds0, Fy, Fx, Mdim, dyds, dxds)
        LHSG, stressErrB = diffEq(s, p, xCoeffs, yCoeffs, cI-finiteDifferenceCI, dgds0, Fy, Fx, Mdim, dyds, dxds)
        derivative = (stressErrF-stressErrB)/(finiteDifferenceCI*2)
        if derivative==0:
            pass
        else:
            # determine next cI:
            cI = cI - stressErr/derivative
        if not np.isfinite(derivative):
            print("something wrong:",stressErr, derivative)
            print(stressErr, stressErrB, stressErrF)
        if iiii>45:
            print(iiii,"--------------------------------------")
        iiii+=1
    print("done rootfinding for stress  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF")
    # print ("\033[A\033[A", end="\n")
    # print("i in ODE function:",int(np.floor(s/step)))
    cI_dumper[int(np.floor(s/step))] = cI
    # print("made it throught once")
    print(cI,"                                   ", end="\r")
    return LHS

def diffEq(s, p, xCoeffs, yCoeffs, cI, dgds0, Fy, Fx, Mdim, dyds, dxds):
    LHS = np.empty(3)

    xcoord0 = coord(0,xCoeffs)
    ycoord0 = coord(0,yCoeffs)

    LHS[0] = dgds0
    LHS[0] = LHS[0] + Fy/Mdim*(p[1]-xcoord0) - Fx/Mdim*(p[2]-ycoord0)
    if s==0:
        # print(LHS[0],dgds0)
        assert(np.isclose(LHS[0],dgds0,rtol=1e-3))
        global lABPrevOuter 
        lABPrevOuter = [0,0]
    LHS[1] = np.cos(np.arctan2(dyds,dxds)+p[0]) # -(g*y +/- (g^2 - y^2 + 1)^(1/2))/(g^2 + 1)
    LHS[2] = np.sin(np.arctan2(dyds,dxds)+p[0])

    if np.sin(np.arctan2(dyds,dxds))==0:
        LHS[1] = LHS[1]*dxds/np.cos(np.arctan2(dyds,dxds))
        LHS[2] = LHS[2]*dxds/np.cos(np.arctan2(dyds,dxds))
    else:
        LHS[1] = LHS[1]*dyds/np.sin(np.arctan2(dyds,dxds))
        LHS[2] = LHS[2]*dyds/np.sin(np.arctan2(dyds,dxds))
    
    # get outer geometry:
    rn = r_n(s, xCoeffs, yCoeffs)
    
    la, lb = l_a_l_b_rootfinding_with_direct_params(lABPrevOuter,rn,cI,False)
    lABPrevOuter = [la, lb]
    # if not (np.isfinite((la))):
    #     print(la,lb)
    # print(lABPrevOuter)
    a = rn-la
    b = rn+lb
    #use for stress
    if np.isfinite(rn):
        stressInner = abs(E*(1-rn/a)*rn*LHS[0])
        stressOuter = abs(E*(1-rn/b)*rn*LHS[0])
        stress      = max(stressInner, stressOuter)
        stressErr   = stress-maxStress
    else:
        stressInner = abs(E*LHS[0]*la)
        stress      = stressInner
        stressErr   = stress-maxStress
    if np.isnan(stressInner):
        stressInner = 0
        stressOuter = 0
        stressErr   = 0
    return LHS, stressErr

n = 2
fullArcLength = 5

E = 27.5*10**6
outPlaneThickness = .375

globalInnerRadiusLimit = 0.75
globalOuterRadiusLimit = 6/2

R0 = 2.2/2
R3 = 5.9/2

R1 = (R0+R3)/2-.0125
R2 = (R0+R3)/2

betaB = 160/3*deg2rad*.8
betaC = 2*160/3*deg2rad
betaD = 160*deg2rad
beta0 = betaD

x0 = R0
y0 = 0

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
cIs  = np.array([.008, .001, .018])

ctrlcIs      = np.array([0,fullArcLength*.5,fullArcLength])
ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

# maxTorque = 4554.5938

dragVector0 = [R0, R1, R2, R3, \
               betaB, betaC, beta0, \
               cIs[0], cIs[1], cIs[2], \
               ctrlcIs[1], \
               ctrlLengths[1], ctrlLengths[2],
               fullArcLength]

geometryDef , smesh = drag_vector_spring(dragVector0)
step = smesh[1]-smesh[0]
cI_dumper = np.zeros(len(smesh))
rn = r_n(smesh,geometryDef[0],geometryDef[1])

dxdxi = d_coord_d_s(smesh, geometryDef[0])
dydxi = d_coord_d_s(smesh, geometryDef[1])
dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
dxds  = dxdxi*dxids
dyds  = dydxi*dxids

res, sF, divergeFlag = deform_spring_by_torque(5000,geometryDef,fullArcLength,ODE_with_rootfinding_for_stress)

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

la = np.empty(len(smesh))
lb = np.empty(len(smesh))
h = np.empty(len(smesh))
start = time.time()
lABPrev = [0, 0]
for i in range(len(smesh)):
    # print(i)
    lAB = l_a_l_b_rootfinding_with_direct_params(lABPrev, rn[i], cI_dumper[i], False)
    # print(AB)
    la[i] = lAB[0]
    lb[i] = lAB[1]
    h[i] = lb[i]+la[i]
    lABPrev = lAB
end=time.time()
ecc = cI_dumper/(outPlaneThickness*h*rn)
print("rootfinding time,", end-start)
# print(smesh)
# print(rn)'
a = rn-la
b = rn+lb



print(cI_dumper)
print(res)