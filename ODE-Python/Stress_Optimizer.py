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

globalCIList = np.empty(globalLen)

def fixed_rk4_extra(fun, y0, xmesh, *args): # (fun, alphaCoeffs, cICoeffs, y0, Fx, Fy, xmesh)
    step = xmesh[1]-xmesh[0]
    auxRes = np.empty(globalLen)
    if hasattr(y0, '__len__'):
        res = np.empty((len(y0), len(xmesh)))
        for i in range(len(xmesh)):
            # print(i)
            if i == 0:
                for j in range(len(y0)):
                    res[j,i]  = y0[j]
                    auxRes[i] = 
            else:
                stepRes, aux = rk4_step(fun, xmesh[i], y0, step, args)
                # print(stepRes)
                for j in range(len(y0)):
                    res[j,i] = stepRes[j]
                auxRes[i] = aux
                y0 = stepRes
    else:
        res = np.empty((len(xmesh)))
        for i in range(len(xmesh)):
            # print(i)
            if i == 0:
                res[i] = y0
            else:
                stepRes = rk4_step(fun, xmesh[i], y0, step, args)
                res[i] = stepRes
                y0 = stepRes
    # print("gonna return")
    return res, cI