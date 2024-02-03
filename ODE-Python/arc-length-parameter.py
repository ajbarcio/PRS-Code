import time
import math
import numpy as np
import scipy as sp
from scipy import integrate
from scipy import optimize as op
import numpy.linalg as lin
import scipy.linalg as slin
import matplotlib.pyplot as plt
import matplotlib.transforms as tfm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from materials import *
from openpyxl import *
from scipy.integrate import solve_ivp as ODE45

deg2rad = np.pi/180
E = 10000

globalRes = 100
globalLen = globalRes+1
globalStep = 1/globalRes
globalMaxIndex = globalLen-1
smod           = 2

p0 = np.array([0, .001])
p1 = np.array([90*deg2rad, .001])
p2 = np.array([-45*deg2rad, .001])
p3 = np.array([-90*deg2rad, .001])

outPlaneThickness = 0.375

def get_coeffs():
    D = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
    P = np.array([p0,p1,p2,p3])
    coeffs = D.dot(P)
    print(coeffs)
    return coeffs

def get_metaphor_coords(s, coeffs):
    U = np.array([s**3, s**2, s, 1])
    metaphorCoords = U.dot(coeffs)
    return metaphorCoords

def alpha_s(s, coeffs):
    # print(get_metaphor_coords(s, coeffs)[0])
    return get_metaphor_coords(s, coeffs)[0]

def genStiffness_s(s, coeffs):
    # print(get_metaphor_coords(s, coeffs)[1])
    return get_metaphor_coords(s, coeffs)[1]

def d_alpha_d_s(s, coeffs):
    dAlphadS = coeffs[0,0]*s**2+coeffs[1,0]*s+coeffs[2,0]
    return dAlphadS

def d_genStiffness_d_s(s, coeffs):
    dGenStiffnessdS = coeffs[0,1]*s**2+coeffs[1,1]*s+coeffs[2,1]
    return dGenStiffnessdS

def r_n(s, coeffs):
    rn = 1/d_alpha_d_s(s, coeffs)
    return rn

def b_a_rootfinding(s, coeffs, printBool):
    rn =  r_n(s, coeffs)
    genStiff = genStiffness_s(s, coeffs)

    c = np.log(1.2)/0.2

    def func(x, rn, genStiff):
        # x0 = rc, x1 = h
        f1 = (x[1]-x[0])/(np.log(x[1]/x[0]))-rn
        # f1 = f1**2
        f2 = outPlaneThickness*(x[1]-x[0])*((x[1]+x[0])/2-rn)*rn-genStiff
        # f2 = f2**2
        # return lin.norm(np.array([f1, f2]))
        return np.array([f1, f2])
    def jac(x, rn, genStiff):
        return np.array([[(x[1]-x[0])/(x[0]*np.log(x[1]/x[0])**2)-1/(np.log(x[1]/x[0])), \
                          (x[1]*np.log(x[1]/x[0])-x[1]+x[0])/(x[1]*np.log(x[1]/x[0])**2)], \
                         [-rn*outPlaneThickness*(x[0]-rn), rn*outPlaneThickness*(x[1]-rn)]])
    err = 1

    c = np.log(1.2)/0.2
    a0 = c*rn
    b0 =rn+np.sqrt(rn**2-4*0.5*(c-c**2/2))

    x0 = np.array([a0,b0])
    x = x0
    if np.isinf(rn):
        x = np.array([float('inf'), genStiff/outPlaneThickness])
    else:
        while err > 10**-6:
            xprev = x
            x = x - np.transpose(lin.inv(jac(x, rn, genStiff)).dot(func(x, rn, genStiff)))
            err = lin.norm(x-xprev)
            # print(err)
    if(printBool):
        print(x0)
        print(x)
    return x

coeffs = get_coeffs()
smesh = np.linspace(0,1,globalLen)

plot0 = np.empty(globalLen)
plot1 = np.empty(globalLen)
plot3 = np.empty(globalLen)
for i in range(len(smesh)):
    plot0[i] = alpha_s(i*globalStep,coeffs)
    plot1[i] = genStiffness_s(i*globalStep,coeffs)
    plot3[i] = r_n(i*globalStep, coeffs)
print("alpha:",plot0)
print("genStiffness:",plot1)
print("rn:",plot3)
plt.figure(0)
plt.plot(plot0, plot1)
plt.figure(1)
plt.plot(smesh, plot3)

plot = np.empty([globalLen, 2])
plotrn = np.empty(globalLen)

for i in range(globalLen):
    plot[i,:] = b_a_rootfinding(i*globalStep,coeffs,False)
    plotrn[i] = r_n(i*globalStep,coeffs)

plt.figure(2)
plt.plot(np.linspace(0,1,globalLen), plot[:,0]) # a
plt.plot(np.linspace(0,1,globalLen), plot[:,1]) # b
plt.plot(np.linspace(0,1,globalLen), plotrn)
plt.show()

# def rc_diffeq(dydx, s, coeffs):
#     pass
# sSpan = [0,1]
# dydx = 
# dy_dx = ODE45(rc_diffeq, sSpan, dydx, args=(F,), dense_output=True, method='LSODA')

def deform_ODE(s, gamma, F):
    Fx = F*np.sin(alpha_s(sMax, coeffs))
    Fy = F*np.cos(alpha_s(sMax, coeffs))
    dgSds = d_genStiffness_d_s(s, coeffs)
    LHS = np.empty(2)
    LHS[0] = gamma[1]
    LHS[1] = (Fx*np.sin(alpha_s(s, coeffs)+gamma[0])-Fy*np.cos(alpha_s(s, coeffs)+gamma[0])-E*dgSds*gamma[1])/(genStiffness_s(s, coeffs)*E)
    return LHS

# sSpan = [0,1]
# sMax = sSpan[-1]
# gamma = np.zeros(2)
# F = 1000
# dfrmation = ODE45(deform_ODE, sSpan, gamma, args=(F,), dense_output=True, method='LSODA')

