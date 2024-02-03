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

fullArcLength = 6
globalRes = 100
globalLen = globalRes+1
globalStep = fullArcLength/globalRes
globalMaxIndex = globalLen-1

alphas = [0,100*deg2rad,190*deg2rad,135*deg2rad]
hs     = [1,0.25,1]
xy0 = [1,0]


si = fullArcLength/3.0
sj = 2.0*fullArcLength/3.0
sk = fullArcLength/2.0
sf = fullArcLength

outPlaneThickness = 0.375

def h_poly():
    Mat = np.array([[0,0,0,0,1], \
                    [sk**4,sk**3,sk**2,sk,1], \
                    [sf**4,sf**3,sf**2,sf,1], \
                    [0,0,0,1,0], \
                    [4*sf**3,3*sf**2,2*sf**1,1,0], \
                    ])
    Targ = np.array([[hs[0]], \
                     [hs[1]], \
                     [hs[2]], \
                     [0], \
                     [0], \
                     ])
    hCoeffs = lin.solve(Mat, Targ)
    print(hCoeffs)
    return hCoeffs

def h_s(s, hCoeffs):
    U = np.array([s**4, s**3, s**2, s, 1])
    h = U.dot(hCoeffs)
    return h[0]

def alpha_poly():
    Mat = np.array([[0,0,0,0,0,1], \
                    [si**5,si**4,si**3,si**2,si,1], \
                    [sj**5,sj**4,sj**3,sj**2,sj,1], \
                    [sf**5,sf**4,sf**3,sf**2,sf,1], \
                    [0,0,0,0,1,0], \
                    [5*sf**4,4*sf**3,3*sf**2,2*sf,1,0]])
    Targ = np.array([[alphas[0]],[alphas[1]],[alphas[2]],[alphas[3]],[0],[0]])
    alphaCoeffs = lin.solve(Mat,Targ)
    print(alphaCoeffs)
    return alphaCoeffs

def alpha(s, alphaCoeffs):
    U = np.array([s**5, s**4, s**3, s**2, s, 1])
    alpha = U.dot(alphaCoeffs)
    return alpha[0]

def d_alpha_d_s(s, alphaCoeffs):
    U = np.array([5*s**4, 4*s**3, 3*s**2, 2*s, 1, 0])
    dAlphadS = U.dot(alphaCoeffs)
    return dAlphadS[0]

def r_n(s, alphaCoeffs):
    rn = abs(1/d_alpha_d_s(s, alphaCoeffs))
    if s == fullArcLength:
        rn = float('inf')
    return rn

def rc_rootfinding(s, alphaCoeffs, hCoeffs, printBool):
    rn =  r_n(s, alphaCoeffs)
    h = h_s(s, hCoeffs)

    def func(x, rn, h):
        f1 = h/((np.log((x+0.5*h)/(x-0.5*h))))-rn
        return f1
    def jac(x, rn, h):
        jac = 4*h**2/((4*x**2-h**2)*np.log((x+0.5*h)/(x-0.5*h))**2)
        return jac 
    err = 1
    x0 = h*1.5
    x = x0
    if np.isinf(rn):
        x = float('inf')
    else:
        ii = 0
        while err > 10**-6 or ii < 3000:
            xprev = x
            x = x - (func(x, rn, h))/jac(x, rn, h)
            if x < h/2:
                x = h/2+0.00001
            err = abs(x-xprev)
            ii+=1
    if(printBool):
        print(x0)
        print(x)
        print(rn)
        print(err)
    return x

alphaCoeffs = alpha_poly()
hCoeffs     = h_poly()
smesh = np.linspace(0,fullArcLength,globalLen)
print(smesh)
print(alpha(0, alphaCoeffs))

plot0 = np.empty(globalLen)
plot1 = np.empty(globalLen)
plot2 = np.empty(globalLen)
plot3 = np.empty(globalLen)
plot4 = np.empty(globalLen)
for i in range(len(smesh)):
    plot0[i] = alpha(i*globalStep,alphaCoeffs)
    plot1[i] = d_alpha_d_s(i*globalStep,alphaCoeffs)
    plot2[i] = h_s(i*globalStep,hCoeffs)
    plot3[i] = r_n(i*globalStep,alphaCoeffs)
    print("starting to rootfind", i)
    plot4[i] = rc_rootfinding(i*globalStep,alphaCoeffs,hCoeffs,False)
    print("finished rootfinding")
print("alpha:",plot0)
print("dalpha/ds:",plot1)
print("rn:",plot3)
plt.figure(0)
plt.plot(smesh, plot0)
plt.plot(smesh, plot1)
plt.plot(smesh, plot2)
plt.plot(smesh, plot3)
plt.plot(smesh, plot4)
plt.show()

dx = np.empty(globalLen)
dy = np.empty(globalLen)
x  = np.empty(globalLen)
y  = np.empty(globalLen)
for i in range(len(smesh)):
    if i==0:
        dx[i]=0
        dy[i]=0
        x[i] = xy0[0]
        y[i] = xy0[1]
    else:
        dx[i]=globalStep*np.cos(alpha(i*globalStep,alphaCoeffs))
        dy[i]=globalStep*np.sin(alpha(i*globalStep,alphaCoeffs))
        x[i] = x[i-1]+dx[i]
        y[i] = y[i-1]+dy[i]
plt.figure(1)
plt.plot(x,y)
plt.show()