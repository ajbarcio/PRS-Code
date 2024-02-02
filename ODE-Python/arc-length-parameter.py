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

startingAngle = 0*deg2rad
intersectOffset = 30
endingAngle = (180-intersectOffset)*deg2rad

stiffFactor = 0.01

p1Angle = (startingAngle+(endingAngle-startingAngle)/3)
p2Angle = (startingAngle+2*(endingAngle-startingAngle)/3)

p0 = np.array([startingAngle, stiffFactor])
p1 = np.array([p1Angle, 0])
p2 = np.array([p2Angle, 0])
p3 = np.array([endingAngle, stiffFactor])

outPlaneThickness = 0.375

def get_coeffs():
    D = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
    P = np.array([p0,p1,p2,p3])
    coeffs = D.dot(P)
    return coeffs

def get_metaphor_coords(s, coeffs):
    U = np.array([s**3, s**2, s, 1])
    metaphorCoords = U.dot(coeffs)
    return metaphorCoords

def alpha_s(s, coeffs):
    return get_metaphor_coords(s, coeffs)[0]

def genStiffness_s(s, coeffs):
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

def h_rc_rootfinding(s, coeffs, printBool):
    rn =  r_n(s, coeffs)
    genStiff = get_metaphor_coords(s, coeffs)[1]
    def func(x, rn, genStiff):
        f1 = x[1]/((np.log((x[0]+0.5*x[1])/(x[0]-0.5*x[1]))))-rn
        # f1 = f1**2
        f2 = outPlaneThickness*x[1]*(x[0]-rn)*rn-genStiff
        # f2 = f2**2
        # return lin.norm(np.array([f1, f2]))
        return np.array([f1, f2])
    def jac(x, rn, genStiff):
        num = np.log((x[0]+0.5*x[1])/(x[0]-0.5*x[1]))*x[1]**2+4*x[1]*x[0]-4*x[0]**2*np.log((x[0]+0.5*x[1])/(x[0]-0.5*x[1]))
        den = np.log((x[0]+0.5*x[1])/(x[0]-0.5*x[1]))**2*(x[1]-2*x[0])*(x[1]+2*x[0])
        return np.array([[4*x[1]**2/((4*x[0]**2-x[1]**2)*np.log((x[0]+0.5*x[1])/(x[0]-0.5*x[1]))**2), \
                          num/den], \
                         [outPlaneThickness*x[1]*rn, outPlaneThickness*(x[0]-rn)*rn]])
    err = 1
    x0 = np.array([1,1])
    x=x0
    while err > 10**-6:
        xprev = x
        x = x - np.transpose(lin.inv(jac(x, rn, genStiff)).dot(func(x, rn, genStiff)))
        err = lin.norm(x-xprev)
        # print(err)
    # res = op.fsolve(func, guess, args = (rn,genStiff))
    if(printBool):
        print(x0)
        print(x)
    return x

coeffs = get_coeffs()
plot = np.empty([globalLen, 2])
plotrn = np.empty(globalLen)

for i in range(globalLen):
    plot[i,:] = h_rc_rootfinding(i*globalStep,coeffs,False)
    plotrn[i] = r_n(i*globalStep,coeffs)
plt.plot(np.linspace(0,1,globalLen), plot[:,0])
plt.plot(np.linspace(0,1,globalLen), plot[:,1])
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

sSpan = [0,1]
sMax = sSpan[-1]
gamma = np.zeros(2)
F = 1000
dfrmation = ODE45(deform_ODE, sSpan, gamma, args=(F,), dense_output=True, method='LSODA')

