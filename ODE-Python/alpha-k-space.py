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
ks     = [.05,.01,.05]
xy0 = [1,0]


si = fullArcLength/3.0
sj = 2.0*fullArcLength/3.0
sk = fullArcLength/2.0
sf = fullArcLength

outPlaneThickness = 0.375

def k_poly():
    Mat = np.array([[0,0,0,0,1], \
                    [sk**4,sk**3,sk**2,sk,1], \
                    [sf**4,sf**3,sf**2,sf,1], \
                    [0,0,0,1,0], \
                    [4*sf**3,3*sf**2,2*sf**1,1,0], \
                    ])
    Targ = np.array([[ks[0]], \
                     [ks[1]], \
                     [ks[2]], \
                     [0], \
                     [0], \
                     ])
    hCoeffs = lin.solve(Mat, Targ)
    print(hCoeffs)
    return hCoeffs

def k_s(s, kCoeffs):
    U = np.array([s**4, s**3, s**2, s, 1])
    k = U.dot(kCoeffs)
    return k[0]

def d_k_d_s(s, kCoeffs):
    U = np.array([4*s**3, 3*s**2, 2*s, 1, 0])
    dKdS = U.dot(kCoeffs)
    return dKdS[0]

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

def a_b_rootfinding(s, alphaCoeffs, kCoeffs, printBool):
    rn =  r_n(s, alphaCoeffs)
    k = k_s(s, kCoeffs)

    def func(x, rn, k):

        # x0 = a, x1 = b

        f1 = (x[1]-x[0])/(np.log(x[1]/x[0]))-rn
        f2 = outPlaneThickness*(x[1]-x[0])*((x[1]+x[0])/2-rn)*rn-k
        
        return np.array([f1, f2])
    def jac(x, rn, k):

        return np.array([[(x[1]-x[0])/(x[0]*np.log(x[1]/x[0])**2)-1/(np.log(x[1]/x[0])), \
                          (x[1]*np.log(x[1]/x[0])-x[1]+x[0])/(x[1]*np.log(x[1]/x[0])**2)], \
                         [-rn*outPlaneThickness*(x[0]-rn), rn*outPlaneThickness*(x[1]-rn)]])
    
    err = 1

    factor = 1.2

    c = (np.log(factor)-np.log(1))/(factor-1)

    a0 = c*rn
    b0 = rn + np.sqrt(rn**2+4*0.5*((c-c**2/2)*rn-k/(outPlaneThickness*rn)))

    x0 = [a0, b0]
    x = x0
    if np.isinf(rn):
        x = float('inf')
    else:
        ii = 0
        while err > 10**-6 or ii < 3000:
            xprev = x
            x = x - np.transpose(lin.inv(jac(x, rn, k)).dot(func(x, rn, k)))
            err = lin.norm(x-xprev)
            ii+=1
    if(printBool):
        print(x0)
        print(x)
        print(rn)
        print(err)
    return x


def create_profiles(alphaCoeffs, kCoeffs, smesh):
    plot0 = np.empty(globalLen)
    plot1 = np.empty(globalLen)
    plot2 = np.empty(globalLen)
    plot3 = np.empty(globalLen)
    plot4 = np.empty((globalLen,2))
    for i in range(len(smesh)):
        plot0[i] = alpha(i*globalStep,alphaCoeffs)
        plot1[i] = d_alpha_d_s(i*globalStep,alphaCoeffs)
        plot2[i] = k_s(i*globalStep,kCoeffs)
        plot3[i] = r_n(i*globalStep,alphaCoeffs)
        # print("starting to rootfind", i)
        plot4[i] = a_b_rootfinding(i*globalStep,alphaCoeffs,kCoeffs,False)
        # print("finished rootfinding")
    print("alpha:",plot0)
    print("dalpha/ds:",plot1)
    print("rn:",plot3)
    print("ab:",plot4)
    plt.figure(0)
    plt.plot(smesh, plot0)
    plt.plot(smesh, plot1)
    plt.plot(smesh, plot2)
    plt.plot(smesh, plot3)
    plt.plot(smesh, plot4)
    plt.show()
    return plot0, plot1, plot2, plot3, plot4

def create_geometry(abProfile, alphaProfile, dAlphaProfile, rnProfile, alphaCoeffs, smesh):
    hProfile  = abProfile[:,1]-abProfile[:,0]
    rcProfile = (abProfile[:,1]+abProfile[:,0])/2
    eProfile  = rcProfile - rnProfile
    for i in range(globalLen):
        if np.isnan(eProfile[i]):
            eProfile[i] = 0
    print(eProfile)
    rnDx = np.empty(globalLen)
    rnDy = np.empty(globalLen)
    rnX  = np.empty(globalLen)
    rnY  = np.empty(globalLen)
    rcX  = np.empty(globalLen)
    rcY  = np.empty(globalLen)
    for i in range(len(smesh)):
        if i==0:
            rnDx[i]=0
            rnDy[i]=0
            rnX[i] = xy0[0]
            rnY[i] = xy0[1]
            rcX[i] = xy0[0]+eProfile[i]*np.sin(alpha(i*globalStep,alphaCoeffs))*d_alpha_d_s(i*globalStep,alphaCoeffs)
            rcY[i] = xy0[1]+eProfile[i]*-np.cos(alpha(i*globalStep,alphaCoeffs))*d_alpha_d_s(i*globalStep,alphaCoeffs)
        else:
            rnDx[i]=globalStep*np.cos(alpha(i*globalStep,alphaCoeffs))
            rnDy[i]=globalStep*np.sin(alpha(i*globalStep,alphaCoeffs))
            rnX[i] = rnX[i-1]+rnDx[i]
            rnY[i] = rnY[i-1]+rnDy[i]
            rcX[i] = rnX[i]+eProfile[i]*np.sin(alpha(i*globalStep,alphaCoeffs))*d_alpha_d_s(i*globalStep,alphaCoeffs)
            rcY[i] = rnY[i]+eProfile[i]*-np.cos(alpha(i*globalStep,alphaCoeffs))*d_alpha_d_s(i*globalStep,alphaCoeffs)
    plt.figure(1)
    plt.plot(rnX,rnY)
    plt.plot(rcX,rcY)
    plt.show()
    return rnX, rnY


def create_initial_spring(geoBool):
    
    alphaCoeffs = alpha_poly()
    kCoeffs     = k_poly()
    smesh = np.linspace(0,fullArcLength,globalLen)
    print(smesh)
    print(alpha(0, alphaCoeffs))

    if geoBool:
           
        [alphaProfile, dAlphaProfile, kProfile, rnProfile, abProfile] = create_profiles(alphaCoeffs, kCoeffs, smesh)
        [rnX, rnY] = create_geometry(abProfile, alphaProfile, dAlphaProfile, rnProfile, alphaCoeffs, smesh)
    
    return rnX, rnY

def deform_ODE(s, gamma, F, ang):
    
    Fx = F*np.sin(ang)
    Fy = F*-np.cos(ang)
    dgSds = d_k_d_s(s, kCoeffs)
    LHS = np.empty(2)
    LHS[0] = gamma[1]
    LHS[1] = (Fx*np.sin(alpha(s, alphaCoeffs)+gamma[0])-Fy*np.cos(alpha(s, alphaCoeffs)+gamma[0])-E*dgSds*gamma[1])/(k_s(s, kCoeffs)*E)
    return LHS




def integral_s(fun, array, mesh):
    #centered Riemman sum
    i = 1
    sum = 0
    while i < len(array):
        step = mesh[i]-mesh[i-1]
        sum += 0.5*(fun(mesh[i], array[i])+fun(mesh[i-1], array[i-1]))*step
        i += 1
    return sum

def integrand_x(s, gamma):
    return np.cos(alpha(s, alphaCoeffs)+gamma)
def integrand_y(s, gamma):
    return np.sin(alpha(s, alphaCoeffs)+gamma)


[rnX, rnY] = create_initial_spring(True)

sMax = fullArcLength
alphaCoeffs = alpha_poly()
kCoeffs = k_poly()
F = 600
dGdS0 = 0

err = 1

g0 = np.array([0, dGdS0])
ang = alpha(sMax, alphaCoeffs)
res = ODE45(deform_ODE, [0, sMax], g0, args=(F,ang), method='LSODA')
gamma = res.y[0,:]
xdis = integral_s(integrand_x, gamma, res.t)
ydis = integral_s(integrand_y, gamma, res.t)
diff0 = np.sqrt(xdis**2+ydis**2)-np.sqrt(rnX**2+rnY**2)[-1]
ang = np.arctan2((rnY[-1]-ydis),(rnX[-1]-xdis))
print(diff0, ang)

dGdS0 = 0.01
dGdS0prev = 0
gain = 1
iii = 0
while abs(err)>1.1e-3:
    g0 = np.array([0, dGdS0])
    res = ODE45(deform_ODE, [0, sMax], g0, args=(F,ang), method='LSODA')
    gamma = res.y[0,:]
    xdis = integral_s(integrand_x, gamma, res.t)
    ydis = integral_s(integrand_y, gamma, res.t)
    diff1 = np.sqrt(xdis**2+ydis**2)-np.sqrt(rnX**2+rnY**2)[-1]
    ang = np.arctan2((rnY[-1]-ydis),(rnX[-1]-xdis))
    grad = (diff1-diff0)/(dGdS0-dGdS0prev)
    err  = abs(diff1)
    dGdS0prev = dGdS0
    dGdS0 = dGdS0-gain*diff1/grad
    diff0 = diff1
    iii +=1
    
    print(diff0, dGdS0, ang, iii)

Fx = F*np.sin(alpha(sMax, alphaCoeffs))
Fy = F*-np.cos(alpha(sMax, alphaCoeffs))
dBeta = np.arccos(integral_s(integrand_x, gamma, res.t)/np.sqrt(rnX**2+rnY**2)[-1])-(alphas[-1]-alphas[0])
Torque = 2*(xdis*Fy-ydis*Fx+E*k_s(sMax, kCoeffs)*res.y[1,-1])
springK = Torque/dBeta

print("dBeta:",dBeta,"Torque:",Torque,"spring k:",springK)


