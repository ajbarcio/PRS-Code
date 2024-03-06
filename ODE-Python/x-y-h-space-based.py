from mpl_toolkits import mplot3d
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from materials import *
from openpyxl import *
from scipy.integrate import solve_ivp as ODE45
from scipy import optimize as op
import os
import random
import math

finiteDifferenceLength = 0.0000001

fullArcLength = 6.2
globalRes = 500 # 2284 for desired acuracy
globalLen = globalRes+1
globalStep = fullArcLength/globalRes
globalMaxIndex = globalLen-1

deg2rad = np.pi/180

E = 10000000
outPlaneThickness = .375

R0 = 2.2
R1 = (2.2+5.9)/2
R2 = (2.2+5.9)/2
R3 = 5.9

betaB = 45*deg2rad
betaC = 90*deg2rad
betaD = 135*deg2rad
beta0 = betaD

x0 = R0
y0 = 0

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
hs  = np.array([.75, .25, .6])

ctrlHs      = np.array([0,fullArcLength*.5,fullArcLength])
ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

def h_poly(hs, ctrlHs):
    Mat = np.array([[0,0,0,0,1], \
                    [ctrlHs[1]**4,ctrlHs[1]**3,ctrlHs[1]**2,ctrlHs[1],1], \
                    [ctrlHs[2]**4,ctrlHs[2]**3,ctrlHs[2]**2,ctrlHs[2],1], \
                    [0,0,0,1,0], \
                    [4*ctrlHs[2]**3,3*ctrlHs[2]**2,2*ctrlHs[2]**1,1,0], \
                    ])
    Targ = np.array([[hs[0]], \
                     [hs[1]], \
                     [hs[2]], \
                     [0], \
                     [0], \
                     ])
    hCoeffs = lin.solve(Mat, Targ)
    return hCoeffs

def xy_poly(pts, ctrlLengths):
    xTarg = np.array([[pts[0,0]], \
                      [pts[1,0]], \
                      [pts[2,0]], \
                      [pts[3,0]], \
                      [1], \
                      [pts[3,0]]]) # /lin.norm(pts[3])
    
    yTarg = np.array([[pts[0,1]], \
                      [pts[1,1]], \
                      [pts[2,1]], \
                      [pts[3,1]], \
                      [0], \
                      [pts[3,1]]])
    
    Mat = np.array([[0,0,0,0,0,1], \
                    [ctrlLengths[1]**5,ctrlLengths[1]**4,ctrlLengths[1]**3,ctrlLengths[1]**2,ctrlLengths[1],1], \
                    [ctrlLengths[2]**5,ctrlLengths[2]**4,ctrlLengths[2]**3,ctrlLengths[2]**2,ctrlLengths[2],1], \
                    [ctrlLengths[3]**5,ctrlLengths[3]**4,ctrlLengths[3]**3,ctrlLengths[3]**2,ctrlLengths[3],1], \
                    [0,0,0,0,1,0], \
                    [5*ctrlLengths[3]**4,4*ctrlLengths[3]**3,3*ctrlLengths[3]**2,2*ctrlLengths[3],1,0]])
    
    xCoeffs = lin.solve(Mat,xTarg)
    yCoeffs = lin.solve(Mat,yTarg)

    return xCoeffs, yCoeffs

def coord(s, coeffs):
    
    if hasattr(s, "__len__"):
        coord = np.empty(len(s))
        i = 0
        for value in s:
            U = np.array([value**5, value**4, value**3, value**2, value, 1])
            coord[i] = U.dot(coeffs)[0]
            i+=1
        return coord
    else:
        U = np.array([s**5, s**4, s**3, s**2, s, 1])
        coord = U.dot(coeffs)
        return coord[0]
    
def d_coord_d_s(s, coeffs):
    if hasattr(s, "__len__"):
        dCds = np.empty(len(s))
        i=0
        for value in s:
            U = np.array([5*value**4, 4*value**3, 3*value**2, 2*value, 1, 0])
            dCds[i] = U.dot(coeffs)[0]
            i+=1
        return dCds
    else:
        U = np.array([5*s**4, 4*s**3, 3*s**2, 2*s, 1, 0])
        dCds = U.dot(coeffs)
        return dCds[0]

def alpha_polyfit(smesh,yCoeffs, xCoeffs ):

    alphaList = np.arctan2(d_coord_d_s(smesh, yCoeffs),d_coord_d_s(smesh, xCoeffs))
    for i in range(len(alphaList)):
        if alphaList[i]<0:
            alphaList[i]=alphaList[i]+2*math.pi

    alphaPoly = np.polyfit(smesh,alphaList,6)
    return alphaPoly


def alpha_xy(s, xCoeffs, yCoeffs):
    smesh = np.linspace(0, fullArcLength, globalLen)
    apoly = alpha_polyfit(smesh, yCoeffs, xCoeffs)
    alpha = apoly[0]*s**6+apoly[1]*s**5+apoly[2]*s**4+apoly[3]*s**3+apoly[4]*s**2+apoly[5]*s**1+apoly[6]*s**0
    
    return alpha

def r_n(s, xCoeffs, yCoeffs):
    smesh = np.linspace(0, fullArcLength, globalLen)
    apoly = alpha_polyfit(smesh, yCoeffs, xCoeffs)
    dAlphadS = 6*apoly[0]*s**5+5*apoly[1]*s**5+4*apoly[2]*s**3+3*apoly[3]*s**2+2*apoly[4]*s**1+1*apoly[5]*s**0+0*apoly[6]*s**0
    rn = abs(1/dAlphadS)
    return rn

def h_s(s, hCoeffs):
    if hasattr(s, "__len__"):
        h = np.empty(len(s))
        i = 0
        for value in s:
            U = np.array([value**4, value**3, value**2, value, 1])
            h[i] = U.dot(hCoeffs)[0]
            i+=1
        return h
    else:
        U = np.array([s**4, s**3, s**2, s, 1])
        h = U.dot(hCoeffs)
        return h[0]

def rc_rootfinding(s, xCoeffs, yCoeffs, hCoeffs, printBool):
    rn =  r_n(s, xCoeffs, yCoeffs)
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
        while err > 10**-6 and ii < 1500:
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

def cI_s(s, xCoeffs, yCoeffs, hCoeffs):
    cI = h_s(s, hCoeffs)*outPlaneThickness* \
         (rc_rootfinding(s, xCoeffs, yCoeffs, hCoeffs, False)-r_n(s, xCoeffs, yCoeffs)) \
         *r_n(s, xCoeffs, yCoeffs)
    return cI

def d_cI_d_s(s, xCoeffs, yCoeffs, hCoeffs):
    dcIds = (cI_s(s+finiteDifferenceLength, xCoeffs, yCoeffs, hCoeffs)-cI_s(s-finiteDifferenceLength, xCoeffs, yCoeffs, hCoeffs))/(2*finiteDifferenceLength)
    return dcIds

def form_spring(pts, hs, ctrlLengths, ctrlHs):
    xCoeffs = xy_poly(pts, ctrlLengths)[0]
    yCoeffs = xy_poly(pts, ctrlLengths)[1]
    hCoeffs = h_poly(hs, ctrlHs)
    return xCoeffs, yCoeffs, hCoeffs

def deform_ODE_Barcio(s, p, *args): #Fx, Fy, geometryDef, dgds0

    Fx          = args[0][0][0][0]
    Fy          = args[0][0][0][1]
    geometryDef = args[0][0][0][2]
    dgds0       = args[0][0][0][3]

    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    hCoeffs = geometryDef[2]

    Tdim = E*cI_s(s, xCoeffs, yCoeffs, hCoeffs)
    yorg = coord(s, yCoeffs)
    xorg = coord(s, xCoeffs)
    dxds = d_coord_d_s(s, xCoeffs)
    dyds = d_coord_d_s(s, yCoeffs)

    alpha = alpha_xy(s, xCoeffs, yCoeffs)

    LHS = np.empty(2)
    LHS[0] = p[1]
    LHS[1] = (Fx*np.sin(alpha+p[0])-Fy*np.cos(alpha+p[0])-E*d_cI_d_s(s, xCoeffs, yCoeffs, hCoeffs)*p[1])/(cI_s(s, xCoeffs, yCoeffs, hCoeffs)*E)
    return LHS

def deform_ODE_Barcio_to_Thomas(s, p, *args):
    
    
    pass

def deform_ODE_Thomas(s, p, *args): #Fx, Fy, geometryDef, dgds0

    # deal with args in the shittiest way possible

    Fx          = args[0][0][0][0]
    Fy          = args[0][0][0][1]
    geometryDef = args[0][0][0][2]
    dgds0       = args[0][0][0][3]

    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    hCoeffs = geometryDef[2]

    # Tdim = EI (dimensionalize bending moment)
    # dxds, dyds are ANALYTICAL derivatives based on x/y polynomials

    Mdim = E*cI_s(s, xCoeffs, yCoeffs, hCoeffs)

    xorg = coord(s, xCoeffs)
    yorg = coord(s, yCoeffs)

    dxds = d_coord_d_s(s, xCoeffs)
    dyds = d_coord_d_s(s, yCoeffs)

    dydx = dyds/dxds

    alpha = alpha_xy(s, xCoeffs, yCoeffs)
    # constants
    LHS = np.empty(3)

    LHS[0] = dgds0 # - Fx/Tdim*yorg + Fy/Tdim*xorg
    # LHS[1] = 0 #-dxds
    # LHS[2] = 0 #-dyds

    # print(LHS[0], LHS[1], LHS[2])
    
    # rateMag1 = lin.norm([p[1], p[2]])
    rateMag = lin.norm([dyds, dxds])
    # rateMag =  np.average([rateMag1, rateMag2])

    LHS[0] = LHS[0] + Fy/Mdim*p[1] - Fx/Mdim*p[2]
    LHS[1] = np.cos(np.arctan2(dyds,dxds)+p[0]) # -(g*y +/- (g^2 - y^2 + 1)^(1/2))/(g^2 + 1)
    LHS[2] = np.sin(np.arctan2(dyds,dxds)+p[0])
    print(LHS[1], LHS[2])
    print("=>")
    LHS[1] = LHS[1]*dxds/np.cos(np.arctan2(dyds,dxds))
    LHS[2] = LHS[2]*dyds/np.sin(np.arctan2(dyds,dxds))
    print(LHS[1], LHS[2])
    print(dxds, dyds)
    print("--------------")

    # print(LHS[0], LHS[1], LHS[2])
    # print(dxds, np.cos(alpha+p[0]), dxds/np.cos(alpha))
    # print(dxds/np.cos(alpha))
    # print("--------------------")

    return LHS

def fixed_rk4(fun, y0, xmesh, *args): # (fun, alphaCoeffs, cICoeffs, y0, Fx, Fy, xmesh)
    step = xmesh[1]-xmesh[0]
    if hasattr(y0, '__len__'):
        res = np.empty((len(y0), len(xmesh)))
        for i in range(len(xmesh)):
            # print(i)
            if i == 0:
                for j in range(len(y0)):
                    res[j,i] = y0[j]
            else:
                stepRes = rk4_step(fun, xmesh[i], y0, step, args)
                # print(stepRes)
                for j in range(len(y0)):
                    res[j,i] = stepRes[j]
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
    return res

def rk4_step(fun, x0, y0, dx, *args): # (fun, alphaCoeffs, cICoeffs, x0, y0, du, Fx, Fy)
    #
    #  Get four sample values of the derivative.
    #
    f1 = fun ( x0,            y0, args)
    f2 = fun ( x0 + dx / 2.0, y0 + dx * f1 / 2.0, args)
    f3 = fun ( x0 + dx / 2.0, y0 + dx * f2 / 2.0, args)
    f4 = fun ( x0 + dx,       y0 + dx * f3, args)
    #
    #  Combine them to estimate the solution gamma at time T1 = T0 + DT.
    #
    y1 = y0 + dx * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0

    return y1

geometryDef = form_spring(pts, hs, ctrlLengths, ctrlHs)
smesh = np.linspace(0, fullArcLength, globalLen)
res = fixed_rk4(deform_ODE_Thomas, np.array([15*deg2rad,x0,y0]), smesh, (0, 0, geometryDef, 0))
# res2 = fixed_rk4(deform_ODE_Barcio, np.array([0,0]), smesh, (0, 0, geometryDef, 0))

xorg = coord(smesh, geometryDef[0])
yorg = coord(smesh, geometryDef[1])

# print(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

plt.figure("thomas method")
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

# plt.figure("barcio method")
# plt.plot(xorg, yorg)
# plt.plot(xorg)

plt.figure(0)
plt.plot(np.transpose(res))
plt.plot(xorg)
plt.plot(yorg)
# for i in range(len(smesh)):
#     print(xorg[i], yorg[i], res[1,i], res[2,i])

plt.show()