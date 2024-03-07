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

def cI_poly(cIs, ctrlcIs):
    Mat = np.array([[0,0,0,0,1], \
                    [ctrlcIs[1]**4,ctrlcIs[1]**3,ctrlcIs[1]**2,ctrlcIs[1],1], \
                    [ctrlcIs[2]**4,ctrlcIs[2]**3,ctrlcIs[2]**2,ctrlcIs[2],1], \
                    [0,0,0,1,0], \
                    [4*ctrlcIs[2]**3,3*ctrlcIs[2]**2,2*ctrlcIs[2]**1,1,0], \
                    ])
    Targ = np.array([[cIs[0]], \
                     [cIs[1]], \
                     [cIs[2]], \
                     [0], \
                     [0], \
                     ])
    cICoeffs = lin.solve(Mat, Targ)
    return cICoeffs

def xy_poly(pts, ctrlLengths):
    xTarg = np.array([[pts[0,0]], \
                      [pts[1,0]], \
                      [pts[2,0]], \
                      [pts[3,0]], \
                      [1], \
                      [pts[3,0]], \
                      [0], \
                      [0]]) # /lin.norm(pts[3])
    
    yTarg = np.array([[pts[0,1]], \
                      [pts[1,1]], \
                      [pts[2,1]], \
                      [pts[3,1]], \
                      [0], \
                      [pts[3,1]], \
                      [0], \
                      [0]])
    
    Mat = np.array([[0,0,0,0,0,0,0,1], \
                    [ctrlLengths[1]**7,ctrlLengths[1]**6,ctrlLengths[1]**5,ctrlLengths[1]**4,ctrlLengths[1]**3,ctrlLengths[1]**2,ctrlLengths[1],1], \
                    [ctrlLengths[2]**7,ctrlLengths[2]**6,ctrlLengths[2]**5,ctrlLengths[2]**4,ctrlLengths[2]**3,ctrlLengths[2]**2,ctrlLengths[2],1], \
                    [ctrlLengths[3]**7,ctrlLengths[3]**6,ctrlLengths[3]**5,ctrlLengths[3]**4,ctrlLengths[3]**3,ctrlLengths[3]**2,ctrlLengths[3],1], \
                    [0,0,0,0,0,0,1,0], \
                    [7*ctrlLengths[3]**6,6*ctrlLengths[3]**5,5*ctrlLengths[3]**4,4*ctrlLengths[3]**3,3*ctrlLengths[3]**2,2*ctrlLengths[3],1,0],
                    [0,0,0,0,0,2,0,0],
                    [42*ctrlLengths[3]**5,30*ctrlLengths[3]*4,20*ctrlLengths[3]**3,12*ctrlLengths[3]**2,6*ctrlLengths[3]**1,2,0,0]])
    
    xCoeffs = lin.solve(Mat,xTarg)
    yCoeffs = lin.solve(Mat,yTarg)

    return xCoeffs, yCoeffs

def coord(s, coeffs):
    
    if hasattr(s, "__len__"):
        coord = np.empty(len(s))
        i = 0
        for value in s:
            U = np.array([value**7, value**6, value**5, value**4, value**3, value**2, value, 1])
            coord[i] = U.dot(coeffs)[0]
            i+=1
        return coord
    else:
        U = np.array([s**7, s**6, s**5, s**4, s**3, s**2, s, 1])
        coord = U.dot(coeffs)
        return coord[0]
    
def d_coord_d_s(s, coeffs):
    if hasattr(s, "__len__"):
        dCds = np.empty(len(s))
        i=0
        for value in s:
            U = np.array([7*value**6, 6*value**5, 5*value**4, 4*value**3, 3*value**2, 2*value, 1, 0])
            dCds[i] = U.dot(coeffs)[0]
            i+=1
        return dCds
    else:
        U = np.array([7*s**6, 6*s**5, 5*s**4, 4*s**3, 3*s**2, 2*s, 1, 0])
        dCds = U.dot(coeffs)
        return dCds[0]

def d2_coord_d_s2(s, coeffs):
    if hasattr(s, "__len__"):
        dCds = np.empty(len(s))
        i=0
        for value in s:
            U = np.array([42*value**5, 30*value**4, 20*value**3, 12*value**2, 6*value**1, 2, 0, 0])
            dCds[i] = U.dot(coeffs)[0]
            i+=1
        return dCds
    else:
        U = np.array([42*s**5, 30*s**4, 20*s**3, 12*s**2, 6*s**1, 2, 0, 0])
        dCds = U.dot(coeffs)
        return dCds[0]

# def alpha_polyfit(smesh,yCoeffs, xCoeffs ):

#     alphaList = np.arctan2(d_coord_d_s(smesh, yCoeffs),d_coord_d_s(smesh, xCoeffs))
#     for i in range(len(alphaList)):
#         if alphaList[i]<0:
#             alphaList[i]=alphaList[i]+2*math.pi

#     alphaPoly = np.polyfit(smesh,alphaList,9)
#     return alphaPoly


def alpha_xy(s, xCoeffs, yCoeffs):
    if hasattr(s, "__len__"):
        alphaList = np.arctan2(d_coord_d_s(smesh, yCoeffs),d_coord_d_s(smesh, xCoeffs))
        for i in range(len(alphaList)):
            if alphaList[i]<0:
                alphaList[i]=alphaList[i]+2*math.pi
        return alphaList
    else:
        alpha = np.arctan2(d_coord_d_s(s, yCoeffs),d_coord_d_s(s, xCoeffs))
        if alpha<0:
            alpha=alpha+2*math.pi
        return alpha

def r_n(s, xCoeffs, yCoeffs):
    d2yds2 = d2_coord_d_s2(s, yCoeffs)
    d2xds2 = d2_coord_d_s2(s, xCoeffs)
    dyds   = d_coord_d_s  (s, yCoeffs)
    dxds   = d_coord_d_s  (s, xCoeffs)
    # dAlphadS = (d2yds2/dxds-d2xds2*dyds/dxds)/(1+dyds**2/dxds)
    if (d2yds2/dxds-d2xds2*dyds/dxds)==0:
        rn = float('inf')
    else:
        rn = abs((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds))
    # if s == fullArcLength:
    #     print("-----------------------------------------")
    #     print(dyds, dxds, d2yds2, d2xds2)
    #     print("-----------------------------------------")
    return rn
    

def cI_s(s, cICoeffs):
    if hasattr(s, "__len__"):
        cI = np.empty(len(s))
        i = 0
        for value in s:
            U = np.array([value**4, value**3, value**2, value, 1])
            cI[i] = U.dot(cICoeffs)[0]
            i+=1
        return cI
    else:
        U = np.array([s**4, s**3, s**2, s, 1])
        cI = U.dot(cICoeffs)
        return cI[0]

# def d_cI_d_s(s, xCoeffs, yCoeffs, hCoeffs):
#     dcIds = (cI_s(s+finiteDifferenceLength, xCoeffs, yCoeffs, hCoeffs)-cI_s(s-finiteDifferenceLength, xCoeffs, yCoeffs, hCoeffs))/(2*finiteDifferenceLength)
#     return dcIds

def d_cI_d_s(s, cICoeffs):
    if hasattr(s, "__len__"):
        dcIcs = np.empty(len(s))
        i = 0
        for value in s:
            U = np.array([4*value**3, 3*value**2, 2*value**1, 1, 0])
            dcIcs[i] = U.dot(cICoeffs)[0]
            i+=1
        return dcIcs
    else:
        U = np.array([4*s**3, 3*s**2, 2*s, 1, 0])
        dcIcs = U.dot(cICoeffs)
        return dcIcs[0]

def rc_rootfinding(s, xCoeffs, yCoeffs, cICoeffs, printBool):
    rn =  r_n(s, xCoeffs, yCoeffs)
    cI = cI_s(s, cICoeffs)
    # h = (cI/(outPlaneThickness*(x-rn)*rn))

    def func(x, rn, h):
        f1 = h/((np.log((x+0.5*h)/(x-0.5*h))))-rn
        return f1
    def jac(x, rn, h):
        jac = 4*h**2/((4*x**2-h**2)*np.log((x+0.5*h)/(x-0.5*h))**2)
        return jac 
    err = 1
    
    # x0 = h*1.5
    x0 = 1000
    x = x0
    ii = 0
    if np.isinf(rn):
        x = float('inf')
    else:
        while err > 10**-6 and ii < 3000:
            xprev = x
            h = (cI/(outPlaneThickness*(x-rn)*rn))
            # if ii==0:
            #     print(x,h)
            x = x - (func(x, rn, h))/jac(x, rn, h)
            if x < h/2:
                x = h/2+0.00001
            err = abs(x-xprev)
            ii+=1
    if(printBool):
        print(x0)
        print(x)
        print(ii)
        print(rn)
        print(err)
    # here x represents rc
    return x

def a_b_rootfinding(s, xCoeffs, yCoeffs, cICoeffs, printBool):
    rn =  r_n(s, xCoeffs, yCoeffs)
    cI = cI_s(s, cICoeffs)
    def func(x, rn, cI):
        # x0 = a, x1 = b
        f1 = (x[1]-x[0])/(np.log(x[1]/x[0]))-rn
        f2 = outPlaneThickness*(x[1]-x[0])*((x[1]+x[0])/2-rn)*rn-cI        
        return np.array([f1, f2])
    def jac(x, rn, cI):
        return np.array([[(x[1]-x[0])/(x[0]*np.log(x[1]/x[0])**2)-1/(np.log(x[1]/x[0])), \
                          (x[1]*np.log(x[1]/x[0])-x[1]+x[0])/(x[1]*np.log(x[1]/x[0])**2)], \
                         [-rn*outPlaneThickness*(x[0]-rn), rn*outPlaneThickness*(x[1]-rn)]])
    err = 1
    # establish expected closesness range for a/b compared to rn to establish initial guess based on linearization about 1
    factor = 1.2
    # approimation: ln(x)=c(x-1), c such that c(x-1) = ln(factor)
    c = (np.log(factor)-np.log(1))/(factor-1)
    a0 = c*rn
    b0 = rn + np.sqrt(rn**2+4*0.5*((c-c**2/2)*rn-cI/(outPlaneThickness*rn)))
    x0 = [a0, b0]
    x = x0
    ii = 0
    if np.isinf(rn):
        x = np.array([float('inf'), float('inf')])
    else:
        while err > 10**-6 or ii < 3000:
            xprev = x
            x = x - np.transpose(lin.inv(jac(x, rn, cI)).dot(func(x, rn, cI)))
            err = lin.norm(x-xprev)
            ii+=1
    if(printBool):
        print(x0)
        print(x)
        print(ii)
        print(rn)
        print(err)
    return x # here x is (a,b)??

# def cI_s(s, xCoeffs, yCoeffs, hCoeffs):
    if r_n(s, xCoeffs, yCoeffs) == float('inf'):
        cI = cI_s(s, hCoeffs)**3*outPlaneThickness/12
    else:
        cI = cI_s(s, hCoeffs)*outPlaneThickness* \
            (rc_rootfinding(s, xCoeffs, yCoeffs, hCoeffs, False)-r_n(s, xCoeffs, yCoeffs)) \
            *r_n(s, xCoeffs, yCoeffs)
    return cI

def form_spring(pts, hs, ctrlLengths, ctrlcIs):
    xCoeffs = xy_poly(pts, ctrlLengths)[0]
    yCoeffs = xy_poly(pts, ctrlLengths)[1]
    cICoeffs = cI_poly(cIs, ctrlcIs)
    return xCoeffs, yCoeffs, cICoeffs

# def deform_ODE_Barcio(s, p, *args): #Fx, Fy, geometryDef, dgds0

    Fx          = args[0][0][0][0]
    Fy          = args[0][0][0][1]
    geometryDef = args[0][0][0][2]
    dgds0       = args[0][0][0][3]

    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    cICoeffs = geometryDef[2]

    Tdim = E*cI_s(s, xCoeffs, yCoeffs, cICoeffs)
    yorg = coord(s, yCoeffs)
    xorg = coord(s, xCoeffs)
    dxds = d_coord_d_s(s, xCoeffs)
    dyds = d_coord_d_s(s, yCoeffs)

    alpha = alpha_xy(s, xCoeffs, yCoeffs)

    LHS = np.empty(2)
    LHS[0] = p[1]
    LHS[1] = (Fx*np.sin(alpha+p[0])-Fy*np.cos(alpha+p[0])-E*d_cI_d_s(s, xCoeffs, yCoeffs, cICoeffs)*p[1])/(cI_s(s, xCoeffs, yCoeffs, cICoeffs)*E)
    return LHS

# def deform_ODE_Barcio_to_Thomas(s, p, *args):
    
    
    pass

def deform_ODE(s, p, *args): #Fx, Fy, geometryDef, dgds0

    # deal with args in the shittiest way possible

    Fx          = args[0][0][0][0]
    Fy          = args[0][0][0][1]
    geometryDef = args[0][0][0][2]
    dgds0       = args[0][0][0][3]

    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    cICoeffs = geometryDef[2]

    # Tdim = EI (dimensionalize bending moment)
    # dxds, dyds are ANALYTICAL derivatives based on x/y polynomials

    Mdim = E*cI_s(s, cICoeffs)

    dxds = d_coord_d_s(s, xCoeffs)
    dyds = d_coord_d_s(s, yCoeffs)
    # constants
    LHS = np.empty(3)

    LHS[0] = dgds0

    LHS[0] = LHS[0] + Fy/Mdim*p[1] - Fx/Mdim*p[2]
    LHS[1] = np.cos(np.arctan2(dyds,dxds)+p[0]) # -(g*y +/- (g^2 - y^2 + 1)^(1/2))/(g^2 + 1)
    LHS[2] = np.sin(np.arctan2(dyds,dxds)+p[0])
    # print(LHS[1], LHS[2])
    # print("=>")
    LHS[1] = LHS[1]*dxds/np.cos(np.arctan2(dyds,dxds))
    LHS[2] = LHS[2]*dyds/np.sin(np.arctan2(dyds,dxds))
    # print(LHS[1], LHS[2])
    # print(dxds, dyds)
    # print("--------------")

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

def deform_spring_by_torque(torqueTarg, geometryDef):
    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    cICoeffs = geometryDef[2]

    xorg = coord(fullArcLength, xCoeffs)
    yorg = coord(fullArcLength, yCoeffs)

    SF = np.array([0, 0, 0])
    err = np.ones(2)
    i = 0
    while lin.norm(err) > 10e-6:
        # establish results of guess
        err, res = int_error_result(SF, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
        # estimate jacobian with finite difference
        J = fd_J(SF, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
        # establish next guess:
        SF = SF - lin.pinv(J).dot(err)
        print("torque deform error:",lin.norm(err))

    return res
        

def fd_J(SF, xCoeffs, yCoeffs, cICoeffs, xorg, yorg):
    # finite difference in Fx
    SFback = SF-np.array([finiteDifferenceForce,0,0])
    SFfrwd = SF+np.array([finiteDifferenceForce,0,0])
    errBack, resG = int_error_result(SFback, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    errFrwd, resG = int_error_result(SFfrwd, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    derrRdFx = (errFrwd[0]-errBack[0])/(2*finiteDifferenceForce)
    derrGdFx = (errFrwd[1]-errBack[1])/(2*finiteDifferenceForce)
    # finite difference in Fy
    SFback = SF-np.array([0,finiteDifferenceForce,0])
    SFfrwd = SF+np.array([0,finiteDifferenceForce,0])
    errBack, resG = int_error_result(SFback, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    errFrwd, resG = int_error_result(SFfrwd, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    derrRdFy = (errFrwd[0]-errBack[0])/(2*finiteDifferenceForce)
    derrGdFy = (errFrwd[1]-errBack[1])/(2*finiteDifferenceForce)
    # finite difference in T
    SFback = SF-np.array([0,0,finiteDifferenceTorque])
    SFfrwd = SF+np.array([0,0,finiteDifferenceTorque])
    errBack, resG = int_error_result(SFback, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    errFrwd, resG = int_error_result(SFfrwd, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    derrRdT = (errFrwd[0]-errBack[0])/(2*finiteDifferenceTorque)
    derrGdT = (errFrwd[1]-errBack[1])/(2*finiteDifferenceTorque)

    J = np.array([[derrRdFx, derrRdFy, derrRdT],[derrGdFx, derrGdFy, derrGdT]])

    return J

def int_error_result(SF, xCoeffs, yCoeffs, cICoeffs, xorg, yorg):

    smesh = np.linspace(0, fullArcLength, globalLen)
    geometryDef = [xCoeffs, yCoeffs, cICoeffs]

    dgds0 = (SF[2]/(n) + yorg*SF[0] - xorg*SF[1])/(E*cI_s(0, cICoeffs))
    res = fixed_rk4(deform_ODE, np.array([0,x0,y0]), smesh, (SF[0], SF[1], geometryDef, dgds0))
    Rinitial = lin.norm([xorg,yorg])
    Rfinal   = lin.norm([res[1,-1],res[2,-1]])
    dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg,xorg))
    # print(np.arctan2(res[2,-1],res[1,-1])/deg2rad, np.arctan2(yorg,xorg)/deg2rad)
    # print(dBeta)
    err = np.array([Rinitial-Rfinal, res[0,-1]-dBeta])

    return err, res

def tune_stiffness(stiffnessTarg, dBetaTarg):
    
    torqueTarg = dBetaTarg*stiffnessTarg

    # treat out of plane thickness, # of arms as constant
    dragVector0 = [R0, R1, R2, R3, \
                  betaB, betaC, beta0, \
                  cIs[0], cIs[1], cIs[2], \
                  ctrlcIs[1], \
                  ctrlLengths[1], ctrlLengths[2],
                  fullArcLength]
    dragVector = dragVector0
    err = 1
    pctErr = 100
    stepSize = 0.000000001
    i = 0
    while pctErr>10e-3:
        # create spring for guess
        
        geometryDef, smesh = drag_vector_spring(dragVector)

        xorg = coord(smesh, geometryDef[0])
        yorg = coord(smesh, geometryDef[1])
        # establish error in stiffness as a result of guess
        res = deform_spring_by_torque(torqueTarg, geometryDef)
        dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
        stiffness = torqueTarg/dBeta
        err = stiffness - stiffnessTarg
        pctErr = err/stiffnessTarg
        if pctErr < 10e-3:
            break
        print("stiffness error:",err)
        # estimate gradient (pain)
        grad = estimate_grad(dragVector, torqueTarg, stiffness, err, yorg, xorg)
        print('grad', grad)
        # print("next step:",lin.norm(stepSize*grad))
        Jinv = 1/grad
        print('jinv', Jinv)
        # determine next set of parameters
        stepSize = 0.0000000001*err
        dragVector = dragVector - stepSize*grad
        print("deltaSize:",lin.norm(stepSize*grad))
        # dragVector = dragVector - Jinv*err
        i+=1
    print("stiffness refinement iterations:", i)
    return stiffness, res, dragVector, dragVector0
    
def estimate_grad(dragVector, torqueTarg, stiffness, err, yorg, xorg):
    grad = np.empty(len(dragVector))
    dragVector1 = dragVector
    # take a different finite difference based on unit of axis
    print(len(dragVector))
    for i in range(len(dragVector)):
        if i < 4 or i>9:
            print("length")
            dragVectorU=dragVector1
            dragVectorU[i]=dragVector1[i]+finiteDifferenceLength
            geometryDefNew, smesh = drag_vector_spring(dragVectorU)
            res = deform_spring_by_torque(torqueTarg, geometryDefNew)
            dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
            stiffnessd = torqueTarg/dBeta
            errf = stiffness - stiffnessd
            grad[i]=(errf-err)/(finiteDifferenceLength)
        if i>3 and i<7:
            print("angle")
            dragVectorU=dragVector1
            dragVectorU[i]=dragVector1[i]+finiteDifferenceAngle
            geometryDefNew, smesh = drag_vector_spring(dragVectorU)
            res = deform_spring_by_torque(torqueTarg, geometryDefNew)
            dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
            stiffnessd = torqueTarg/dBeta
            errf = stiffness - stiffnessd
            grad[i]=(errf-err)/(finiteDifferenceAngle)
        if i>6 and i<10:
            print("Ic")
            dragVectorU=dragVector1
            dragVectorU[i]=dragVector1[i]+finiteDifferenceLength
            geometryDefNew, smesh = drag_vector_spring(dragVectorU)
            res = deform_spring_by_torque(torqueTarg, geometryDefNew)
            dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
            stiffnessd = torqueTarg/dBeta
            errf = stiffness - stiffnessd
            grad[i]=(errf-err)/(finiteDifferenceCI)
    return grad
    
def drag_vector_spring(dragVector):
    x0 = dragVector[0]
    y0 = 0
    pts = np.array([[x0, y0],[dragVector[1]*np.cos(dragVector[4]),dragVector[1]*np.sin(dragVector[4])], \
                    [dragVector[2]*np.cos(dragVector[5]),dragVector[2]*np.sin(dragVector[5])], \
                    [dragVector[3]*np.cos(dragVector[6]),dragVector[3]*np.sin(dragVector[6])]])
    cIs  = np.array([dragVector[7], dragVector[8], dragVector[9]])

    ctrlcIs      = np.array([0,dragVector[10],dragVector[13]])
    ctrlLengths = np.array([0,dragVector[11],dragVector[12],dragVector[13]])

    geometryDef = form_spring(pts, cIs, ctrlLengths, ctrlcIs)
    smesh = np.linspace(0, fullArcLength, globalLen)

    return geometryDef, smesh

deg2rad = np.pi/180
n = 2
finiteDifferenceLength = 0.0001
finiteDifferenceAngle  = .1*deg2rad
finiteDifferenceForce  = 0.5
finiteDifferenceTorque = 1
finiteDifferenceCI     = 0.00001

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


R0 = 2.2/2
R3 = 5.9/2*.9

R1 = (R0+R3)/2+.25
R2 = (R0+R3)/2-.25

betaB = 150/3*deg2rad*.5
betaC = 2*150/3*deg2rad
betaD = 150*deg2rad
beta0 = betaD

x0 = R0
y0 = 0

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
cIs  = np.array([.01, .0005, .01])

ctrlcIs      = np.array([0,fullArcLength*.6,fullArcLength])
ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

maxTorque = 13541.64
maxDBeta  = 0.087266463

stiffness, res, dragVector, dragVector0 = tune_stiffness(maxTorque/maxDBeta, maxDBeta)
geometryDef, smesh = drag_vector_spring(dragVector)

print(dragVector-dragVector0)
print(dragVector0)
print(dragVector)

xorg = coord(smesh, geometryDef[0])
yorg = coord(smesh, geometryDef[1])

R0 = dragVector[0]
R1 = dragVector[1]
R2 = dragVector[2]
R3 = dragVector[3]
fullArcLength = dragVector[-1]

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
rc = np.empty(len(smesh))
rn = np.empty(len(smesh))
Ic = np.empty(len(smesh))
for i in range(len(smesh)):
    rc[i] = rc_rootfinding(smesh[i], geometryDef[0], geometryDef[1], geometryDef[2], False)
    rn[i] = r_n(smesh[i], geometryDef[0], geometryDef[1])
    Ic[i] = cI_s(smesh[i], geometryDef[2])
e = rc-rn
for i in range(len(e)):
    if np.isnan(e[i]):
        e[i] = 0
# print(e)
# print(rn)
print(rc[5])
rcx = xorg-e*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
rcy = yorg+e*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
# print("sin",np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1])))
# print("cos",np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1])))
plt.plot(rcx,rcy)

plt.figure(0)
plt.plot(np.transpose(res), label='results')
plt.plot(xorg)
plt.plot(yorg)
plt.legend()
# for i in range(len(smesh)):
#     print(xorg[i], yorg[i], res[1,i], res[2,i])
a = np.empty(len(smesh))
b = np.empty(len(smesh))
h = np.empty(len(smesh))
for i in range(len(smesh)):
    AB = a_b_rootfinding(smesh[i], geometryDef[0], geometryDef[1], geometryDef[2], False)
    # print(AB)
    a[i] = AB[0]
    b[i] = AB[1]
    if np.isinf(a[i]) or np.isinf(b[i]):
        h[i] = np.cbrt(12*cI_s(smesh[i],geometryDef[2])/outPlaneThickness)
    else:
        h[i] = b[i]-a[i]

e = Ic/(outPlaneThickness*h*rn)
nb = h/2+e
na = h/2-e

xb = -nb*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
xa = -na*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
yb = nb*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
ya = na*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

plt.figure("geometry results")
plt.plot(xorg+xb,yorg+yb)
plt.plot(xorg-xa,yorg-ya)
plt.plot(-(xorg+xb),-(yorg+yb))
plt.plot(-(xorg-xa),-(yorg-ya))

plt.figure(99)
plt.plot(smesh, cI_s(smesh, geometryDef[2]))
plt.plot(smesh, h)
plt.plot(smesh,a)
plt.plot(smesh,b)
plt.plot(smesh,rn)

# for i in range(len(smesh)):
#     print(coord(smesh[i],geometryDef[0]), coord(smesh[i],geometryDef[1]))

plt.show()