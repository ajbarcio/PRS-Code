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
from copy import deepcopy as dc

def cI_poly(cIs, ctrlcIs):
    Mat = np.array([[0,0,0,0,0,1], \
                    [ctrlcIs[1]**5,ctrlcIs[1]**4,ctrlcIs[1]**3,ctrlcIs[1]**2,ctrlcIs[1],1], \
                    [ctrlcIs[2]**5,ctrlcIs[2]**4,ctrlcIs[2]**3,ctrlcIs[2]**2,ctrlcIs[2],1], \
                    [0,0,0,0,1,0], \
                    [5*ctrlcIs[2]**4,4*ctrlcIs[2]**3,3*ctrlcIs[2]**2,2*ctrlcIs[2]**1,1,0], \
                    [5*ctrlcIs[1]**4,4*ctrlcIs[1]**3,3*ctrlcIs[1]**2,2*ctrlcIs[1]**1,1,0]
                    ])
    Targ = np.array([[cIs[0]], \
                     [cIs[1]], \
                     [cIs[2]], \
                     [0], \
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
        alphaList = np.arctan2(d_coord_d_s(s, yCoeffs),d_coord_d_s(s, xCoeffs))
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
    if (d2yds2/dxds-d2xds2*dyds/dxds**2)==0:
        rn = float('inf')
    else:
        rn = abs((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
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
            U = np.array([value**5, value**4, value**3, value**2, value, 1])
            cI[i] = U.dot(cICoeffs)[0]
            i+=1
        return abs(cI)
    else:
        U = np.array([s**5, s**4, s**3, s**2, s, 1])
        cI = U.dot(cICoeffs)
        return abs(cI[0])

# def d_cI_d_s(s, xCoeffs, yCoeffs, hCoeffs):
#     dcIds = (cI_s(s+finiteDifferenceLength, xCoeffs, yCoeffs, hCoeffs)-cI_s(s-finiteDifferenceLength, xCoeffs, yCoeffs, hCoeffs))/(2*finiteDifferenceLength)
#     return dcIds

def d_cI_d_s(s, cICoeffs):
    if hasattr(s, "__len__"):
        dcIcs = np.empty(len(s))
        i = 0
        for value in s:
            U = np.array([5*value**4, 4*value**3, 3*value**2, 2*value**1, 1, 0])
            dcIcs[i] = U.dot(cICoeffs)[0]
            i+=1
        return dcIcs
    else:
        U = np.array([5*s**4, 4*s**3, 3*s**2, 2*s, 1, 0])
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

def form_spring(pts, cIs, ctrlLengths, ctrlcIs):
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

    SF = np.array([0, 0, torqueTarg])
    # print("before forward integration", SF[2])
    err = np.ones(2)
    i = 0
    while lin.norm(err) > 10e-6 and i<500:
        # establish results of guess
        errPrev = err
        err, res = int_error_result(SF, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
        # estimate jacobian with finite difference
        J = fd_J(SF, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
        # establish next guess:
        SF = SF - lin.pinv(J).dot(err)
        # print("torque deform error:",lin.norm(err))
        i+=1
        if lin.norm(err)>lin.norm(errPrev):
            print("torque deform diverging", i)
    # print("torque iterations: ",i)
    # print("after forward integration", SF[2])
    return res
        

def fd_J(SF, xCoeffs, yCoeffs, cICoeffs, xorg, yorg):
    # finite difference in Fx
    SFback = dc(SF-np.array([finiteDifferenceForce,0,0]))
    SFfrwd = dc(SF+np.array([finiteDifferenceForce,0,0]))
    errBack, resG = int_error_result(SFback, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    errFrwd, resG = int_error_result(SFfrwd, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    derrRdFx = (errFrwd[0]-errBack[0])/(2*finiteDifferenceForce)
    derrGdFx = (errFrwd[1]-errBack[1])/(2*finiteDifferenceForce)
    # finite difference in Fy
    SFback = dc(SF-np.array([0,finiteDifferenceForce,0]))
    SFfrwd = dc(SF+np.array([0,finiteDifferenceForce,0]))
    errBack, resG = int_error_result(SFback, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    errFrwd, resG = int_error_result(SFfrwd, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    derrRdFy = (errFrwd[0]-errBack[0])/(2*finiteDifferenceForce)
    derrGdFy = (errFrwd[1]-errBack[1])/(2*finiteDifferenceForce)
    # finite difference in T
    # SFback = SF-np.array([0,0,finiteDifferenceTorque])
    # SFfrwd = SF+np.array([0,0,finiteDifferenceTorque])
    # errBack, resG = int_error_result(SFback, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    # errFrwd, resG = int_error_result(SFfrwd, xCoeffs, yCoeffs, cICoeffs, xorg, yorg)
    # derrRdT = (errFrwd[0]-errBack[0])/(2*finiteDifferenceTorque)
    # derrGdT = (errFrwd[1]-errBack[1])/(2*finiteDifferenceTorque)

    J = np.array([[derrRdFx, derrRdFy, 0],[derrGdFx, derrGdFy, 0]])

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


def tune_stiffness(stiffnessTarg, dBetaTarg, dragVector, discludeVector):
    
    torqueTarg = dBetaTarg*stiffnessTarg
    print("initial guess check,", dragVector)
    # treat out of plane thickness, # of arms as constant
    # #PASS DRAG VECTOR
    dragVector0=dc(dragVector)
    # dragVector0 = [R0, R1, R2, R3, \
    #               betaB, betaC, beta0, \
    #               cIs[0], cIs[1], cIs[2], \
    #               ctrlcIs[1], \
    #               ctrlLengths[1], ctrlLengths[2],
    #               fullArcLength]
    # # print(dragVector0)
    # dragVector = dc(dragVector0)
    err = 1
    relErr = 1
    convErr = 1
    stepSizeCoeff = 1
    j = 0
    resetCounter = 0
    resetStatus  = 0
    while relErr>10e-6 and convErr>10e-6:
        print("stiffness refinement iteration:",j)
        # create spring for guess
        
        geometryDef, smesh = drag_vector_spring(dragVector)
        # print(geometryDef[2],dragVector[7])
        xorg = coord(smesh, geometryDef[0])
        yorg = coord(smesh, geometryDef[1])
        # establish error in stiffness as a result of guess
        relErrPrev = relErr
        res = deform_spring_by_torque(torqueTarg, geometryDef)
        dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
        stiffness = torqueTarg/dBeta
        err = abs(stiffness - stiffnessTarg)
        relErr = abs(err/stiffnessTarg)
        convErr = abs(relErr-relErrPrev)
        print("stiffness error:",err)
        print("stiffness relative error:",relErr)
        print("chenge in rel error", convErr)
        # estimate gradient (pain)
        # print(dragVector0)
        if abs(relErr) > abs(relErrPrev):
            dragVector = dc(dragVectorPrev)
            convErr = 0
        else:
            grad = estimate_grad(dragVector, torqueTarg, stiffness, err, yorg, xorg)
            Jinv = 1/grad
            for i in range(len(grad)):
                if discludeVector[i]==False:
                    grad[i]=0
                    Jinv[i]=0
                #### WTF IS GOING ON HERE
            # print(dragVector0)
            # print('grad', grad)
            # print("next step:",lin.norm(stepSize*grad))
            # Jinv = 1/grad
            # print('jinv', Jinv)
            # determine next set of parameters
            stepSize = stepSizeCoeff*lin.norm(grad)
            dragVectorPrev = dc(dragVector)
            # dragVector = dragVector + stepSize*grad
            # print(Jinv)
            dragVector = dragVector + err*Jinv*stepSize
            # print(dragVector-dragVectorPrev)
            for i in range(len(grad)):
                if discludeVector[i]==False:
                    assert(dragVector0[i]==dragVector[i])
        if violates_bounds(dragVector):
            print("reset---------------------------------------------------reset")
            # print(dragVector)
            # print(dragVector0)
            dragVector = dc(dragVectorPrev)
            # print(dragVector0)
            # print(dragVector)
            if resetStatus==0:
                if j==0:
                    stepSizeCoeff = stepSizeCoeff*0.1
                else:
                    stepSizeCoeff = stepSizeCoeff*0.1
            else:
                stepSizeCoeff = stepSizeCoeff*0.1
            relErr = 1
            resetCounter+=1
            resetStatus = 1
            print("this is the ",resetCounter,"th reset")
        else:
            assert(not violates_bounds(dragVector))
            print("reduction:", stepSizeCoeff)
            # stepSizeCoeff = stepSizeCoeff*1.1
            j+=1
            resetCounter = 0
            resetStatus  = 0
    assert(not violates_bounds(dragVector))
    for i in range(len(dragVector)):
                if discludeVector[i]==False:
                    # if not dragVector0[i]==dragVector[i]:
                    #     print(dragVector0[i], dragVector[i])
                    assert(dragVector0[i]==dragVector[i])
                    # print(i)
    print("stiffness refinement iterations:", j)
    return stiffness, res, dragVector, dragVector0

def violates_bounds(dragVector):
    truth0 = dragVector[0] < globalInnerRadiusLimit or dragVector[0] > globalOuterRadiusLimit
    truth1 = dragVector[1] < globalInnerRadiusLimit or dragVector[1] > globalOuterRadiusLimit or dragVector[1] < dragVector[0] or dragVector[1] > dragVector[3]
    truth2 = dragVector[2] < globalInnerRadiusLimit or dragVector[2] > globalOuterRadiusLimit or dragVector[2] < dragVector[0] or dragVector[2] > dragVector[3]
    truth3 = dragVector[3] < globalInnerRadiusLimit or dragVector[3] > globalOuterRadiusLimit

    truth4 = False # dragVector[4] > 90*deg2rad
    truth5 = dragVector[5] < dragVector[4] or dragVector[5] > dragVector[6]
    truth6 = False # dragVector[6] > 170*deg2rad

    truth7 = dragVector[7] < dragVector[8] or dragVector[7] > .05
    truth8 = dragVector[8] < 0 or dragVector[8] > .025
    truth9 = dragVector[9] < dragVector[8] or dragVector[9] > .05

    truth  = truth0 or truth1 or truth2 or truth3 or truth4 or truth5 or truth6 or truth7 or truth8 or truth9 
    return truth    

def estimate_grad(dragVector, torqueTarg, stiffness, err, yorg, xorg):
    grad = np.empty(len(dragVector))
    dragVectorBase = dc(dragVector)
    # print("base after assigned", dragVectorBase)
    # take a different finite difference based on unit of axis
    # print(len(dragVector))
    # print("doing this to prove its only called once")
    for i in range(len(dragVector)):
        # print(i)

        # DO THIS LATER

        # if (certain range)
            #difference = finite difference length
        # if eeeee

        # dragVectorFD[i] = dragVectorFD+difference

        if i < 4 or i>9:
            # print("length")
            # print("base", dragVectorBase)
            dragVectorFD=dc(dragVectorBase)
            dragVectorFD[i]=dragVectorFD[i]+finiteDifferenceLength
            geometryDefNew, smesh = drag_vector_spring(dragVectorFD)
            res = deform_spring_by_torque(torqueTarg, geometryDefNew)
            dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
            stiffnessd = torqueTarg/dBeta
            errf = stiffness - stiffnessd
            grad[i]=(errf-err)/(finiteDifferenceLength)
        if i>3 and i<7:
            # print("angle")
            # print("base", dragVectorBase)
            dragVectorFD=dc(dragVectorBase)
            dragVectorFD[i]=dragVectorFD[i]+finiteDifferenceAngle
            geometryDefNew, smesh = drag_vector_spring(dragVectorFD)
            res = deform_spring_by_torque(torqueTarg, geometryDefNew)
            dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
            stiffnessd = torqueTarg/dBeta
            errf = stiffness - stiffnessd
            grad[i]=(errf-err)/(finiteDifferenceAngle)
        if i>6 and i<10:
            # print("Ic")
            # print("base", dragVectorBase)
            dragVectorFD=dc(dragVectorBase)
            dragVectorFD[i]=dragVectorFD[i]+finiteDifferenceCI
            geometryDefNew, smesh = drag_vector_spring(dragVectorFD)
            res = deform_spring_by_torque(torqueTarg, geometryDefNew)
            dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
            stiffnessd = torqueTarg/dBeta
            errf = stiffness - stiffnessd
            grad[i]=(errf-err)/(finiteDifferenceCI)

    grad = grad/lin.norm(grad)

    return grad
    
def drag_vector_spring(dragVectorArg):
    x0 = dragVectorArg[0]
    y0 = 0
    pts_dv = np.array([[x0, y0],[dragVectorArg[1]*np.cos(dragVectorArg[4]),dragVectorArg[1]*np.sin(dragVectorArg[4])], \
                    [dragVectorArg[2]*np.cos(dragVectorArg[5]),dragVectorArg[2]*np.sin(dragVectorArg[5])], \
                    [dragVectorArg[3]*np.cos(dragVectorArg[6]),dragVectorArg[3]*np.sin(dragVectorArg[6])]])
    cIs_dv  = np.array([dragVectorArg[7], dragVectorArg[8], dragVectorArg[9]])
    # print(cIs_dv)

    ctrlcIs_dv      = np.array([0,dragVectorArg[10],dragVectorArg[13]])
    ctrlLengths_dv = np.array([0,dragVectorArg[11],dragVectorArg[12],dragVectorArg[13]])

    geometryDef = form_spring(pts_dv, cIs_dv, ctrlLengths_dv, ctrlcIs_dv)
    smesh = np.linspace(0, fullArcLength, globalLen)

    return geometryDef, smesh

deg2rad = np.pi/180
n = 2
finiteDifferenceLength = 0.0001
finiteDifferenceAngle  = .1*deg2rad
finiteDifferenceForce  = 0.5
finiteDifferenceTorque = 1
finiteDifferenceCI     = 0.00000001

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
betaD = 150*deg2rad
beta0 = betaD

x0 = R0
y0 = 0

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
cIs  = np.array([.008, .001, .008])

ctrlcIs      = np.array([0,fullArcLength*.6,fullArcLength])
ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

maxTorque = 13541.64
maxDBeta  = 0.087266463

dragVector0 = [R0, R1, R2, R3, \
               betaB, betaC, beta0, \
               cIs[0], cIs[1], cIs[2], \
               ctrlcIs[1], \
               ctrlLengths[1], ctrlLengths[2],
               fullArcLength]