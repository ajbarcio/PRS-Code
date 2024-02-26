from mpl_toolkits import mplot3d
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from materials import *
from openpyxl import *
from scipy.integrate import solve_ivp as ODE45
from scipy import optimize as op
import os

def cI_poly(cIs, sk, sf):
    Mat = np.array([[0,0,0,0,1], \
                    [sk**4,sk**3,sk**2,sk,1], \
                    [sf**4,sf**3,sf**2,sf,1], \
                    [0,0,0,1,0], \
                    [4*sf**3,3*sf**2,2*sf**1,1,0], \
                    ])
    Targ = np.array([[cIs[0]], \
                     [cIs[1]], \
                     [cIs[2]], \
                     [0], \
                     [0], \
                     ])
    cICoeffs = lin.solve(Mat, Targ)
    # print(hCoeffs)
    return cICoeffs

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

def d_cI_d_s(s, cICoeffs):
    U = np.array([4*s**3, 3*s**2, 2*s, 1, 0])
    dcIdS = U.dot(cICoeffs)
    return dcIdS[0]

def alpha_poly(alphas, si, sj, sf):
    Mat = np.array([[0,0,0,0,0,1], \
                    [si**5,si**4,si**3,si**2,si,1], \
                    [sj**5,sj**4,sj**3,sj**2,sj,1], \
                    [sf**5,sf**4,sf**3,sf**2,sf,1], \
                    [0,0,0,0,1,0], \
                    [5*sf**4,4*sf**3,3*sf**2,2*sf,1,0]])
    Targ = np.array([[alphas[0]],[alphas[1]],[alphas[2]],[alphas[3]],[0],[0]])
    alphaCoeffs = lin.solve(Mat,Targ)
    # print(alphaCoeffs)
    return alphaCoeffs

def alpha(s, alphaCoeffs):

    if hasattr(s, "__len__"):
        alpha = np.empty(len(s))
        i = 0
        for value in s:
            U = np.array([value**5, value**4, value**3, value**2, value, 1])
            alpha[i] = U.dot(alphaCoeffs)[0]
            i+=1
        return alpha
    else:
        U = np.array([s**5, s**4, s**3, s**2, s, 1])
        alpha = U.dot(alphaCoeffs)
        return alpha[0]

def d_alpha_d_s(s, alphaCoeffs):
    U = np.array([5*s**4, 4*s**3, 3*s**2, 2*s, 1, 0])
    dAlphadS = U.dot(alphaCoeffs)
    return dAlphadS[0]

def r_n(s, alphaCoeffs):
    rn = abs(1/d_alpha_d_s(s, alphaCoeffs))
    # this may not always be valid
    if s == fullArcLength:
        rn = float('inf')
    return rn
## Only for geometry finding AFTER desired stiffness acquired

def a_b_rootfinding(s, alphaCoeffs, cICoeffs, printBool):
    rn =  r_n(s, alphaCoeffs)
    k = cI_s(s, cICoeffs)
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
    # establish expected closesness range for a/b compared to rn to establish initial guess based on linearization about 1
    factor = 1.2
    # approimation: ln(x)=c(x-1), c such that c(x-1) = ln(factor)
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

def create_profiles(alphaCoeffs, cICoeffs, smesh, geoBool):
    plotAlpha = np.empty(globalLen)
    plotCurvature = np.empty(globalLen)
    plotCI = np.empty(globalLen)
    plotRn = np.empty(globalLen)
    plotAB = np.empty((globalLen,2))
    plotH = np.empty(globalLen)
    for i in range(len(smesh)):
        plotAlpha[i] = alpha(i*globalStep,alphaCoeffs)
        plotCurvature[i] = d_alpha_d_s(i*globalStep,alphaCoeffs)
        plotCI[i] = cI_s(i*globalStep,cICoeffs)
        plotRn[i] = r_n(i*globalStep,alphaCoeffs)
        if geoBool:
        # print("starting to rootfind", i)
            plotAB[i] = a_b_rootfinding(i*globalStep,alphaCoeffs,cICoeffs,False)
        # print("finished rootfinding")
    # print("alpha:",plot0)
    # print("dalpha/ds:",plot1)
    # print("rn:",plot3)
    for i in range(len(smesh)):
        plotH[i] = plotAB[i][1]-plotAB[i][0] 
    if geoBool:
        plt.figure(0)
        plt.plot(smesh, plotAlpha)
        plt.plot(smesh, plotCurvature)
        plt.plot(smesh, plotCI)
        plt.plot(smesh, plotRn)
        plt.plot(smesh, plotAB)
        plt.plot(smesh, plotH)
        plt.show()
    if geoBool:
        return plotAlpha, plotCurvature, plotCI, plotRn, plotAB, plotH
    else:
        return plotAlpha, plotCurvature, plotCI, plotRn

def mesh_original_geometry(alphaCoeffs, smesh):
    rnDx = np.empty(len(smesh))
    rnDy = np.empty(len(smesh))
    rnX  = np.empty(len(smesh))
    rnY  = np.empty(len(smesh))
    for i in range(len(smesh)):
        if i==0:
            rnDx[i]=0
            rnDy[i]=0
            rnX[i] = xy0[0]
            rnY[i] = xy0[1]
        else:
            step = smesh[i]-smesh[i-1]
            rnDx[i]=step*np.cos(alpha(i*step,alphaCoeffs))
            rnDy[i]=step*np.sin(alpha(i*step,alphaCoeffs))
            rnX[i] = rnX[i-1]+rnDx[i]
            rnY[i] = rnY[i-1]+rnDy[i]
    return rnX, rnY

def mesh_deformed_geometry(alphaCoeffs, gamma, smesh):
    rnDx = np.empty(len(smesh))
    rnDy = np.empty(len(smesh))
    rnX  = np.empty(len(smesh))
    rnY  = np.empty(len(smesh))
    for i in range(len(smesh)):
        if i==0:
            rnDx[i]=0
            rnDy[i]=0
            rnX[i] = xy0[0]
            rnY[i] = xy0[1]
        else:
            step = smesh[i]-smesh[i-1]
            rnDx[i]=step*np.cos(alpha(i*step,alphaCoeffs)+gamma[i])
            rnDy[i]=step*np.sin(alpha(i*step,alphaCoeffs)+gamma[i])
            rnX[i] = rnX[i-1]+rnDx[i]
            rnY[i] = rnY[i-1]+rnDy[i]

    return rnX, rnY
    

def create_initial_spring(alphas, cIs, si, sj, sk, sf):
    
    alphaCoeffs = alpha_poly(alphas, si, sj, sf)
    cICoeffs     = cI_poly(cIs, sk, sf)
    smesh = np.linspace(0,fullArcLength,globalLen)
    
    return alphaCoeffs, cICoeffs, smesh

def deform_ODE(s, alphaCoeffs, cICoeffs, gamma, Fx, Fy):
    
    LHS = np.empty(2)
    LHS[0] = gamma[1]
    LHS[1] = (Fx*np.sin(alpha(s, alphaCoeffs)+gamma[0])-Fy*np.cos(alpha(s, alphaCoeffs)+gamma[0])-E*d_cI_d_s(s, cICoeffs)*gamma[1])/(cI_s(s, cICoeffs)*E)
    return LHS

def rk4_step(fun, alphaCoeffs, cICoeffs, u0, gamma0, du, Fx, Fy):
    #
    #  Get four sample values of the derivative.
    #
    f1 = fun ( u0,            alphaCoeffs, cICoeffs, gamma0, Fx, Fy)
    f2 = fun ( u0 + du / 2.0, alphaCoeffs, cICoeffs, gamma0 + du * f1 / 2.0, Fx, Fy)
    f3 = fun ( u0 + du / 2.0, alphaCoeffs, cICoeffs, gamma0 + du * f2 / 2.0, Fx, Fy)
    f4 = fun ( u0 + du,       alphaCoeffs, cICoeffs, gamma0 + du * f3, Fx, Fy)
    #
    #  Combine them to estimate the solution gamma at time T1 = T0 + DT.
    #
    gamma1 = gamma0 + du * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0

    return gamma1

def fixed_rk4(fun, alphaCoeffs, cICoeffs, gamma0, Fx, Fy, smesh):
    res = np.empty((len(gamma0), len(smesh)))
    step = smesh[1]-smesh[0]
    for i in range(len(smesh)):
        if i == 0:
            res[0,i] = gamma0[0]
            res[1,i] = gamma0[1]
        else:
            stepRes = rk4_step(fun, alphaCoeffs, cICoeffs, smesh[i-1], gamma0, step, Fx, Fy)
            res[0,i] = stepRes[0]
            res[1,i] = stepRes[1]
            gamma0 = stepRes        
    return res

def integral_s(fun, array, mesh):
    #right-hand Riemman sum
    i = 1
    sum = 0
    while i < len(array):
        step = mesh[i]-mesh[i-1]
        sum += (fun(mesh[i], array[i]))*step
        i += 1
    return sum

def integrand_x(s, gamma):
    return np.cos(alpha(s, alphaCoeffs)+gamma)
def integrand_y(s, gamma):
    return np.sin(alpha(s, alphaCoeffs)+gamma)
# attempting 0 force


def forward_integration_result(fun, alphaCoeffs, cICoeffs, gamma0, Fx, Fy, smesh, plotBool, itr, resolution):
    res = fixed_rk4(fun, alphaCoeffs, cICoeffs, gamma0, Fx, Fy, smesh)
    gamma = res[0,:]
    [rnXd, rnYd] = mesh_deformed_geometry(alphaCoeffs, gamma, smesh)
    [rnX, rnY]   = mesh_original_geometry(alphaCoeffs, smesh)

    dgdsL = res[1,-1]

    xdis = rnXd[-1]
    ydis = rnYd[-1]

    xorg = rnX[-1]
    yorg = rnY[-1]

    rdis   = np.sqrt(xdis**2+ydis**2)
    rorg   = np.sqrt(xorg**2+yorg**2)

    dBeta = (np.arctan2(ydis,xdis)-np.arctan2(yorg,xorg))

    rErr = rdis-rorg
    gammaErr = gamma[-1]-dBeta

    if plotBool:
        
        cmap = plt.get_cmap('cool')
        color = cmap(itr/float(resolution))
        plt.figure(0)
        plt.plot(rnXd,rnYd, c=color, linewidth=3.0)
        plt.figure(1)
        cmap = plt.get_cmap('winter')
        color = cmap(itr/float(resolution))
        if itr/float(resolution) == 1:
            color = 'orange'
        plt.plot(smesh, gamma/deg2rad, c = color)
        cmap = plt.get_cmap('autumn')
        color = cmap(itr/float(resolution))
        plt.plot(smesh, res[1,:]/deg2rad, c = color)
        if itr == 0:
            cmap = plt.get_cmap('cool')
            plt.figure(0)
            color = cmap(1.0)
            plt.plot(rnX,rnY, c=color, linewidth=3.0)
            theta = np.linspace(0, 2*np.pi, 100)
            outerCircleX = rorg*np.cos(theta)
            outerCircleY = rorg*np.sin(theta)
            plt.plot(outerCircleX,outerCircleY)

    return rErr, gammaErr, dBeta, xdis, ydis, dgdsL, res

def deform_spring_byDeflection(alphaCoeffs, cICoeffs, torqueTarg):
    n=2

def deform_spring_byTorque(alphaCoeffs, cICoeffs, torqueTarg, torqueStepRes, plotBool):  
    n = 2
    itr = 0
    smesh = np.linspace(0, fullArcLength, globalLen)
    hJ = 1 # Lbs
    hT = torqueTarg/torqueStepRes # lbs-in

    F  = [0, 0, hT]
    dF = [0, 0, 0]
    gains = np.array([1, 1, 50])

    Fx = F[0]
    Fy = F[1]
    torque = F[2]
    fevalCounter = 0
    e = np.array([1, 1, 1])

    # reduce error in all 3 directions simultaneously:
    while lin.norm(e) > 10e-6:

        F = F+dF
        Fx = F[0]
        Fy = F[1]
        torque = F[2]
        [rnX, rnY]   = mesh_original_geometry(alphaCoeffs, smesh)
        dgds0 = (torque/(n) + rnY[0]*Fx - rnX[0]*Fy)/(E*cI_s(0, cICoeffs))
        torqueCheck = n*(Fy*rnX[0]-Fx*rnY[0]+E*cI_s(0, cICoeffs)*dgds0)
        # print("torquecheck error:",torqueCheck-torque)
        gamma0 = np.array([0, dgds0])
        xorg = rnX[-1]
        yorg = rnY[-1]
        rorg = np.sqrt(xorg**2+yorg**2)
        betaorg = np.arctan2(yorg,xorg)

        # above lines should be moved into forward_integration_result

        cErr, gammaErr, dBeta, xdis, ydis, dgdsL, resG = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx, Fy, smesh, plotBool, itr, torqueStepRes)
        tErr = n*(rnX[0]*Fy-rnY[0]*Fx+E*cI_s(0, cICoeffs)*dgds0)-torqueTarg
        fevalCounter+=1

        e = np.array([cErr, gammaErr, tErr]) # RESULT
        print(e)

        # finite difference in Fx
        Fx = Fx+hJ
        Fy = F[1]
        J11, J21, dBeta, xdisG, ydisG, dgdsLG, resG = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx, Fy, smesh, False, itr, torqueStepRes)
        J31 = n*(rnX[0]*Fy-rnY[0]*Fx+E*cI_s(0, cICoeffs)*dgds0)-torqueTarg
        fevalCounter+=1
        # finite difference in Fy
        Fx = F[0]
        Fy = Fy+hJ
        J12, J22, dBeta, xdisG, ydisG, dgdsLG, resG = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx, Fy, smesh, False, itr, torqueStepRes)
        J32 = n*(rnX[0]*Fy-rnY[0]*Fx+E*cI_s(0, cICoeffs)*dgds0)-torqueTarg
        fevalCounter+=1
        # finite difference in T
        Fx = F[0]
        Fy =  F[1]
        torque = torque+hJ
        dgds0 = (torque/(n) + rnY[0]*Fx - rnX[0]*Fy)/(E*cI_s(0, cICoeffs))
        J13, J23, dBeta, xdisG, ydisG, dgdsLG, resG = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx, Fy, smesh, False, itr, torqueStepRes)
        J33 = n*(rnX[0]*Fy-rnY[0]*Fx+E*cI_s(0, cICoeffs)*dgds0)-torqueTarg
        fevalCounter+=1
        # approximate jacobian with small finite difference:
        J = np.array([[J11-e[0],J12-e[0],J13-e[0]],[J21-e[1],J22-e[1],J23-e[1]],[J31-e[2],J32-e[2],J33-e[2]]])*1/(hJ)
        # print(J)
        # print(lin.pinv(J))
        # dF = np.array(lin.inv(J).dot(-e))
        # print("orig dF", dF)
        dF = np.array(lin.inv(J).dot(-e).dot(gains*(np.identity(3))))
        print("new dF", dF)

    # print solid mechanics properties of result
    if plotBool:
        rnXd, rnYd = mesh_deformed_geometry(alphaCoeffs, res[0,:], smesh)
        plt.figure("Moment, Shear")
        Fs = -Fx*np.sin(alpha(smesh, alphaCoeffs)+res[0,:])+Fy*np.cos(alpha(smesh, alphaCoeffs)+res[0,:])
        plt.plot(smesh, Fs, label = "shear force N") # shear FORCE along the beam
        plt.plot(smesh, (res[1,:]*E*(cI_s(smesh, cICoeffs))), label="bending M") # internal bending moment along the beam
        Fsx = -Fs*np.sin(alpha(smesh, alphaCoeffs)+res[0,:])
        Fsy =  Fs*np.cos(alpha(smesh, alphaCoeffs)+res[0,:])
        plt.plot(smesh, rnXd*Fy-rnYd*Fx, label="M due to force?") # shear MOMENT along the beam
        Tt = (res[1,:]*E*(cI_s(smesh, cICoeffs)))+rnXd*Fy-rnYd*Fx #vs Fsy and Fsx
        plt.plot(smesh, Tt, label="total torque reaction") # total torque load
        plt.axhline(torqueTarg/2) #torque reaction target
        plt.legend(loc='best')

    print("globalRes:", globalRes)
    print('feval per torque iteration:', fevalCounter/torqueStepRes)
    print('total feval', fevalCounter)
    # print("original torque:", torqueTarg)
    initialTorque = n*(Fy*rnX[0]-Fx*rnY[0]+E*cI_s(0, cICoeffs)*dgds0)
    finalTorque = n*(Fy*xdis-Fx*ydis+E*cI_s(fullArcLength, cICoeffs)*dgdsL)
    print("Torque Integration Error:", finalTorque-torqueTarg)
    print("initial Torque Error:", initialTorque-torqueTarg)

    return dBeta, cErr, gammaErr


deg2rad = np.pi/180
# E = 27.5*10**6
E = 1000000

fullArcLength = 6
globalRes = 100
globalLen = globalRes+1
globalStep = fullArcLength/globalRes
globalMaxIndex = globalLen-1

alphas  = [0,100*deg2rad,190*deg2rad,135*deg2rad]
# alphas = [0,0,0,0]
# alphas  = [0,30*deg2rad,60*deg2rad,90*deg2rad]
# alphas  = [0,30*deg2rad,-30*deg2rad,0]
cIs     = [.1,.05,.1]
# alphas = [0,45*deg2rad,135*deg2rad,180*deg2rad]
# cIs     = [.05,.05,.05]

xy0    = [1,0]

si = fullArcLength/3.0
sj = 2.0*fullArcLength/3.0
sk = fullArcLength/2.0
sf = fullArcLength

outPlaneThickness = 0.375

[alphaCoeffs, cICoeffs, smesh] = create_initial_spring(alphas, cIs, si, sj, sk, sf)
globalRes = 25
globalLen = globalRes+1
[dBeta, rErr, gammaErr] = deform_spring_byTorque(alphaCoeffs, cICoeffs, 13541.641, 50, True)
# for i in range(6):
#     globalRes = 2*globalRes
#     globalLen = globalRes+1
#     plotBool = False
#     if i == 5:
#         plotBool = True
#     [dBeta, rErr, gammaErr] = deform_spring_byTorque(alphaCoeffs, cICoeffs, 13541.641, 50, plotBool)

print(dBeta*1/deg2rad)
print(rErr, gammaErr)
plt.show()

#  [alph, curvature, cI, rn, ab, h] = create_profiles(alphaCoeffs, cICoeffs, smesh, True)
# dBetaTarg = 10*deg2rad
# Fguess    = 1000
# [torque, Fnew, springK] = deform_spring_byAngle(alphaCoeffs, cICoeffs, dBetaTarg, Fguess)
# print("torque",torque,"force:",Fnew,"K:",springK)
# approxK = eval_stiffness(alphaCoeffs, cICoeffs, 5000, 10, True)
# print(approxK)
# create_profiles(alphaCoeffs, cICoeffs, smesh, True)

#  finalParameters = tune_stiffness(5000,3000)





