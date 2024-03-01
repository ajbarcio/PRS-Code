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
            # print("im trying")
        # print("finished rootfinding")
    # print("alpha:",plot0)
    # print("dalpha/ds:",plot1)
    # print("rn:",plot3)
    for i in range(len(smesh)):
        plotH[i] = plotAB[i][1]-plotAB[i][0]
        if np.isnan(plotH[i]):
            plotH[i] = np.cbrt(12*cI_s(smesh[i], cICoeffs)/outPlaneThickness)
    if geoBool:
        plt.figure(98)
        plt.plot(smesh, plotAlpha)
        plt.plot(smesh, plotCurvature)
        plt.plot(smesh, plotCI)
        plt.plot(smesh, plotRn)
        plt.plot(smesh, plotAB)
        plt.plot(smesh, plotH)
        # plt.show()
    if geoBool:
        return plotAlpha, plotCurvature, plotCI, plotRn, plotAB, plotH
    else:
        return plotAlpha, plotCurvature, plotCI, plotRn

def mesh_original_geometry(alphaCoeffs, smesh):
    # Riemann Sum

    rnDx = np.empty(len(smesh))
    rnDy = np.empty(len(smesh))
    rnXG  = np.empty(len(smesh))
    rnYG  = np.empty(len(smesh))
    for i in range(len(smesh)):
        if i==0:
            rnDx[i]=0
            rnDy[i]=0
            rnXG[i] = xy0[0]
            rnYG[i] = xy0[1]
        else:
            step = smesh[i]-smesh[i-1]
            rnDx[i]=step*np.cos(alpha(i*step,alphaCoeffs))
            rnDy[i]=step*np.sin(alpha(i*step,alphaCoeffs))
            rnXG[i] = rnXG[i-1]+rnDx[i]
            rnYG[i] = rnYG[i-1]+rnDy[i]

    # RK
    
    rnX = fixed_rk4(xdis_integrand, xy0[0], smesh, alphaCoeffs)
    rnY = fixed_rk4(ydis_integrand, xy0[1], smesh, alphaCoeffs)

    # print(rnX[-1]-rnXG[-1])
    # print(rnY[-1]-rnYG[-1])

    return rnX, rnY

def mesh_deformed_geometry(alphaCoeffs, gamma, smesh):

    step = fullArcLength/globalRes

    # rnDx = np.empty(len(smesh))
    # rnDy = np.empty(len(smesh))
    # rnX  = np.empty(len(smesh))
    # rnY  = np.empty(len(smesh))
    # for i in range(len(smesh)):
    #     if i==0:
    #         rnDx[i]=0
    #         rnDy[i]=0
    #         rnX[i] = xy0[0]
    #         rnY[i] = xy0[1]
    #     else:
    #         step = smesh[i]-smesh[i-1]
    #         rnDx[i]=step*np.cos(alpha(i*step,alphaCoeffs)+gamma[i])
    #         rnDy[i]=step*np.sin(alpha(i*step,alphaCoeffs)+gamma[i])
    #         rnX[i] = rnX[i-1]+rnDx[i]
    #         rnY[i] = rnY[i-1]+rnDy[i]

    rnX = fixed_rk4(xdis_integrand, xy0[0], smesh, alphaCoeffs, gamma, step)
    rnY = fixed_rk4(ydis_integrand, xy0[1], smesh, alphaCoeffs, gamma, step)

    return rnX, rnY
    

def create_initial_spring(alphas, cIs, si, sj, sk, sf):
    
    alphaCoeffs = alpha_poly(alphas, si, sj, sf)
    cICoeffs     = cI_poly(cIs, sk, sf)
    smesh = np.linspace(0,fullArcLength,globalLen)
    
    return alphaCoeffs, cICoeffs, smesh

def xdis_integrand(s, y, *args): # alphaCoeffs, gamma, step
    # print(args[0][0])
    alphaCoeffs = args[0][0][0]
    
    # print("without gamma", alpha(s, alphaCoeffs))
    if len(args[0][0])>1:
        step = args[0][0][2]
        i = int(s/step-1)
        gamma   = args[0][0][1]
        argument    = alpha(s, alphaCoeffs) + gamma[i]
        argument = argument
        # print("with gamma", argument)
    else:
        argument    = alpha(s, alphaCoeffs)   
        argument = argument
    return np.cos(argument)

def ydis_integrand(s, y, *args):

    alphaCoeffs = args[0][0][0]
    # print("without gamma", alpha(s, alphaCoeffs))
    if len(args[0][0])>1:
        step = args[0][0][2]
        i = int(s/step-1)
        gamma   = args[0][0][1]
        argument    = alpha(s, alphaCoeffs) + gamma[i]
        argument = argument
        # print("with gamma", argument)
    else:
        argument    = alpha(s, alphaCoeffs) 
        argument = argument 
    return np.sin(argument)

def deform_ODE(s, gamma, *args): # alphaCoeffs, cICoeffs, Fx, Fy
    
    # print(args[0][0][0])

    alphaCoeffs = args[0][0][0]
    cICoeffs =    args[0][0][1]
    Fx =          args[0][0][2]
    Fy =          args[0][0][3]

    LHS = np.empty(2)
    LHS[0] = gamma[1]
    LHS[1] = (Fx*np.sin(alpha(s, alphaCoeffs)+gamma[0])-Fy*np.cos(alpha(s, alphaCoeffs)+gamma[0])-E*d_cI_d_s(s, cICoeffs)*gamma[1])/(cI_s(s, cICoeffs)*E)
    return LHS

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

def fixed_rk4(fun, y0, xmesh, *args): # (fun, alphaCoeffs, cICoeffs, y0, Fx, Fy, xmesh)
    step = xmesh[1]-xmesh[0]
    if hasattr(y0, '__len__'):
        res = np.empty((len(y0), len(xmesh)))
        for i in range(len(xmesh)):
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
            if i == 0:
                res[i] = y0
            else:
                stepRes = rk4_step(fun, xmesh[i], y0, step, args)
                res[i] = stepRes
                y0 = stepRes
            
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
    res = fixed_rk4(fun, gamma0, smesh, alphaCoeffs, cICoeffs, Fx, Fy) # fun, y0, xmesh, *args
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
        color = cmap(itr/552.0)
        plt.figure("geometry")
        plt.plot(rnXd,rnYd, c=color, linewidth=3.0)
        plt.figure("results")
        cmap = plt.get_cmap('winter')
        color = cmap(itr/552.0)
        if itr/float(resolution) == 1:
            color = 'orange'
        plt.figure("results")
        plt.plot(smesh, gamma/deg2rad, c = color)
        cmap = plt.get_cmap('autumn')
        color = cmap(itr/552.0)
        plt.plot(smesh, res[1,:]/deg2rad, c = color)
        if itr == 0:
            print("this should be original geometry")
            cmap = plt.get_cmap('cool')
            plt.figure("geometry")
            color = 'red'
            plt.plot(rnX,rnY, c=color, linewidth=4.0)
            plt.plot(-rnX,-rnY, c=color, linewidth=4.0)
            theta = np.linspace(0, 2*np.pi, 100)
            outerCircleX = rorg*np.cos(theta)
            outerCircleY = rorg*np.sin(theta)
            plt.plot(outerCircleX,outerCircleY)
        if itr >= 5000:
            plt.figure("geometry")
            color = 'green'
            plt.plot(rnXd,rnYd, c=color, linewidth=4.0)
            plt.plot(-rnXd,-rnYd, c=color, linewidth=4.0)
            theta = np.linspace(0, 2*np.pi, 100)
            outerCircleX = rorg*np.cos(theta)
            outerCircleY = rorg*np.sin(theta)
            plt.plot(outerCircleX,outerCircleY)
            plt.axis('equal')

            

    return rErr, gammaErr, dBeta, xdis, ydis, dgdsL, res

def deform_spring_byTorque(alphaCoeffs, cICoeffs, torqueTarg, torqueStepRes, plotBool):  
    Tgain = 0.9999

    divergingFlag = 0
    
    n = 2
    itr = 0
    smesh = np.linspace(0, fullArcLength, globalLen)
    hJ = 1 # Lbs
    hT = torqueTarg/torqueStepRes # lbs-in

    F = [0, 0, hT]

    Fx = F[0]
    Fy = F[1]
    torque = F[2]
    fevalCounter = 0

    # test error for small amount of pure bending:
    while abs(torque-torqueTarg)>10e-10:
        print("torque increasing step")
        [rnX, rnY]   = mesh_original_geometry(alphaCoeffs, smesh)
        # print(rnX)
        # print(rnY)
        dgds0 = (torque/(n) + rnY[0]*Fx - rnX[0]*Fy)/(E*cI_s(0, cICoeffs))
        torqueCheck = n*(Fy*rnX[0]-Fx*rnY[0]+E*cI_s(0, cICoeffs)*dgds0)
        # print("torquecheck error:",torqueCheck-torque)
        gamma0 = np.array([0, dgds0])
        xorg = rnX[-1]
        yorg = rnY[-1]
        rorg = np.sqrt(xorg**2+yorg**2)
        betaorg = np.arctan2(yorg,xorg)

        # above lines should be moved into forward_integration_result

        cErr, gammaErr, dBeta, xdis, ydis, dgdsL, resG = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx, Fy, smesh, plotBool, fevalCounter, torqueStepRes)
        fevalCounter+=1

        e = np.array([cErr, gammaErr]) # RESULT
        eprev = e

        # print(e)

        # finite difference in Fx
        Fx1 = Fx+hJ
        Fy1 = F[1]
        J11, J21, dBetaG, xdisG, ydisG, dgdsLG, resG = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx1, Fy1, smesh, False, fevalCounter, torqueStepRes)
        fevalCounter+=1
        # finite difference in Fy
        Fx2 = F[0]
        Fy2 = Fy+hJ
        J12, J22, dBetaG, xdisG, ydisG, dgdsLG, resG = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx2, Fy2, smesh, False, fevalCounter, torqueStepRes)
        fevalCounter+=1
        # approximate jacobian with small finite difference:
        J = np.array([[J11-e[0],J12-e[0]],[J21-e[1],J22-e[1]]])*1/(hJ)
        # print(J)
        dF = np.array(lin.pinv(J).dot(-e))
        dF = np.append(dF, 0)
        # print(dF)
        errorDecreaseCount = 0

        while lin.norm(e)>10**-10:
            print("error decreasing step")
            errGain = 1
            F = F+dF
            # see how good your jacobian approximation was
            Fx = F[0]
            Fy = F[1]
            torque = F[2]
            dgds0 = (torque/(n) + rnY[0]*Fx - rnX[0]*Fy)/(E*cI_s(0, cICoeffs))
            gamma0 = np.array([0, dgds0])
            
            cErr, gammaErr, dBeta, xdis, ydis, dgdsL, res = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx, Fy, smesh, plotBool, fevalCounter, torqueStepRes)
            fevalCounter+=1

            e = np.array([cErr, gammaErr])
            if lin.norm(e)>lin.norm(eprev):
                divergingFlag = 1
                print("diverging")
            # if divergingFlag:
            #     break

            # finite difference in Fx
            Fx1 = Fx+hJ
            Fy1 = F[1]
            J11, J21, dBetaG, xdisG, ydisG, dgdsLG, resG = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx1, Fy1, smesh, False, fevalCounter, torqueStepRes)
            fevalCounter+=1
            # finite difference in Fy
            Fx2 = F[0]
            Fy2 = Fy+hJ
            J12, J22, dBetaG, xdisG, ydisG, dgdsLG, resG = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx2, Fy2, smesh, False, fevalCounter, torqueStepRes)
            fevalCounter+=1
            # approximate jacobian with small finite difference:
            J = np.array([[J11-e[0],J12-e[0]],[J21-e[1],J22-e[1]]])*1/(hJ)
            # print(J)
            dF = np.array(lin.pinv(J).dot(-e))*errGain
            eprev = e
            dF = np.append(dF, 0)
            # print(dF)
            # print(lin.norm(e))
            errorDecreaseCount += 1
            itr+=1

        # if i >= torqueStepRes-1:
        #     print("FINAL planar force vector:", F, "i", i)
        # else:

        if abs(torque-torqueTarg)>10e-10:
            F[2] = F[2]+(torqueTarg-torque)*Tgain
            torque = F[2]
        else:
            break
        itr+=1
        
        # print("Final planar force vector:", F)
    print("total feval", fevalCounter)
    # print("final torque try:", torque)
    # print("ideal dgds0             ", (torqueTarg/(n) + rnY[0]*Fx - rnX[0]*Fy)/(E*cI_s(0, cICoeffs)))
    # print("final dgds0 (calculated)", (torque/(n) + rnY[0]*Fx - rnX[0]*Fy)/(E*cI_s(0, cICoeffs)))
    # print("final dgds0 (called)    ", dgds0)
    # print("torque based on final dgds0: (calculated)", n*(Fy*rnX[0]-Fx*rnY[0]+E*cI_s(0, cICoeffs)*(torque/(n) + rnY[0]*Fx - rnX[0]*Fy)/(E*cI_s(0, cICoeffs))))
    # print("torque based on final dgds0: (called)    ", n*(Fy*rnX[0]-Fx*rnY[0]+E*cI_s(0, cICoeffs)*dgds0))
    # print("torque based on ideal dgds0:             ", n*(Fy*rnX[0]-Fx*rnY[0]+E*cI_s(0, cICoeffs)*(torque/(n) + rnY[0]*Fx - rnX[0]*Fy)/(E*cI_s(0, cICoeffs))))
    # print("torque target                            ", torqueTarg)

    dgds0 = (torque/(n) + rnY[0]*Fx - rnX[0]*Fy)/(E*cI_s(0, cICoeffs))
    gamma0 = np.array([0, dgds0])
    print(F)
    cErrG, gammaErrG, dBetaG, xdis, ydis, dgdsL, res = forward_integration_result(deform_ODE, alphaCoeffs, cICoeffs, gamma0, Fx, Fy, smesh, True, 100000, torqueStepRes)
    initialTorque = n*(Fy*rnX[0]-Fx*rnY[0]+E*cI_s(0, cICoeffs)*dgds0)
    finalTorque = n*(Fy*xdis-Fx*ydis+E*cI_s(fullArcLength, cICoeffs)*dgdsL)

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
        Tt = (res[1,:]*E*(cI_s(smesh, cICoeffs)))+rnXd*Fy-rnYd*Fx 
        plt.plot(smesh, Tt, label="total torque reaction") # total torque load
        TtAlt = (res[1,:]*E*(cI_s(smesh, cICoeffs)))+rnX*Fy-rnY*Fx # ???
        # plt.plot(smesh, TtAlt, label="total torque reaction (but I fucked with it)") # total torque load
        plt.axhline(torqueTarg/2,xmin = 0, xmax = fullArcLength) #torque reaction target
        plt.legend(loc='best')   

        plt.figure("geometry")


    # print("Torque Integration Error:", finalTorque-torqueTarg)
    # print("initial Torque Error:", initialTorque-torqueTarg)
    print(globalRes, finalTorque-torqueTarg)
    # print(cErr, gammaErr)
    # print(cErrG, gammaErrG)
    return dBeta, cErr, gammaErr, res

def tune_stiffness(desiredStiffness):
    # constants:
    outerRadius = 5.5
    
    # first guess:
    fullArcLength = 6
    innerRadius = 0.5  

    xy0      = np.array([innerRadius,0])
    alphas   = np.array([0,90*deg2rad,180*deg2rad,162*deg2rad])
    cIs      = np.array([.01,.005,.01])
    s_checks = np.array([0.25,0.5,0.333])

    si = s_checks[0]
    sj = s_checks[1]
    sk = s_checks[2]
    sf = fullArcLength
    
    outPlaneThickness = 0.375

    init_perturb = 0.01

    parameterVector = np.concatenate((alphas[1:-1],cIs,s_checks,outPlaneThickness))
    err = 1
    while err>0.01:
        springK = check_stiffness(parameterVector)


def check_stiffness(parameterVector, torque):
    alphas =   parameterVector[0:2]
    cIs    =   parameterVector[3:5]
    s_checks = parameterVector[6:8]

    si = s_checks[0]
    sj = s_checks[1]
    sk = s_checks[2]
    sf = fullArcLength

    outPlaneThickness = parameterVector[-1]

    [alphaCoeffs, cICoeffs, smesh] = create_initial_spring(alphas, cIs, si, sj, sk, sf)
    [dBeta, cErr, gammaErr, res] = deform_spring_byTorque(alphaCoeffs, cICoeffs, torque, 50, True)
    springK = torque/dBeta

    return springK


deg2rad = np.pi/180
E = 27.5*10**6
# E = 1000000

fullArcLength = 6
globalRes = 100
globalLen = globalRes+1
globalStep = fullArcLength/globalRes
globalMaxIndex = globalLen-1

xy0    = np.array([.25,0])

alphas  = [0,90*deg2rad,180*deg2rad,162*deg2rad]
# alphas = [0,0,0,0]
# alphas  = [0,15*deg2rad,30*deg2rad,45*deg2rad]
# alphas  = [0,30*deg2rad,-30*deg2rad,0]
cIs     = np.array([.01,.005,.01])
# alphas = [0,45*deg2rad,135*deg2rad,180*deg2rad]
# cIs     = [.05,.005,.05]

si = fullArcLength/4.0
sj = fullArcLength/2.0
sk = fullArcLength*1/3.0
sf = fullArcLength

outPlaneThickness = 0.25
globalRes = 100 # 2284 for desired acuracy
globalLen = globalRes+1
globalStep = fullArcLength/globalRes
globalMaxIndex = globalLen-1
[alphaCoeffs, cICoeffs, smesh] = create_initial_spring(alphas, cIs, si, sj, sk, sf)

[dBeta, cErr, gammaErr, res] = deform_spring_byTorque(alphaCoeffs, cICoeffs, 13541.641, 50, True)
improvements = 0
for i in range(improvements):
    globalRes = 2*globalRes
    globalLen = globalRes+1
    globalStep = fullArcLength/globalRes
    globalMaxIndex = globalLen-1
    plotBool = False
    if i == improvements-1:
        plotBool = True
    [dBeta, cErr, gammaErr, res] = deform_spring_byTorque(alphaCoeffs, cICoeffs, 13541.641, 50, plotBool)

print(dBeta*1/deg2rad)
print(cErr, gammaErr)
print(lin.norm([cErr, gammaErr]))
# plt.show()
meshAlpha, meshCurvature, meshCI, meshRn, meshAB, meshH = create_profiles(alphaCoeffs, cICoeffs, smesh, True)
# print(meshAB)

meshA = np.empty(len(smesh))
meshB = np.empty(len(smesh))

for i in range(len(smesh)):
        meshA[i] = meshAB[i][0]
        meshB[i] = meshAB[i][1]
# print(meshA)
# print(meshB)
# print(meshH)
sigmaA = (E*(1-meshRn/meshA)*meshRn*res[1,:])
# print(sigmaA)
sigmaB = (E*(1-meshRn/meshB)*meshRn*res[1,:])
# print(sigmaB)
for i in range(len(smesh)):
    if np.isnan(sigmaA[i]) and np.isnan(sigmaB[i]):
        sigmaA[i] = -E*res[1,i]*0.5*meshH[i]
        sigmaB[i] =  E*res[1,i]*0.5*meshH[i]
maxSigma = max(np.max(abs(sigmaA)), np.max(abs(sigmaB)))
print(maxSigma)
yield_stress=309700
print(yield_stress)

plt.figure(99)
plt.plot(smesh,sigmaA)
plt.plot(smesh,sigmaB)
plt.plot(smesh,meshRn)
plt.plot(smesh,meshH*10**5)
plt.plot(smesh,-meshH*10**5)
plt.axhline(yield_stress, xmin = 0, xmax = fullArcLength)
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