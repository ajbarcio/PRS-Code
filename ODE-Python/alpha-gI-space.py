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
    U = np.array([s**4, s**3, s**2, s, 1])
    k = U.dot(cICoeffs)
    return k[0]

def d_cI_d_s(s, cICoeffs):
    U = np.array([4*s**3, 3*s**2, 2*s, 1, 0])
    dKdS = U.dot(cICoeffs)
    return dKdS[0]

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

def deform_ODE(s, alphaCoeffs, cICoeffs, gamma, F, ang):
    
    Fx = F*-np.sin(ang)
    Fy = F*np.cos(ang)
    LHS = np.empty(2)
    LHS[0] = gamma[1]
    LHS[1] = (Fx*np.sin(alpha(s, alphaCoeffs)+gamma[0])-Fy*np.cos(alpha(s, alphaCoeffs)+gamma[0])-E*d_cI_d_s(s, cICoeffs)*gamma[1])/(cI_s(s, cICoeffs)*E)
    return LHS

def rk4_step(fun, alphaCoeffs, cICoeffs, u0, gamma0, du, F, ang):
    #
    #  Get four sample values of the derivative.
    #
    f1 = fun ( u0,            alphaCoeffs, cICoeffs, gamma0, F, ang )
    f2 = fun ( u0 + du / 2.0, alphaCoeffs, cICoeffs, gamma0 + du * f1 / 2.0, F, ang )
    f3 = fun ( u0 + du / 2.0, alphaCoeffs, cICoeffs, gamma0 + du * f2 / 2.0, F, ang )
    f4 = fun ( u0 + du,       alphaCoeffs, cICoeffs, gamma0 + du * f3, F, ang )
    #
    #  Combine them to estimate the solution gamma at time T1 = T0 + DT.
    #
    gamma1 = gamma0 + du * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0

    return gamma1

def fixed_rk4(fun, alphaCoeffs, cICoeffs, gamma0, F, ang, smesh):
    res = np.empty((len(gamma0), len(smesh)))
    step = smesh[1]-smesh[0]
    for i in range(len(smesh)):
        if i == 0:
            res[0,i] = gamma0[0]
            res[1,i] = gamma0[1]
        else:
            stepRes = rk4_step(fun, alphaCoeffs, cICoeffs, smesh[i-1], gamma0, step, F, ang)
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


def deform_spring_byAngle(alphaCoeffs, cICoeffs, dBetaTarg, Fguess):
    print("first guess")
    [dBeta, torque, springK, gamma, iii, finalGuess] = deform_spring_byForce(alphaCoeffs, cICoeffs, Fguess, Fguess, True)
    plt.show()
    err = (abs(dBeta)-dBetaTarg)*1/deg2rad
    errPrev = err

    factor = 1.01

    Fprev = Fguess
    Fnew  = Fprev*factor
    print("ANGLE DEFORM RESULT: angle error:",err,"next force guess", Fnew)
    iiii = 0
    print("second guess; arbitrary step")
    while abs(err)>10e-6:
        print("gonna try a force")
        [dBeta, torque, springK, gamma, iii, finalGuess] = deform_spring_byForce(alphaCoeffs, cICoeffs, Fnew, Fnew, True)
        print("tried a force")
        plt.show()
        err = (abs(dBeta)-dBetaTarg)*1/deg2rad
        jac = (err-errPrev)/(Fnew-Fprev)
        Fprev = Fnew
        errPrev = err
        Fnew = Fprev - errPrev/jac
        print("ANGLE DEFORM RESULT: angle error:",err,"next force guess", Fnew)
        iiii+=1

    return torque, Fnew, springK
def deform_spring_byForce_bisect(alphaCoeffs, cICoeffs, F, Fmax, plotBool):
    sMax = fullArcLength
    dgdsLeft = 0
    dgdsRight = 0.01

    angLeft = alpha(sMax, alphaCoeffs)
    angRight = ang+0.5*deg2rad

    exitFlag = 0
    while exitFlag == 0:
        gamma0 = [0, dgdsLeft]
        ang = angLeft
        leftRes = fixed_rk4(deform_ODE, alphaCoeffs, cICoeffs, gamma0, F, ang, smesh)
        gamma0 = [0, dgdsRight]
        ang = angRight
        

    gamma0 = [0, dgds0]
    ang = alpha(sMax, alphaCoeffs)
    smesh = np.linspace(0, fullArcLength, globalLen)


    pass
def deform_spring_byForce(alphaCoeffs, cICoeffs, F, Fmax, plotBool):
    print("trying a force:", F)
    sMax = fullArcLength
    dgds0 = 0

    gamma0 = [0, dgds0]
    ang = alpha(sMax, alphaCoeffs)
    angOrig = ang
    smesh = np.linspace(0, fullArcLength, globalLen)
    print("integrating forward - first iteration")
    res = fixed_rk4(deform_ODE, alphaCoeffs, cICoeffs, gamma0, F, ang, smesh)
    print("done")
    gamma = res[0,:]

    Fx = F*-np.sin(ang)
    Fy = F*np.cos(ang)

    [rnXd, rnYd] = mesh_deformed_geometry(alphaCoeffs, gamma, smesh)
    [rnX, rnY]   = mesh_original_geometry(alphaCoeffs, smesh)
    
    if plotBool:
        # plt.figure(F)
        cmap = plt.get_cmap('viridis')
        color = cmap(float(0))
        plt.plot(rnX,rnY, c=color, linewidth=2.0)
    
    xdis = rnXd[-1]
    ydis = rnYd[-1]

    xorg = rnX[-1]
    yorg = rnY[-1]

    rdis   = np.sqrt(xdis**2+ydis**2)
    rorg = np.sqrt(xorg**2+yorg**2)
    # print("orig difference:",rdis,rorg)
    diff0 = rdis-rorg
    err = abs(diff0)

    dgds0Prev = dgds0
    dgds0 = 0.001 # is there a better way to estimate this initial step??

    angPrev = ang
    ang = ang+0.5*deg2rad

    errPrev = abs(diff0)
    gain = 0.005 # also fuck this

    iii = 1
    rprev = rorg
    restartFlag = 0
    while abs(err)>1.1e-12:
        gamma0 = np.array([0, dgds0])
        res = fixed_rk4(deform_ODE, alphaCoeffs, cICoeffs, gamma0, F, ang, smesh)
        gamma = res[0,:]
        [rnXd, rnYd] = mesh_deformed_geometry(alphaCoeffs, gamma, smesh)
        [rnX, rnY]   = mesh_original_geometry(alphaCoeffs, smesh)

        xdis = rnXd[-1]
        ydis = rnYd[-1]

        xorg = rnX[-1]
        yorg = rnY[-1]

        rdis   = np.sqrt(xdis**2+ydis**2)
        rorg = np.sqrt(xorg**2+yorg**2)
        rprev = rorg

        diff1 = rdis-rorg
        err  = abs(diff1)
        jacInv = np.array([(err-errPrev)/(dgds0-dgds0Prev), (err-errPrev)/(ang-angPrev)])
        jacInv = jacInv/lin.norm(jacInv)**2
        if lin.norm(jacInv) == 0:
            break

        Fx = F*-np.sin(ang)
        Fy = F*np.cos(ang)
        dgds0Prev = dgds0
        diff0 = diff1
        errPrev = abs(diff1)
        angPrev = ang
        # print(err*jac[0]*1/(lin.norm(jac)**2),err*jac[1]*1/(lin.norm(jac)**2))
        dgds0 = dgds0-err*jacInv[0] # *1/(lin.norm(jacInv)**2)
        ang   = ang  -err*jacInv[1] # *1/(lin.norm(jacInv)**2)
        if dgds0-dgds0Prev == 0:
            print("no change in dgds")
        if ang-angPrev == 0:
            print("no change in ang")
        iii +=1

        if iii>500:
            print("restarting!!!")
            restartFlag += 1
            gamma0 = [0, 0]
            ang = alpha(sMax, alphaCoeffs)
            smesh = np.linspace(0, fullArcLength, globalLen)
            print("integrating forward - first iteration")
            res = fixed_rk4(deform_ODE, alphaCoeffs, cICoeffs, gamma0, F, ang, smesh)
            print("done")
            gamma = res[0,:]

            Fx = F*-np.sin(ang)
            Fy = F*np.cos(ang)

            [rnXd, rnYd] = mesh_deformed_geometry(alphaCoeffs, gamma, smesh)
            [rnX, rnY]   = mesh_original_geometry(alphaCoeffs, smesh)
            
            if plotBool:
                # plt.figure(F)
                cmap = plt.get_cmap('viridis')
                color = cmap(float(0))
                plt.plot(rnX,rnY, c=color, linewidth=2.0)
            
            xdis = rnXd[-1]
            ydis = rnYd[-1]

            xorg = rnX[-1]
            yorg = rnY[-1]

            rdis   = np.sqrt(xdis**2+ydis**2)
            rorg = np.sqrt(xorg**2+yorg**2)
            # print("orig difference:",rdis,rorg)
            diff0 = rdis-rorg
            err = abs(diff0)

            dgds0Prev = dgds0
            dgds0 = 0.1*restartFlag # is there a better way to estimate this initial step??

            angPrev = ang
            ang = ang+1*deg2rad*restartFlag

            errPrev = abs(diff0)
            gain = 0.005 # also fuck this

            iii = 1
            rprev = rorg

    """ fig = plt.figure()
    ax = plt.axes(projection='3d')
    for i in range(len(dGdS0_plot)):
        cmap = plt.get_cmap('viridis')
        color = cmap(float(i/len(dGdS0_plot)))
        ax.plot3D(dGdS0_plot[i:i+5+1],ang_plot[i:i+5+1],z_plot[i:i+5+1], c=color)
    plt.show()

    # plt.figure(F)
    plt.plot(res.t, rorg_vec)
    print(res.t)
    print(rorg_vec) """
    print(iii)
    if plotBool:
        cmap = plt.get_cmap('viridis')
        color = cmap(abs(float(F))/float(Fmax))
        plt.plot(rnXd,rnYd, c=color, linewidth=7.0)

    [rnXd, rnYd] = mesh_deformed_geometry(alphaCoeffs, gamma, smesh)
    [rnX, rnY]   = mesh_original_geometry(alphaCoeffs, smesh)

    xdis = rnXd[-1]
    ydis = rnYd[-1]

    xorg = rnX[-1]
    yorg = rnY[-1]

    rdis   = np.sqrt(xdis**2+ydis**2)
    rorg = np.sqrt(xorg**2+yorg**2)

    dBeta = -(np.arctan2(ydis,xdis)-np.arctan2(yorg,xorg))
    Torque = 2*(xdis*Fy-ydis*Fx+E*cI_s(sMax, cICoeffs)*res[1,-1])
    print("FORCE DEFORM RESULT:", dBeta*1/deg2rad, "derivative IC", dgds0, "angle range",ang-angOrig)
    if not dBeta == 0:
        springK = Torque/dBeta*deg2rad
    else: 
        springK = 'unknown'
    if plotBool:
        theta = np.linspace(0, 2*np.pi, 100)
        outerCircleX = rorg*np.cos(theta)
        outerCircleY = rorg*np.sin(theta)
        plt.plot(outerCircleX,outerCircleY)

    finalGuess = {"angle": ang, "dgds0": dgds0}
    
    return dBeta, Torque, springK, gamma, iii, finalGuess

def eval_stiffness(alphaCoeffs, cICoeffs, Fmax, Fres, plotBool):

    Fres = 5
    Fmax = 50000
    Fvec = np.linspace(-Fmax, Fmax, 2*Fres+1)

    dBetaPlot = np.empty(len(Fvec))
    TPlot =     np.empty(len(Fvec))
    kPlot =     np.empty(len(Fvec))

    ploti = 0
    recalli = len(Fvec)+2
    for F in Fvec:
        [dBeta, Torque, springK, gamma, iii, finalGuess]  = deform_spring_byForce(alphaCoeffs, cICoeffs, F, Fmax, plotBool)
        # print("iterations:", iii, "dBeta:",dBeta*1/deg2rad,"degrees","Torque:",Torque)
        dBetaPlot[ploti] = dBeta*1/deg2rad
        TPlot[ploti] =Torque
        if springK == 'unknown':
            kPlot[ploti] = kPlot[ploti-1]
            recalli = ploti
        else:
            kPlot[ploti] = springK
        if ploti == recalli+1:
            kPlot[recalli] = (kPlot[ploti-2]+kPlot[ploti])/2
        ploti+=1
        # print("final angle:", finalGuess["angle"], "final dgds0", finalGuess["dgds0"])
    approxK = abs((TPlot[-1]-TPlot[0])/(dBetaPlot[-1]-dBetaPlot[0]))
    print("overall K:", approxK)
    

    if plotBool:
        plt.figure('defl vs T')
        plt.plot(TPlot,dBetaPlot)
        plt.xlabel('Torque (lb-in)')
        plt.ylabel('displacement (deg)')
        plt.figure('k vs defl')
        plt.plot(dBetaPlot,kPlot)
        plt.show()

    return approxK

def tune_stiffness(goalK, Feval, **kwargs):

    def stiffnessObjectiveFun(dragVector, goalK, Feval):
        [alphas, cIs, si, sj, sk ,sf] = generate_params(dragVector)
        [alphaCoeffs, cICoeffs, smesh] = create_initial_spring(alphas, cIs, si, sj, sk, sf)

        currentK = eval_stiffness(alphaCoeffs, cICoeffs, Feval, 5, False)
        dist = currentK-goalK
        print("dist:",dist)
        print("params:",dragVector)
        return currentK-goalK


    def generate_params(dragVector):

        alphas = [dragVector[0],dragVector[1],dragVector[2],dragVector[3]]
        cIs    = [dragVector[6],dragVector[7],dragVector[8]]
        si     = dragVector[10]*dragVector[4]
        sj     = dragVector[10]*dragVector[5]
        sk     = dragVector[10]*dragVector[9]
        sf     = dragVector[10]

        return alphas, cIs, si, sj, sk, sf
    
   
    dragVectorParams = {'alpha0':       0,            \
                        'alpha1':       100*deg2rad,  \
                        'alpha2':       190*deg2rad,  \
                        'alpha3':       135*deg2rad,  \
                        'alpha1factor': 1/3.0,        \
                        'alpha2factor': 2/3.0,        \
                        'cI0':          .1,           \
                        'cI1':          .05,          \
                        'cI2':          .1,           \
                        'cIfactor':     0.5,          \
                        'fullLength':   6             }
    
    for key, value in kwargs.items():
        dragVectorParams[key] = value

    dragVector = np.array(list(dragVectorParams.values()))

    angleRange = 30*deg2rad
    propRange = 1/5.0
    lenRange  = 0.375

    bnds = ((0,0), (dragVector[1]-angleRange,dragVector[1]+angleRange), \
                   (dragVector[2]-angleRange,dragVector[2]+angleRange),  \
                   (dragVector[3]-angleRange,dragVector[3]+angleRange),  \
                   (dragVector[4]-propRange,dragVector[4]+propRange),  \
                   (dragVector[5]-propRange,dragVector[5]+propRange),  \
                   (0,None),  \
                   (0,None),  \
                   (0,None),  \
                   (dragVector[9]-propRange,dragVector[9]+propRange),  \
                   (dragVector[10]-lenRange,dragVector[10]+lenRange))

    res = op.minimize(stiffnessObjectiveFun, dragVector, args=(goalK,Feval), bounds=bnds, method='Nelder-Mead')
    # err = stiffnessObjectiveFun(dragVector, goalK, Feval)
    # err0 = err
    # dragVector0 = dragVector
    # dragVector = dragVector*1.01
    # dragVector[1] = dragVector0[1]
    # dragVector[2] = dragVector0[2]
    # dragVector[3] = dragVector0[3]
    # jac = np.empty(len(dragVector0))
    # gain = 1.5
    # while err>1:
    #     err = stiffnessObjectiveFun(dragVector, goalK, Feval)
    #     for i in range(len(dragVector0)):
    #         if i == 0 or i == 1 or i == 2 or i == 3:
    #             jac[i] = 1
    #         else:
    #             jac[i] = (err-err0)/(dragVector[i]-dragVector0[i])
    #             dragVector0[i] = dragVector[i]
    #         dragVector[i] = dragVector[i] - gain*err * jac[i]/(lin.norm(jac)**2)
        
    #     err0 = err

    # finalParameters = dragVector

    if res.success:
        print("good shit")
    finalParameters = res.x
    
    [alphas, cIs, si, sj, sk ,sf] = generate_params(finalParameters)
    [alphaCoeffs, cICoeffs, smesh] = create_initial_spring(alphas, cIs, si, sj, sk, sf)
    finalK = eval_stiffness(alphaCoeffs, cICoeffs, Feval, 5, True)
    print("final stiffness:",finalK)

    return finalParameters

deg2rad = np.pi/180
# E = 27.5*10**6
E = 1000000

fullArcLength = 6.2
globalRes = 150
globalLen = globalRes+1
globalStep = fullArcLength/globalRes
globalMaxIndex = globalLen-1

alphas  = [0,100*deg2rad,190*deg2rad,135*deg2rad]
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
#  [alph, curvature, cI, rn, ab, h] = create_profiles(alphaCoeffs, cICoeffs, smesh, True)
dBetaTarg = 5*deg2rad
Fguess    = 1000
[torque, Fnew, springK] = deform_spring_byAngle(alphaCoeffs, cICoeffs, dBetaTarg, Fguess)
print("torque",torque,"force:",Fnew,"K:",springK)
# approxK = eval_stiffness(alphaCoeffs, cICoeffs, 5000, 10, True)
# print(approxK)
create_profiles(alphaCoeffs, cICoeffs, smesh, True)

#  finalParameters = tune_stiffness(5000,3000)





