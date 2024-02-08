from mpl_toolkits import mplot3d
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from materials import *
from openpyxl import *
from scipy.integrate import solve_ivp as ODE45
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

def deform_ODE(s, gamma, F, ang):
    
    Fx = F*-np.sin(ang)
    Fy = F*np.cos(ang)
    LHS = np.empty(2)
    LHS[0] = gamma[1]
    LHS[1] = (Fx*np.sin(alpha(s, alphaCoeffs)+gamma[0])-Fy*np.cos(alpha(s, alphaCoeffs)+gamma[0])-E*d_cI_d_s(s, cICoeffs)*gamma[1])/(cI_s(s, cICoeffs)*E)
    return LHS

def rk4_step(fun, u0, gamma0, du, F, ang):
    #
    #  Get four sample values of the derivative.
    #
    f1 = fun ( u0,            gamma0, F, ang )
    f2 = fun ( u0 + du / 2.0, gamma0 + du * f1 / 2.0, F, ang )
    f3 = fun ( u0 + du / 2.0, gamma0 + du * f2 / 2.0, F, ang )
    f4 = fun ( u0 + du,       gamma0 + du * f3, F, ang )
    #
    #  Combine them to estimate the solution gamma at time T1 = T0 + DT.
    #
    gamma1 = gamma0 + du * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0

    return gamma1

def fixed_rk4(fun, gamma0, F, ang, smesh):
    res = np.empty((len(gamma0), len(smesh)))
    step = smesh[1]-smesh[0]
    for i in range(len(smesh)):
        if i == 0:
            res[0,i] = gamma0[0]
            res[1,i] = gamma0[1]
        else:
            stepRes = rk4_step(fun, smesh[i-1], gamma0, step, F, ang)
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

def deform_spring(alphaCoeffs, cICoeffs, F, plotBool):

    sMax = fullArcLength
    dcIds0 = 0

    # err = 1

    # dcIds_plot=[]
    # ang_plot=[]
    # z_plot=[]

    gamma0 = [0, 0]
    ang = alpha(sMax, alphaCoeffs)
    smesh = np.linspace(0, fullArcLength, globalLen)

    res = fixed_rk4(deform_ODE, gamma0, F, ang, smesh)

    gamma = res[0,:]
        # res = ODE45(deform_ODE, [0, sMax], gamma0, args=(F,ang), method='RK45', max_step=globalStep, atol=1, rtol=1)

    Fx = F*-np.sin(ang)
    Fy = F*np.cos(ang)

    [rnXd, rnYd] = mesh_deformed_geometry(alphaCoeffs, gamma, smesh)
    [rnX, rnY]   = mesh_original_geometry(alphaCoeffs, smesh)
    
    if plotBool:
        plt.figure(F)
        cmap = plt.get_cmap('viridis')
        color = cmap(float(0))
        plt.plot(rnXd,rnYd, c=color)
    
    xdis = rnXd[-1]
    ydis = rnYd[-1]

    xorg = rnX[-1]
    yorg = rnY[-1]

    rdis   = np.sqrt(xdis**2+ydis**2)
    rorg = np.sqrt(xorg**2+yorg**2)
    # print("orig difference:",rdis,rorg)
    diff0 = rdis-rorg
    err = abs(diff0)

    # dcIds_plot.append(dcIds0)
    # ang_plot.append(ang)
    # z_plot.append(abs(diff0))

    dcIdsPrev = dcIds0
    dcIds0 = 0.001 # is there a better way to estimate this initial step??

    angPrev = ang
    ang = ang+0.5*deg2rad

    errPrev = abs(diff0)
    gain = 0.005 # also fuck this

    iii = 1
    rprev = rorg
    while abs(err)>1.1e-12:
        # use the current shooting variables (dGdS and ang) to solve
        gamma0 = np.array([0, dcIds0])
    
        res = fixed_rk4(deform_ODE, gamma0, F, ang, smesh)

        gamma = res[0,:]
        [rnXd, rnYd] = mesh_deformed_geometry(alphaCoeffs, gamma, smesh)
        [rnX, rnY]   = mesh_original_geometry(alphaCoeffs, smesh)
        # print(res.t[-1])
        
        # if plotBool:
        #     color = cmap(float(iii)/float(27))
        #     plt.plot(rnXd,rnYd, c=color, label=iii)

        xdis = rnXd[-1]
        ydis = rnYd[-1]

        xorg = rnX[-1]
        yorg = rnY[-1]

        rdis   = np.sqrt(xdis**2+ydis**2)
        rorg = np.sqrt(xorg**2+yorg**2)
        
        # if rorg == rprev:
        #     print("difference:",rdis,rorg)

        # stepPrev = step
        rprev = rorg

        diff1 = rdis-rorg
        err  = abs(diff1)
        # dcIds_plot.append(dcIds0)
        # ang_plot.append(ang)
        # z_plot.append(err)

        # estimate jacobian from prev. 2 guesses
        jac = np.array([(err-errPrev)/(dcIds0-dcIdsPrev), (err-errPrev)/(ang-angPrev)])
        # these two lines probably no longer necessary
        if lin.norm(jac) == 0:
            break

        Fx = F*-np.sin(ang)
        Fy = F*np.cos(ang)

        # save guess for next time
        dcIdsPrev = dcIds0
        diff0 = diff1
        errPrev = abs(diff1)
        angPrev = ang

        #shitty newtons method to determine next guess

        dcIds0 = dcIds0-err*jac[0]*1/(lin.norm(jac)**2)
        ang   = ang  -err*jac[1]*1/(lin.norm(jac)**2)
        iii +=1
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # for i in range(len(dGdS0_plot)):
    #     cmap = plt.get_cmap('viridis')
    #     color = cmap(float(i/len(dGdS0_plot)))
    #     ax.plot3D(dGdS0_plot[i:i+5+1],ang_plot[i:i+5+1],z_plot[i:i+5+1], c=color)
    # plt.show()

    # # plt.figure(F)
    # plt.plot(res.t, rorg_vec)
    # print(res.t)
    # print(rorg_vec)
        

    if plotBool:
        cmap = plt.get_cmap('viridis')
        color = cmap(float(1))
        plt.plot(rnXd,rnYd, c=color)

    dBeta = -(np.arctan2(ydis,xdis)-np.arctan2(yorg,xorg))
    Torque = 2*(xdis*Fy-ydis*Fx+E*cI_s(sMax, cICoeffs)*res[1,-1])
    if not dBeta == 0:
        springK = Torque/dBeta*deg2rad
    else: 
        springK = 'unknown'
    if plotBool:
        theta = np.linspace(0, 2*np.pi, 100)
        outerCircleX = rorg*np.cos(theta)
        outerCircleY = rorg*np.sin(theta)
        plt.plot(outerCircleX,outerCircleY)
    
    return dBeta, Torque, springK, gamma, iii

def eval_stiffness(alphaCoeffs, cICoeffs, Fmax, Fres, plotBool):

    Fres = 5
    Fmax = 5000
    Fvec = np.linspace(-Fmax, Fmax, 2*Fres+1)

    dBetaPlot = np.empty(len(Fvec))
    TPlot =     np.empty(len(Fvec))
    kPlot =     np.empty(len(Fvec))

    ploti = 0
    recalli = len(Fvec)+2
    for F in Fvec:
        [dBeta, Torque, springK, gamma, iii]  = deform_spring(alphaCoeffs, cICoeffs, F, plotBool)
        print("iterations:", iii, "dBeta:",dBeta*1/deg2rad,"degrees","Torque:",Torque,"spring k:",springK)
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
    approxK = (TPlot[-1]-TPlot[0])/(dBetaPlot[-1]-dBetaPlot[0])

    if plotBool:
        plt.figure('defl vs T')
        plt.plot(TPlot,dBetaPlot)
        plt.xlabel('Torque (lb-in)')
        plt.ylabel('displacement (deg)')
        plt.figure('k vs defl')
        plt.plot(dBetaPlot,kPlot)
        plt.show()

    return approxK

def optimize_stiffness(goalK, s_events_initial, ):


deg2rad = np.pi/180
E = 27.5*10**6

fullArcLength = 6.2
globalRes = 150
globalLen = globalRes+1
globalStep = fullArcLength/globalRes
globalMaxIndex = globalLen-1

alphas = [0,100*deg2rad,190*deg2rad,135*deg2rad]
cIs     = [.1,.05,.1]
# alphas = [0,45*deg2rad,135*deg2rad,180*deg2rad]
# cIs     = [.05,.05,.05]

xy0    = [1,0]

si = fullArcLength/3.0
sj = 2.0*fullArcLength/3.0
sk = fullArcLength/2.0
sf = fullArcLength

outPlaneThickness = 0.375

# [alph, curvature, cI, rn, ab, h] = create_profiles(alphaCoeffs, cICoeffs, smesh, True)
[alphaCoeffs, cICoeffs, smesh] = create_initial_spring(alphas, cIs, si, sj, sk, sf)
approxK = eval_stiffness(alphaCoeffs, cICoeffs, 5000, 10, True)
print(approxK)





