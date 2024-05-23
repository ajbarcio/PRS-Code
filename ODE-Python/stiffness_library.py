from mpl_toolkits import mplot3d
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp as ODE45
from scipy import optimize as op
import os
import random
import math
from copy import deepcopy as dc
from StatProfiler import SSProfile
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import sympy as sp

import subprocess as sb
import win32com.client
import pythoncom

deg2rad = np.pi/180

def make_SW_part(Spring):
    startSW()
    sw = connectToSW()
    newFile(sw, "newPart.SLDPRT")

def Ic_multiPoly(IcPts, IcArcLens):

    # Determine how many parabolas you need
    numPolys = (len(IcPts)-1)*2
    numSegments = int(numPolys/2)

    # Determine the x, y locations of your control points
    dys = np.empty(numPolys)
    dxs = np.empty(numPolys)

    for i in range((numSegments)):
        dys[i] = (IcPts[i+1]-IcPts[i])/2
        dxs[i] = (IcArcLens[i+1]-IcArcLens[i])/2

    ctrlX = np.empty(numPolys+1)
    ctrlY = np.empty(numPolys+1)

    for i in range(numSegments):
        ctrlX[2*i]   = IcArcLens[i]
        ctrlX[2*i+1] = IcArcLens[i]+dxs[i]
        ctrlY[2*i]   = IcPts[i]
        ctrlY[2*i+1] = IcPts[i]+dys[i]
    ctrlX[-1] = IcArcLens[-1]
    ctrlY[-1] = IcPts[-1]

    # Solve for coeffs of each parabola
    coeffs = np.empty((numPolys,3))
    for i in range(numPolys):
        if not i % 2:
            Targ = np.array([[ctrlY[i]],[ctrlY[i+1]],[0]])
            Mat  = np.array([[ctrlX[i]**2,ctrlX[i],1],
                            [ctrlX[i+1]**2,ctrlX[i+1],1],
                            [2*ctrlX[i],1,0]])
        else:
            Targ = np.array([[ctrlY[i]],[ctrlY[i+1]],[0]])
            Mat  = np.array([[ctrlX[i]**2,ctrlX[i],1],
                            [ctrlX[i+1]**2,ctrlX[i+1],1],
                            [2*ctrlX[i+1],1,0]]) ## right hand point is
        # print(Targ)
        # print(Mat)
        # print(lin.solve(Mat,Targ))
        coeffs[i,:] = lin.solve(Mat,Targ).T

    return coeffs, ctrlX
    # for i in range(numPolys):
    #     Targ = np.array([[IcPts[i]],[],[]])

def Ic_spline(IcPts, IcArcLens): # Deprecated
    inputPts = np.hstack((np.atleast_2d(IcArcLens).T,np.atleast_2d(IcPts).T))
    dy1 = (inputPts[0,1]-inputPts[1,1])/2
    dy2 = (inputPts[2,1]-inputPts[1,1])/2

    # dx1 = (inputPts[1,0]-inputPts[0,0])/4
    # dx2 = (inputPts[2,0]-inputPts[1,0])/4
    dx1   = np.sqrt(dy1**2+((inputPts[1,0]-inputPts[0,0])/4)**2)
    dx2   = np.sqrt(dy2**2+((inputPts[2,0]-inputPts[1,0])/4)**2)
    print(dx1)
    P0  = inputPts[0,:]
    P4  = inputPts[1,:]
    P8  = inputPts[2,:]

    P1  = P0+np.array([dx1,0])
    P3  = P4-np.array([dx1,0])
    P5  = P4+np.array([dx2,0])
    P7  = P8-np.array([dx2,0])
    P2  = P0+np.array([2*dx1,-dy1])
    P6  = P4+np.array([2*dx2, dy2])

    ctrlPts = np.vstack((P0,P1,P2,P3,P4,P5,P6,P7,P8))
    splPts  = np.array([P0[0],P2[0],P4[0],P6[0],P8[0]])
    print(ctrlPts)
    numSplines = len(IcPts)+1
    degSpline  = 2
    solns = [None]*numSplines
    diffs = [None]*numSplines
    print(len(solns))
    for i in range(numSplines):
        print(i)
        index = i*(degSpline)
        pts = ctrlPts[index:index+degSpline+1,:]
        q   = np.array([[pts[0,0],pts[1,0],pts[2,0]],[pts[0,1],pts[1,1],pts[2,1]],[1,1,1]])

        print("q:", q)
        Q_bar = (lin.det(q)*lin.inv(q)).T
        Q_bar = sp.nsimplify(Q_bar,tolerance=1e-15,rational=True)
        # Q_bar = sp.simplify(Q_bar)
        # Q_bar[0,0] = 0
        u = np.atleast_2d(Q_bar[:,0])
        v = np.atleast_2d(Q_bar[:,1])
        w = np.atleast_2d(Q_bar[:,2])
        Q = 2*(u*w.T+w*u.T)-v*v.T
        # Q = sp.nsimplify(Q,tolerance=1e-10,rational=True)
        # Q = sp.simplify(Q)
        print("Q:", Q)
        print("Q_bar:", Q_bar)
        x = sp.Symbol('x')
        y = sp.Symbol('y')
        eqn = np.array([x,y,1]).dot(Q).dot(np.array([x,y,1]))
        print(eqn)
        soln = sp.solvers.solve(eqn,y)
        solns[i] = soln[0]
        diffs[i] = sp.diff(soln[0],x)

    return solns, diffs, splPts

def Ic_poly(IcPts, IcArcLens): # Deprecated

    # for the Ic curve, we want:
    # flat at either end, and at the minimum
    # So: for n points we need:
    #     n at y(s)
    #     n at y'(s)
    #     n at y''(s)

    # so for n points we need nDerivs*n constraints
    # CURRENTLY ONLY USING FIRST DERIVATIVE CONSTRAINT
    # DEGREE = 2*3-1 = 5th degree polynomial, +1 6th degree polynomial (7 coefficients)

    nDerivs = 2
    nConstraints = nDerivs*len(IcPts) # degree + 1
    extraConstraint = False
    # Always keep it an even function???
    if not np.mod(nConstraints,2):
        nConstraints+=2
        extraConstraint=True

    Mat = np.empty((nConstraints,nConstraints))
    Targ = np.empty((nConstraints,1))
    for i in range(len(IcArcLens)):
        Targ[nDerivs*i]=IcPts[i]
        Targ[nDerivs*i+1]=0

        # Targ[nDerivs*i+2]=0# (-.001)**(i+1)
        # Targ[nDerivs*i+3]=0# (-.001)
        for j in range(nConstraints):
            Mat[nDerivs*i,j]   = IcArcLens[i]**(nConstraints-1-j)
            Mat[nDerivs*i+1,j] = (nConstraints-1-j)*IcArcLens[i]**max(nConstraints-2-j,0) ### FIRST DERIVATIVE
            # Mat[nDerivs*i+2,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)*IcArcLens[i]**max(nConstraints-3-j,0)
            # Mat[nDerivs*i+3,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)*max(nConstraints-3-j,0)*IcArcLens[i]**max(nConstraints-4-j,0)
    if extraConstraint:
        Targ[-1] = 0 # but how to tell what we're constraining?
        Targ[-2] = 0
        for j in range(nConstraints):
            Mat[-1,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)*max(nConstraints-3-j,0)*IcArcLens[0]**max(nConstraints-4-j,0) # fuck it lets just try the third derivativeone shall we?
            Mat[-2,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)*max(nConstraints-3-j,0)*IcArcLens[-1]**max(nConstraints-4-j,0) # on the first and last one shall we?
    # print(Targ)
    cICoeffs = lin.solve(Mat, Targ)
    return cICoeffs

def xy_poly(pts, XYArcLens):

    # if we have n points we need:
    # n constraints on first derivative
    # 2 constraints on first derivative (angle at start and end)
    # 2 constraints on second derivative (straightness at start and end)
    # n + 4 constraints

    # initialize matrix sizes
    nConstraints = len(pts)+4
    Mat = np.empty((nConstraints,nConstraints))
    YTarg = np.empty((nConstraints,1))
    XTarg = np.empty((nConstraints,1))
    # FIRST ASSIGN 0th DERIVATIVE CONSTRAINTS (c(s)=p)
    for i in range(len(XYArcLens)):
        # target correct x-y value
        XTarg[i]=pts[i,0]
        YTarg[i]=pts[i,1]
        # at associated s value

        for j in range(nConstraints):
            Mat[i,j]   = XYArcLens[i]**(nConstraints-1-j)
    # NOW ASSIGN FIRST AND SECOND DERIVATIVE CONSTRAINTS
    for i in range(-1,-5,-1):
        # target correct index (first or last)
        index = (i % 2)*(len(XYArcLens)-1)
        # print(index)
        if i < -2:
            # with correct first derivative (radial direction)
            XTarg[i] = pts[index,0]/lin.norm(pts[index,:])
            YTarg[i] = pts[index,1]/lin.norm(pts[index,:])
            # at correct s value
            for j in range(nConstraints):
                Mat[i,j] = (nConstraints-1-j)*XYArcLens[index]**max(nConstraints-2-j,0) ## FIRST TWO ARE FIRST DERIVATIVES
        else:
            # and with zero second derivative at both points
            XTarg[i] = 0
            YTarg[i] = 0
            for j in range(nConstraints):
                Mat[i,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)*XYArcLens[index]**max(nConstraints-3-j,0) ## LAST TWO ARE SECOND DERIVATIVES
    # print(XTarg)
    # print(YTarg)
    # print(Mat)
    XCoeffs = lin.solve(Mat, XTarg)
    YCoeffs = lin.solve(Mat, YTarg)
    # print(XCoeffs, YCoeffs)
    return XCoeffs, YCoeffs

def PPoly_Eval(x, coeffs, deriv=0, ranges=0):
    # FIRST DECISION: Are there multiple polynomials?
    # print(ranges)
    if hasattr(ranges, "__len__"):
        # SECOND DECISION: Are there multiple s points to evaluate?
        if hasattr(x, "__len__"):
            y = np.empty(len(x))
            i = 0
            for value in x:
                U = np.empty(coeffs.shape[1])
                for j in range(coeffs.shape[1]):
                    preCoeff = 1
                    for k in range(deriv):
                        preCoeff = preCoeff*max(coeffs.shape[1]-j-(k+1),0)
                    U[j] = preCoeff*value**max((coeffs.shape[1]-j-1-deriv),0)
                for l in range(len(ranges)-1):
                    if value >= ranges[l] and value < ranges[l+1]:
                        index = l
                        break
                    index = len(ranges)-2
                y[i] = U.dot(coeffs[index,:])
                i+=1
            return (y)
        else:
            U = np.empty(coeffs.shape[1])
            for j in range(coeffs.shape[1]):
                preCoeff = 1
                for k in range(deriv):
                    preCoeff = preCoeff*max(coeffs.shape[1]-j-(k+1),0)
                U[j] = preCoeff*x**max((coeffs.shape[1]-j-1-deriv),0)
            for l in range(len(ranges)-1):
                if x >= ranges[l] and x < ranges[l+1]:
                    index = l
                    break
                index = len(ranges)-2
            y = U.dot(coeffs[index,:])
            return (y)
    else:
        if hasattr(x, "__len__"):
            y = np.empty(len(x))
            i = 0
            for value in x:
                U = np.empty(len(coeffs))
                for j in range(len(coeffs)):
                    preCoeff = 1
                    for k in range(deriv):
                        preCoeff = preCoeff*max(len(coeffs)-j-(k+1),0)
                    U[j] = preCoeff*value**max((len(coeffs)-j-1-deriv),0)
                y[i] = U.dot(coeffs)
                i+=1
            return (y)
        else:
            U = np.empty(len(coeffs))
            for j in range(len(coeffs)):
                preCoeff = 1
                for k in range(deriv):
                    preCoeff = preCoeff*max(len(coeffs)-j-(k+1),0)
                U[j] = preCoeff*x**max((len(coeffs)-j-1-deriv),0)
            y = U.dot(coeffs)
            return (y[0])

def d_xi_d_s(ximesh, XCoeffs, YCoeffs):
    dxdxi = PPoly_Eval(ximesh, XCoeffs, deriv=1)
    dydxi = PPoly_Eval(ximesh, YCoeffs, deriv=1)
    dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
    return dxids

class Deform_Wrapper:
    def __init__(self, function):
        self.all_output = []
        self.all_moments = []
        self.function = function

    def __call__(self, torqueTarg, ODE, torqueResolution=10, SF = np.array([0,0,0]), breakBool=False):
        self.all_moments.append(torqueTarg)
        SSProfile("BVP").tic()
        res, SF, divergeFlag, i  = self.function(torqueTarg, ODE, torqueResolution, SF, breakBool)
        SSProfile("BVP").toc()
        self.all_output.append(SF)
        return res, SF, divergeFlag, i

class Spring:
    def __init__(self,
                 # whole bunch of default values:
                 E                 = 27.5*10**6,
                 designStress      = 270000,
                 n                 = 2,
                 fullArcLength     = 4.8,
                 outPlaneThickness = 0.375,
                 radii             = np.array([1.1,2.025,2.025,2.95]),
                 betaAngles        = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
                 IcPts             = np.array([0.008, 0.001, 0.008]),
                 IcArcLens         = np.array([0.5]),                              # IcPts - 2
                 XYArcLens         = np.array([0.333,0.667]),             # radii - 2
                 resolution        = 200
                 ):

        # stick all the arguments in the object
        self.E = E
        self.designStress = designStress
        self.n = n
        self.t = outPlaneThickness
        self.fullArcLength = fullArcLength
        self.radii = radii
        self.betaAngles = betaAngles

        # create control points for polynomials
        self.pts = np.empty((len(radii),2))
        for i in range(len(radii)):
            self.pts[i,:] = [self.radii[i]*np.cos(self.betaAngles[i]),self.radii[i]*np.sin(self.betaAngles[i])]

        self.IcPts = IcPts

        self.IcArcLens = np.empty(len(IcPts))
        for i in range(len(IcPts)):
            if i==0:
                self.IcArcLens[i] = 0
            elif i==len(IcPts)-1:
                self.IcArcLens[i] = self.fullArcLength
            else:
                self.IcArcLens[i] = self.fullArcLength*IcArcLens[i-1]
        # print(self.IcArcLens)

        self.XYArcLens = np.empty(len(radii))
        for i in range(len(radii)):
            if i==0:
                self.XYArcLens[i] = 0
            elif i==len(radii)-1:
                self.XYArcLens[i] = self.fullArcLength
            else:
                self.XYArcLens[i] = self.fullArcLength*XYArcLens[i-1]
        # print(self.XYArcLens)

        self.parameterVector = np.concatenate((self.radii, self.betaAngles[1:],
                                              self.IcPts, self.IcArcLens[1:-1],
                                              self.XYArcLens[1:-1], np.atleast_1d(self.fullArcLength)))

        # assign some constants for approximating derivatives
        self.finiteDifferenceLength = 0.001
        self.finiteDifferenceAngle  = .1*deg2rad
        self.finiteDifferenceForce  = 0.1
        self.finiteDifferenceTorque = 0.5
        self.finiteDifferenceCI     = 20

        # assign some constants for meshing the spring
        self.res = resolution
        self.len = self.res+1
        self.step = fullArcLength/self.res
        self.endIndex = self.len-1

        # create uniform mesh for spring
        self.ximesh = np.linspace(0,self.fullArcLength,self.len)

        self.geometry_coeffs()
        self.dxids = d_xi_d_s(self.ximesh, self.XCoeffs, self.YCoeffs)
        # this takes way too many lines but fuck it IDC, better to see the math
        # Find coordinates of frames:
        # at base of spring
        self.x0 = PPoly_Eval(0,self.XCoeffs)
        self.y0 = PPoly_Eval(0,self.YCoeffs)
        # at end of spring
        self.xL = PPoly_Eval(self.fullArcLength,self.XCoeffs)
        self.yL = PPoly_Eval(self.fullArcLength,self.YCoeffs)
        # at end of spring, interms of frame at base of spring
        self.momentArmX = self.xL - self.x0
        self.momentArmY = self.yL - self.y0

        # inherit deform wrapper
        self.wrapped_torque_deform = Deform_Wrapper(self.deform_by_torque)

    def geometry_coeffs(self):
        # sticks the coefficients in the object
        self.IcCoeffs, self.domains = Ic_multiPoly(self.IcPts, self.IcArcLens)
        self.XCoeffs, self.YCoeffs  = xy_poly(self.pts, self.XYArcLens)
        # compiles all the coeffs into one variable
        self.geometryDef = [self.XCoeffs, self.YCoeffs, self.IcCoeffs]
        return self.geometryDef

    def smart_initial_load_guess(self, torqueTarg, ODE):
        # this method makes an informed guess about the forcing at full
        # loading condition based off the magnitude of the torque
        print("getting initial guess")
        guessWrapper = Deform_Wrapper(self.deform_by_torque)

        SF = np.array([0,0,torqueTarg*0.05])
        guessWrapper(SF[2],ODE,SF=SF)
        SF = np.array([0,0,torqueTarg*0.15])
        guessWrapper(SF[2],ODE,SF=SF)

        torques = np.array([torqueTarg*0.1,torqueTarg*0.2])
        angles = np.empty(len(guessWrapper.all_output))
        magnitudes = np.empty(len(guessWrapper.all_output))
        for i in range(len(guessWrapper.all_output)):
            angles[i] = np.arctan2(guessWrapper.all_output[i][1],guessWrapper.all_output[i][0])
            magnitudes[i] = lin.norm(guessWrapper.all_output[i][0:2])
        anglePoly = np.polyfit(torques,angles,1)
        magPoly = np.polyfit(torques,magnitudes,1)
        angleGuess  = np.polyval(anglePoly, torqueTarg)
        magGuess  = np.polyval(magPoly, torqueTarg)
        SFGuess = np.array([magGuess*np.cos(angleGuess), magGuess*np.sin(angleGuess), torqueTarg])
        print("done getting initial guess")
        return(SFGuess)

    def deform_by_torque_slowRamp(self, torqueTarg, ODE, torqueResolution=15):
        self.rampWrapper = Deform_Wrapper(self.deform_by_torque)
        # this method goes straight to full loading condition lmao
        # TODO: MAKE THIS AUTOMATIC
        SFStart = np.array([0,0,0])
        # the err is a two vector, so make it arbitrarily high to enter loop
        err = np.ones(2)*float('inf')
        # 0th iteration
        i = 0
        j = 0
        k = 0
        divergeFlag = 0
        # limit to 100 iterations to converge
        n = len(err)
        stepTorque = 0
        while stepTorque <= torqueTarg:
            print(stepTorque, j, i)
            SFStart[2] = stepTorque
            res, SF, divergeFlag, i = self.rampWrapper(stepTorque,ODE,SF=SFStart,breakBool=True)
            if divergeFlag:
                if stepTorque:
                    stepTorque-=torqueTarg/torqueResolution
                torqueResolution *= 2
            else:
                stepTorque+=torqueTarg/torqueResolution
                SFStart = SF
                SFStart[2] += torqueTarg/torqueResolution
            j+=1
            k+=i
        return res, SF, divergeFlag, k

    def deform_by_torque(self,torqueTarg,ODE,torqueResolution=10,SF=np.array([0,0,0]),breakBool=False):
        # this method goes straight to full loading condition lmao

        # the err is a two vector, so make it arbitrarily high to enter loop
        err = np.ones(2)*float('inf')
        # 0th iteration
        i = 0
        divergeFlag = 0
        # limit to 100 iterations to converge
        while i <100:
            errPrev = err
            # determine boundary condition compliance, estimate jacobian
            err, self.res = self.forward_integration(ODE,SF,torqueTarg)
            J = self.estimate_jacobian_fd(ODE,SF,torqueTarg)
            # freak the fuck out if it didnt work
            # print(err)
            if lin.norm(err)>lin.norm(errPrev):
                print("torque deform diverging", i)
                # print(err, errPrev)
                # print(SF)
                print(J)
                divergeFlag = 1
                if breakBool:
                    break
            # make a new guess if it did
            elif lin.norm(err,2) > 10e-6:
                SF = SF-lin.inv(J).dot(err)
                # print(lin.inv(J))
            else:
                break
            i+=1
        self.solnSF = SF
        return self.res, self.solnSF, divergeFlag, i

    def forward_integration(self, ODE, SF, torqueTarg):
        # set up initial condition
        dgds0 = (SF[2]/self.n + self.momentArmY*SF[0] - self.momentArmX*SF[1])/(self.E*PPoly_Eval(0, self.IcCoeffs, ranges=self.domains))
        # print(dgds0)
        # print(self.IcCoeffs)
        # print("Ic(0):", PPoly_Eval(0, self.IcCoeffs, ranges=self.domains))
        # perform forward integration
        self.res  = fixed_rk4(ODE, np.array([0,self.x0,self.y0]), self.ximesh, (SF[0], SF[1], dgds0))
        # calcualte error values
        Rinitial = lin.norm([self.xL,self.yL])
        # print(self.xL, self.yL)
        Rfinal   = lin.norm([self.res[1,-1],self.res[2,-1]])
        # print("radius before integration:",Rinitial)
        # print("radius assumed by integration:",Rfinal)
        # print(res[1,-2], res[2,-1])
        dBeta    = np.arctan2(self.res[2,-1],self.res[1,-1])-self.betaAngles[-1]
        # print("deflection:",dBeta)
        # Err = diff. in radius, diff between gamma(L) and beta(L)
        err = np.array([Rinitial-Rfinal, self.res[0,-1]-dBeta, SF[2]-torqueTarg])
        return err, self.res

    def estimate_jacobian_fd(self, ODE, SF, torqueTarg, n=3):
        Jac = np.empty((3,3))
        for i in range(n):
            finiteDifference = np.zeros(len(SF))
            if i < 2:
                finiteDifference[i] = self.finiteDifferenceForce
                # step = 1
            else:
                finiteDifference[i] = self.finiteDifferenceTorque
                # step = 0.1
            errBack, resG = self.forward_integration(ODE, SF-finiteDifference, torqueTarg)
            errForw, resG = self.forward_integration(ODE, SF+finiteDifference, torqueTarg)
            Jac[0,i]     = (errForw[0]-errBack[0])/(2*lin.norm(finiteDifference))
            Jac[1,i]     = (errForw[1]-errBack[1])/(2*lin.norm(finiteDifference))
            Jac[2,i]     = (errForw[2]-errBack[2])/(2*lin.norm(finiteDifference))
        Jac = Jac[0:n,0:n]
        # print(Jac)
        return Jac

    def deform_ODE(self, xi, p, *args):
        # print(args)
        Fx    = args[0][0][0][0]
        Fy    = args[0][0][0][1]
        dgds0 = args[0][0][0][2]

        Mdim = self.E*PPoly_Eval(xi, self.IcCoeffs, ranges=self.domains)
        dxdxi = PPoly_Eval(xi, self.XCoeffs, deriv=1)
        dydxi = PPoly_Eval(xi, self.YCoeffs, deriv=1)

        dxids = d_xi_d_s(xi, self.XCoeffs, self.YCoeffs)
        # print(dxids)
        dxds = dxids*dxdxi
        dyds = dxids*dydxi

        LHS = np.empty(3)
        # print(dgds0)
        LHS[0] = (dgds0 + Fy/Mdim*(p[1]-self.x0) - Fx/Mdim*(p[2]-self.y0))/dxids
        if xi==0:
            # print(LHS[0],dgds0)
            assert(np.isclose(LHS[0],dgds0,rtol=1e-5))
        LHS[1] = np.cos(np.arctan2(dyds,dxds)+p[0])/dxids
        LHS[2] = np.sin(np.arctan2(dyds,dxds)+p[0])/dxids

        # if np.sin(np.arctan2(dydxi,dxdxi))==0:
        #     LHS = LHS*dxdxi/np.cos(np.arctan2(dydxi,dxdxi))
        #     # print(LHS)
        # else:
        #     LHS = LHS*dydxi/np.sin(np.arctan2(dydxi,dxdxi))
            # print(LHS)
        return LHS

    def l_a_l_b_rootfinding(self, s, lABPrev, printBool=0):
        # get the relevant values at a given point
        rn = r_n(s, self.XCoeffs, self.YCoeffs)
        Ic = PPoly_Eval(s, self.IcCoeffs, ranges = self.domains)
        # define the system to be solved
        def func(x, rn, cI):
            f1 = (x[0]+x[1])/(np.log((rn+x[1])/(rn-x[0])))-rn
            # f2 = (x[0]+x[1])*(outPlaneThickness*(x[1]-(x[1]+x[0])/2)*rn)-cI
            f2 = self.t*rn*(x[0]+x[1])*(x[1]/2-x[0]/2)-cI
            return np.array([f1, f2])
        def jac(x, rn, cI):
            return np.array([[1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn-x[0])), \
                            1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn+x[1]))], \
                            [-rn*self.t*x[0], rn*self.t*x[1]]])
        # some error checking an escapes for possible non-convergent cases
        if not np.isinf(rn):
            if lABPrev[0]==lABPrev[1]:
                x0 = [lABPrev[0], lABPrev[1]+0.001]
            else:
                x0 = lABPrev
            err = 1
        # excape for locally straight beam
        else:
            err = 0
            l = np.cbrt(12*Ic/self.t)/2
            x0 = [l, l]
        x = x0

        iii = 0
        # solve the problem
        while err > 10**-6 and iii <500:
            # print("entered while")
            xprev = x
            x = x - np.transpose(lin.inv(jac(x, rn, Ic)).dot(func(x, rn, Ic)))
            # print(x)
            err = lin.norm(x-xprev)
            # err = lin.norm(func(x, rn, cI))
            iii+=1
        # escape for beam so uncurved that its _basically_ straight:
        if(lin.norm(x)>lin.norm(lABPrev)*100):
            l = np.cbrt(12*Ic/self.t)/2
            x = [l,l]
        if(printBool):
            print(x0)
            print(x)
            print(iii)
            print(rn)
            print(err)
            if iii > 2999:
                print("--------------------DID NOT CONVERGE-------------------------")

        return x    # here x is [la, lb]

    def generate_surfaces(self):
        lABPrev = [0,0]

        self.la = np.empty(len(self.ximesh))
        self.lb = np.empty(len(self.ximesh))
        self.h = np.empty(len(self.ximesh))
        self.Ic = PPoly_Eval(self.ximesh, self.IcCoeffs, ranges=self.domains)
        self.rn = r_n(self.ximesh, self.XCoeffs, self.YCoeffs)
        for i in range(len(self.ximesh)):
            SSProfile("lAB rootfinding").tic()
            lAB = self.l_a_l_b_rootfinding(self.ximesh[i], lABPrev)
            SSProfile("lAB rootfinding").toc()
            self.la[i] = lAB[0]
            self.lb[i] = lAB[1]
            self.h[i] = self.lb[i]+self.la[i]
            lABPrev = lAB
        self.ecc = self.Ic/(self.t*self.h*self.rn)
        self.a = self.rn-self.la
        self.b = self.rn+self.la

        self.alpha = alpha_xy(self.ximesh, self.XCoeffs, self.YCoeffs)
        # generate xy paths for surfaces
        self.undeformedBSurface = self.undeformedNeutralSurface+np.hstack((np.atleast_2d(-self.lb*np.sin(self.alpha)).T, np.atleast_2d(self.lb*np.cos(self.alpha)).T))
        self.undeformedASurface = self.undeformedNeutralSurface-np.hstack((np.atleast_2d(-self.la*np.sin(self.alpha)).T, np.atleast_2d(self.la*np.cos(self.alpha)).T))

        return self.undeformedASurface, self.undeformedBSurface

    def spring_geometry(self, plotBool=1, deformBool=1):

        ## Things to make for the undeformed state:

        # Generate neutral radius path and give it nicely formatted class variables
        self.undeformedNeutralSurface = np.hstack((np.atleast_2d(PPoly_Eval(self.ximesh, self.XCoeffs)).T, np.atleast_2d(PPoly_Eval(self.ximesh, self.YCoeffs)).T))
        # Generate outer and inner surfaces
        self.generate_surfaces()
        # generate centroidal surface
        self.undeformedCentroidalSurface = self.undeformedNeutralSurface+np.hstack((np.atleast_2d(self.ecc*np.sin(self.alpha)).T, np.atleast_2d(self.ecc*np.cos(self.alpha)).T))
        if plotBool:
            # plot shit
            plt.figure(1)
            plt.plot(self.undeformedNeutralSurface[:,0],self.undeformedNeutralSurface[:,1])
            plt.plot(self.undeformedCentroidalSurface[:,0],self.undeformedCentroidalSurface[:,1])
            plt.plot(self.undeformedASurface[:,0],self.undeformedASurface[:,1])
            plt.plot(self.undeformedBSurface[:,0],self.undeformedBSurface[:,1])

        ## Things to make for the deformed state:

        if deformBool:
            # generate neutral surface after deformation (and give nice format)
            self.deformedNeutralSurface = np.hstack((np.atleast_2d(self.res[1,:]).T, np.atleast_2d(self.res[2,:]).T))
            if plotBool:
                plt.plot(self.deformedNeutralSurface[:,0],self.deformedNeutralSurface[:,1])
        if plotBool:
            plt.axis("equal")
            plt.show()

def alpha_xy(s, xCoeffs, yCoeffs):
    if hasattr(s, "__len__"):
        alphaList = np.arctan2(PPoly_Eval(s, yCoeffs, deriv=1),PPoly_Eval(s, xCoeffs, deriv=1))
        for i in range(len(alphaList)):
            if alphaList[i]<0:
                alphaList[i]=alphaList[i]+2*math.pi
        return alphaList
    else:
        alpha = np.arctan2(PPoly_Eval(s, yCoeffs, deriv=1),PPoly_Eval(s, xCoeffs, deriv=1))
        if alpha<0:
            alpha=alpha+2*math.pi
        return alpha

def r_n(s, xCoeffs, yCoeffs):
    d2yds2 = PPoly_Eval(s, yCoeffs, deriv=2)
    d2xds2 = PPoly_Eval(s, xCoeffs, deriv=2)
    dyds   = PPoly_Eval(s, yCoeffs, deriv=1)
    dxds   = PPoly_Eval(s, xCoeffs, deriv=1)
    # dAlphadS = (d2yds2/dxds-d2xds2*dyds/dxds)/(1+dyds**2/dxds)
    if hasattr(s, "__len__"):
        rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
        for i in range(len(rn)):
            if rn[i] == float('nan'):
                rn[i] = float('inf')
            if abs((d2yds2[i]/dxds[i]-d2xds2[i]*dyds[i]/dxds[i]**2)) <= 10**-13:
                rn[i] = float('inf')*np.sign(rn[i-1])
    else:
        if abs(d2yds2/dxds-d2xds2*dyds/dxds**2)<=10**-13:
            rn = float('inf')
        else:
            rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
    return rn

def d_rn_d_s(s, xCoeffs, yCoeffs):
    d3yds3 = PPoly_Eval(s, yCoeffs, deriv=3)
    d3xds3 = PPoly_Eval(s, xCoeffs, deriv=3)
    d2yds2 = PPoly_Eval(s, yCoeffs, deriv=2)
    d2xds2 = PPoly_Eval(s, xCoeffs, deriv=2)
    dyds   = PPoly_Eval(s, yCoeffs, deriv=1)
    dxds   = PPoly_Eval(s, yCoeffs, deriv=1)

    denominator = (d2yds2*dxds-d2xds2*dyds)**2
    numerator   = ( dxds**3*d3yds3 - dyds**3*d3xds3 +
                    2*dyds*d2xds2**2*dxds - 2*dyds*d2yds2**2*dxds - dxds**2*d3xds3*dyds -
                    2*dxds**2*d2yds2*d2xds2 + 2*dyds**2*d2yds2*d2xds2 +
                    dyds**2*d3yds3*dxds )

    drnds = -numerator/denominator

    return drnds

def deform_ODE(s, p, *args): #Fx, Fy, geometryDef, dgds0 THIS IS THE NON_OBJECT ORIENTED VERSION

    # deal with args in the shittiest way possible

    # p = [gamma, x_dis, y_dis]

    Fx          = args[0][0][0][0]
    Fy          = args[0][0][0][1]
    geometryDef = args[0][0][0][2]
    dgds0       = args[0][0][0][3]

    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    cICoeffs = geometryDef[2]

    # Mdim = EI (dimensionalize bending moment)
    # dxds, dyds are ANALYTICAL derivatives based on x/y polynomials

    Mdim = E*PPoly_Eval(s, cICoeffs)

    dxds = d_coord_d_s(s, xCoeffs)
    dyds = d_coord_d_s(s, yCoeffs)
    # constants
    LHS = np.empty(3)

    xcoord0 = coord(0,xCoeffs)
    ycoord0 = coord(0,yCoeffs)

    LHS[0] = dgds0
    # LHS[0] = LHS[0] + Fy/Mdim*(p[1]-xcoord0) - Fx/Mdim*(p[2]-ycoord0)
    LHS[0] = LHS[0] + Fy/Mdim*(p[1]-xcoord0) - Fx/Mdim*(p[2]-ycoord0)
    if s==0:
        # print(LHS[0],dgds0)
        assert(np.isclose(LHS[0],dgds0,rtol=1e-3))
    LHS[1] = np.cos(np.arctan2(dyds,dxds)+p[0]) # -(g*y +/- (g^2 - y^2 + 1)^(1/2))/(g^2 + 1)
    LHS[2] = np.sin(np.arctan2(dyds,dxds)+p[0])

    if np.sin(np.arctan2(dyds,dxds))==0:
        LHS[0] = LHS[0]*dxds/np.cos(np.arctan2(dyds,dxds))
        LHS[1] = LHS[1]*dxds/np.cos(np.arctan2(dyds,dxds))
        LHS[2] = LHS[2]*dxds/np.cos(np.arctan2(dyds,dxds))
    else:
        LHS[0] = LHS[0]*dyds/np.sin(np.arctan2(dyds,dxds))
        LHS[1] = LHS[1]*dyds/np.sin(np.arctan2(dyds,dxds))
        LHS[2] = LHS[2]*dyds/np.sin(np.arctan2(dyds,dxds))
    # print(LHS[1], LHS[2])
    # print(dxds, dyds)
    # print("--------------")

    # print(LHS[0], LHS[1], LHS[2])
    # print(dxds, np.cos(alpha+p[0]), dxds/np.cos(alpha))
    # print(dxds/np.cos(alpha))
    # print("--------------------")

    return LHS

def d_rn_d_s_numerical(s, xCoeffs, yCoeffs, eps=1e-5):
    rnf = r_n(s+eps, xCoeffs, yCoeffs)
    rnb = r_n(s-eps, xCoeffs, yCoeffs)
    derivative = (rnf-rnb)/(2*eps)
    return derivative

def ndiff(s, f, eps=1e-5):
    rnf = f(s+eps)
    rnb = f(s-eps)
    derivative = (rnf-rnb)/(2*eps)
    return derivative

# def xy_poly(pts, ctrlLengths):

#     normalizeFinalAngle = lin.norm(pts[3,:])

#     xTarg = np.array([[pts[0,0]], # x(0)       = P0_x
#                       [pts[1,0]], # x(s1)      = P1_x
#                       [pts[2,0]], # x(s2)      = P2_x
#                       [pts[3,0]], # x(s3)      = P3_x
#                       [1],        # dxdxi(0)   = 1 (flat, dydx=0)
#                       [pts[3,0]/normalizeFinalAngle], # dxdxi(s3)  = P3_x (radial, dydx = P3_y/P3_x)
#                       [0],        # d2xdxi2(0) = 0    (d2ydx2(0)=0, straight)
#                       [0]])       # d2xdxi2(0) = 0    (d2ydx2(L)=0, straight)

#     yTarg = np.array([[pts[0,1]], # y(0)       = P0_y
#                       [pts[1,1]], # y(s1)      = P1_y
#                       [pts[2,1]], # y(s2)      = P2_y
#                       [pts[3,1]], # y(s3)      = P3_y
#                       [0],        # dydxi(0)   = 0 (flat, dydx=0)
#                       [pts[3,1]/normalizeFinalAngle], # dydxi(s3)  = P3_y (radial, dydx = P3_y/P3_x)
#                       [0],        # d2ydxi2(0) = 0    (d2ydx2(0)=0, straight)
#                       [0]])       # d2ydxi2(0) = 0    (d2ydx2(L)=0, straight)

#     Mat = np.array([[0,0,0,0,0,0,0,1],                                                                                                              # coord(0)
#                     [ctrlLengths[1]**7,ctrlLengths[1]**6,ctrlLengths[1]**5,ctrlLengths[1]**4,ctrlLengths[1]**3,ctrlLengths[1]**2,ctrlLengths[1],1], # coord(s1)
#                     [ctrlLengths[2]**7,ctrlLengths[2]**6,ctrlLengths[2]**5,ctrlLengths[2]**4,ctrlLengths[2]**3,ctrlLengths[2]**2,ctrlLengths[2],1], # coord(s2)
#                     [ctrlLengths[3]**7,ctrlLengths[3]**6,ctrlLengths[3]**5,ctrlLengths[3]**4,ctrlLengths[3]**3,ctrlLengths[3]**2,ctrlLengths[3],1], # coord(s3)
#                     [0,0,0,0,0,0,1,0],                                                                                                              # dcoorddxi(0)
#                     [7*ctrlLengths[3]**6,6*ctrlLengths[3]**5,5*ctrlLengths[3]**4,4*ctrlLengths[3]**3,3*ctrlLengths[3]**2,2*ctrlLengths[3],1,0],     # dcoorddxi(s3)
#                     [0,0,0,0,0,2,0,0],                                                                                                              # d2coorddxi2(0)
#                     [42*ctrlLengths[3]**5,30*ctrlLengths[3]**4,20*ctrlLengths[3]**3,12*ctrlLengths[3]**2,6*ctrlLengths[3]**1,2,0,0]])                # d2coorddxi2(s3)
#     # print(Mat)
#     xCoeffs = lin.solve(Mat,xTarg)
#     yCoeffs = lin.solve(Mat,yTarg)

#     assert(np.allclose(np.dot(Mat,xCoeffs),xTarg))
#     assert(np.allclose(np.dot(Mat,yCoeffs),yTarg))

#     # print("These need to match:",np.dot(Mat,xCoeffs),xTarg)

#     return xCoeffs, yCoeffs

# def coord(s, coeffs):

#     if hasattr(s, "__len__"):
#         coord = np.empty(len(s))
#         i = 0
#         for value in s:
#             U = np.array([value**7, value**6, value**5, value**4, value**3, value**2, value, 1])
#             coord[i] = U.dot(coeffs)[0]
#             i+=1
#         return coord
#     else:
#         U = np.array([s**7, s**6, s**5, s**4, s**3, s**2, s, 1])
#         coord = U.dot(coeffs)
#         return coord[0]

# def d_coord_d_s(s, coeffs):
#     if hasattr(s, "__len__"):
#         dCds = np.empty(len(s))
#         i=0
#         for value in s:
#             U = np.array([7*value**6, 6*value**5, 5*value**4, 4*value**3, 3*value**2, 2*value, 1, 0])
#             dCds[i] = U.dot(coeffs)[0]
#             i+=1
#         return dCds
#     else:
#         U = np.array([7*s**6, 6*s**5, 5*s**4, 4*s**3, 3*s**2, 2*s, 1, 0])
#         dCds = U.dot(coeffs)
#         return dCds[0]

# def d2_coord_d_s2(s, coeffs):
#     if hasattr(s, "__len__"):
#         dCds = np.empty(len(s))
#         i=0
#         for value in s:
#             U = np.array([42*value**5, 30*value**4, 20*value**3, 12*value**2, 6*value**1, 2, 0, 0])
#             dCds[i] = U.dot(coeffs)[0]
#             i+=1
#         return dCds
#     else:
#         U = np.array([42*s**5, 30*s**4, 20*s**3, 12*s**2, 6*s**1, 2, 0, 0])
#         # print(U)
#         dCds = U.dot(coeffs)
#         return dCds[0]

# def d3_coord_d_s3(s, coeffs):
    if hasattr(s, "__len__"):
        dCds = np.empty(len(s))
        i=0
        for value in s:
            U = np.array([210*value**4, 120*value**3, 60*value**2, 24*value**1, 6, 0, 0, 0])
            dCds[i] = U.dot(coeffs)[0]
            i+=1
        return dCds
    else:
        U = np.array([210*s**4, 120*s**3, 60*s**2, 24*s**1, 6, 0, 0, 0])
        dCds = U.dot(coeffs)
        return dCds[0]

# def alpha_polyfit(smesh,yCoeffs, xCoeffs ):

#     alphaList = np.arctan2(d_coord_d_s(smesh, yCoeffs),d_coord_d_s(smesh, xCoeffs))
#     for i in range(len(alphaList)):
#         if alphaList[i]<0:
#             alphaList[i]=alphaList[i]+2*math.pi

#     alphaPoly = np.polyfit(smesh,alphaList,9)
#     return alphaPoly


# def alpha_xy(s, xCoeffs, yCoeffs): DEPRECATED
#     if hasattr(s, "__len__"):
#         alphaList = np.arctan2(d_coord_d_s(s, yCoeffs),d_coord_d_s(s, xCoeffs))
#         for i in range(len(alphaList)):
#             if alphaList[i]<0:
#                 alphaList[i]=alphaList[i]+2*math.pi
#         return alphaList
#     else:
#         alpha = np.arctan2(d_coord_d_s(s, yCoeffs),d_coord_d_s(s, xCoeffs))
#         if alpha<0:
#             alpha=alpha+2*math.pi
#         return alpha

def d_alpha_d_s(s, xCoeffs, yCoeffs):
    d2yds2 = d2_coord_d_s2(s, yCoeffs)
    d2xds2 = d2_coord_d_s2(s, xCoeffs)
    dyds   = d_coord_d_s  (s, yCoeffs)
    dxds   = d_coord_d_s  (s, xCoeffs)
    dads = ((d2yds2/dxds-d2xds2*dyds/dxds**2)/(1+dyds**2/dxds**2))
    return dads

def d_2_alpha_d_s_2(s, xCoeffs, yCoeffs):
    d2ads2 = -d_alpha_d_s(s, xCoeffs, yCoeffs)**2*d_rn_d_s(s, xCoeffs, yCoeffs)
    return d2ads2

# def r_n(s, xCoeffs, yCoeffs): DEPRECATED
#     d2yds2 = d2_coord_d_s2(s, yCoeffs)
#     d2xds2 = d2_coord_d_s2(s, xCoeffs)
#     dyds   = d_coord_d_s  (s, yCoeffs)
#     dxds   = d_coord_d_s  (s, xCoeffs)
#     # dAlphadS = (d2yds2/dxds-d2xds2*dyds/dxds)/(1+dyds**2/dxds)
#     if hasattr(s, "__len__"):
#         rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#         for i in range(len(rn)):
#             if abs((d2yds2[i]/dxds[i]-d2xds2[i]*dyds[i]/dxds[i]**2)) <= 10**-13:
#                 rn[i] = float('inf')*np.sign(rn[i-1])
#     else:
#         if abs(d2yds2/dxds-d2xds2*dyds/dxds**2)<=10**-13:
#             rn = float('inf')
#         else:
#             rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#     return rn

def d_rn_d_s(s, xCoeffs, yCoeffs):
    d3yds3 = d3_coord_d_s3(s, yCoeffs)
    d3xds3 = d3_coord_d_s3(s, xCoeffs)
    d2yds2 = d2_coord_d_s2(s, yCoeffs)
    d2xds2 = d2_coord_d_s2(s, xCoeffs)
    dyds   = d_coord_d_s  (s, yCoeffs)
    dxds   = d_coord_d_s  (s, xCoeffs)

    denominator = (d2yds2*dxds-d2xds2*dyds)**2
    numerator   = ( dxds**3*d3yds3 - dyds**3*d3xds3 +
                    2*dyds*d2xds2**2*dxds - 2*dyds*d2yds2**2*dxds - dxds**2*d3xds3*dyds -
                    2*dxds**2*d2yds2*d2xds2 + 2*dyds**2*d2yds2*d2xds2 +
                    dyds**2*d3yds3*dxds )

    drnds = -numerator/denominator

    return drnds



# def d_cI_d_s(s, xCoeffs, yCoeffs, hCoeffs):
#     dcIds = (cI_s(s+finiteDifferenceLength, xCoeffs, yCoeffs, hCoeffs)-cI_s(s-finiteDifferenceLength, xCoeffs, yCoeffs, hCoeffs))/(2*finiteDifferenceLength)
#     return dcIds

def d_cI_d_s(s, cICoeffs):
    if hasattr(s, "__len__"):
        dcIcs = np.empty(len(s))
        i = 0
        for value in s:
            U = np.array([7*value**6, 6*value**5, 5*value**4, 4*value**3, 3*value**2, 2*value**1, 1, 0])
            dcIcs[i] = U.dot(cICoeffs)[0]
            i+=1
        return dcIcs
    else:
        U = np.array([7*s**6, 6*s**5, 5*s**4, 4*s**3, 3*s**2, 2*s, 1, 0])
        dcIcs = U.dot(cICoeffs)
        return dcIcs[0]


def rc_rootfinding(s, xCoeffs, yCoeffs, cICoeffs, printBool):
    rn =  r_n(s, xCoeffs, yCoeffs)
    cI = PPoly_Eval(s, cICoeffs)
    # h = (cI/(outPlaneThickness*(x-rn)*rn))

    def func(x, rn, h):
        f1 = h/((np.log((x+0.5*h)/(x-0.5*h))))-rn
        return f1
    def jac(x, rn, h):
        jac = 4*h**2/((4*x**2-h**2)*np.log((x+0.5*h)/(x-0.5*h))**2)
        return jac
    err = 1

    # x0 = h*1.5
    x0 = rn*1.1
    x = x0
    ii = 0
    if np.isinf(rn):
        x = float('inf')
    else:
        while err > 10**-6 and ii < 3000:
            xprev = x
            h = (cI/(self.t*(x-rn)*rn))
            print('h,',h)
            # if ii==0:
            #     print(x,h)
            x = x - (func(x, rn, h))/jac(x, rn, h)
            print("next x")
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
    cI = PPoly_Eval(s, cICoeffs)
    def func(x, rn, cI):
        # x0 = a, x1 = b
        f1 = (x[1]-x[0])/(np.log(x[1]/x[0]))-rn
        f2 = self.t*(x[1]-x[0])*((x[1]+x[0])/2-rn)*rn-cI
        return np.array([f1, f2])
    def jac(x, rn, cI):
        return np.array([[(x[1]-x[0])/(x[0]*np.log(x[1]/x[0])**2)-1/(np.log(x[1]/x[0])), \
                          (x[1]*np.log(x[1]/x[0])-x[1]+x[0])/(x[1]*np.log(x[1]/x[0])**2)], \
                         [-rn*self.t*(x[0]-rn), rn*self.t*(x[1]-rn)]])
    err = 1
    # establish expected closesness range for a/b compared to rn to establish initial guess based on linearization about 1
    factor = 1.2
    # approimation: ln(x)=c(x-1), c such that c(x-1) = ln(factor)
    c = (np.log(factor)-np.log(1))/(factor-1)
    a0 = c*rn
    b0 = rn + np.sqrt(rn**2+4*0.5*((c-c**2/2)*rn-cI/(self.t*rn)))
    x0 = [a0, b0]
    x = x0
    ii = 0
    if np.isinf(rn):
        x = np.array([float('inf'), float('inf')])
    else:
        while err > 10**-6 and ii < 3000:
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

def h_e_rootfinding(s, xCoeffs, yCoeffs, cICoeffs, printbool):
    rn = r_n(s, xCoeffs, yCoeffs)
    cI = PPoly_Eval(s, cICoeffs)
    # x = [e h]
    def func(x, rn, cI):
        f1 = (x[1])/(np.log((rn+x[0]+x[1]/2)/(rn+x[0]-x[1]/2)))-rn
        f2 = cI/(self.t*x[0]*rn)-x[1]
        return np.array([f1, f2])
    def jac(x, rn, cI):
        return np.array([[1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn-x[0])), \
                          1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn+x[1]))], \
                         [cI/(2*rn*self.t*(x[1]-(x[1]+x[0])/2)**2)-1, -cI/(2*rn*self.t*(x[1]-(x[1]+x[0])/2)**2)-1]])


def l_a_l_b_rootfinding(s, lABPrev, xCoeffs, yCoeffs, cICoeffs, printBool):
    rn = r_n(s, xCoeffs, yCoeffs)
    cI = PPoly_Eval(s, cICoeffs)
    def func(x, rn, cI):
        f1 = (x[0]+x[1])/(np.log((rn+x[1])/(rn-x[0])))-rn
        # f2 = (x[0]+x[1])*(outPlaneThickness*(x[1]-(x[1]+x[0])/2)*rn)-cI
        f2 = self.t*rn*(x[0]+x[1])*(x[1]/2-x[0]/2)-cI
        return np.array([f1, f2])
    def jac(x, rn, cI):
        return np.array([[1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn-x[0])), \
                          1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn+x[1]))], \
                         [-rn*self.t*x[0], rn*self.t*x[1]]])
    # print(rn)
    if not np.isinf(rn):
        if lABPrev[0]==lABPrev[1]:
            x0 = [lABPrev[0], lABPrev[1]+0.001]
        else:
            x0 = lABPrev
        err = 1
    else:
        err = 0
        l = np.cbrt(12*cI/self.t)/2
        x0 = [l, l]
        # print("entered")
    x = x0
    # print("x0:",x0)
    iii = 0
    while err > 10**-6 and iii <500:
        # print("entered while")
        xprev = x
        x = x - np.transpose(lin.inv(jac(x, rn, cI)).dot(func(x, rn, cI)))
        # print(x)
        err = lin.norm(x-xprev)
        # err = lin.norm(func(x, rn, cI))
        iii+=1
    if(lin.norm(x)>lin.norm(lABPrev)*100):
        l = np.cbrt(12*cI/self.t)/2
        x = [l,l]
    if(printBool):
        print(x0)
        print(x)
        print(iii)
        print(rn)
        print(err)
        if iii > 2999:
            print("--------------------DID NOT CONVERGE-------------------------")

    # print("DONE WITH S = ",s)
    return x    # here x is [la, lb]



# def cI_s(s, xCoeffs, yCoeffs, hCoeffs):
    if r_n(s, xCoeffs, yCoeffs) == float('inf'):
        cI = PPoly_Eval(s, hCoeffs)**3*self.t/12
    else:
        cI = PPoly_Eval(s, hCoeffs)*self.t* \
            (rc_rootfinding(s, xCoeffs, yCoeffs, hCoeffs, False)-r_n(s, xCoeffs, yCoeffs)) \
            *r_n(s, xCoeffs, yCoeffs)
    return cI

def form_spring(pts, cIs, ctrlLengths, ctrlcIs):
    xCoeffs = xy_poly(pts, ctrlLengths)[0]
    yCoeffs = xy_poly(pts, ctrlLengths)[1]
    cICoeffs = Ic_poly(cIs, ctrlcIs)
    return xCoeffs, yCoeffs, cICoeffs

# def deform_ODE_Barcio(s, p, *args): #Fx, Fy, geometryDef, dgds0

    Fx          = args[0][0][0][0]
    Fy          = args[0][0][0][1]
    geometryDef = args[0][0][0][2]
    dgds0       = args[0][0][0][3]

    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    cICoeffs = geometryDef[2]

    Tdim = E*PPoly_Eval(s, xCoeffs, yCoeffs, cICoeffs)
    yorg = coord(s, yCoeffs)
    xorg = coord(s, xCoeffs)
    dxds = d_coord_d_s(s, xCoeffs)
    dyds = d_coord_d_s(s, yCoeffs)

    alpha = alpha_xy(s, xCoeffs, yCoeffs)

    LHS = np.empty(2)
    LHS[0] = p[1]
    LHS[1] = (Fx*np.sin(alpha+p[0])-Fy*np.cos(alpha+p[0])-E*d_cI_d_s(s, xCoeffs, yCoeffs, cICoeffs)*p[1])/(PPoly_Eval(s, xCoeffs, yCoeffs, cICoeffs)*E)
    return LHS

# def deform_ODE_Barcio_to_Thomas(s, p, *args):

    pass

def deform_with_stress_ODE(s, p, *args):
    Fx          = args[0][0][0][0]
    Fy          = args[0][0][0][1]
    geometryDef = args[0][0][0][2]
    dgds0       = args[0][0][0][3]
    stressTarg  = args[0][0][0][4]

    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    cICoeffs = geometryDef[2]
    rn     = r_n(s, geometryDef[0], geometryDef[1])
    cI     = PPoly_Eval(s, geometryDef[2])
    cIPrev = cI
    Mdim = E*PPoly_Eval(s, cICoeffs)

    dxds = d_coord_d_s(s, xCoeffs)
    dyds = d_coord_d_s(s, yCoeffs)
    stressError = 1
    lABPrev = [0,0]
    while stressError>10^-6:
        LHS = np.empty(3)

        LHS[0] = dgds0

        LHS[0] = LHS[0] + Fy/Mdim*p[1] - Fx/Mdim*p[2]
        LHS[1] = np.cos(np.arctan2(dyds,dxds)+p[0]) # -(g*y +/- (g^2 - y^2 + 1)^(1/2))/(g^2 + 1)
        LHS[2] = np.sin(np.arctan2(dyds,dxds)+p[0])
        # print(LHS[1], LHS[2])
        # print("=>")
        LHS[1] = LHS[1]*dxds/np.cos(np.arctan2(dyds,dxds))
        LHS[2] = LHS[2]*dyds/np.sin(np.arctan2(dyds,dxds))

        lAB = l_a_l_b_rootfinding(s, lABPrev, geometryDef[0], geometryDef[1], cI, False)
        lABPrev = lAB
        stressMax = np.max([E*(1-rn/(rn-lAB[0]))*rn*LHS[0], E*(1-rn/(rn+lAB[1]))*rn*LHS[0]])
    return LHS, stressMax

def deform_ODE(s, p, *args): #Fx, Fy, geometryDef, dgds0

    # deal with args in the shittiest way possible

    # p = [gamma, x_dis, y_dis]

    Fx          = args[0][0][0][0]
    Fy          = args[0][0][0][1]
    geometryDef = args[0][0][0][2]
    dgds0       = args[0][0][0][3]

    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    cICoeffs = geometryDef[2]

    # Mdim = EI (dimensionalize bending moment)
    # dxds, dyds are ANALYTICAL derivatives based on x/y polynomials

    Mdim = E*PPoly_Eval(s, cICoeffs)

    dxds = d_coord_d_s(s, xCoeffs)
    dyds = d_coord_d_s(s, yCoeffs)
    # constants
    LHS = np.empty(3)

    xcoord0 = coord(0,xCoeffs)
    ycoord0 = coord(0,yCoeffs)

    LHS[0] = dgds0
    # LHS[0] = LHS[0] + Fy/Mdim*(p[1]-xcoord0) - Fx/Mdim*(p[2]-ycoord0)
    LHS[0] = LHS[0] + Fy/Mdim*(p[1]-xcoord0) - Fx/Mdim*(p[2]-ycoord0)
    if s==0:
        # print(LHS[0],dgds0)
        assert(np.isclose(LHS[0],dgds0,rtol=1e-3))

    LHS[1] = np.cos(np.arctan2(dyds,dxds)+p[0]) # -(g*y +/- (g^2 - y^2 + 1)^(1/2))/(g^2 + 1)
    LHS[2] = np.sin(np.arctan2(dyds,dxds)+p[0])

    if np.sin(np.arctan2(dyds,dxds))==0:
        LHS[0] = LHS[0]*dxds/np.cos(np.arctan2(dyds,dxds))
        LHS[1] = LHS[1]*dxds/np.cos(np.arctan2(dyds,dxds))
        LHS[2] = LHS[2]*dxds/np.cos(np.arctan2(dyds,dxds))
    else:
        LHS[0] = LHS[0]*dyds/np.sin(np.arctan2(dyds,dxds))
        LHS[1] = LHS[1]*dyds/np.sin(np.arctan2(dyds,dxds))
        LHS[2] = LHS[2]*dyds/np.sin(np.arctan2(dyds,dxds))
    # print(LHS[1], LHS[2])
    # print(dxds, dyds)
    # print("--------------")

    # print(LHS[0], LHS[1], LHS[2])
    # print(dxds, np.cos(alpha+p[0]), dxds/np.cos(alpha))
    # print(dxds/np.cos(alpha))
    # print("--------------------")

    return LHS

def d_rn_d_s_numerical(s, xCoeffs, yCoeffs, eps=1e-5):
    rnf = r_n(s+eps, xCoeffs, yCoeffs)
    rnb = r_n(s-eps, xCoeffs, yCoeffs)
    derivative = (rnf-rnb)/(2*eps)
    return derivative

def ndiff(s, f, eps=1e-5):
    rnf = f(s+eps)
    rnb = f(s-eps)
    derivative = (rnf-rnb)/(2*eps)
    return derivative

# ndiff(s, lambda s: r_n(s, geometryDef[0], geometryDef[1]), eps=1e-6)
class geo_ODE_wrapper:
    def __init__(self):
        self.all_inputs_ever = []
        self.all_p_inputs_ever = []
        self.all_outputs_ever = []

        pass
    def __call__(self, s, p , geometryDef):
        self.all_inputs_ever.append(s)
        self.all_p_inputs_ever.append(np.array(p))
        ret = self.geo_ODE(s,p,geometryDef)
        self.all_outputs_ever.append(ret)
        return ret

    def geo_ODE(self, s, p, geometryDef):
        geometryDef = geometryDef[0][0]
        # print(geometryDef)

        rn     = r_n(s, geometryDef[0], geometryDef[1])
        drnds  = d_rn_d_s(s, geometryDef[0], geometryDef[1])
        # drnds  = d_rn_d_s_numerical(s, geometryDef[0], geometryDef[1])
        cI     = PPoly_Eval(s, geometryDef[2])
        dcIds  = d_cI_d_s(s, geometryDef[2])
        dads   = d_alpha_d_s(s, geometryDef[0], geometryDef[1])
        d2ads2 = d_2_alpha_d_s_2(s, geometryDef[0], geometryDef[1])
        # print(rn, s)
        geoFunc1 = (dcIds*dads+cI*d2ads2)/self.t
        # geoFunc2 = rn*(-drnds*((p[0]+p[1])/rn**2+(-p[0]-p[1])/((rn+p[1])*(rn-p[0]))))
        geoFunc2 = d2ads2*(p[0]+p[1])*(p[1]-p[0]-p[0]*p[1]*dads)
        # print(geoFunc1, geoFunc2)
        geoFuncs = np.array([[geoFunc1], [geoFunc2]])
        states   = np.array([[-p[0], p[1]], [p[0]*dads*(1+p[1]*dads), -p[1]*dads*(1-p[0]*dads)]])# [p[0]*p[1]*rn-p[0]*rn**2, p[0]*p[1]*rn-p[1]*rn**2]])
        # print("states matrix:",states)
        if np.all(np.isfinite(states)):
            if lin.matrix_rank(states) < 2:
                LHS = np.array([float("nan"), float("nan")])
                return LHS
            else:
                LHS = lin.inv(states).dot(geoFuncs)
        else:
                LHS = lin.inv(states).dot(geoFuncs)
        # fuckYou = False
        # for i in [0,1]:
        #     for j in [0,1]:
        #         if np.isnan(states[i,j]):
        #             fuckYou = True
        # if not fuckYou:
        #     if lin.matrix_rank(states) < 2:
        #         LHS = np.array([float("nan"), float("nan")])
        #         return LHS
        #     else:
        #         LHS = lin.inv(states).dot(geoFuncs)
        # else:
        #     LHS = lin.inv(states).dot(geoFuncs)
        # print("state rate of change:", np.array([LHS[0][0], LHS[1][0]]))
        return np.array([LHS[0][0], LHS[1][0]])

geo_ODE = geo_ODE_wrapper()

def fixed_rk4(fun, y0, xmesh, *args): # (fun, alphaCoeffs, cICoeffs, y0, Fx, Fy, xmesh)
    step = xmesh[1]-xmesh[0]
    y0 = np.atleast_1d(y0)
    if hasattr(y0, '__len__'):
        res = np.empty((len(y0), len(xmesh)))
        # print(range(len(xmesh))[-1])
        for i in range(len(xmesh)):
            # print("i in rk:",i)
            if i == 0:
                # stepRes = rk4_step(fun, xmesh[i], y0, step, args)
                for j in range(len(y0)):
                #     res[j,i] = stepRes[j]
                # y0 = stepRes
                    res[j,i] = y0[j]
            else:
                stepRes = rk4_step(fun, xmesh[i-1], y0, step, args)
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
                stepRes = rk4_step(fun, xmesh[i-1], y0, step, args)
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

def spring_deform_eval(SF, geometryDef, totalSpringLength, deformFunc):
    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    cICoeffs = geometryDef[2]

    err = np.ones(2)*100
    i=0
    divergeFlag = 0
    stepSize = 1
    finiteDifferenceForce = 0.1
    finiteDifferenceTorque = 0.5
    while lin.norm(err,2) > 10e-6 and i<100:
        # establish results of guess
        errPrev = err
        err, res = int_error_result(deformFunc, SF, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, coord(totalSpringLength, xCoeffs), coord(totalSpringLength, yCoeffs))
        dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(coord(totalSpringLength, yCoeffs),coord(totalSpringLength, xCoeffs)))
        J = fd_J(deformFunc, SF, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, coord(totalSpringLength, xCoeffs), coord(totalSpringLength, yCoeffs), finiteDifferenceForce)

        if lin.norm(err)>lin.norm(errPrev):
            print("torque deform diverging", i)
            print(err, errPrev)
            print(SF)
            print(J)
            divergeFlag = 1
            # finiteDifferenceForce = finiteDifferenceForce*1.1
        else:
            SF = SF - lin.pinv(J).dot(err)*stepSize
        i+=1
    return res, SF, divergeFlag

def deform_initial_guess(torqueEps, torqueTarg, geometryDef, totalSpringLength, deformFunc):
    SF = np.array([0, 0, torqueEps])
    res1, SF1, divergeFlag = spring_deform_eval(SF, geometryDef, totalSpringLength, deformFunc)
    ang1 = np.arctan2(SF1[1],SF1[0])
    mag1 = lin.norm([SF1[0],SF1[1]])
    SF = np.array([0, 0, torqueEps*torqueTarg*0.1])
    res2, SF2, divergeFlag = spring_deform_eval(SF, geometryDef, totalSpringLength, deformFunc)
    ang2 = np.arctan2(SF2[1],SF2[0])
    mag2 = lin.norm([SF2[0],SF2[1]])

    print("ang1, mag1",ang1, mag1)
    print("ang2, mag2",ang2, mag2)

    ang_slope = (ang2-ang1)/(SF2[2]-SF1[2])
    mag_slope = (mag2-mag1)/(SF2[2]-SF1[2])

    ang_intercept = ang1-ang_slope*SF1[2]
    mag_intercept = mag1-mag_slope*SF1[2]

    ang_coeffs = np.array([ang_slope, ang_intercept])
    mag_coeffs = np.array([mag_slope, mag_intercept])

    plt.figure("first two guesses")
    plt.plot(res1[1,:],res1[2,:])
    plt.plot(res2[1,:],res2[2,:])
    theta = np.linspace(0, 2*np.pi, 100)
    R3 = lin.norm([coord(totalSpringLength,geometryDef[0]),coord(totalSpringLength,geometryDef[1])])
    R0 = lin.norm([coord(0,geometryDef[0]),coord(0,geometryDef[1])])
    outerCircleX = R3*np.cos(theta)
    outerCircleY = R3*np.sin(theta)
    innerCircleX = R0*np.cos(theta)
    innerCircleY = R0*np.sin(theta)
    plt.plot(outerCircleX,outerCircleY)
    plt.plot(innerCircleX,innerCircleY)
    plt.axis('equal')

    return ang_coeffs, mag_coeffs

def deform_spring_by_torque2(torqueTarg, geometryDef, totalSpringLength, deformFunc):
    # Shitty estimation for

    torqueEps = 1
    ang_coeffs, mag_coeffs = deform_initial_guess(torqueEps, torqueTarg, geometryDef, totalSpringLength, deformFunc)
    ang_guess = ang_coeffs.dot([torqueTarg,1])
    mag_guess = mag_coeffs.dot([torqueTarg,1])
    Fx = mag_guess*np.cos(ang_guess)
    Fy = mag_guess*np.sin(ang_guess)

    print("ang_guess, mag_guess", ang_guess, mag_guess)
    checkang = np.arctan2(Fy,Fx)
    checkmag = lin.norm([Fy,Fx])
    print("check ang, check mag",checkang,checkmag)

    SF = np.array([Fx, Fy, torqueTarg])
    print("initial guess:",SF)
    res, SF, divergeFlag = spring_deform_eval(SF, geometryDef, totalSpringLength, deformFunc)
    return res, SF, divergeFlag

def deform_spring_by_torque(torqueTarg, geometryDef, totalSpringLength, deformFunc):
    # print("called deform spring")
    xCoeffs = geometryDef[0]
    yCoeffs = geometryDef[1]
    cICoeffs = geometryDef[2]

    overallXDistance = coord(totalSpringLength, xCoeffs)-coord(0,xCoeffs)
    overallYDistance = coord(totalSpringLength, yCoeffs)-coord(0,yCoeffs)

    SF = np.array([0, 0, torqueTarg])
    # SF = Fx Fy T
    # print("before forward integration", SF[2])
    err = np.ones(2)*100
    i = 0
    divergeFlag = 0
    stepSize = 0.5
    finiteDifferenceForce  = 10
    finiteDifferenceTorque = 1
    while lin.norm(err,2) > 10e-6 and i<100:
        # establish results of guess
        errPrev = err
        err, res = int_error_result(deformFunc, SF, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, coord(totalSpringLength, xCoeffs), coord(totalSpringLength, yCoeffs))
        # print("called int err                                 ", end="\r")
        # print("err result:",err)
        # plt.figure("geometry results")
        # plt.plot(res[1,:],res[2,:])
        dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(coord(totalSpringLength, yCoeffs),coord(totalSpringLength, xCoeffs)))
        # print(i, dBeta/deg2rad, SF)
        # estimate jacobian with finite difference
        # print(SF)
        # print(err)
        # print(res)
        J = fd_J(deformFunc, SF, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, coord(totalSpringLength, xCoeffs), coord(totalSpringLength, yCoeffs), finiteDifferenceForce)

        if lin.norm(err)>lin.norm(errPrev):
            print("torque deform diverging", i)
            print(err, errPrev)
            print(SF)
            print(J)
            divergeFlag = 1
            finiteDifferenceForce = finiteDifferenceForce*1.1
        else:
        # print(J)
        # establish next guess:
        # print(SF)
            SF = SF - lin.pinv(J).dot(err)*stepSize
        # print("torque deform error:",lin.norm(err))
        i+=1
        # print(i)
        # print("error norm",lin.norm(err,2))
        # print("dBeta error:",res[0,-1]-dBeta)

    # print("torque iterations: ",i)
    # print("after forward integration", SF[2])
    # print(SF)
    # print("Inside full arc length:", totalSpringLength)
    return res, SF, divergeFlag


def fd_J(deformFunc, SF, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, xorg, yorg, ffForce):
    # finite difference in Fx
    # print(ffForce)

    SFback = dc(SF-np.array([ffForce,0,0]))
    # print(SFback)
    SFfrwd = dc(SF+np.array([ffForce,0,0]))

    errBack, resG = int_error_result(deformFunc, SFback, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, xorg, yorg)
    errFrwd, resG = int_error_result(deformFunc, SFfrwd, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, xorg, yorg)
    derrRdFx = (errFrwd[0]-errBack[0])/(2*ffForce)
    derrGdFx = (errFrwd[1]-errBack[1])/(2*ffForce)
    # finite difference in Fy
    SFback = dc(SF-np.array([0,ffForce,0]))
    SFfrwd = dc(SF+np.array([0,ffForce,0]))
    errBack, resG = int_error_result(deformFunc, SFback, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, xorg, yorg)
    errFrwd, resG = int_error_result(deformFunc, SFfrwd, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, xorg, yorg)
    derrRdFy = (errFrwd[0]-errBack[0])/(2*ffForce)
    derrGdFy = (errFrwd[1]-errBack[1])/(2*ffForce)
    # finite difference in T
    SFback = SF-np.array([0,0,finiteDifferenceTorque])
    SFfrwd = SF+np.array([0,0,finiteDifferenceTorque])
    errBack, resG = int_error_result(deformFunc, SFback, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, xorg, yorg)
    errFrwd, resG = int_error_result(deformFunc, SFfrwd, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, xorg, yorg)
    derrRdT = (errFrwd[0]-errBack[0])/(2*finiteDifferenceTorque)
    derrGdT = (errFrwd[1]-errBack[1])/(2*finiteDifferenceTorque)

    J = np.array([[derrRdFx, derrRdFy, 0],[derrGdFx, derrGdFy, 0]])

    return J

def int_error_result(deformFunc, SF, xCoeffs, yCoeffs, cICoeffs, totalSpringLength, outerXCoord, outerYCoord):

    smesh = np.linspace(0, totalSpringLength, globalLen)
    geometryDef = [xCoeffs, yCoeffs, cICoeffs]

    xorg = coord(0,xCoeffs)
    yorg = coord(0,yCoeffs)

    dgds0 = (SF[2]/(n) + (outerYCoord-yorg)*SF[0] - (outerXCoord-xorg)*SF[1])/(E*PPoly_Eval(0, cICoeffs))
    # print("dgds0",dgds0)
    res = fixed_rk4(deformFunc, np.array([0,x0,y0]), smesh, (SF[0], SF[1], geometryDef, dgds0))
    # print("called fixed rk4", end="\r")
    Rinitial = lin.norm([outerXCoord,outerYCoord])
    Rfinal   = lin.norm([res[1,-1],res[2,-1]])

    dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(coord(totalSpringLength,yCoeffs),coord(totalSpringLength,xCoeffs))) # res[2,:] = y(s), res[1,:] = x(s)
    # if np.isnan(dBeta):
    #     # print("dBeta IS NANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
    # if np.isnan(res[0,-1]):
        # print("gamma(L) IS NANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
    # print("dBeta interior:",dBeta/deg2rad)
    # print("full arc length inside inside",totalSpringLength)
    # print(np.arctan2(res[2,-1],res[1,-1])/deg2rad, np.arctan2(yorg,xorg)/deg2rad)
    # print(dBeta)

    err = np.array([Rinitial-Rfinal, res[0,-1]-dBeta])
    # print("err:",err)

    return err, res


def tune_stiffness(stiffnessTarg, dBetaTarg, dragVector, discludeVector):

    torqueTarg = dBetaTarg*stiffnessTarg
    print("initial guess check,", dragVector)
    # treat out of plane thickness, # of arms as constant
    # #PASS DRAG VECTOR
    dragVector0=dc(dragVector)
    dragVectorPrev = dragVector0
    err = 1
    relErr = 1000000000
    convErr = 1
    stepSizeCoeff = 1*10**-7
    j = 0
    resetCounter = 0
    resetStatus  = 0
    procReset = 0
    while relErr>10e-6 and convErr>10e-6:
        SSProfile("Stiffness Tuning").tic()
        print("stiffness refinement iteration:",j)
        # create spring for guess

        geometryDef, smesh = drag_vector_spring(dragVector)
        # print(geometryDef[2],dragVector[7])
        xorg = coord(smesh, geometryDef[0])
        yorg = coord(smesh, geometryDef[1])
        # establish error in stiffness as a result of guess
        relErrPrev = relErr
        totalSpringLength = smesh[-1]
        res, SFG, f = deform_spring_by_torque(torqueTarg, geometryDef, totalSpringLength, deform_ODE)

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
        if f:
            procReset=1
            print("f")
        elif abs(relErr) > abs(relErrPrev):
            dragVector = dc(dragVectorPrev)
            print("error increase", j)
            procReset = 1
        else:
            grad = estimate_grad(dragVector, torqueTarg, stiffness, err, yorg, xorg)
            print(grad)
            Jinv = 1/grad
            for i in range(len(grad)):
                if discludeVector[i]==False:
                    grad[i]=0
                    Jinv[i]=0
                else:
                    grad[i]=grad[i]*discludeVector[i]
                    Jinv[i]=Jinv[i]*discludeVector[i]
                #### WTF IS GOING ON HERE

            stepSize = stepSizeCoeff*lin.norm(grad)
            dragVectorPrev = dc(dragVector)
            # dragVector = dragVector + stepSize*grad
            # print(Jinv)
            dragVector = dragVector + err*Jinv*stepSize
            # print(dragVector-dragVectorPrev)
            for i in range(len(grad)):
                if discludeVector[i]==False:
                    assert(dragVector0[i]==dragVector[i])
        # if procReset == 0:
        #     procReset = violates_bounds(dragVector)
        #     print("violate bounds reset")
        if procReset:
            print("reset---------------------------------------------------reset")
            dragVector = dc(dragVectorPrev)
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
            procReset = 0
            print("this is the ",resetCounter,"th reset")
        else:
            # assert(not violates_bounds(dragVector))
            print("reduction:", stepSizeCoeff)
            # stepSizeCoeff = stepSizeCoeff*1.1
            j+=1
            resetCounter = 0
            resetStatus  = 0
        SSProfile("Stiffness Tuning").toc()
    # assert(not violates_bounds(dragVector))
    for i in range(len(dragVector)):
                if discludeVector[i]==False:
                    # if not dragVector0[i]==dragVector[i]:
                    #     print(dragVector0[i], dragVector[i])
                    assert(dragVector0[i]==dragVector[i])
                    # print(i)
    print("stiffness refinement iterations:", j)
    if convErr <10e-6:
        print("this is as close as you can get with starting guess and gains")
    if relErr <10e-6:
        print("wow we actually found a workable solution")

    return stiffness, res, dragVector, dragVector0

def stiffness_stress_error(stiffnessTarg, torqueTarg, allowStress, dragVector):
    geometryDef, smesh = drag_vector_spring(dragVector)
    # print(geometryDef[2],dragVector[7])
    xorg = coord(smesh, geometryDef[0])
    yorg = coord(smesh, geometryDef[1])
    # establish error in stiffness and max stress as a result of guess
    totalSpringLength = (dragVector[-1])
    res, SF, f = deform_spring_by_torque(torqueTarg, geometryDef, totalSpringLength)
    la, lb   = outer_geometry(smesh, geometryDef, False)

    rn = r_n(smesh, geometryDef[0], geometryDef[1])
    dgds0 = (SF[2]/(n) + yorg*SF[0] - xorg*SF[1])/(E*PPoly_Eval(0, geometryDef[2]))
    Fx  = SF[0]
    Fy  = SF[1]
    Mdim = E*PPoly_Eval(smesh, geometryDef[2])
    dgds = dgds0 + Fy/Mdim*res[1,:] - Fx/Mdim*res[2,:]


    innerStress = (E*(1-rn/(rn-la))*rn*dgds)
    outerStress = (E*(1-rn/(rn-lb))*rn*dgds)
    for i in range(len(smesh)):
        if not (np.isfinite(innerStress[i]) or np.isfinite(outerStress[i])):
            innerStress[i] = E*dgds[i]*la[i]
            outerStress[i] = innerStress[i]
    maxInnerStress = np.max(innerStress)
    maxOuterStress = np.max(outerStress)
    maxStress      = np.max([maxInnerStress, maxOuterStress])

    dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
    stiffness = torqueTarg/dBeta

    errStiffness = abs(stiffness - stiffnessTarg)
    errStress    = abs(maxStress - allowStress)
    return errStiffness, errStress, stiffness, maxStress, res, f

def estimate_grad_stress(dragVector, torqueTarg, stiffnessTarg, stiffness, allowStress, relErr, yorg, xorg):
    grad = np.empty(len(dragVector))
    dragVectorBase = dc(dragVector)
    for i in range(len(dragVector)):
        if i < 4 or i>9:
            # print("length")
            # print("base", dragVectorBase)
            dragVectorFD=dc(dragVectorBase)
            dragVectorFD[i]=dragVectorFD[i]+finiteDifferenceLength
            geometryDefNew, smesh = drag_vector_spring(dragVectorFD)
            totalSpringLength = dragVectorFD[-1]
            errStiffness, errStress, stiffness, maxStress, res, f = stiffness_stress_error(stiffnessTarg, torqueTarg, allowStress, dragVector)
            relErrf = lin.norm(np.array([(errStiffness/stiffnessTarg), (errStress/allowStress)]))
            grad[i]=(relErrf-relErr)/(finiteDifferenceLength)
        if i>3 and i<7:
            # print("angle")
            # print("base", dragVectorBase)
            dragVectorFD=dc(dragVectorBase)
            dragVectorFD[i]=dragVectorFD[i]+finiteDifferenceAngle
            totalSpringLength = dragVectorFD[-1]
            geometryDefNew, smesh = drag_vector_spring(dragVectorFD)
            totalSpringLength = dragVectorFD[-1]
            errStiffness, errStress, stiffness, maxStress, res, f = stiffness_stress_error(stiffnessTarg, torqueTarg, allowStress, dragVector)
            relErrf = lin.norm(np.array([(errStiffness/stiffnessTarg), (errStress/allowStress)]))
            grad[i]=(relErrf-relErr)/(finiteDifferenceAngle)
        if i>6 and i<10:
            # print("Ic")
            # print("base", dragVectorBase)
            dragVectorFD=dc(dragVectorBase)
            dragVectorFD[i]=dragVectorFD[i]+finiteDifferenceCI
            geometryDefNew, smesh = drag_vector_spring(dragVectorFD)
            totalSpringLength = dragVectorFD[-1]
            errStiffness, errStress, stiffness, maxStress, res, f = stiffness_stress_error(stiffnessTarg, torqueTarg, allowStress, dragVector)
            relErrf = lin.norm(np.array([(errStiffness/stiffnessTarg), (errStress/allowStress)]))
            grad[i]=(relErrf-relErr)/(finiteDifferenceCI)

    grad = grad/lin.norm(grad)

    return grad

def stress_stiffness_tuning(stiffnessTarg, dBetaTarg, allowStress, dragVector, discludeVector):
    torqueTarg = dBetaTarg*stiffnessTarg
    print("initial guess check,", dragVector)
    # treat out of plane thickness, # of arms as constant
    # #PASS DRAG VECTOR
    dragVector0=dc(dragVector)
    #set up iteration variables
    err = 1
    relErr = 1
    convErr = 1
    # in general case, set  this to 1
    stepSizeCoeff = 1*10**-7
    j = 0
    resetCounter = 0
    resetStatus  = 0
    procReset = 0
    while relErr>10e-6 and convErr>10e-6: # quit when you either reach the target or can't go any further
        SSProfile("Stress/Stiffness Tuning").tic()
        print("stiffness refinement iteration:",j)
        geometryDef, smesh = drag_vector_spring(dragVector)
        xorg = coord(smesh, geometryDef[0])
        yorg = coord(smesh, geometryDef[1])
        relErrPrev = relErr
        errStiffness, errStress, stiffness, maxStress, res, f = stiffness_stress_error(stiffnessTarg, torqueTarg, allowStress, dragVector)
        relErr = lin.norm(np.array([(errStiffness/stiffnessTarg), (errStress/allowStress)]))
        convErr = abs(relErr-relErrPrev)
        print("stiffness error:",errStiffness)
        print("stress error:", errStress)
        print("stiffness relative error:",relErr)
        print("chenge in rel error", convErr)
        # estimate gradient (pain)
        if f:
            procReset=1
        if abs(relErr) > abs(relErrPrev):
            dragVector = dc(dragVectorPrev)
            procReset = 1
        else:
            grad = estimate_grad_stress(dragVector, torqueTarg, stiffnessTarg, stiffness, allowStress, relErr, yorg, xorg)
            Jinv = 1/grad
            for i in range(len(grad)):
                if discludeVector[i]==False:
                    grad[i]=0
                    Jinv[i]=0
                else:
                    grad[i]=grad[i]*discludeVector[i]
                    Jinv[i]=Jinv[i]*discludeVector[i]
            # determine next set of parameters
            stepSize = stepSizeCoeff*lin.norm(grad)
            dragVectorPrev = dc(dragVector)
            dragVector = dragVector + err*Jinv*stepSize
            for i in range(len(grad)):
                if discludeVector[i]==False:
                    assert(dragVector0[i]==dragVector[i])
        # if procReset == 0:
        #     procReset = violates_bounds(dragVector)
        if procReset:
            print("reset---------------------------------------------------reset")
            dragVector = dc(dragVectorPrev)
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
            procReset=0
            print("this is the ",resetCounter,"th reset")
        else:
            # assert(not violates_bounds(dragVector))
            print("reduction:", stepSizeCoeff)
            # stepSizeCoeff = stepSizeCoeff*1.1
            j+=1
            resetCounter = 0
            resetStatus  = 0
        SSProfile("Stress/Stiffness Tuning").toc()
    # assert(not violates_bounds(dragVector))
    for i in range(len(dragVector)):
                if discludeVector[i]==False:
                    # if not dragVector0[i]==dragVector[i]:
                    #     print(dragVector0[i], dragVector[i])
                    assert(dragVector0[i]==dragVector[i])
                    # print(i)
    print("stiffness refinement iterations:", j)
    if convErr <10e-6:
        print("this is as close as you can get with starting guess and gains")
    if relErr <10e-6:
        print("wow we actually found a workable solution")

    return stiffness, maxStress, res, dragVector, dragVector0

def outer_geometry(smesh, geometryDef, printBool):
    # IF YOU'RE USING THIS TO GET THE OUTER GEOMETRY AT A SINGLE POINT, PASS

    if hasattr(smesh,"__len__"):
        la = np.empty(len(smesh))
        lb = np.empty(len(smesh))
        hLALB = np.empty(len(smesh))
        lABPrev = [0, 0]
        for i in range(len(smesh)):
            lAB = l_a_l_b_rootfinding(smesh[i], lABPrev, geometryDef[0], geometryDef[1], geometryDef[2], printBool)
            # print(AB)
            la[i] = lAB[0]
            lb[i] = lAB[1]
            hLALB[i] = lb[i]+la[i]
            lABPrev = lAB
    else:
        lABPrev = [0,0]
        lAB = l_a_l_b_rootfinding(smesh, lABPrev, geometryDef[0], geometryDef[1], geometryDef[2], printBool)
        la = lAB[0]
        lb = lAB[1]
    return la, lb

def violates_bounds(dragVector):
    truth0 = dragVector[0] < globalInnerRadiusLimit or dragVector[0] > globalOuterRadiusLimit
    truth1 = dragVector[1] < globalInnerRadiusLimit or dragVector[1] > globalOuterRadiusLimit or dragVector[1] < dragVector[0] or dragVector[1] > dragVector[3]
    truth2 = dragVector[2] < globalInnerRadiusLimit or dragVector[2] > globalOuterRadiusLimit or dragVector[2] < dragVector[0] or dragVector[2] > dragVector[3]
    truth3 = dragVector[3] < globalInnerRadiusLimit or dragVector[3] > globalOuterRadiusLimit

    truth4 = False # dragVector[4] > 90*deg2rad
    truth5 = dragVector[5] < dragVector[4] or dragVector[5] > dragVector[6]
    truth6 = False # dragVector[6] > 170*deg2rad

    truth7 = dragVector[7] < dragVector[8] #or dragVector[7] > .05
    truth8 = dragVector[8] < 0 #or dragVector[8] > .025
    truth9 = dragVector[9] < dragVector[8] #or dragVector[9] > .05
    truths = [truth0,truth1,truth2,truth3,truth4,truth5,truth6,truth7,truth8,truth9]
    for value, index in zip(truths, [0,1,2,3,4,5,6,7,8,9]):
        if not value:
            code=index
            print(code)
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
            totalSpringLength = dragVectorFD[-1]
            res, SFG, f = deform_spring_by_torque(torqueTarg, geometryDefNew, totalSpringLength,deform_ODE)
            dBeta = (np.arctan2(res[2,-1],res[1,-1])-np.arctan2(yorg[-1],xorg[-1]))
            stiffnessd = torqueTarg/dBeta
            errf = stiffness - stiffnessd
            grad[i]=(errf-err)/(finiteDifferenceLength)
        if i>3 and i<7:
            # print("angle")
            # print("base", dragVectorBase)
            dragVectorFD=dc(dragVectorBase)
            dragVectorFD[i]=dragVectorFD[i]+finiteDifferenceAngle
            totalSpringLength = dragVectorFD[-1]
            geometryDefNew, smesh = drag_vector_spring(dragVectorFD)
            totalSpringLength = dragVectorFD[-1]
            res, SFG, f = deform_spring_by_torque(torqueTarg, geometryDefNew, totalSpringLength,deform_ODE)
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
            totalSpringLength = dragVectorFD[-1]
            res, SFG, f = deform_spring_by_torque(torqueTarg, geometryDefNew,totalSpringLength,deform_ODE)
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
    smesh = np.linspace(0, dragVectorArg[-1], globalLen)

    return geometryDef, smesh

def colorline(
    x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0),
        linewidth=3, alpha=1.0):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc


def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

def startSW():
    ## Starts Solidworks
    SW_PROCESS_NAME = r'C:/Program Files/SOLIDWORKS Corp/SOLIDWORKS/SLDWORKS.exe'
    sb.Popen(SW_PROCESS_NAME)

def shutSW():
    ## Kills Solidworks
    sb.call('Taskkill /IM SLDWORKS.exe /F')

def connectToSW():
    ## With Solidworks window open, connects to application
    sw = win32com.client.Dispatch("SLDWORKS.Application")
    return sw

def openFile(sw, Path):
    ## With connection established (sw), opens part, assembly, or drawing file
    f = sw.getopendocspec(Path)
    model = sw.opendoc7(f)
    return model

def newFile(sw, Path):
    template = "C:\ProgramData\SolidWorks\SOLIDWORKS 2014\templates\Part.prtdot"
    model = sw.NewDocument(template, 0,0,0)
    model.SaveAs(Path, 0, 2)
    return model

def updatePrt(model):
    ## Rebuilds the active part, assembly, or drawing (model)
    model.EditRebuild3


deg2rad = np.pi/180
n = 2
finiteDifferenceLength = 0.001
finiteDifferenceAngle  = .1*deg2rad
ffForce  = .5
finiteDifferenceTorque = 0.1
finiteDifferenceCI     = 20

straights = []

fullArcLength = 4.8
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

lABPrevOuter = [0,0]