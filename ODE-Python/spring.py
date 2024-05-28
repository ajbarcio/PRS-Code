import numpy as np
import math
import numpy.linalg as lin
import matplotlib.pyplot as plt
from StatProfiler import SSProfile

from utils import numerical_fixed_mesh_diff, fixed_rk4, rk4_step, colorline
from polynomials import Ic_multiPoly, xy_poly, PPoly_Eval

deg2rad = np.pi/180

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
        # print("getting initial guess")
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
        # print("done getting initial guess")
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
            elif lin.norm(err,2) > 10e-10:
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
        self.dBeta    = np.arctan2(self.res[2,-1],self.res[1,-1])-self.betaAngles[-1]
        # print("deflection:",dBeta)
        # Err = diff. in radius, diff between gamma(L) and beta(L)
        err = np.array([Rinitial-Rfinal, self.res[0,-1]-self.dBeta, SF[2]-torqueTarg])
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
        # Deal with args in the most stupid way possible
        Fx    = args[0][0][0][0]
        Fy    = args[0][0][0][1]
        dgds0 = args[0][0][0][2]
        # prepare moment dimensionalizer
        Mdim = self.E*PPoly_Eval(xi, self.IcCoeffs, ranges=self.domains)
        # prepare xi->s space transform and all derivatives
        dxdxi = PPoly_Eval(xi, self.XCoeffs, deriv=1)
        dydxi = PPoly_Eval(xi, self.YCoeffs, deriv=1)

        dxids = d_xi_d_s(xi, self.XCoeffs, self.YCoeffs)
        # print(dxids)
        dxds = dxids*dxdxi
        dyds = dxids*dydxi

        ### DO THE DIFFEQ ###

        LHS = np.empty(3)
        LHS[0] = (dgds0 + Fy/Mdim*(p[1]-self.x0) - Fx/Mdim*(p[2]-self.y0))
        LHS[1] = np.cos(np.arctan2(dyds,dxds)+p[0])
        LHS[2] = np.sin(np.arctan2(dyds,dxds)+p[0])
        # check to make sure the math makes sense
        if xi==0:
            # print(LHS[0],dgds0)
            assert(np.isclose(LHS[0],dgds0,rtol=1e-5))
        # transfer from s space to xi space
        LHS = LHS/dxids
        return LHS

    def l_a_l_b_rootfinding(self, s, lABPrev, printBool=0):
        # get the relevant values at a given point
        rn = r_n(s, self.XCoeffs, self.YCoeffs)
        Ic = PPoly_Eval(s, self.IcCoeffs, ranges = self.domains)
        # define the system to be solved
        def func(x, rn, Ic):
            f1 = (x[0]+x[1])/(np.log((rn+x[1])/(rn-x[0])))-rn
            # f2 = (x[0]+x[1])*(outPlaneThickness*(x[1]-(x[1]+x[0])/2)*rn)-cI
            f2 = self.t*rn*(x[0]+x[1])*(x[1]/2-x[0]/2)-Ic
            return np.array([f1, f2])
        def jac(x, rn, Ic):
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
        # this threshold is arbitrary
        if(lin.norm(x)>lin.norm(lABPrev)*100):
            l = np.cbrt(12*Ic/self.t)/2
            x = [l,l]
        if(printBool):
            print(x0)
            print(x)
            print(iii)
            print(rn)
            print(err)
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

    def calculate_stresses(self):
        dgdxi = numerical_fixed_mesh_diff(self.res[0,:], self.ximesh)
        self.dgds = dgdxi*self.dxids
        self.innerSurfaceStress = np.empty(len(self.ximesh))
        self.outerSurfaceStress = np.empty(len(self.ximesh))
        for i in range(len(self.innerSurfaceStress)):
            if not np.isinf(self.rn[i]):
                self.innerSurfaceStress[i] = abs(self.E*(1-self.rn[i]/self.a[i])*self.rn[i]*self.dgds[i])
                self.outerSurfaceStress[i] = abs(self.E*(1-self.rn[i]/self.b[i])*self.rn[i]*self.dgds[i])
            else:
                self.innerSurfaceStress[i] = self.E*self.dgds[i]*0.5*self.h[i]

        self.normalizedInnerStress =self.innerSurfaceStress/self.designStress
        self.normalizedOuterStress =self.outerSurfaceStress/self.designStress

        self.maxStresses = np.empty(len(self.normalizedInnerStress))
        for i in range(len(self.normalizedInnerStress)):
            if self.normalizedInnerStress[i] > self.normalizedOuterStress[i]:
                self.maxStresses[i] = self.normalizedInnerStress[i]
            else:
                self.maxStresses[i] = self.normalizedOuterStress[i]
        # print(self.maxStresses)
        self.maxStress = np.nanmax([self.innerSurfaceStress, self.outerSurfaceStress])
        print("max stress:", self.maxStress)
        print("des stress:", self.designStress)
        return self.maxStress, self.maxStresses

    def spring_geometry(self, plotBool=1, deformBool=1):

        ## Things to make for the undeformed state:

        # Generate neutral radius path and give it nicely formatted class variables
        self.undeformedNeutralSurface = np.hstack((np.atleast_2d(PPoly_Eval(self.ximesh, self.XCoeffs)).T, np.atleast_2d(PPoly_Eval(self.ximesh, self.YCoeffs)).T))
        # Generate outer and inner surfaces
        self.generate_surfaces() # A and B surface come from here
        # generate centroidal surface
        self.undeformedCentroidalSurface = self.undeformedNeutralSurface+np.hstack((np.atleast_2d(self.ecc*np.sin(self.alpha)).T, np.atleast_2d(self.ecc*np.cos(self.alpha)).T))

        if plotBool:
            plt.figure(1)
            if not deformBool:
                # only plot neutral and centroidal surface if undeformed
                plt.plot(self.undeformedNeutralSurface[:,0],self.undeformedNeutralSurface[:,1])
                plt.plot(self.undeformedCentroidalSurface[:,0],self.undeformedCentroidalSurface[:,1])
            # in any case plot the inner and outer surfaces
            plt.plot(self.undeformedASurface[:,0],self.undeformedASurface[:,1],"--b")
            plt.plot(self.undeformedBSurface[:,0],self.undeformedBSurface[:,1],"--b")
        if deformBool:
            # generate neutral surface after deformation (and create nicely formatted class variables)
            # (some very silly python array handling happens here)
            self.deformedNeutralSurface = np.hstack((np.atleast_2d(self.res[1,:]).T, np.atleast_2d(self.res[2,:]).T))

            self.deformedBSurface = self.deformedNeutralSurface+np.hstack((np.atleast_2d(-self.lb*np.sin(self.alpha+self.res[0,:])).T,np.atleast_2d(self.lb*np.cos(self.alpha+self.res[0,:])).T))
            self.deformedASurface = self.deformedNeutralSurface-np.hstack((np.atleast_2d(-self.la*np.sin(self.alpha+self.res[0,:])).T,np.atleast_2d(self.la*np.cos(self.alpha+self.res[0,:])).T))
            # calculate stress in a deformed beam
            self.calculate_stresses()
            # flag if overstressed
            if self.maxStress > self.designStress:
                # print("oops you broke lmao")
                print("design stress exceeded")
            if plotBool:
                # plot the neutral surface with max stress at each point
                colorline(self.deformedNeutralSurface[:,0],self.deformedNeutralSurface[:,1],self.maxStresses,cmap=plt.get_cmap('rainbow'))
                # plot inner and outer surfaces
                plt.plot(self.deformedASurface[:,0],self.deformedASurface[:,1],"-k")
                plt.plot(self.deformedBSurface[:,0],self.deformedBSurface[:,1],"-k")
                # plot scale
                colorline(np.ones(101)*(np.nanmax(self.deformedASurface[:,0])+.5),np.linspace(0,2,101),np.linspace(0,1,101),cmap=plt.get_cmap('rainbow'),linewidth=10)
        if plotBool:
            outerCircle = plt.Circle([0,0],self.radii[3],color ="k",fill=False)
            innerCircle = plt.Circle([0,0],self.radii[0],color ="k",fill=False)
            fig = plt.gcf()
            ax = fig.gca()
            ax.add_patch(outerCircle)
            ax.add_patch(innerCircle)
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

def d_xi_d_s(ximesh, XCoeffs, YCoeffs):
    dxdxi = PPoly_Eval(ximesh, XCoeffs, deriv=1)
    dydxi = PPoly_Eval(ximesh, YCoeffs, deriv=1)
    dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
    return dxids