import numpy as np
import scipy as scp
import scipy.stats as stat
import scipy.special as spec
import numpy.linalg as lin
import matplotlib.pyplot as plt
import os

from typing import Optional

from modules.StatProfiler import SSProfile
from modules.utils import fixed_rk4, numerical_fixed_mesh_diff, colorline, identify_quadrant, deg2rad
from modules.PATHDEF import Path, RadiallyEndedPolynomial
from modules.CRSCDEF import Crsc, Piecewise_Ic_Control

import modules.materials as materials
from modules.materials import Material, Titanium5, TestMaterial, Maraging300Steel

import json

class Spring:
    def __init__(self, path: Path, crsc: Crsc, material,
                 resolution = 200,
                 torqueCapacity = 3000,
                 name       = None):
        # self.parameters = {key: value for key, value in locals().items() if not key.startswith('__') and key != 'self'}
        self.name = name
        # Parameters and data structures from passed objects
        # Path accessor from path definition
        self.path     = path
        # Geometry accessor
        self.crsc     = crsc
        # Material info
        self.material = material
        # Design intent
        self.torqueCapacity = torqueCapacity
        # Geometry parameterization
        self.fullArcLength = path.arcLen
        # thickness arg
        self.t = crsc.t

        # Mesh information
        self.resl = resolution
        self.len = self.resl+1
        self.step = self.fullArcLength/self.resl
        self.endIndex = self.len-1

        self.ximesh = np.linspace(0,self.fullArcLength,self.len)

        # finite difference information
        self.finiteDifferenceLength = 0.125
        # self.path.arcLen/(2*self.resl)
        self.finiteDifferenceAngle  = 1*deg2rad
        self.finiteDifferenceForce  = 5
        self.finiteDifferenceTorque = 10
        self.finiteDifferenceIc     = 0.0001
        self.finiteDifferenceFactor = 0.001

        self.differenceThreshold = 0.0006

        # Default deforming mode:
        self.deformMode = self.deform_by_torque

        # fuck your memory, extract commonly used attributes from passed
        # definition objects:
        try:
            # material properties:
            self.E = self.material.E
            self.designStress = self.material.yieldStress*.9

            # path properties:
            self.x0 = self.path.startPoint[0]
            self.y0 = self.path.startPoint[1]
            self.xL = self.path.endPoint[0]
            self.yL = self.path.endPoint[1]

            self.momentArmX = self.path.momentArm[0]
            self.momentArmY = self.path.momentArm[1]

            self.n = self.path.n
            if self.path.innerRadius is not None:
                self.innerRadius = self.path.innerRadius
                self.outerRadius = self.path.outerRadius
            else:
                self.innerRadius = lin.norm(self.path.startPoint)
                self.outerRadius = lin.norm(self.path.endPoint)

        except self.Error as e:
            print("path, cross section, or material definitions are missing \
                  \n standard attributes")

        self.parameters = self.get_parameters()

    def get_parameters(self):
        self.parameters = {"path": self.path,
                           "crsc": self.crsc,
                           "material": self.material,
                           "resolution": self.resl,
                           "torqueCapacity": self.torqueCapacity,
                           "name": self.name
                           }
        return self.parameters

    @classmethod
    def from_param(cls, paramFile):
        
        # Access all three paramfiles
        path_suffix = "_path.param"
        crsc_suffix = "_crsc.param"
        sprg_suffix = "_sprg.param"
        
        # Load spring params to dict
        with open(paramFile+sprg_suffix, 'r') as file:
            params = json.load(file)
        # Extract names (strings) of path, crsc, and material objects
        pathType = params.pop("path", None)
        crscType = params.pop("crsc", None)
        matlType = params.pop("material", None)
        # nameNew  = params.pop("name", params["name"]+"1")

        # Look through the globals() env to steal the content of the vars associated with the strings
        # print(globals())
        pathClass = globals()[pathType]
        crscClass = globals()[crscType]
        matlObjct = globals()[matlType]

        # Create path and crsc objects with appropriate params to load to spring object
        pathFromParams = pathClass.from_param(paramFile+path_suffix)
        crscFromParams = crscClass.from_param(paramFile+crsc_suffix, pathFromParams)
        
        # Create a spring with appropriate parameters
        return cls(path=pathFromParams, crsc=crscFromParams, material=matlObjct, **params)

    class Error(Exception):
        pass

    def deform_withTension_ODE(self, xi, deforms, *args):

        M, Fx, Fy, dgds0 = args[0]
        gamma, x, y = deforms

        alpha = self.path.get_alpha(xi)
        Mdim = self.E*self.crsc.get_Ic(xi)
        Area = self.crsc.get_Area(xi)
        Fax = Fx*np.cos(alpha+gamma) + Fy*np.sin(alpha+gamma)

        dxids = self.path.get_dxi_n(xi)

        LHS = np.empty(3)
        # Bending equation remains the same
        LHS[0] = (dgds0 - Fy/Mdim*(deforms[1] -self.x0) +
                                                   Fx/Mdim*(deforms[2]-self.y0))
        # Coordinate equations change
        # Combination of off-axis components and tensive components:
        LHS[1] = np.cos(alpha+gamma)*(1+Fax/(self.E*Area))
        LHS[2] = np.sin(alpha+gamma)*(1+Fax/(self.E*Area))

        # check to make sure the math makes sense (I.C. satisfied)
        if xi==0:
            assert(np.isclose(LHS[0],dgds0,rtol=1e-5))

        # transfer back from s space to xi space
        LHS = LHS/dxids
        return LHS

    def deform_ODE(self, xi, deforms, *args):
        # Deal with args in a way that is hopefully better
        M, Fx, Fy, dgds0 = args[0]
        gamma, x, y = deforms
        # Prepare moment dimensionalizer
        Mdim = self.E*self.crsc.get_Ic(xi)

        dxids = self.path.get_dxi_n(xi)

        LHS = np.empty(3)
        # LHS[0] = (M/self.n - (y-self.y0)*Fx + (x-self.x0)*Fy)/ \
        #               (self.E*self.crsc.get_Ic(xi))
        LHS[0] = (dgds0 - Fy/Mdim*(deforms[1]-self.x0) +
                                                   Fx/Mdim*(deforms[2]-self.y0))
        LHS[1] = np.cos(self.path.get_alpha(xi)+gamma)
        LHS[2] = np.sin(self.path.get_alpha(xi)+gamma)
        
        # check to make sure the math makes sense (I.C. satisfied)
        if xi==0:
            assert(np.isclose(LHS[0],dgds0,rtol=1e-5))

        # transfer back from s space to xi space
        LHS = LHS/dxids
        return LHS

    def forward_integration(self, ODE, SF, torqueTarg, startingIndex=0):

        """
        THIS FUNCTION:
            Evaluates a spring's deflection under a certain loading condition
            by integrating forward across a specified ODE (there is currently
            only one supported)
        """
        Ic0 = self.crsc.get_Ic(0)

        # set up initial condition
        self.dgds0 = (SF[2]/self.n - self.momentArmY*SF[0] + self.momentArmX*SF[1])/ \
                      (self.E*Ic0)
        # self.dgds0 = (SF[2] + self.momentArmY*SF[0] - self.momentArmX*SF[1])/ \
        #               (self.E*Ic0)
        # perform forward integration
        # arguments: ODE, Initial Conditions, Mesh, args
        self.res  = fixed_rk4(ODE, np.array([0,self.x0,self.y0]), self.ximesh,
                                            # (SF[2],
                                            (SF[2], SF[0], SF[1], self.dgds0))
        # calcualte error values (difference in pre-and post-iteration radii)
        #                        (difference in final gamma and final beta)
        Rinitial = lin.norm([self.xL,self.yL])
        Rfinal   = lin.norm([self.res[1,-1],self.res[2,-1]])
        # Calculate dBeta using arccos
        ptInitial = np.array([self.xL, self.yL])
        ptFinal   = np.array([self.res[1,-1],self.res[2,-1]])

        ptInitialNormal = np.array([self.xL, self.yL])/ \
                                          lin.norm(np.array([self.xL, self.yL]))
        ptFinalNormal   = np.array([self.res[1,-1],self.res[2,-1]])/ \
                             lin.norm(np.array([self.res[1,-1],self.res[2,-1]]))

        quadInitial = identify_quadrant(ptInitialNormal)
        quadFinal   = identify_quadrant(ptFinalNormal)

        if quadFinal < quadInitial:
            quad = quadFinal
        else:
            quad = quadInitial

        if quad%2:
            axis=np.array([1,0])
        else:
            axis=np.array([0,1])
        if quad>=3:
            axis=axis*-1

        Angle1 = np.arccos(min(ptInitialNormal.dot(axis),1.0))
        Angle2 = np.arccos(min(ptFinalNormal.dot(axis),1.0))

        # cosAngle = min(ptFinal.dot(ptInitial),1.0)
        self.dBeta    = Angle2-Angle1
        self.dist     = lin.norm(ptFinal-ptInitial)
        # Err = diff. in radius, diff between gamma(L) and beta(L), distance
        #       from target torque (just to make matrix square)
        err = np.array([Rinitial-Rfinal, (self.res[0,-1])-(self.dBeta),
                                                              SF[2]-torqueTarg])
        return err, self.res

    def BC_jacobian(self, ODE, SF, torqueTarg, n=3):

        """
        THIS FUNCTION:
            estimates the jacobian of the forward integration with respect to
            the spatial force vector using a finite difference method
        """
        # print("JACOBIAN SHIT HAPPENS HERE")
        # Jacobian is square because square matrix good
        Jac = np.empty((3,3))
        # Use these values for finite differences (TUNE)
        self.finiteDifferenceTorque=0.1*torqueTarg
        self.finiteDifferenceForce = self.finiteDifferenceTorque/lin.norm(self.path.momentArm)
        # assign each row at a time
        for i in range(n):
            finiteDifference = np.zeros(len(SF))
            # decide based on whether you're using a force or torque which
            # difference to use
            if i < 2:
                finiteDifference[i] = self.finiteDifferenceForce
            else:
                finiteDifference[i] = self.finiteDifferenceTorque
            # evaluate the ODE at a forwards and backwards step
            errBack, resG = self.forward_integration(ODE, SF-finiteDifference, torqueTarg)
            errForw, resG = self.forward_integration(ODE, SF+finiteDifference, torqueTarg)
            # assign row of jacobian
            Jac[0,i]     = (errForw[0]-errBack[0])/(2*lin.norm(finiteDifference))
            Jac[1,i]     = (errForw[1]-errBack[1])/(2*lin.norm(finiteDifference))
            Jac[2,i]     = (errForw[2]-errBack[2])/(2*lin.norm(finiteDifference))
        # this line included in case you ever want to use a partial jacobian
        # Jac = Jac[0:n,0:n]

        return Jac

    def deform_by_torque(self, torqueTarg, ODE,
                         SF=np.array([0,0,0]),
                         breakBool=False):
        # print("called")
        """
        return self.res, self.solnSF, divergeFlag, i

        THIS FUNCTION:
            Solves for a deformation given a certain torque loading condition
        """

        print("Used Direct Torque Method")

        # the err is a two vector, so make it arbitrarily high to enter loop
        err = np.ones(2)*float('inf')
        # 0th iteration
        i = 0
        divergeFlag = 0
        # limit to 100 iterations to converge

        gain=1
        while i < 100:
            # print(i)
            errPrev = err
            # determine boundary condition compliance, estimate jacobian
            err, self.res = self.forward_integration(ODE,SF,torqueTarg)
            # print(lin.norm(err))
            J = self.BC_jacobian(ODE,SF,torqueTarg,n=3)
            # print(J)
            # freak out if it didnt work
            if np.any((SF > 10*10**5)):
                print("I frew up")
                break
            if lin.norm(err)>lin.norm(errPrev):
                # print information on what is happening
                print("torque deform diverging", i)
                # print(J)
                # print(err, errPrev)
                # print(SF)
                divergeFlag = 1
                # If break bool is true, break if you diverge even once
                # usually for debug purposes, I just let it run, and see if
                # it can find its way back to a convergent solution
                if breakBool:
                    break
            # regardless, make a new guess if it did work
            if lin.norm(err,2) > 1e-8:
                # (according to a newton's method)
                # print(SF)
                # print(J)
                # print(lin.pinv(J))
                # print("err:", err)
                # print(lin.norm(err))
                # print(SF)
                SF = SF-lin.pinv(J).dot(err)*gain
                # print(SF)
            # and if the error is small enough to converge, break out
            else:
                break
            i+=1
        # store the final solution in the object
        self.solnSF = SF
        self.solnerr = lin.norm(err)
        self.solnerrVec = err
        return self.res, self.solnSF, divergeFlag, i

    def detailed_deform_regression(self, torqueTarg, ODE, resl, degree):
        steps = np.linspace(0,1*degree,resl+1)*torqueTarg
        print(steps)
        # mags   = []
        # angles = []
        Fxs   = []
        Fys = []
        solnSF = np.array([0,0,steps[0]])
        for step in  steps:
            solnSF[2]=step
            print("trying", step, "inlb")
            print("initial guess:", solnSF)
            gres, solnSF, divergeFlag, i = self.deform_by_torque(step,
                                                                 ODE,
                                                                 SF=solnSF,
                                                                 breakBool=False)
            print("close?", self.solnerr)
            Fxs.append(solnSF[0])
            Fys.append(solnSF[1])
            # mags.append(lin.norm(solnSF))
            # angles.append(np.arctan2(solnSF[1],solnSF[0]))
            # recreateSolnSF = [lin.norm(solnSF)*np.cos(np.arctan2(solnSF[1],solnSF[0])),
            #                     lin.norm(solnSF)*np.sin(np.arctan2(solnSF[1],solnSF[0]))]
            # print("next initial guess should be (real answer):", solnSF)
            # print("next initial guess should be:", recreateSolnSF)
            # self.plot_deform(showBool=False)
            # print(self.solnerr)

        FxRegress = stat.linregress(steps, Fxs)
        FyRegress = stat.linregress(steps[1:], Fys[1:])

        # Plotting for debug
        plt.figure("Fx regression")
        plt.scatter(steps, Fxs)
        plt.plot(steps, FxRegress.slope*np.array(steps)+FxRegress.intercept)
        plt.figure("Fy regression")
        plt.scatter(steps, Fys)
        plt.plot(steps, FyRegress.slope*np.array(steps)+FyRegress.intercept)
        plt.figure("regression graphics")
        self.plot_deform(showBool=False)

    def deform_by_torque_predict_forces(self, torqueTarg, ODE, breakBool=False):
        # print(torqueTarg, ODE)
        # print(type(torqueTarg))

        print("Used Predict Forces Method")

        defBreakBool = breakBool
        steps = [torqueTarg*0.05, torqueTarg*0.2]
        Fxs   = []
        Fys = []
        solnSF = np.array([0,0,steps[0]])
        for step in  steps:
            # print("step:", step)
            solnSF[2] = step
            # print(solnSF)
            gres, solnSF, divergeFlag, i = self.deform_by_torque(step,
                                                                 ODE,
                                                                 SF=solnSF,
                                                                 breakBool=defBreakBool)
            Fxs.append(solnSF[0])
            Fys.append(solnSF[1])
            # print(solnSF)
            # self.plot_deform(showBool=False)

        FxRegress = stat.linregress(steps, Fxs)
        FyRegress = stat.linregress(steps, Fys)

        # upper bound (hopefully)
        FxPredict = FxRegress.slope*(torqueTarg)+FxRegress.intercept
        FyPredict = FyRegress.slope*(torqueTarg)+FyRegress.intercept

        # Plotting for debug
        # plt.figure("magnitudes regression")
        # stepsPlot = steps+[torqueTarg]
        # magsPlot = Fxs+[FxPredict]
        # angsPlot = Fys+[FyPredict]
        # plt.scatter(stepsPlot, magsPlot)
        # plt.plot(stepsPlot, FxRegress.slope*np.array(stepsPlot)+FxRegress.intercept)
        # plt.figure("angles regression")
        # plt.scatter(stepsPlot, angsPlot)
        # plt.plot(stepsPlot, FyRegress.slope*np.array(stepsPlot)+FyRegress.intercept)

        SFPredict = np.array([FxPredict,
                              FyPredict,
                              torqueTarg])
        # print("Predicted Guess:", SFPredict)

        res, solnSF, divergeFlag, i = self.deform_by_torque(torqueTarg, ODE,
                                                            SF=SFPredict,
                                                            breakBool=False) # Not sure about this one

        # print("Solution       :", solnSF)
        # print(self.solnerrVec)
        return res, solnSF, divergeFlag, i

    def deform_by_torque_predict_angle(self, torqueTarg, ODE, breakBool=False):

        print("Used Predict Force Angle Method")

        defBreakBool = breakBool
        steps = [torqueTarg*0.05, torqueTarg*0.2]
        mags   = []
        angles = []
        for step in  steps:
            gres, solnSF, divergeFlag, i = self.deform_by_torque(step,
                                                                 ODE,
                                                                 breakBool=defBreakBool)
            mags.append(lin.norm(solnSF))
            angles.append(np.arctan2(solnSF[1],solnSF[0]))
            # self.plot_deform(showBool=False)

        magRegress = stat.linregress(steps, mags)
        angRegress = stat.linregress(steps, angles)

        # upper bound
        magPredict = magRegress.slope*(torqueTarg)+magRegress.intercept
        angPredict = angRegress.slope*(torqueTarg)+angRegress.intercept

        # Plotting for debug
        # plt.figure("magnitudes regression")
        # print("seriously what the fuck", steps)
        # stepsPlot = steps+[torqueTarg]
        # print("man what the fuck", stepsPlot)
        # magsPlot = mags+[magPredict]
        # angsPlot = angles+[angPredict]
        # plt.scatter(stepsPlot, magsPlot)
        # plt.plot(stepsPlot, magRegress.slope*np.array(stepsPlot)+magRegress.intercept)
        # plt.figure("angles regression")
        # plt.scatter(stepsPlot, angsPlot)
        # plt.plot(stepsPlot, angRegress.slope*np.array(stepsPlot)+angRegress.intercept)

        # print(magPredict, angPredict)
        SFPredict = np.array([magPredict*np.cos(angPredict),
                            magPredict*np.sin(angPredict),
                            torqueTarg])
        # print(SFPredict)

        res, solnSF, divergeFlag, i = self.deform_by_torque(torqueTarg, ODE,
                                                            SF=SFPredict,
                                                            breakBool=True) # Not sure about this one

        return res, solnSF, divergeFlag, i

    def calculate_stresses(self):

        """
        THIS FUNCTION:
            Calcualtes the maximum stress at any point along the beam
        """

        # calculate gamma derivative for stress calcs
        dgdxi = numerical_fixed_mesh_diff(self.res[0,:], self.ximesh)
        dxids = self.path.get_dxi_n(self.ximesh)
        self.dgds = dgdxi*dxids

        rn = self.path.get_rn(self.ximesh)

        self.a = rn - self.la
        self.b = rn + self.lb

        # prepare stress arrays
        self.innerSurfaceStress = np.empty(len(self.ximesh))
        self.outerSurfaceStress = np.empty(len(self.ximesh))
        # calculate stress differently depending on whether or not beam is
        # locally straight
        for i in range(len(self.innerSurfaceStress)):
            if not np.isinf(rn[i]):
                self.innerSurfaceStress[i] =  \
                    abs(self.E*(1-rn[i]/self.a[i])*rn[i]*self.dgds[i])
                self.outerSurfaceStress[i] =  \
                    abs(self.E*(1-rn[i]/self.b[i])*rn[i]*self.dgds[i])
            else:
                self.innerSurfaceStress[i] = self.E*self.dgds[i]*0.5*self.h[i]

        # create arrays normalized to the allowable stress (for plotting)
        self.normalizedInnerStress =self.innerSurfaceStress/self.designStress
        self.normalizedOuterStress =self.outerSurfaceStress/self.designStress

        # record only the max stress between inner and outer at any point
        self.maxStresses = np.empty(len(self.normalizedInnerStress))
        for i in range(len(self.normalizedInnerStress)):
            if self.normalizedInnerStress[i] > self.normalizedOuterStress[i]:
                self.maxStresses[i] = self.normalizedInnerStress[i]
            else:
                self.maxStresses[i] = self.normalizedOuterStress[i]

        # record the maximum stress along whole beam
        self.maxStress = np.nanmax([self.innerSurfaceStress, self.outerSurfaceStress])
        if self.maxStress>self.designStress:
            print("~~~~~~~~~~~ DESIGN STRESS EXCEEDED ~~~~~~~~~~~~")

        return self.maxStress, self.maxStresses

    def plot_spring(self, showBool=False, trans=1, targetAxes=None):

        self.A, self.B =    self.crsc.get_outer_geometry(self.resl)
        # print(self.A)
        # print(self.B)
        self.Sn        =    self.path.get_neutralSurface(self.resl)
        try:
            self.Sc    = self.path.get_centroidalSurface(self.resl)
        except:
            pass

        # plot the principal leg
        if targetAxes is None:
            plt.figure("Graphic Results")
            ax = plt.gca()
        else:
            ax = targetAxes

        ax.plot(self.A[:,0],self.A[:,1],"k", alpha=trans)
        ax.plot(self.B[:,0],self.B[:,1],"k", alpha=trans)
        ax.plot(self.Sn[:,0],self.Sn[:,1],"--b",label="netural", alpha=trans)
        if "self.Sc" in locals():
            ax.plot(self.Sc[:,0],self.Sc[:,1],"--r",label="centroidal", alpha=trans)
        if hasattr(self.path, "pts"):
            ax.plot(self.path.pts[:,0],self.path.pts[:,1])
        ax.axis("equal")
        ax.legend()

        # plot the other legs
        ang = 2*np.pi/self.n
        for j in np.linspace(1,self.n-1,self.n-1):
            th = ang*(j)
            R = np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
            transformedA = R.dot(self.A.T).T
            transformedB = R.dot(self.B.T).T
            transformedSn = R.dot(self.Sn.T).T
            if "self.Sc" in locals():
                transformedSc = R.dot(self.Sc.T).T
                ax.plot(transformedSc[:,0],transformedSc[:,1],"--r", alpha=trans)
            ax.plot(transformedA[:,0],transformedA[:,1],"k", alpha=trans)
            ax.plot(transformedB[:,0],transformedB[:,1],"k", alpha=trans)
            ax.plot(transformedSn[:,0],transformedSn[:,1],"--b", alpha=trans)


        # plot geometry of inner and outer rotor of spring
        outerCircle = plt.Circle([0,0],self.outerRadius,color ="k",fill=False)
        innerCircle = plt.Circle([0,0],self.innerRadius,color ="k",fill=False)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_patch(outerCircle)
        ax.add_patch(innerCircle)

        if showBool:
            plt.show()

        return ax

    def plot_deform(self, showBool=True, targetAxes=None):

        ax = self.plot_spring(showBool=False, trans=0.25, targetAxes=targetAxes)

        alpha = self.path.get_alpha(self.ximesh)
        self.la, self.lb = self.crsc.get_neutralDistances(self.resl)
        self.h = self.crsc.get_Thk(self.ximesh)

        if hasattr(self, 'res'):
            self.deformedNeutralSurface = np.hstack(
                                               (np.atleast_2d(self.res[1,:]).T,
                                                np.atleast_2d(self.res[2,:]).T))
            self.deformedBSurface = self.deformedNeutralSurface+np.hstack(
                    (np.atleast_2d(-self.lb*np.sin(alpha+self.res[0,:])).T,
                     np.atleast_2d(self.lb*np.cos(alpha+self.res[0,:])).T))
            self.deformedASurface = self.deformedNeutralSurface-np.hstack(
                    (np.atleast_2d(-self.la*np.sin(alpha+self.res[0,:])).T,
                     np.atleast_2d(self.la*np.cos(alpha+self.res[0,:])).T))

            self.calculate_stresses()
            colorline(self.deformedBSurface[:,0],self.deformedBSurface[:,1],
                      self.normalizedOuterStress,cmap=plt.get_cmap('rainbow'))
            colorline(self.deformedASurface[:,0],self.deformedASurface[:,1],
                      self.normalizedInnerStress,cmap=plt.get_cmap('rainbow'))

            spot = self.outerRadius*np.sin(45*deg2rad)
            if targetAxes is None:
                ax.text(spot+.25, -(spot+.25), f"deformation: {self.dBeta/deg2rad:.2f}")
                ax.text(spot+.25, -(spot+.375), f"root thk: {self.h[0]:.2f}")
                ax.text(spot+.25, -(spot+.625), f"tip thk: {self.h[-1]:.2f}")
                # ax.text(spot+.25, -(spot+.875), f"min thk: {self.minThk:.2f}")
                # ax.ylim((-(spot+.875)*1.2, self.outerRadius*1.2))

            if showBool:
                plt.show()

        else:
            raise AttributeError("yo wait you didnt fucking twisht it yet")

    def export_surfaces(self):
        if hasattr(self, 'A'):
            A = np.hstack([self.A, np.zeros_like(self.A)])
            B = np.hstack([self.B, np.zeros_like(self.B)])
            fileExtStrs  = [".csv", ".txt", ".sldcrv"]
            surfaces     = ["_A", "_B"]
            surfacesData = [A,B]
            currDir = os.getcwd()
            paths = []
            for ext in fileExtStrs:
                i = 0
                for surf in surfaces:
                    path = os.path.join(
                         os.path.relpath(currDir),"surfaces",self.name+surf+ext)
                    np.savetxt(path, surfacesData[i],fmt='%f',delimiter=',')
                    paths.append(path)
            return paths
        else:
            print("Too early to call; please generate inner and outer surfaces before export")
        return float('nan')

    def export_parameters(self):

        def convert_values(data):
            for key, value in data.items():
                if isinstance(value, np.ndarray):
                    data[key] = value.tolist()  # Convert NumPy array to list
                if isinstance(value, Path) or isinstance(value, Crsc):
                    data[key] = str(value.__class__.__name__)
                if isinstance(value, Material):
                    data[key] = value.name
            return data

        pathParameters = convert_values(self.path.get_parameters())
        crscParameters = convert_values(self.crsc.get_parameters())
        sprgParameters = convert_values(self.get_parameters())

        lists = [pathParameters, crscParameters, sprgParameters]
        names = ['_path', '_crsc', '_sprg']
        currDir = os.getcwd()
        ext = '.param'
        paths = []
        for i in range(len(lists)):
            path = os.path.join(os.path.relpath(currDir),"springs",self.name+names[i]+ext)
            with open(path, 'w') as fp:
                json.dump(lists[i], fp)
            paths.append(path)
        return paths

    def load_values(self, paramFile):
        with open(paramFile, 'r') as parameterFile:
            parameters = json.load(parameterFile)
        for key, value in parameters:
            setattr(self, key, value)
        # IN PROGRESS

class Optimized_Spring(Spring):
    def __init__(self, pathDef: Path, material: materials.Material,
                 threshold = 5,
                 transition = 4,
                 t = 0.375,
                 resolution=200,
                 torqueCapacity=3000,
                 stiffness=52000,
                 useFatigue=False,
                 name = None):
        # Name the spring
        self.name = name
        # "inherit" the path definition
        self.path = pathDef
        # Design intent
        self.torqueCapacity   = torqueCapacity
        self.desStiffness = stiffness

        # Assign Material
        self.material = material
        # Material info
        self.E = material.E
        self.yieldStress = material.yieldStress
        # Stress constraints
        if useFatigue:
            material.reset_designStress()
        self.designStress = material.designStress

        # Geometry parameterization
        self.fullArcLength = self.path.arcLen
        self.t = t

        # Mesh information
        self.resl = resolution
        self.len = self.resl+1
        self.step = self.fullArcLength/self.resl
        self.endIndex = self.len-1

        self.ximesh = np.linspace(0,self.fullArcLength,self.len)

        # Extract properties from passed objects
        try:
            # path properties:
            self.x0 = self.path.startPoint[0]
            self.y0 = self.path.startPoint[1]
            self.xL = self.path.endPoint[0]
            self.yL = self.path.endPoint[1]

            self.momentArmX = self.path.momentArm[0]
            self.momentArmY = self.path.momentArm[1]

            self.n = self.path.n
            if self.path.innerRadius is not None:
                self.innerRadius = self.path.innerRadius
                self.outerRadius = self.path.outerRadius
            else:
                self.innerRadius = lin.norm(self.path.startPoint)
                self.outerRadius = lin.norm(self.path.endPoint)

        except self.Error as e:
            print("path definition is missing standard attributes")

        # Weightings for straight vs curved behavior
        self.weightingsPrepared = 0

        # default values
        self.threshold  = threshold # (in)
        self.transition = transition # Bigger transition is sharper

        self.stressData = []
        self.xiData = []
        self.weightData = []

    def prepare_weightings(self):
        self.lambda1 = lambda x : self.threshold**self.transition/(x**self.transition+self.threshold**self.transition)
        self.lambda2 = lambda x : 1-self.lambda1(x)

        # self.lambda2 = lambda x : 1/np.pi*(np.arctan(transition*np.pi*(x-threshold))
        #                                    -np.arctan(transition*np.pi*(-threshold)))
        # self.lambda1 = lambda x : -self.lambda2(x) + 1

        self.switching = lambda x : (-1)**x/2+1/2

        self.weightingsPrepared = 1

    def tune_weightings(self):
        err = 1000
        i = 0
        gain = 1
        while err > 10 and i < 100:
            #
            # print("ONE LOOP")
            self.stressData = []
            self.weightingsPrepared=0
            threshold = self.threshold
            transition = self.transition # load current
            self.forward_integration(self.weighted_ODE, np.array([0,0,self.torqueCapacity]), self.torqueCapacity)
            stressArray = np.array(self.stressData)
            errVec = stressArray-np.ones_like(stressArray)*self.designStress
            err = lin.norm(errVec)
            # Finite difference:
            # Forward
            self.weightingsPrepared=0
            self.threshold = threshold+0.01
            self.transition = transition
            self.forward_integration(self.weighted_ODE, np.array([0,0,self.torqueCapacity]), self.torqueCapacity)
            stressArray = np.array(self.stressData)
            errVec = stressArray-np.ones_like(stressArray)*self.designStress
            errf = lin.norm(errVec)
            # Backward
            self.weightingsPrepared=0
            self.threshold = threshold-0.01
            self.transition = transition
            self.forward_integration(self.weighted_ODE, np.array([0,0,self.torqueCapacity]), self.torqueCapacity)
            stressArray = np.array(self.stressData)
            errVec = stressArray-np.ones_like(stressArray)*self.designStress
            errb = lin.norm(errVec)

            derrdthresh = (errf-errb)/0.02

            # Forward
            self.weightingsPrepared=0
            self.threshold = threshold
            self.transition = transition+0.01
            self.forward_integration(self.weighted_ODE, np.array([0,0,self.torqueCapacity]), self.torqueCapacity)
            stressArray = np.array(self.stressData)
            errVec = stressArray-np.ones_like(stressArray)*self.designStress
            errf = lin.norm(errVec)
            # Backward
            self.weightingsPrepared=0
            self.threshold = threshold
            self.transition = transition-0.01
            self.forward_integration(self.weighted_ODE, np.array([0,0,self.torqueCapacity]), self.torqueCapacity)
            stressArray = np.array(self.stressData)
            errVec = stressArray-np.ones_like(stressArray)*self.designStress
            errb = lin.norm(errVec)

            derrdtrans = (errf-errb)/0.02

            self.threshold = threshold - err/derrdthresh*gain
            self.transition = transition - err/derrdthresh*gain

            print(self.threshold, self.transition, err)
        print(i)
        print(err)
        self.threshold = threshold
        self.transition = transition
        return self.threshold, self.transition

    def weighted_ODE(self, xi, states, loads):
        # print("----------------------------------")
        # print(states)
        # Repackage states
        gamma, x, y, la, lb = states
        # Repackage loads (IC)
        # print(loads)
        Fx, Fy, dgds0 = loads
        # Get distortion derivative
        dxids = self.path.get_dxi_n(xi)

        alpha = self.path.get_alpha(xi)
        # treating path as known (nominal)
        rn = self.path.get_rn(xi)
        drnds = self.path.get_drn(xi)

        # Prepare moment dimensionalizer
        if not lb==la:
            Ic = self.t*rn/2*(lb**2-la**2)
        else:
            Ic = 1/12*self.t*(2*la)**3
        Mdim = self.E*Ic

        # T0 = dgds0*Mdim

        LHS = np.empty(5)
        # Solve deformation ODE (one step)
        LHS[0] = (dgds0 + Fy/Mdim*(x-self.x0) - Fx/Mdim*(y-self.y0))
        LHS[1] = np.cos(alpha+gamma)
        LHS[2] = np.sin(alpha+gamma)

        # check to make sure the math makes sense (I.C. satisfied)
        if xi==0:
            try:
                assert(np.isclose(LHS[0],dgds0,rtol=1e-5))
            except:
                print(loads[0])
                print(LHS[0], dgds0)
                assert(False)

        # Translate LHS of deformation ODE for use in geo ODE
        dgammads = LHS[0]
        dxds     = LHS[1]
        dyds     = LHS[2]

        # Prepare forcing function F* and its derivative
        FStar    = Fy*(x-self.x0)-Fx*(y-self.y0)
        dFStards = Fy*dxds - Fx*dyds

        # Prepare weighting matrix W
        if not self.weightingsPrepared:
            self.prepare_weightings()
        lam1 = self.lambda1(abs(rn))
        lam2 = self.lambda2(abs(rn))
        self.weightData.append(lam1)

        W = np.array([[lam1, 0, 0, 0],[0, lam1, 0, 0],[0, 0, lam2, 0],[0, 0, 0, lam2]])
        # print(rn, lam1, lam2)

        # Determine stress case
        innerStressFactor = abs(la/(rn-la))
        outerStressFactor = abs(lb/(rn+lb))
        dominantStressFactor = max(innerStressFactor, outerStressFactor)

        if dominantStressFactor==innerStressFactor:
            dDen  = rn-la
            dSide = la
            phi   = 0
            # print("a case")
        elif dominantStressFactor==outerStressFactor:
            dDen  = rn+lb
            dSide = lb
            phi   = 1
            # print("b case")
        else:
            print("OH MY GOD ~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            assert(False)

        if np.isfinite(rn):
            stress = abs(dgammads*self.E*rn*dominantStressFactor)
        else:
            # stress = 1/(dgds0*la)
            stress = abs(dgammads*self.E*la)
        # print("stress:", stress)

        # print("error:", stress-self.designStress)
        self.stressData.append(stress)
        self.xiData.append(xi)

        # Just give us nice, short symbols to use
        # sigma = self.designStress
        # TEST: Try with actual updating stress
        sigma = stress
        t = self.t

        # Prepare full system:
        stressConstrDenominator = t*(dDen*sigma - dgds0*dSide)**2
        stressConstrNumerator   = sigma*dSide*(FStar*drnds-dFStards*dDen)+dSide**2*dFStards*dgds0

        if np.isinf(rn):
            # print("should be finite:", stressConstrDenominator)
            QFrac = FStar/(stressConstrDenominator)
            if FStar==0:
                stressConstrNumerator = 0 + dSide**2*dFStards*dgds0
        else:
            QFrac = FStar*sigma*rn/(stressConstrDenominator)

        Q = np.array([[  la/(rn*(rn-la)),                 -lb/(rn*(rn+lb))               ], \
                      [-(la+(QFrac)*self.switching(phi)),  lb+(QFrac)*self.switching(phi)]])

        P = np.array([[1,0],[-1,1]])
        QP = np.vstack((Q,P))

        qStar = np.array([[drnds*(1/(rn-la)-1/(rn+lb)-(la+lb)/(rn**2))],
                          [stressConstrNumerator/stressConstrDenominator]])
        pStar = np.array([[3/4*dFStards/(la*t*sigma)],[0]])
        qStarPStar = np.vstack((qStar,pStar))
        # print(drnds, rn)
        # print(dFStards, FStar)
        # print(qStarPStar)

        # Solve geometry portion of ODE (one step)
        if rn==float('inf'):
            # print(states)
            # print(la, lb, QFrac)
            # print(Q, P)
            # print(qStar, pStar)
            # print("I'm at either end of the spring")
            pass
        qdot = lin.pinv(W.dot(QP)).dot(W).dot(qStarPStar)
        # print(qdot)
        LHS[3] = qdot[0]
        LHS[4] = qdot[1]

        # print(LHS)

        # transfer back from s space to xi space
        LHS = LHS/dxids

        return LHS

    def forward_integration(self, ODE, SF, torqueTarg):
        """
        THIS FUNCTION:
            Evaluates a spring's deflection under a certain loading condition
            by integrating forward across a specified ODE (there is currently
            only one supported)
        """
        ## Assume that rn0 = infinity
        rn0 = self.path.get_rn(0)
        # print(rn0)
        if not rn0==float('inf'):
            print("Not supported for optimization")
            return 0
        else:
            # print(self.designStress)
            h0 = np.sqrt(6*SF[2]/self.n/(self.t*self.designStress))
            Ic0 = self.t*(h0)**3/12
            print(h0)
            la0=h0/2
            lb0=la0
            # print("init thickness:", la0)
            print(Ic0)

        # set up initial condition
        self.dgds0 = (SF[2]/self.n + self.momentArmY*SF[0] - self.momentArmX*SF[1])/ \
                      (self.E*Ic0)
        # perform forward integration
        # arguments: ODE, Initial Conditions, Mesh, args
        self.res  = fixed_rk4(ODE, np.array([0,self.x0,self.y0,la0,lb0]), self.ximesh, np.array([SF[0], SF[1], self.dgds0]))
        # calcualte error values (difference in pre-and post-iteration radii)
        #                        (difference in final gamma and final beta)
        Rinitial = lin.norm([self.xL,self.yL])
        Rfinal   = lin.norm([self.res[1,-1],self.res[2,-1]])
        # Calculate dBeta using arccos
        ptInitial = np.array([self.xL, self.yL])/ \
                                          lin.norm(np.array([self.xL, self.yL]))
        ptFinal   = np.array([self.res[1,-1],self.res[2,-1]])/ \
                             lin.norm(np.array([self.res[1,-1],self.res[2,-1]]))

        quadInitial = identify_quadrant(ptInitial)
        quadFinal   = identify_quadrant(ptFinal)

        if quadFinal < quadInitial:
            quad = quadFinal
        else:
            quad = quadInitial

        if quad%2:
            axis=np.array([1,0])
        else:
            axis=np.array([0,1])
        if quad>=3:
            axis=axis*-1

        Angle1 = np.arccos(min(ptInitial.dot(axis),1.0))
        Angle2 = np.arccos(min(ptFinal.dot(axis),1.0))

        # cosAngle = min(ptFinal.dot(ptInitial),1.0)
        self.dBeta    = Angle2-Angle1
        # Err = diff. in radius, diff between gamma(L) and beta(L), distance
        #       from target torque (just to make matrix square)
        err = np.array([Rinitial-Rfinal, (self.res[0,-1])-(self.dBeta),
                                                              SF[2]-torqueTarg])
        return err, self.res


def determineFastestSolver(spring: Spring, torqueGain=1):
    testDeformLevel = spring.torqueCapacity*torqueGain

    testFailed = np.zeros(3)
    testTime = np.zeros(3)
    testModes = [spring.deform_by_torque, spring.deform_by_torque_predict_forces, spring.deform_by_torque_predict_angle]
    testNames = ["nonlinear", "predict forces", "predict angles"]

    for j in range(len(testNames)):
        # Test each mode of solving for full-torque deformation:
        SSProfile(testNames[j]).tic()
        res, solnSF, divergeFlag, i = testModes[j](testDeformLevel,spring.deform_ODE)
        err, res = spring.forward_integration(spring.deform_ODE, solnSF, testDeformLevel)
        # print(err)
        # print(testNames[j])
        SSProfile(testNames[j]).toc()
        if not np.isclose(lin.norm(err),0):
            testFailed[j]=1
            testTime[j]=1000
        else:
            testTime[j]=SSProfile(testNames[j]).agg

    fastest = testTime.argmin()
    spring.deformMode = testModes[fastest]

    return str(testModes[fastest])