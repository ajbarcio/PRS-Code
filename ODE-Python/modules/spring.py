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
from modules.PATHDEF import Path
from modules.CRSCDEF import Crsc

import modules.materials as materials

class Optimized_Spring:
    def __init__(self, pathDef, material: materials.Material, 
                 resolution=200, 
                 torqueCapacity=3000, 
                 stiffness=52000, 
                 useFatigue=False,
                 name = None):
        # Name the spring
        self.name = name
        # "inherit" the path definition
        self.path=pathDef
        # Assign Material
        self.material = material
        # Design intent
        self.fullTorque   = torqueCapacity
        self.desStiffness = stiffness
        # Material
        self.E = material.E
        self.yieldStress = material.yieldStress
        # Stress constraints
        if useFatigue:
            material.reset_designStress()
        self.designStress = material.designStress

        self.fullLength = pathDef

class Spring:
    def __init__(self, path: Path, crsc: Crsc, material,
                 resolution = 200,
                 torqueCapacity = 3000,
                 name       = None):

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
        self.resl = resolution
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

    class Error(Exception):
        pass

    def deform_ODE(self, xi, deforms, *args):
        # Deal with args in a way that is hopefully better

        Fx, Fy, dgds0 = args[0]
        # Prepare moment dimensionalizer
        Mdim = self.E*self.crsc.get_Ic(xi)

        dxids = self.path.get_dxi_n(xi)

        dxds = self.path.get_dxdy_n(xi, 'x')*dxids
        dyds = self.path.get_dxdy_n(xi, 'y')*dxids

        LHS = np.empty(3)

        LHS[0] = (dgds0 + Fy/Mdim*(deforms[1]-self.x0) -
                                                   Fx/Mdim*(deforms[2]-self.y0))
        LHS[1] = np.cos(np.arctan2(dyds,dxds)+deforms[0])
        LHS[2] = np.sin(np.arctan2(dyds,dxds)+deforms[0])
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
        self.dgds0 = (SF[2]/self.n + self.momentArmY*SF[0] - self.momentArmX*SF[1])/ \
                      (self.E*Ic0)
        # perform forward integration
        # arguments: ODE, Initial Conditions, Mesh, args
        self.res  = fixed_rk4(ODE, np.array([0,self.x0,self.y0]), self.ximesh,
                                            # (SF[2], 
                                            (SF[0], SF[1], self.dgds0))
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

    def BC_jacobian(self, ODE, SF, torqueTarg, n=3):

        """
        THIS FUNCTION:
            estimates the jacobian of the forward integration with respect to
            the spatial force vector using a finite difference method
        """
        # print("JACOBIAN SHIT HAPPENS HERE")
        # Jacobian is square because square matrix good
        Jac = np.empty((3,3))
        # assign each row at a time
        self.finiteDifferenceTorque=0.1*torqueTarg
        self.finiteDifferenceForce = self.finiteDifferenceTorque/lin.norm(self.path.momentArm)
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

        """
        return self.res, self.solnSF, divergeFlag, i

        THIS FUNCTION:
            Solves for a deformation given a certain torque loading condition
        """

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

    def plot_spring(self, showBool=False, trans=1):

        self.A, self.B =    self.crsc.get_outer_geometry(self.resl)
        # print(self.A)
        # print(self.B)
        self.Sn        =    self.path.get_neutralSurface(self.resl)
        try:
            self.Sc    = self.path.get_centroidalSurface(self.resl)
        except:
            pass

        # plot the principal leg
        plt.figure("Graphic Results")
        plt.plot(self.A[:,0],self.A[:,1],"k", alpha=trans)
        plt.plot(self.B[:,0],self.B[:,1],"k", alpha=trans)
        plt.plot(self.Sn[:,0],self.Sn[:,1],"--b",label="netural", alpha=trans)
        if "self.Sc" in locals():
            plt.plot(self.Sc[:,0],self.Sc[:,1],"--r",label="centroidal", alpha=trans)
        # plt.plot(self.path.pts[:,0],self.path.pts[:,1])
        plt.axis("equal")
        plt.legend()

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
                plt.plot(transformedSc[:,0],transformedSc[:,1],"--r", alpha=trans)
            plt.plot(transformedA[:,0],transformedA[:,1],"k", alpha=trans)
            plt.plot(transformedB[:,0],transformedB[:,1],"k", alpha=trans)
            plt.plot(transformedSn[:,0],transformedSn[:,1],"--b", alpha=trans)
            

        # plot geometry of inner and outer rotor of spring
        outerCircle = plt.Circle([0,0],self.outerRadius,color ="k",fill=False)
        innerCircle = plt.Circle([0,0],self.innerRadius,color ="k",fill=False)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_patch(outerCircle)
        ax.add_patch(innerCircle)

        if showBool:
            plt.show()

    def plot_deform(self, showBool=True):

        self.plot_spring(showBool=False, trans=0.25)

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

            plt.text(spot+.25, -(spot+.25), f"deformation: {self.dBeta/deg2rad:.2f}")

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
            for ext in fileExtStrs:
                i = 0
                for surf in surfaces:
                    path = os.path.join(
                         os.path.relpath(currDir),"surfaces",self.name+surf+ext)
                    np.savetxt(path, surfacesData[i],fmt='%f',delimiter=',')
            #         i+=1
            # if platform.system()=='Windows':
            #     # Use the windows format for path names (ewwwwww)
            #     np.savetxt(".\\surfaces\\"+self.name+"_A.csv", A,
            #             fmt='%f',delimiter=',')
            #     np.savetxt(".\\surfaces\\"+self.name+"_Bcsv", B,
            #             fmt='%f',delimiter=',').
            #     np.savetxt(".\\surfaces\\"+self.name+"_A.txt", A,
            #             fmt='%f',delimiter=',')
            #     np.savetxt(".\\surfaces\\"+self.name+"_B.txt", B,
            #             fmt='%f',delimiter=',')
            #     np.savetxt(".\\surfaces\\"+self.name+"_A.sldcrv", A,
            #             fmt='%f',delimiter=',')
            #     np.savetxt(".\\surfaces\\"+self.name+"_B.sldcrv", B,
            #             fmt='%f',delimiter=',')
            # elif platform.system()=='Linux':
            #     # Use the linux format for path names (yay!)
            #     np.savetxt("./surfaces/"+self.name+"_A.csv", A,
            #             fmt='%f',delimiter=',')
            #     np.savetxt("./surfaces/"+self.name+"_B.csv", B,
            #             fmt='%f',delimiter=',')
            #     np.savetxt("./surfaces/"+self.name+"_A.txt", A,
            #             fmt='%f',delimiter=',')
            #     np.savetxt("./surfaces/"+self.name+"_B.txt", B,
            #             fmt='%f',delimiter=',')
            #     np.savetxt("./surfaces/"+self.name+"_A.sldcrv", A,
            #             fmt='%f',delimiter=',')
            #     np.savetxt("./surfaces/"+self.name+"_B.sldcrv", B,
            #             fmt='%f',delimiter=',')
        else:
            print("Too early to call; please generate inner and outer surfaces before export")

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