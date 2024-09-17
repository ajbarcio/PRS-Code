import numpy as np
import scipy as scp
import scipy.stats as stat
import numpy.linalg as lin
import matplotlib.pyplot as plt
import os

from materials import Maraging300Steel
from utils import fixed_rk4, numerical_fixed_mesh_diff, colorline

deg2rad = np.pi/180

class Spring:
    def __init__(self, geoDef, material,
                 resolution = 200,
                 name       = None):

        self.name=name
        # Parameters and data structures from passed objects
        # Path accessor from path definition
        self.path     = geoDef.path
        # Geometry accessor
        self.crsc     = geoDef
        # Material info
        self.material = material

        self.fullArcLength = self.crsc.fullParamLength
        self.resl = resolution

        # thickness arg
        self.t = self.crsc.t

        # Mesh information
        self.resl = resolution
        self.len = self.resl+1
        self.step = self.fullArcLength/self.resl
        self.endIndex = self.len-1

        self.ximesh = np.linspace(0,self.fullArcLength,self.len)

        ### SILLY THINGS FOR DEBUGGING GEO ODE ###
        self.singularityCounter  = 0
        self.zeroSubCounter = 0
        # -1 has a local minima for error correction
        self.geoFeedbackGain     = -1
        self.geoDerivativeGain   = -0.1
        self.corr_prev = np.array([0,0])
        self.rnErr = []
        self.IcErr = []
        self.thickXi = []

        # finite difference information
        self.finiteDifferenceLength = self.path.fullParamLength/(2*self.resl)
        self.finiteDifferenceAngle  = .1*deg2rad
        self.finiteDifferenceForce  = 0.1
        self.finiteDifferenceTorque = 0.5
        self.finiteDifferenceIc     = 0.00001
        self.finiteDifferenceFactor = 0.001

        # fuck your memory, extract commonly used attributes from passed
        # definition objects:
        try:
            # material properties:
            self.E = self.material.E
            self.designStress = self.material.yieldStress*.9

            # path properties:
            self.x0 = self.path.x0
            self.y0 = self.path.y0
            self.xL = self.path.xL
            self.yL = self.path.yL

            self.momentArmX = self.path.momentArmX
            self.momentArmY = self.path.momentArmY

            self.n = self.path.n

            self.innerRadius = self.path.innerRadius
            self.outerRadius = self.path.outerRadius

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

        LHS[0] = (dgds0 + Fy/Mdim*(deforms[1]-self.x0) -  \
                                                   Fx/Mdim*(deforms[2]-self.y0))
        LHS[1] = np.cos(np.arctan2(dyds,dxds)+deforms[0])
        LHS[2] = np.sin(np.arctan2(dyds,dxds)+deforms[0])
        # check to make sure the math makes sense (I.C. satisfied)
        if xi==0:
            assert(np.isclose(LHS[0],dgds0,rtol=1e-5))

        # transfer back from s space to xi space
        LHS = LHS/dxids
        return LHS

    def geo_ODE(self, xi, q):

        # treating path as known (nominal)
        # treating IC   as known (preliminary)
        rn = self.path.get_rn(xi)
        Ic = self.crsc.get_Ic(xi)

        dIcdxi = self.crsc.get_dIc(xi)
        drndxi = self.path.get_drn(xi)
        # # transform to arc length space
        dxids = self.path.get_dxi_n(xi)
        dIcds = dIcdxi*dxids
        drnds = drndxi*dxids

        # Version without transform

        # dIcds = dIcdxi
        # drnds = drndxi
        def jac(x, rn): # IC doesn't happen to occur in the derivatives
            return np.array([[1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn-x[0])), \
                            1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn+x[1]))], \
                            [-rn*self.t*x[0], rn*self.t*x[1]]])
        # Adjust with feedback gain
        # (Single step try first)
        rnCheck = (q[0]+q[1])/(np.log((rn+q[1])/(rn-q[0])))
        IcCheck = -self.t*rn*(q[0]+q[1])*(q[0]-q[1])/2
        if not xi == 0:
            err  = [rnCheck-rn, IcCheck-Ic]
            self.rnErr.append(err[0])
            self.IcErr.append(err[1])
            self.thickXi.append(xi)
            corr = lin.inv(jac(q, rn)).dot(np.array(err))
            # print(err)
            # print(corr)
            # geoFeedbackGain set to a negative number
            q += corr*self.geoFeedbackGain

        # q = prime     states (la, lb)
        # p = secondary states (Ic, rn)
        # Q = assembled prime     quasi-state vectors
        # P = assembled secondary quasi-state vectors

        Q = np.array([[-q[0],            q[1]],
                      [1/(rn-q[0])-1/rn, 1/(rn+q[1])-1/rn]])
        P = np.array([[1/(rn*self.t),    -Ic/(rn**2*self.t)],
                      [0,                1/(rn-q[0])-1/(rn+q[1])-(q[1]+q[0])/rn**2]])
        p_dot = np.array([[dIcds],[drnds]])

        # This is the diffeq step
        try:
            q_dot = lin.inv(Q).dot(P).dot(p_dot)
        # There is a chance the matrix is singular, in that case:
        except lin.LinAlgError:
            # Catch exception
            print("singular exception at", xi)
            # Take pseudoinverse instead
            q_dot = lin.pinv(Q).dot(P).dot(p_dot)
            self.singularityCounter+=1
        # If we have a bogus answer
        if np.invert(np.isfinite(q_dot)).any(): # or lin.norm(q_dot) > 1:
            # Just substitute that particular step with a numerical
            # approximation of the derivative
            print("bogus answer", xi, lin.norm(q_dot))
            # xi_next      = xi+self.finiteDifferenceLength
            # lalb_current = self.crsc.get_lalb(xi)
            # lalb_next    = self.crsc.get_lalb(xi_next)
            # q1_dot = numerical_fixed_mesh_diff(np.array([lalb_current[0], lalb_next[0]]), np.array([xi, xi_next]))
            # q2_dot = numerical_fixed_mesh_diff(np.array([lalb_current[1], lalb_next[1]]), np.array([xi, xi_next]))
            # q_dot = np.array([q1_dot[0], q2_dot[0]])
            self.zeroSubCounter+=1
            
            # Just call it zero 5head
            q_dot = np.array([0, 0])
            # print("q_dot", q_dot)
            # print("p_dot", p_dot)

        # Back to xi space
        q_dot = q_dot/dxids
        return q_dot.flatten()

    def forward_integration(self, ODE, SF, torqueTarg):


        """
        THIS FUNCTION:
            Evaluates a spring's deflection under a certain loading condition
            by integrating forward across a specified ODE (there is currently
            only one supported)
        """

        # set up initial condition
        dgds0 = (SF[2]/self.n + self.momentArmY*SF[0] - self.momentArmX*SF[1])/ \
                      (self.E*self.crsc.get_Ic(0))

        # perform forward integration
        self.res  = fixed_rk4(ODE, np.array([0,self.x0,self.y0]), self.ximesh,
                                                          (SF[0], SF[1], dgds0))
        # calcualte error values (difference in pre-and post-iteration radii)
        #                        (difference in final gamma and final beta)
        Rinitial = lin.norm([self.xL,self.yL])
        Rfinal   = lin.norm([self.res[1,-1],self.res[2,-1]])
        # Calculate dBeta using arccos
        ptInitial = np.array([self.xL, self.yL])/ \
                                          lin.norm(np.array([self.xL, self.yL]))
        ptFinal   = np.array([self.res[1,-1],self.res[2,-1]])/ \
                             lin.norm(np.array([self.res[1,-1],self.res[2,-1]]))

        cosAngle = min(ptFinal.dot(ptInitial),1.0)
        self.dBeta    = np.arccos(cosAngle)*np.sign(torqueTarg)
        # Err = diff. in radius, diff between gamma(L) and beta(L), distance
        #       from target torque (just to make matrix square)
        err = np.array([Rinitial-Rfinal, abs(self.res[0,-1])-abs(self.dBeta),
                                                              SF[2]-torqueTarg])
        return err, self.res

    def BC_jacobian(self, ODE, SF, torqueTarg, n=3):

        """
        THIS FUNCTION:
            estimates the jacobian of the forward integration with respect to
            the spatial force vector using a finite difference method
        """

        # Jacobian is square because square matrix good
        Jac = np.empty((n,n))
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
        Jac = Jac[0:n,0:n]

        return Jac

    def deform_by_torque(self, torqueTarg, ODE,
                         SF=np.array([0,0,0]),
                         breakBool=False):

        """
        THIS FUNCTION:
            Solves for a deformation given a certain torque loading condition
        """

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
            J = self.BC_jacobian(ODE,SF,torqueTarg)
            # print(J)
            # freak out if it didnt work
            if lin.norm(err)>lin.norm(errPrev):
                # print information on what is happening
                print("torque deform diverging", i)
                print(J)
                print(err, errPrev)
                print(SF)
                divergeFlag = 1
                # If break bool is true, break if you diverge even once
                # usually for debug purposes, I just let it run, and see if
                # it can find its way back to a convergent solution
                if breakBool:
                    break
            # make a new guess if it did work
            elif lin.norm(err,2) > 10e-10:
                # (according to a newton's method)
                SF = SF-lin.pinv(J).dot(err)
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
            print("trying", step, "inlb")
            print("initial guess:", solnSF)
            gres, solnSF, divergeFlag, i = self.deform_by_torque(step, 
                                                                 ODE,
                                                                 SF=solnSF,
                                                                 breakBool=False)
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

    def deform_by_torque_predict_forces(self, torqueTarg, ODE, breakBool=False):
        defBreakBool = breakBool
        steps = [torqueTarg*0.05, torqueTarg*0.2]
        Fxs   = []
        Fys = []
        for step in  steps:
            gres, solnSF, divergeFlag, i = self.deform_by_torque(step, 
                                                                 ODE,
                                                                 breakBool=defBreakBool)
            Fxs.append(solnSF[0])
            Fys.append(solnSF[1])
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
        print("Predicted Guess:", SFPredict)

        res, solnSF, divergeFlag, i = self.deform_by_torque(torqueTarg, ODE, 
                                                            SF=SFPredict,
                                                            breakBool=True) # Not sure about this one

        print("Solution       :", solnSF)
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
            self.plot_deform(showBool=False)
        
        magRegress = stat.linregress(steps, mags)
        angRegress = stat.linregress(steps, angles)

        # upper bound
        magPredict = magRegress.slope*(torqueTarg)+magRegress.intercept
        angPredict = angRegress.slope*(torqueTarg)+angRegress.intercept

        # Plotting for debug
        plt.figure("magnitudes regression")
        print("seriously what the fuck", steps)
        stepsPlot = steps+[torqueTarg]
        print("man what the fuck", stepsPlot)
        magsPlot = mags+[magPredict]
        angsPlot = angles+[angPredict]
        plt.scatter(stepsPlot, magsPlot)
        plt.plot(stepsPlot, magRegress.slope*np.array(stepsPlot)+magRegress.intercept)
        plt.figure("angles regression")
        plt.scatter(stepsPlot, angsPlot)
        plt.plot(stepsPlot, angRegress.slope*np.array(stepsPlot)+angRegress.intercept)

        print(magPredict, angPredict)
        SFPredict = np.array([magPredict*np.cos(angPredict), 
                            magPredict*np.sin(angPredict),
                            torqueTarg])
        print(SFPredict)

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

    def plot_spring_ODE(self, geometry, showBool=True, trans=1):

            self.A, self.B =    self.crsc.get_outer_geometry_ODE(self.resl, geometry)
            # print(self.A[-5:])
            # print("--")
            # print(self.B[-5:])
            self.Sn        =    self.path.get_neutralSurface(self.resl)
            self.Sc        = self.path.get_centroidalSurface(self.resl)

            # plot the principal leg
            plt.figure("Graphic Results")
            plt.plot(self.A[:,0],self.A[:,1],"--g", alpha=trans)
            plt.plot(self.B[:,0],self.B[:,1],"--g", alpha=trans)
            plt.plot(self.Sn[:,0],self.Sn[:,1],"--b",label="netural", alpha=trans)
            plt.plot(self.Sc[:,0],self.Sc[:,1],"--r",label="centroidal", alpha=trans)
            plt.plot(self.path.pts[:,0],self.path.pts[:,1])
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
                transformedSc = R.dot(self.Sc.T).T
                plt.plot(transformedA[:,0],transformedA[:,1],"k", alpha=trans)
                plt.plot(transformedB[:,0],transformedB[:,1],"k", alpha=trans)
                plt.plot(transformedSn[:,0],transformedSn[:,1],"--b", alpha=trans)
                plt.plot(transformedSc[:,0],transformedSc[:,1],"--r", alpha=trans)

            # plot geometry of inner and outer rotor of spring
            outerCircle = plt.Circle([0,0],self.outerRadius,color ="k",fill=False)
            innerCircle = plt.Circle([0,0],self.innerRadius,color ="k",fill=True)
            fig = plt.gcf()
            ax = fig.gca()
            ax.add_patch(outerCircle)
            ax.add_patch(innerCircle)

            if showBool:
                plt.show()

    def plot_spring(self, showBool=True, trans=1):

        self.A, self.B =    self.crsc.get_outer_geometry(self.resl)
        # print(self.A[-5:])
        # print("--")
        # print(self.B[-5:])
        self.Sn        =    self.path.get_neutralSurface(self.resl)
        self.Sc        = self.path.get_centroidalSurface(self.resl)

        # plot the principal leg
        plt.figure("Graphic Results")
        plt.plot(self.A[:,0],self.A[:,1],"k", alpha=trans)
        plt.plot(self.B[:,0],self.B[:,1],"k", alpha=trans)
        plt.plot(self.Sn[:,0],self.Sn[:,1],"--b",label="netural", alpha=trans)
        plt.plot(self.Sc[:,0],self.Sc[:,1],"--r",label="centroidal", alpha=trans)
        plt.plot(self.path.pts[:,0],self.path.pts[:,1])
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
            transformedSc = R.dot(self.Sc.T).T
            plt.plot(transformedA[:,0],transformedA[:,1],"k", alpha=trans)
            plt.plot(transformedB[:,0],transformedB[:,1],"k", alpha=trans)
            plt.plot(transformedSn[:,0],transformedSn[:,1],"--b", alpha=trans)
            plt.plot(transformedSc[:,0],transformedSc[:,1],"--r", alpha=trans)

        # plot geometry of inner and outer rotor of spring
        outerCircle = plt.Circle([0,0],self.outerRadius,color ="k",fill=False)
        innerCircle = plt.Circle([0,0],self.innerRadius,color ="k",fill=True)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_patch(outerCircle)
        ax.add_patch(innerCircle)

        if showBool:
            plt.show()

    def plot_deform(self, showBool=True):

        self.plot_spring(showBool=False, trans=0.25)

        alpha = self.path.get_alpha(self.ximesh)
        self.la, self.lb = self.crsc.get_neturalDistances(self.resl)
        self.h = self.crsc.h

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
                    i+=1
        else:
            print("Too early to call; please generate inner and outer surfaces before export")
