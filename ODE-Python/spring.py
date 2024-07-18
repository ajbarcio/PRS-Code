# TODO: Consistent formatting

import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from StatProfiler import SSProfile
import scipy

from copy import deepcopy as dc

from utils import numerical_fixed_mesh_diff, fixed_rk4, rk4_step, colorline
from polynomials import Ic_multiPoly, xy_poly, PPoly_Eval
from spring_utils import alpha_xy, r_n, d_rn, d_xi_d_s, d_alpha, d_2_alpha
from materials import * # This module ONLY includes materials objects so I think this is safe

deg2rad = np.pi/180

class Deform_Wrapper:
    def __init__(self, function):
        # Log all the outputs, input moments
        self.all_output = []
        self.all_moments = []
        # pass the function to be wrapped
        self.function = function

    def __call__(self, torqueTarg, ODE, torqueResolution=10, SF = np.array([0,0,0]), breakBool=False):
        # stick the torque targe in the moments list
        self.all_moments.append(torqueTarg)
        # time how long the deformation takes to solve
        SSProfile("BVP").tic()
        res, SF, divergeFlag, i  = self.function(torqueTarg, ODE, torqueResolution, SF, breakBool)
        SSProfile("BVP").toc()
        # log the solution
        self.all_output.append(SF)
        # return the quantities that the wrapped function would
        return res, SF, divergeFlag, i

class Spring:
    def __init__(self, material, name=None,
                 # whole bunch of default values:
                 n                   = 2,
                 fullParamLength     = 6,
                 outPlaneThickness   = 0.375,
                 radii               = np.array([1.1,2.025,2.025,2.95]),
                 betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
                 IcPts               = np.array([0.008, 0.001, 0.008]),
                 IcParamLens         = np.array([0.5]),                              # IcPts - 2
                 XYParamLens         = np.array([0.333,0.667]),             # radii - 2
                 resolution          = 200
                 ):
        # improperly initialize the spring (unable to start with good guess for arc length)
        E = material.E
        designStress = material.yieldStress*0.8
        self.init(E, designStress, n, fullParamLength, outPlaneThickness, radii,
                  betaAngles, IcPts, IcParamLens, XYParamLens, resolution)
        # iterate spring initialization until xi mesh approximates s mesh as well as possible
        conv=1
        err = 1
        i = 0
        # iterate until you converge on a value
        while conv>10e-6:
            SSProfile("reinit").tic()
            errPrev = err
            # measure real arc length each time
            correctedFullLength = self.measure_length()
            # initialize the spring with the remeasured legnth
            self.init(E, designStress, n, correctedFullLength, ### WE ONLY CHANGE THIS
                    outPlaneThickness, radii,
                    betaAngles, IcPts, IcParamLens, XYParamLens, resolution)
            err = lin.norm(self.dxids-np.ones(len(self.dxids)))
            # err = self.fullParamLength-self.measure_length() # These two error values seem equivalent
            conv = abs(err-errPrev)
            i+=1
            SSProfile("reinit").toc()
        # if the spring was given a name, save its input parameters in a filej
        self.name = name
        if not name==None:
            self.generate_param_file()
        # np.savetxt("springs/")

    def generate_param_file(self):

        labels = np.array(["n", "full_length", "Out_Plane_Thk"])
        for i in range(len(self.radii)):
            labels = np.append(labels, "radius_"+str(i))
        for i in range(len(self.betaAngles)):
            labels = np.append(labels, "beta_ang_"+str(i))
        for i in range(len(self.IcPts)):
            labels = np.append(labels, "ctrl_Ic_"+str(i))
        for i in range(len(self.IcFactors)):
            labels = np.append(labels, "ctrl_Ic_factor_"+str(i))
        for i in range(len(self.XYFactors)):
            labels = np.append(labels, "ctrl_XY_factor_"+str(i))

        saveVariables = np.hstack([self.n, self.fullParamLength, self.t,
                                   self.radii,self.betaAngles,self.IcPts,
                                   self.IcFactors,self.XYFactors])
        out = np.hstack([np.atleast_2d(labels).T, np.atleast_2d(saveVariables).T])
        filepath = "springs\\"+self.name
        print(out)
        np.savetxt(filepath, out, "%s", ",")

    def init(self,
             # whole bunch of default values:
             E                   = 27.5*10**6,
             designStress        = 270000,
             n                   = 2,
             fullParamLength     = 6,
             outPlaneThickness   = 0.375,
             radii               = np.array([1.1,2.025,2.025,2.95]),
             betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
             IcPts               = np.array([0.008, 0.001, 0.008]),
             IcParamLens           = np.array([0.5]),                              # IcPts - 2
             XYParamLens           = np.array([0.333,0.667]),             # radii - 2
             resolution          = 200):

        """
        THIS FUNCTION:
            Takes in default or custom parameters and
            defines the quantities and function needed for analysis
        """

        ############################ PARAMETER KEY ############################
        # E                 - youngs modulus                                  #
        # designStress      - max allowable stress, currently 80% yield       #
        # n                 - # of flexures                                   #
        # fullParamLength   - maximum value of internal parameter             #
        # outPlaneThickness - axial thickness of spring (in)                  #
        # radii             - radii of control points of netural surface (in) #
        # betaAngles        - angle of control points of neutral surface (rad)#
        # IcPts             - curved second moment of area of beam at points  #
        #                     along beam (in^4)                               #
        # IcParamLens       - internal parameter values at which Ic control   #
        #                     points apply                                    #
        # IcParamLens       - internal parameter values at which x-y          #
        #                     (radii-angle) control points apply              #
        # resolution        - number of points in mesh                        #
        ########################################################################

        # stick all the arguments in the object
        self.E = E
        self.designStress = designStress
        self.n = n
        self.t = outPlaneThickness
        self.fullParamLength = fullParamLength
        self.radii = radii
        self.betaAngles = betaAngles

        self.IcFactors = IcParamLens
        self.XYFactors = XYParamLens

        # create control points for polynomials

        # x-y control points
        self.pts = np.empty((len(radii),2))
        for i in range(len(radii)):
            self.pts[i,:] = [self.radii[i]*np.cos(self.betaAngles[i]),self.radii[i]*np.sin(self.betaAngles[i])]

        # Ic control points are input directly
        self.IcPts = IcPts

        # convert IcParamLens, input as proportions of the full parameter length
        # to parameter lengths
        self.IcParamLens = np.empty(len(IcPts))
        for i in range(len(IcPts)):
            if i==0:
                self.IcParamLens[i] = 0
            elif i==len(IcPts)-1:
                self.IcParamLens[i] = self.fullParamLength
            else:
                self.IcParamLens[i] = self.fullParamLength*IcParamLens[i-1]

        # Do the same with the parameter length controls on x-y points
        self.XYParamLens = np.empty(len(radii))
        for i in range(len(radii)):
            if i==0:
                self.XYParamLens[i] = 0
            elif i==len(radii)-1:
                self.XYParamLens[i] = self.fullParamLength
            else:
                self.XYParamLens[i] = self.fullParamLength*XYParamLens[i-1]

        # create a vector of all mutable parameters
        self.parameterVector = np.concatenate((self.radii, self.betaAngles,
                                              self.IcPts, self.IcFactors,
                                              self.XYFactors))
        # assign some constants for approximating derivatives
        # NOTE THAT NOT ALL OF THESE ARE USED YET
        self.finiteDifferenceLength = 0.001
        self.finiteDifferenceAngle  = .1*deg2rad
        self.finiteDifferenceForce  = 0.1
        self.finiteDifferenceTorque = 0.5
        self.finiteDifferenceIc     = 0.00001
        self.finiteDifferenceFactor = 0.001
        # create a vector of finite differences for mutable parameters
        self.finiteDifferenceVector = np.concatenate((np.ones(len(self.radii))*self.finiteDifferenceLength,
                                                      np.ones(len(self.betaAngles))*self.finiteDifferenceAngle,
                                                      np.ones(len(self.IcPts))*self.finiteDifferenceIc,
                                                      np.ones(len(self.IcFactors))*self.finiteDifferenceFactor,
                                                      np.ones(len(self.XYFactors))*self.finiteDifferenceFactor,
                                                      ))

        # assign some constants for meshing the spring
        self.resl = resolution
        self.len = self.resl+1
        self.step = fullParamLength/self.resl
        self.endIndex = self.len-1

        # create uniform mesh for spring
        # 'xi' is the internal parameter
        self.ximesh = np.linspace(0,self.fullParamLength,self.len)

        # generate the coefficients for x-y and Ic polynomials
        self.geometry_coeffs()
        # find the derivative to convert between the internal parameter mesh
        # and an arc length based mesh
        self.dxids = d_xi_d_s(self.ximesh, self.XCoeffs, self.YCoeffs)

        # this takes way too many lines but its better to see the math written
        # out:

        # Find coordinates of frames:
        # at base of spring
        self.x0 = PPoly_Eval(0,self.XCoeffs)
        self.y0 = PPoly_Eval(0,self.YCoeffs)
        # at end of spring
        self.xL = PPoly_Eval(self.fullParamLength,self.XCoeffs)
        self.yL = PPoly_Eval(self.fullParamLength,self.YCoeffs)
        # at end of spring, interms of frame at base of spring
        self.momentArmX = self.xL - self.x0
        self.momentArmY = self.yL - self.y0

        # inherit deform wrapper
        self.wrapped_torque_deform = Deform_Wrapper(self.deform_by_torque)

    def tune_stiffness(self, targetStiffness, targetTorque, discludeVector):
        err = 1
        errPrev = 1
        convergenceLimit = 0.01
        successFlag = 1
        try:
            while err > convergenceLimit:
                errPrev = err
                # make a new spring from current parameter vector
                self.redefine_spring_from_parameterVector()
                # evaluate the stiffness of the spring
                res, SF, divergeFlag, i = self.deform_by_torque_smartGuess(targetTorque, self.deform_ODE)
                if divergeFlag:
                    successFlag = 0
                    break
                stiffness = targetTorque/(self.dBeta/deg2rad)
                print("guess stiffness:", stiffness)
                # calcualte the distance from the target stiffness
                err = (stiffness-targetStiffness)/targetStiffness
                print("err:", err)
                if err>errPrev:
                    successFlag = 0
                    break
                # calcualte the jacobian
                J = self.stiffness_jacobian(self.deform_ODE, targetStiffness, targetTorque, err)
                for i in range(len(J)):
                    J[i]*=discludeVector[i]
                Jinv = lin.pinv(np.atleast_2d(J).T)
                # print(Jinv)
                # update parameter vector guess
                # print("before:", self.parameterVector)
                self.parameterVector = (self.parameterVector - Jinv*err)[0]
                # print("after:", self.parameterVector)
            if successFlag:
                print("successful convergence to ",convergenceLimit*100,"%")
            print("stiffness (lbf/deg)", stiffness)
            print("max stress:", self.maxStress)
            print("des stress:", self.designStress)
        except KeyboardInterrupt:
                # If you get bored and kill it early:
                # calculate and prepare for plotting the full results
                print(self.n)
                self.full_results()
                # print out some results
                print("torque:", targetTorque)
                print("angular deformation:", self.dBeta/deg2rad)
                print("stiffness (lbf/deg)", targetTorque/(self.dBeta/deg2rad))
                print("ideal stiffness:   ", targetTorque/(5))
                print("max stress:", self.maxStress)
                print("des stress:", self.designStress)
                # save parameter information to file
                filepath = "springs\\"+self.name
                np.savetxt(filepath, self.parameterVector, "%f", ",")
                # save surface geometry to file
                A = np.hstack((self.undeformedASurface,np.atleast_2d(np.zeros(len(self.undeformedASurface))).T))
                B = np.hstack((self.undeformedBSurface,np.atleast_2d(np.zeros(len(self.undeformedBSurface))).T))
                np.savetxt("surfaces\\premature_A_surface.txt", A, delimiter=",", fmt='%f')
                np.savetxt("surfaces\\premature_B_surface.txt", B, delimiter=",", fmt='%f')

    def stiffness_jacobian(self, ODE, targetStiffness, targetTorque, currErr):
        print("calculating J")
        initialParameterVector = dc(self.parameterVector)
        J = np.empty(len(self.parameterVector))
        for i in range(len(self.parameterVector)):
            # print("parameter derivative", i, "/15", end="\r")
            # Forward difference
            self.parameterVector[i] += self.finiteDifferenceVector[i]
            self.redefine_spring_from_parameterVector()
            self.deform_by_torque_smartGuess(targetTorque, self.deform_ODE)
            stiffness = targetTorque/self.dBeta
            # calcualte the distance from the target stiffness
            errf = (stiffness-targetStiffness)/targetStiffness
            # reset parameter vector
            self.parameterVector = dc(initialParameterVector)
            # Backwards difference
            self.parameterVector[i] -= self.finiteDifferenceVector[i]
            self.redefine_spring_from_parameterVector()
            self.deform_by_torque_smartGuess(targetTorque, self.deform_ODE)
            stiffness = targetTorque/self.dBeta
            # calcualte the distance from the target stiffness
            errb = (stiffness-targetStiffness)/targetStiffness
            self.parameterVector = dc(initialParameterVector)
            J[i] = (errf-errb)/(2*self.finiteDifferenceVector[i])
        return J

    def sensitivity_study(self, index, factor, resolution, testTorque):

        """
        THIS FUNCTION:
            Runs a sensitivity study across a single design variable, to expose
            how a single variable locally affects max stress and deflection
        """

        # initialize a wrapper to log results
        sensitivityTorqueWrapper = Deform_Wrapper(self.deform_by_torque)

        # store a copy of the original "parameters"
        oldParameterVector = dc(self.parameterVector)
        # print(oldParameterVector)

        # establish the design variable range over which to evaluate
        # deformations
        start = (1-factor)*self.parameterVector[index]
        end   = (1+factor)*self.parameterVector[index]
        # uniform mesh over the range
        values = np.linspace(start, end, resolution+1)
        # prepare arrays to store performance characteristics
        stresses    = np.empty(len(values))
        deflections = np.empty(len(values))
        # over the range
        for i in range(len(values)):\
            # time the study
            SSProfile("sensitivity study").tic()
            # set the parameter in question to the changed value
            self.parameterVector[index] = values[i]
            print("trying:",values[i])
            # re-initialize spring with appropriate design variables
            self.redefine_spring_from_parameterVector()
            # USE SMART GUESS METHOD (see deform_unit_test for documentation)
            SFGuess = self.smart_initial_load_guess(testTorque,self.deform_ODE)
            sensitivityTorqueWrapper(testTorque,self.deform_ODE,SF=SFGuess)
            self.generate_undeformed_surfaces()
            self.calculate_stresses()
            # record stress and deflection with these design variables
            stresses[i] = self.maxStress
            deflections[i] = self.dBeta/deg2rad
            SSProfile("sensitivity study").toc()
        # combine all the results into one array structure (for output to csv)
        output = np.array(sensitivityTorqueWrapper.all_output)
        sensitivityResults = np.hstack((output, np.atleast_2d(deflections).T,
                                        np.atleast_2d(stresses).T,
                                        np.atleast_2d(values).T))
        # reset the spring to its original state
        self.parameterVector = dc(oldParameterVector)
        self.redefine_spring_from_parameterVector()

        return sensitivityResults

    def redefine_spring_from_parameterVector(self):

        """
        THIS FUNCTION:
            regenerates a spring when the parameter vector is changed
        """

        ### This is done in a very dumb way, TODO to make this into a dict and
        #   use for loops

        # record all the new "parameters"
        newRadii      = self.parameterVector[0:len(self.radii)]
        startIndex = len(newRadii)
        newBetaAngles = self.parameterVector[startIndex:
                                             startIndex+len(self.betaAngles)]
        startIndex = len(np.concatenate([newRadii,newBetaAngles]))
        newIcPts = self.parameterVector[startIndex:
                                        startIndex+len(self.IcPts)]
        startIndex = len(np.concatenate([newRadii,newBetaAngles,newIcPts]))
        newIcFactors = self.parameterVector[startIndex:
                                              startIndex+len(self.IcFactors)]
        startIndex = len(np.concatenate([newRadii,newBetaAngles,
                                        newIcPts,newIcFactors]))
        newXYFactors = self.parameterVector[startIndex:
                                              startIndex+len(self.XYFactors)]
        startIndex = len(np.concatenate([newRadii,newBetaAngles,
                                        newIcPts,newIcFactors,newXYFactors]))
        # make sure you created the right length of vector
        assert(startIndex==len(self.parameterVector))


        # re-initialize spring like in __init__
        self.init(E                 = self.E,
                  designStress      = self.designStress,
                  n                 = self.n,
                  fullParamLength   = self.fullParamLength,
                  outPlaneThickness = self.t,
                  radii             = newRadii,
                  betaAngles        = newBetaAngles,
                  IcPts             = newIcPts,
                  IcParamLens       = newIcFactors,
                  XYParamLens       = newXYFactors,
                  resolution        = self.resl)
        correctedFullLength = self.measure_length()
        self.init(E                 = self.E,
                  designStress      = self.designStress,
                  n                 = self.n,
                  fullParamLength   = correctedFullLength,
                  outPlaneThickness = self.t,
                  radii             = self.radii,
                  betaAngles        = self.betaAngles,
                  IcPts             = self.IcPts,
                  IcParamLens       = self.IcFactors,
                  XYParamLens       = self.XYFactors,
                  resolution        = self.resl)

    def geometry_coeffs(self):

        """
        THIS FUNCTION:
            Produces coefficients for the polynomials that define the
            geometry of the spring. Currently this is the x-y coordinates of the
            neutral surface, as well as the second moment of area
        """
        # sticks the coefficients in the object
        self.IcCoeffs, self.domains = Ic_multiPoly(self.IcPts, self.IcParamLens)
        self.XCoeffs, self.YCoeffs  = xy_poly(self.pts, self.XYParamLens)
        # compiles all the coeffs into one variable
        self.geometryDef = [self.XCoeffs, self.YCoeffs, self.IcCoeffs]
        return self.geometryDef

    def smart_initial_load_guess(self, torqueTarg, ODE):

        """
        THIS FUNCTION:
            Takes in a target torque and uses linear extrapolation to find a
            good guess for the spatial force vector solution at full loading
            torque by evaluating the deflection of the spring at two smaller
            torque levels, which should be safe to solve for without fear of
            diverging for most springs
        """

        # initialize guess wrapper for use in linear extrapolation
        guessWrapper = Deform_Wrapper(self.deform_by_torque)


        # evaluate spring at two ("arbitrarily") low-torque points that should
        # converge nicely
        SF = np.array([0,0,torqueTarg*0.05])
        guessWrapper(SF[2],ODE,SF=SF)
        SF = np.array([0,0,torqueTarg*0.15])
        guessWrapper(SF[2],ODE,SF=SF)

        # record angles and magnitudes of [Fx,Fy] vector at each torque level
        torques = np.array([torqueTarg*0.1,torqueTarg*0.2])
        angles = np.empty(len(guessWrapper.all_output))
        magnitudes = np.empty(len(guessWrapper.all_output))
        for i in range(len(guessWrapper.all_output)):
            angles[i] = np.arctan2(guessWrapper.all_output[i][1],guessWrapper.all_output[i][0])
            magnitudes[i] = lin.norm(guessWrapper.all_output[i][0:2])
        # fit a first-order regression to angles and magnitudes
        anglePoly = np.polyfit(torques,angles,1)
        magPoly = np.polyfit(torques,magnitudes,1)
        # evaluate that regression at the true torque target
        angleGuess  = np.polyval(anglePoly, torqueTarg)
        magGuess  = np.polyval(magPoly, torqueTarg)
        # create a spacial force vector guess from that magnitude and angle
        SFGuess = np.array([magGuess*np.cos(angleGuess), magGuess*np.sin(angleGuess), torqueTarg])

        return SFGuess

    def deform_by_torque_smartGuess(self, torqueTarg, ODE):
        SFGuess = self.smart_initial_load_guess(torqueTarg,self.deform_ODE)
        res, SF, divergeFlag, i = self.deform_by_torque(torqueTarg,self.deform_ODE,SF=SFGuess)
        return res, SF, divergeFlag, i

    def deform_by_torque_slowRamp(self, torqueTarg, ODE, torqueResolution=15):

        """
        THIS FUNCTION:
            Slowly increases the torque loading at a fixed rate, using the
            spatial force vector solution at each torque level as the initial
            guess vector for the next torque level
        """
        # initialize wrapper to record results
        self.rampWrapper = Deform_Wrapper(self.deform_by_torque)
        # start at 0
        SFStart = np.array([0,0,0])
        # the err is a two vector, so make it arbitrarily high to enter loop
        err = np.ones(2)*float('inf')
        # 0th iteration
        i = 0
        j = 0
        k = 0
        divergeFlag = 0
        stepTorque = 0
        # until you reach full torque
        while stepTorque <= torqueTarg:
            SFStart[2] = stepTorque
            print("torque level:", stepTorque)
            print(SFStart)
            # defrom using previous SF Force guess at curent torque step
            res, SF, divergeFlag, i = self.rampWrapper(stepTorque,ODE,SF=SFStart,breakBool=True)
            if divergeFlag:
                # if it diverges, and it isn't the first step:
                if stepTorque:
                    # decrease the torque to the previous level
                    stepTorque -= torqueTarg/torqueResolution
                # double the torque resolution (halve the step size)
                torqueResolution *= 2
            # if it's successful:
            else:
                # increase the torque by one step
                stepTorque+=torqueTarg/torqueResolution
                # use the previous guess
                SFStart = SF
            # used for debugging/evaluation purposes: j tracks number of torque
            # steps, k tracks number of forward integrations
            j += 1
            k += i
        return res, SF, divergeFlag, k

    def deform_by_torque(self,torqueTarg,ODE,torqueResolution=10,SF=np.array([0,0,0]),breakBool=False):

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
            J = self.boundary_condition_jacobian(ODE,SF,torqueTarg)
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
        return self.res, self.solnSF, divergeFlag, i

    def forward_integration(self, ODE, SF, torqueTarg):

        """
        THIS FUNCTION:
            Evaluates a spring's deflection under a certain loading condition
            by integrating forward across a specified ODE (there is currently
            only one supported)
        """

        # set up initial condition
        dgds0 = (SF[2]/self.n + self.momentArmY*SF[0] - self.momentArmX*SF[1])/(self.E*PPoly_Eval(0, self.IcCoeffs, ranges=self.domains))

        # perform forward integration
        self.res  = fixed_rk4(ODE, np.array([0,self.x0,self.y0]), self.ximesh, (SF[0], SF[1], dgds0))
        # calcualte error values (difference in pre-and post-iteration radii)
        #                        (difference in final gamma and final beta)
        Rinitial = lin.norm([self.xL,self.yL])
        Rfinal   = lin.norm([self.res[1,-1],self.res[2,-1]])
        # Calculate dBeta using arccos
        ptInitial = self.pts[-1]/lin.norm(self.pts[-1])
        ptFinal   = np.array([self.res[1,-1],self.res[2,-1]])/ \
                             lin.norm(np.array([self.res[1,-1],self.res[2,-1]]))
        # print(ptFinal.dot(ptInitial))
        cosAngle = min(ptFinal.dot(ptInitial),1.0)
        self.dBeta    = np.arccos(cosAngle)*np.sign(torqueTarg)
        # Err = diff. in radius, diff between gamma(L) and beta(L), distance
        #       from target torque (just to make matrix square)
        err = np.array([Rinitial-Rfinal, self.res[0,-1]-self.dBeta, SF[2]-torqueTarg])
        return err, self.res

    def boundary_condition_jacobian(self, ODE, SF, torqueTarg, n=3):

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

    def deform_ODE(self, xi, p, *args):

        """
        THIS FUNCTION:
            Evaluates a single step of the ODE which governs the deformation of
            a curved beam
        """

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
        # dxds and dyds are used in diffeq
        dxds = dxids*dxdxi
        dyds = dxids*dydxi

        ### DO THE DIFFEQ ###

        LHS = np.empty(3)

        LHS[0] = (dgds0 + Fy/Mdim*(p[1]-self.x0) - Fx/Mdim*(p[2]-self.y0))
        LHS[1] = np.cos(np.arctan2(dyds,dxds)+p[0])
        LHS[2] = np.sin(np.arctan2(dyds,dxds)+p[0])
        # check to make sure the math makes sense (I.C. satisfied)
        if xi==0:
            assert(np.isclose(LHS[0],dgds0,rtol=1e-5))

        # transfer back from s space to xi space
        LHS = LHS/dxids
        return LHS

    def deform_ODE_tension(self, xi, p, *args):

        """
        THIS FUNCTION:
            Evaluates a single step of the ODE which governs the deformation of
            a curved beam
        """

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
        # dxds and dyds are used in diffeq
        dxds = dxids*dxdxi
        dyds = dxids*dydxi

        ### DO THE DIFFEQ ###

        LHS = np.empty(3)

        LHS[0] = (dgds0 + Fy/Mdim*(p[1]-self.x0) - Fx/Mdim*(p[2]-self.y0))
        LHS[1] = np.cos(np.arctan2(dyds,dxds)+p[0])
        LHS[2] = np.sin(np.arctan2(dyds,dxds)+p[0])
        # check to make sure the math makes sense (I.C. satisfied)
        if xi==0:
            assert(np.isclose(LHS[0],dgds0,rtol=1e-5))

        # transfer back from s space to xi space
        LHS = LHS/dxids
        return LHS

    def geo_ODE(self, xi, p, *args):
        rn     = r_n(xi, self.XCoeffs, self.YCoeffs)
        drnds  = d_rn(xi, self.XCoeffs, self.YCoeffs)
        # drnds  = d_rn_d_s_numerical(s, geometryDef[0], geometryDef[1])
        cI     = PPoly_Eval(xi, self.IcCoeffs, deriv=0)
        dcIds  = PPoly_Eval(xi, self.IcCoeffs, deriv=1)
        ## OOF
        dads   = d_alpha(xi, self.XCoeffs, self.YCoeffs)
        d2ads2 = d_2_alpha(xi, self.XCoeffs, self.YCoeffs)
        geoFunc1 = (dcIds*dads+cI*d2ads2)/self.t
        # geoFunc2 = rn*(-drnds*((p[0]+p[1])/rn**2+(-p[0]-p[1])/((rn+p[1])*(rn-p[0]))))
        geoFunc2 = d2ads2*(p[0]+p[1])*(p[1]-p[0]-p[0]*p[1]*dads)
        print(p[0], p[1])
        print(dads, d2ads2)
        print(geoFunc2)
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
        ########### I CANNOR REMEMBNER WHAT ANY OF THIS DOES
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

    def l_a_l_b_rootfinding(self, s, lABPrev, printBool=0):

        """
        THIS FUNCTION:
            Uses newton's method to solve the nonlinear system of equations
            which defines the thickness of a bent beam with a given netural
            radius and second moment of area
        """

        # get the relevant values at a given point
        rn = r_n(s, self.XCoeffs, self.YCoeffs)
        Ic = PPoly_Eval(s, self.IcCoeffs, ranges = self.domains)
        # define the system to be solved
        def func(x, rn, Ic):
            f1 = (x[0]+x[1])/(np.log((rn+x[1])/(rn-x[0])))-rn
            f2 = self.t*rn*(x[0]+x[1])*(x[1]/2-x[0]/2)-Ic
            return np.array([f1, f2])
        def jac(x, rn, Ic): # IC doesn't happen to occur in the derivatives
            return np.array([[1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn-x[0])), \
                            1/(np.log((rn+x[1])/(rn-x[0])))-(x[0]+x[1])/(((np.log((rn+x[1])/(rn-x[0])))**2)*(rn+x[1]))], \
                            [-rn*self.t*x[0], rn*self.t*x[1]]])
        # some error checking and escapes for possible non-convergent cases
        if not(np.isinf(rn) or rn > 10e10): # TODO: change this to some large finite threshold
            # check for if the beam was locally straight on the last iteration
            if lABPrev[0]==lABPrev[1]:
                # perturb the initial guess a bit to avoid convergence issues
                x0 = [lABPrev[0], lABPrev[1]+0.001]
            else:
                # THIS IS THE "NORMAL" CASE (in which the previous step was curved)
                x0 = lABPrev
            err = 1
        # escape for locally straight beam
        else:
            # set the error to 0 to not enter the root finder
            err = 0
            # use the straight beam definition of I to find the thickness
            l = np.cbrt(12*Ic/self.t)/2
            x0 = [l, l]
        x = x0

        self.iii = 0
        # solve the problem (do newtons method to rootfind)
        while err > 10**-6 and self.iii <500:
            # newtons method
            xprev = x
            x = x - np.transpose(lin.inv(jac(x, rn, Ic)).dot(func(x, rn, Ic)))
            # error to track convergence
            err = lin.norm(x-xprev)
            self.iii+=1
        # escape if convergence goes towards beam so uncurved that
        # its _basically_ straight (this results in very high thicknesses):
        # this threshold is arbitrary

        # TODO: Delete this
        # if(lin.norm(x)>lin.norm(lABPrev)*100):
        #     # use straight beam definition
        #     l = np.cbrt(12*Ic/self.t)/2
        #     x = [l,l]
        # boolean value used to print out debug info
        if(printBool):
            print(x0)
            print(x)
            print(self.iii)
            print(rn)
            print(err)
        return x    # here x is [la, lb]

    def generate_undeformed_surfaces(self):

        """
        THIS FUNCTION:
            generates the inner and outer surfaces of the spring (in x-y
            coordinates) as well as the centroidal surface
            in the undeformed state
        """

        # Generate neutral radius path and give it nicely formatted class variable
        self.undeformedNeutralSurface = np.hstack((np.atleast_2d(PPoly_Eval(self.ximesh, self.XCoeffs)).T,
                                                   np.atleast_2d(PPoly_Eval(self.ximesh, self.YCoeffs)).T))

        # prepare a whole bunch of output arrays
        lABPrev = [0,0]
        self.la = np.empty(len(self.ximesh))
        self.lb = np.empty(len(self.ximesh))
        self.h = np.empty(len(self.ximesh))
        # create meshed Ic and rn arrays for use in calculating inner/outer
        # surface
        self.Ic = PPoly_Eval(self.ximesh, self.IcCoeffs, ranges=self.domains)
        self.rn = r_n(self.ximesh, self.XCoeffs, self.YCoeffs)
        self.rootfindingSteps = 0
        # perform la/lb rootfinding for each step in the xi mesh
        for i in range(len(self.ximesh)):
            SSProfile("lAB rootfinding").tic()
            lAB = self.l_a_l_b_rootfinding(self.ximesh[i], lABPrev)
            self.rootfindingSteps += self.iii
            SSProfile("lAB rootfinding").toc()
            self.la[i] = lAB[0]
            self.lb[i] = lAB[1]
            # calculate overall thickness
            self.h[i] = self.lb[i]+self.la[i]
            lABPrev = lAB
        # calculate eccentricity as well as radius of curvature for each surface
        self.ecc = self.Ic/(self.t*self.h*self.rn)
        self.a = self.rn-self.la
        self.b = self.rn+self.la
        # Area  = out of plane thickness * in plane thickness
        self.A = self.h*self.t
        # meshed tangent angle
        self.alpha = alpha_xy(self.ximesh, self.XCoeffs, self.YCoeffs)
        # generate xy paths for surfaces
        self.undeformedBSurface = self.undeformedNeutralSurface + \
                                  np.hstack((np.atleast_2d(-self.lb*np.sin(self.alpha)).T,
                                             np.atleast_2d(self.lb*np.cos(self.alpha)).T))
        self.undeformedASurface = self.undeformedNeutralSurface - \
                                np.hstack((np.atleast_2d(-self.la*np.sin(self.alpha)).T,
                                           np.atleast_2d(self.la*np.cos(self.alpha)).T))

        # generate centroidal surface
        self.undeformedCentroidalSurface = self.undeformedNeutralSurface+np.hstack((np.atleast_2d(self.ecc*np.sin(self.alpha)).T, np.atleast_2d(self.ecc*np.cos(self.alpha)).T))

        return self.undeformedASurface, self.undeformedBSurface

    def calculate_stresses(self):

        """
        THIS FUNCTION:
            Calcualtes the maximum stress at any point along the beam
        """

        # calculate gamma derivative for stress calcs
        dgdxi = numerical_fixed_mesh_diff(self.res[0,:], self.ximesh)
        self.dgds = dgdxi*self.dxids

        # prepare stress arrays
        self.innerSurfaceStress = np.empty(len(self.ximesh))
        self.outerSurfaceStress = np.empty(len(self.ximesh))
        # calculate stress differently depending on whether or not beam is
        # locally straight
        for i in range(len(self.innerSurfaceStress)):
            if not np.isinf(self.rn[i]):
                self.innerSurfaceStress[i] =  \
                    abs(self.E*(1-self.rn[i]/self.a[i])*self.rn[i]*self.dgds[i])
                self.outerSurfaceStress[i] =  \
                    abs(self.E*(1-self.rn[i]/self.b[i])*self.rn[i]*self.dgds[i])
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

        return self.maxStress, self.maxStresses

    def full_results(self, plotBool=1, deformBool=1):

        """
        THIS FUNCTION:
            Ensures all final results (e.g. stresses) have been calculated and
            either plots or does not plot a beam in an undeformed or deformed
            state
        """

        # Generate outer and inner surfaces
        self.generate_undeformed_surfaces() # A and B surface come from here

        plottedLines = []

        if plotBool:
            plt.figure("Graphic Results")
            if not deformBool:
                # only plot neutral and centroidal surface if undeformed
                plt.plot(self.undeformedNeutralSurface[:,0],self.undeformedNeutralSurface[:,1])
                plottedLines.append(self.undeformedNeutralSurface)
                plt.plot(self.undeformedCentroidalSurface[:,0],self.undeformedCentroidalSurface[:,1])
                plottedLines.append(self.undeformedCentroidalSurface)
            # in any case plot the inner and outer surfaces
            plt.plot(self.undeformedASurface[:,0],self.undeformedASurface[:,1],"--b")
            plt.plot(self.undeformedBSurface[:,0],self.undeformedBSurface[:,1],"--b")
        if deformBool:
            # generate neutral surface after deformation (and create nicely formatted class variables)
            # (some very silly python array handling happens here)
            self.deformedNeutralSurface = np.hstack((np.atleast_2d(self.res[1,:]).T, np.atleast_2d(self.res[2,:]).T))
            # generate deformed outer surfaces
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
                plottedLines.append(self.deformedNeutralSurface)
                # plot inner and outer surfaces
                plt.plot(self.deformedASurface[:,0],self.deformedASurface[:,1],"-k")
                plottedLines.append(self.deformedASurface)
                plt.plot(self.deformedBSurface[:,0],self.deformedBSurface[:,1],"-k")
                plottedLines.append(self.deformedBSurface)
                # plot scale
                colorline(np.ones(101)*(np.nanmax(self.deformedASurface[:,0])+.5),np.linspace(0,2,101),np.linspace(0,1,101),cmap=plt.get_cmap('rainbow'),linewidth=10)

        # extra "visual sugar" included here, at end of function, for organizational purposes
        if plotBool:
            # plot geometry of inner and outer rotor of spring
            outerCircle = plt.Circle([0,0],self.radii[-1],color ="k",fill=False)
            innerCircle = plt.Circle([0,0],self.radii[0],color ="k",fill=False)
            fig = plt.gcf()
            ax = fig.gca()
            ax.add_patch(outerCircle)
            ax.add_patch(innerCircle)
            # plot the rest of the arms of the spring (rotate calculated geometry)
            ang = 2*np.pi/self.n
            for i in range(len(plottedLines)):
                for j in np.linspace(1,self.n-1,self.n-1):
                    th = ang*(j)
                    R = np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
                    transformedLine = R.dot(plottedLines[i].T)
                    transformedLine = transformedLine.T
                    plt.plot(transformedLine[:,0],transformedLine[:,1],"-k")
            plt.plot(self.pts[:,0],self.pts[:,1])
            # show it
            plt.axis("equal")

    def measure_length(self):

        """
        THIS FUNCTION:
            measures the length of an initialized spring
        """

        # calcualte derivative of x and y with respect to xi (internal parameter)
        self.dxdxi = PPoly_Eval(self.ximesh, self.XCoeffs, deriv = 1)
        self.dydxi = PPoly_Eval(self.ximesh, self.YCoeffs, deriv = 1)
        integrand = np.sqrt(self.dxdxi**2+self.dydxi**2)
        # create a mesh in s according by integrating along the parameter to find
        # arc length
        self.smesh = np.zeros(len(self.ximesh))
        self.smesh[1:len(self.smesh)] = scipy.integrate.cumulative_trapezoid(integrand, self.ximesh, self.ximesh[1]-self.ximesh[0])
        # print(self.smesh)
        assert(self.smesh[0]==0)
        # return the overall length of the spring
        return self.smesh[-1]

def generate_default_spring(name="default_spring"):
    # Fit the Gen 2 actuator form factor
    R0 = 2.0/2
    R3 = 5.0/2

    # All these radii, angles, and Ic come from a previous version of this code and
    # can be considered as arbitrary starting points

    R1 = (R0+R3)/2+.26+.05
    R2 = (R0+R3)/2-.25

    fullAngle = 165

    beta1 = fullAngle/3*deg2rad*.5+.02

    beta2 = 2*fullAngle/3*deg2rad*0.9
    beta0 = fullAngle*deg2rad

    # generate all the input arrays for adjustible parameters
    # (treat material properties, fullParamLength guess, out of plane thickness,
    #  and resolution as fixed)

    inputRadii      = np.array([R0,R1,R2,R3])
    inputBetaAngles = np.array([0,beta1,beta2,beta0])

    Ics = np.array([0.008*.75, 0.00025*.85, 0.00025*1.15, 0.008*1.15])
    IcLens=np.array([0.4*1.2, 0.667*0.8])

    XYParamLens = np.array([0.333,0.667])

    # Generate the spring:

    defaultSpring = Spring(Maraging300Steel(), n = 2, radii=inputRadii,
                                betaAngles=inputBetaAngles,
                                IcPts=Ics,
                                IcParamLens=IcLens, XYParamLens=XYParamLens,
                                name=name)

    A = np.hstack((defaultSpring.undeformedASurface,np.atleast_2d(np.zeros(len(defaultSpring.undeformedASurface))).T))
    B = np.hstack((defaultSpring.undeformedBSurface,np.atleast_2d(np.zeros(len(defaultSpring.undeformedBSurface))).T))
    np.savetxt("surfaces\\"+defaultSpring.name+"_A_surface.txt", A, delimiter=",", fmt='%f')
    np.savetxt("surfaces\\"+defaultSpring.name+"_B_surface.txt", B, delimiter=",", fmt='%f')

    return defaultSpring

def generate_simple_spring():
    # Fit the Gen 2 actuator form factor
    R0 = 2.0/2
    R3 = 5.0/2

    # All these radii, angles, and Ic come from a previous version of this code and
    # can be considered as arbitrary starting points

    R1 =   (R0+R3)/2

    # generate all the input arrays for adjustible parameters
    # (treat material properties, fullParamLength guess, out of plane thickness,
    #  and resolution as fixed)

    inputRadii      = np.array([R0,R1,R3])
    inputBetaAngles = np.array([0,10,15])*deg2rad

    Ics = np.array([0.008*.85, 0.00025*1.15, 0.008*1.15])
    IcLens=np.array([0.5])

    XYParamLens = np.array([0.5])

    # Generate the spring:

    curvedSpring = Spring(Maraging300Steel(), n = 2, radii=inputRadii,
                                betaAngles=inputBetaAngles,
                                IcPts=Ics,
                                IcParamLens=IcLens, XYParamLens=XYParamLens,
                                name="simple_spring")

    return curvedSpring