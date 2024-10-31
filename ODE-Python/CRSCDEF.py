import numpy as np
import numpy.linalg as lin

from utils import PPoly_Eval
from StatProfiler import SSProfile

deg2rad = np.pi/180

"""
##############################################
# STANDARD INTERFACE FOR ALL CRSCDEF OBJECTS #
##############################################
# Methods to include:
#
#   get_outer_geometry(resolution):
#           resolution as int
#
#   get_Ic(coord):
#           coord as scalar or vector float
#
#   get_Thk(coord, #hasPrev#):
#           hasPrev as boolean/float if prev exists
#
# Attributes to include:
#
#   self.t
#
#   self.fullParamLength (SHOULD BE "INHERITED" FROM PREVIOUSLY DEFINED PATHDEF)
#
#
"""

class Constant_Ic():
    def __init__(self, pathDef, t, h0):
        self.path = pathDef
        self.t    = t
        self.h0 = h0
        self.la = h0/2
        self.lb = self.la
        self.Ic0     = 1/12.0*self.t*self.h0**3

        self.fullParamLength = self.path.fullParamLength
        self.returnValue = 1
    
    def l_a_l_b_rootfinding(self, coord, lABPrev, printBool=0):

        """
        THIS FUNCTION:
            Uses newton's method to solve the nonlinear system of equations
            which defines the thickness of a bent beam with a given netural
            radius and second moment of area
        """

        # get the relevant values at a given point
        rn = self.path.get_rn(coord)
        Ic = self.get_Ic(coord)
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
        if not(np.isinf(rn) or abs(rn) > 10e10): # TODO: change this to some large finite threshold
            # check for if the beam was locally straight on the last iteration
            # print(lABPrev)
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
##### Standard Interface #####

    def get_neturalDistances(self, resolution):
        # shouldn't need this for now
        return self.returnValue

    def get_outer_geometry_ODE(self, resolution, lalb):
        # shouldn't need this for now
        return self.returnValue

    def get_outer_geometry(self, resolution):
        undeformedNeutralSurface = self.path.get_neutralSurface(resolution)
        ximesh = np.linspace(0,self.fullParamLength,resolution+1)

        Ic0 = self.get_Ic(0)

        hPrev = np.cbrt(12*Ic0/self.t)
        lABPrev = np.array([hPrev/2, hPrev/2])
        self.la = np.empty(len(ximesh))
        self.lb = np.empty(len(ximesh))
        self.h = np.empty(len(ximesh))

        # print(rn[-5:])
        # perform la/lb rootfinding for each step in the xi mesh
        for i in range(len(ximesh)):
            SSProfile("lAB rootfinding").tic()
            lAB = self.l_a_l_b_rootfinding(ximesh[i], lABPrev)
            SSProfile("lAB rootfinding").toc()
            self.la[i] = lAB[0]
            self.lb[i] = lAB[1]
            # calculate overall thickness
            self.h[i] = self.lb[i]+self.la[i]
            lABPrev = lAB
        # print(h[-5:])
        alpha = self.path.get_alpha(ximesh)
        # generate xy paths for surfaces
        self.undeformedBSurface = undeformedNeutralSurface + \
                                  np.hstack((np.atleast_2d(-self.lb*np.sin(alpha)).T,
                                             np.atleast_2d(self.lb*np.cos(alpha)).T))
        self.undeformedASurface = undeformedNeutralSurface - \
                                np.hstack((np.atleast_2d(-self.la*np.sin(alpha)).T,
                                           np.atleast_2d(self.la*np.cos(alpha)).T))

        return self.undeformedASurface, self.undeformedBSurface

    def get_Thk(self, coord, hasPrev=False):
        lalb = self.get_lalb
        h = lalb[0]+lalb[1]
        return h

    def get_lalb(self, coord, hasPrev=False):
        coord = np.atleast_1d(coord)
        if not hasPrev:
            hPrev = np.cbrt(12*self.get_Ic(coord[0])/self.t)
            lABPrev = np.array([hPrev/2, hPrev/2])
        else:
            lABPrev = hasPrev

        la = np.empty_like(coord)
        lb = np.empty_like(coord)
        h  = np.empty_like(coord)
        i = 0
        for value in coord:
            lAB = self.l_a_l_b_rootfinding(value, lABPrev)
            lABPrev = lAB
            la[i] = lAB[0]
            lb[i] = lAB[1]
            i+=1
        lalb = np.array([la, lb])
        return lalb
    
    def get_Ic(self, coord):
        if hasattr(coord, "__len__"):
            return np.ones_like(coord)*self.Ic0
        else:
            return self.Ic0
    
    def get_dIc(self, coord):
        return 0

class Empty():
    def __init__(self,
                 pathDef, t):

        self.path = pathDef
        self.fullParamLength = self.path.fullParamLength
        self.returnValue = 0
        self.t = t

        self.returnValue=0

#### STANDARD INTERFACE ####

    def get_neturalDistances(self, resolution):
        return self.returnValue

    def get_outer_geometry_ODE(self, resolution, lalb):
        return self.returnValue

    def get_outer_geometry(self, resolution):

        return self.returnValue

    def get_Thk(self, coord, hasPrev=False):
        return self.returnValue

    def get_lalb(self, coord, hasPrev=False):
        return self.returnValue

    def get_Ic(self, coord):
        return self.returnValue
    
    def get_dIc(self, coord):
        return self.returnValue

class Piecewise_Ic_Control():
    def __init__(self,
                 pathDef,
                 outPlaneThickness = 0.375,
                 IcPts             = np.array([0.008, 0.001, 0.008]),
                 IcParamLens       = np.array([0.5])):

        self.path = pathDef

        self.fullParamLength = self.path.fullParamLength

        self.t         = outPlaneThickness
        self.IcFactors = IcParamLens
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

        self.IcCoeffs, self.domains = self.Ic_multiPoly()

    def Ic_multiPoly(self):

        # Determine how many parabolas you need
        numPolys = (len(self.IcPts)-1)*2
        numSegments = int(numPolys/2)

        # Determine the x, y locations of your control points
        dys = np.empty(numPolys)
        dxs = np.empty(numPolys)

        for i in range((numSegments)):
            dys[i] = (self.IcPts[i+1]-self.IcPts[i])/2
            dxs[i] = (self.IcParamLens[i+1]-self.IcParamLens[i])/2

        ctrlX = np.empty(numPolys+1)
        ctrlY = np.empty(numPolys+1)

        for i in range(numSegments):
            ctrlX[2*i]   = self.IcParamLens[i]
            ctrlX[2*i+1] = self.IcParamLens[i]+dxs[i]
            ctrlY[2*i]   = self.IcPts[i]
            ctrlY[2*i+1] = self.IcPts[i]+dys[i]
        ctrlX[-1] = self.IcParamLens[-1]
        ctrlY[-1] = self.IcPts[-1]

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

    def l_a_l_b_rootfinding(self, coord, lABPrev, printBool=0):

        """
        THIS FUNCTION:
            Uses newton's method to solve the nonlinear system of equations
            which defines the thickness of a bent beam with a given netural
            radius and second moment of area
        """

        # get the relevant values at a given point
        rn = self.path.get_rn(coord)
        Ic = self.get_Ic(coord)
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
        if not(np.isinf(rn) or abs(rn) > 10e10): # TODO: change this to some large finite threshold
            # check for if the beam was locally straight on the last iteration
            # print(lABPrev)
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

#### STANDARD INTERFACE ####

    def get_neturalDistances(self, resolution):
        self.get_outer_geometry(resolution)
        return self.la, self.lb

    def get_outer_geometry_ODE(self, resolution, lalb):

        undeformedNeutralSurface = self.path.get_neutralSurface(resolution)
        ximesh = np.linspace(0,self.fullParamLength,resolution+1)

        y0 = self.l_a_l_b_rootfinding(0, np.array([0,0]))

        Ic0 = self.get_Ic(0)

        hPrev = np.cbrt(12*Ic0/self.t)
        lABPrev = np.array([hPrev/2, hPrev/2])
        self.la = np.empty(len(ximesh))
        self.lb = np.empty(len(ximesh))
        self.h = np.empty(len(ximesh))

        self.la = lalb[0,:]
        self.lb = lalb[1,:]
        self.h  = self.la+self.lb

        alpha = self.path.get_alpha(ximesh)
        # generate xy paths for surfaces
        self.undeformedBSurface = undeformedNeutralSurface + \
                                  np.hstack((np.atleast_2d(-self.lb*np.sin(alpha)).T,
                                             np.atleast_2d(self.lb*np.cos(alpha)).T))
        self.undeformedASurface = undeformedNeutralSurface - \
                                np.hstack((np.atleast_2d(-self.la*np.sin(alpha)).T,
                                           np.atleast_2d(self.la*np.cos(alpha)).T))

        return self.undeformedASurface, self.undeformedBSurface


    def get_outer_geometry(self, resolution):

        undeformedNeutralSurface = self.path.get_neutralSurface(resolution)
        ximesh = np.linspace(0,self.fullParamLength,resolution+1)

        Ic0 = self.get_Ic(0)

        hPrev = np.cbrt(12*Ic0/self.t)
        lABPrev = np.array([hPrev/2, hPrev/2])
        self.la = np.empty(len(ximesh))
        self.lb = np.empty(len(ximesh))
        self.h = np.empty(len(ximesh))

        # print(rn[-5:])
        # perform la/lb rootfinding for each step in the xi mesh
        for i in range(len(ximesh)):
            SSProfile("lAB rootfinding").tic()
            lAB = self.l_a_l_b_rootfinding(ximesh[i], lABPrev)
            SSProfile("lAB rootfinding").toc()
            self.la[i] = lAB[0]
            self.lb[i] = lAB[1]
            # calculate overall thickness
            self.h[i] = self.lb[i]+self.la[i]
            lABPrev = lAB
        # print(h[-5:])
        alpha = self.path.get_alpha(ximesh)
        # generate xy paths for surfaces
        self.undeformedBSurface = undeformedNeutralSurface + \
                                  np.hstack((np.atleast_2d(-self.lb*np.sin(alpha)).T,
                                             np.atleast_2d(self.lb*np.cos(alpha)).T))
        self.undeformedASurface = undeformedNeutralSurface - \
                                np.hstack((np.atleast_2d(-self.la*np.sin(alpha)).T,
                                           np.atleast_2d(self.la*np.cos(alpha)).T))

        return self.undeformedASurface, self.undeformedBSurface

    def get_Thk(self, coord, hasPrev=False):
        # I = th^3/12
        # h = cbrt(12I/t)
        coord = np.atleast_1d(coord)
        if not hasPrev:
            hPrev = np.cbrt(12*self.get_Ic(coord[0])/self.t)
            lABPrev = np.array([hPrev/2, hPrev/2])
        else:
            lABPrev = hasPrev

        la = np.empty_like(coord)
        lb = np.empty_like(coord)
        h  = np.empty_like(coord)
        i = 0
        for value in coord:
            lAB = self.l_a_l_b_rootfinding(value, lABPrev)
            lABPrev = lAB
            la[i] = lAB[0]
            lb[i] = lAB[1]
            i+=1
        h = la+lb
        if len(h)==1:
            return h[0]
        else:
            return h

    def get_lalb(self, coord, hasPrev=False):
        coord = np.atleast_1d(coord)
        if not hasPrev:
            hPrev = np.cbrt(12*self.get_Ic(coord[0])/self.t)
            lABPrev = np.array([hPrev/2, hPrev/2])
        else:
            lABPrev = hasPrev

        la = np.empty_like(coord)
        lb = np.empty_like(coord)
        h  = np.empty_like(coord)
        i = 0
        for value in coord:
            lAB = self.l_a_l_b_rootfinding(value, lABPrev)
            lABPrev = lAB
            la[i] = lAB[0]
            lb[i] = lAB[1]
            i+=1
        lalb = np.array([la, lb])
        return lalb

    def get_Ic(self, coord):
        out = PPoly_Eval(coord, self.IcCoeffs, ranges=self.domains)
        return out
    
    def get_dIc(self, coord):
        out = PPoly_Eval(coord, self.IcCoeffs, ranges=self.domains, deriv=1)
        return out