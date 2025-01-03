import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
import scipy
import sympy as sp
import warnings
warnings.filterwarnings("ignore")

from modules.utils import PPoly_Eval, numerical_fixed_mesh_diff, circle_intersection, deg2rad
from modules.StatProfiler import SSProfile

from abc import ABC, abstractmethod
from typing import Optional

class Path(ABC):
    def __init__(self, n: int, arcLen: float,
                 startPoint: tuple = None, endPoint: tuple = None):
        # Absorb universal parameters
        self.n = n
        self.arcLen=arcLen
        # Handle whether start/end points defined for this path type
        optionalPoints = ['startPoint', 'endPoint']
        Mlocals = locals()
        points   = {name: Mlocals[name] for name in optionalPoints}
        handlers = {name: getattr(self, f'calculate_{name}') for name in optionalPoints}
        for name, value in points.items():
            if value is not None:
                setattr(self, name, value)
            else:
                setattr(self, name, handlers[name]())
        self.momentArm = tuple(x - y for x, y in zip(self.endPoint, self.startPoint))

        # Check for existence of either CENTROIDAL or NEUTRAL radius definitions
        if not (self.get_xy_n or self.get_xy_c):
            raise NotImplementedError("Subclasses must implement either method_a or method_b.")
        if not (self.get_dxdy_n or self.get_dxdy_c):
            raise NotImplementedError("Subclasses must implement either method_a or method_b.")
        
    def calculate_startPoint(self):
        """Optional method, can be overridden by subclasses if applicable."""
        raise NotImplementedError("This method is not implemented for this path type.")

    def calculate_endPoint(self):
        """Optional method, can be overridden by subclasses if applicable."""
        raise NotImplementedError("This method is not implemented for this path type.")

    def get_xy_n(self, point, dim):
        """Optional method, can be overridden by subclasses if applicable.
           Only used if this path type defines the netural radius"""
        raise NotImplementedError("This method is not implemented for this path type. You are using a path type which defines the CENTROIDAL radius.")

    def get_dxdy_n(self, point, dim):
        """Optional method, can be overridden by subclasses if applicable.
           Only used if this path type defines the netural radius"""
        raise NotImplementedError("This method is not implemented for this path type. You are using a path type which defines the CENTROIDAL radius.")

    def get_xy_c(self, point, dim):
        """Optional method, can be overridden by subclasses if applicable.
           Only used if this path type defines the centroidal radius"""
        raise NotImplementedError("This method is not implemented for this path type. You are using a path type which defines the NEUTRAL radius.")

    def get_dxdy_c(self, point, dim):
        """Optional method, can be overridden by subclasses if applicable.
           Only used if this path type defines the centroidal radius"""
        raise NotImplementedError("This method is not implemented for this path type. You are using a path type which defines the NEUTRAL radius.")

    def get_neutralSurface(self, resl):
        ximesh = np.linspace(0,self.arcLen,resl+1)
        x = self.get_xy_n(ximesh, 'x')
        y = self.get_xy_n(ximesh, 'y')
        neutralSurface = np.hstack((np.atleast_2d(x).T,
                                    np.atleast_2d(y).T))
        return neutralSurface

    def get_centroidalSurface(self, resl):
        ximesh = np.linspace(0,self.arcLen,resl)
        x = self.get_xy_c(ximesh, 'x')
        y = self.get_xy_c(ximesh, 'x')
        centroidalSurface = np.hstack((np.atleast_2d(x).T,
                                       np.atleast_2d(y).T))
        return centroidalSurface

    @abstractmethod
    def get_dxi_n(self, point):
        pass

    @abstractmethod
    def get_rn(self, point):
        pass

    @abstractmethod
    def get_drn(self, point):
        pass
    
    @abstractmethod
    def get_alpha(self, point):
        pass

    @abstractmethod
    def get_dalpha(self, point):
        pass

    @abstractmethod
    def get_d2alpha(self, point):
        pass

class TestPath(Path):
    def __init__(self, n: int, arcLen: float,
                 startPoint: tuple = None, endPoint: tuple = None,
                 initialRadius: float = None, finalRadius: float = None):
        if (initialRadius is not None) and (finalRadius is not None):
            self.initialRadius = initialRadius
            self.finalRadius   = finalRadius
        elif (initialRadius is not None) ^ (finalRadius is not None):
            raise ValueError("Gotta define both initial AND final radius")
        else:
            self.initialRadius = lin.norm(startPoint)
            self.finalRadius   = lin.norm(endPoint)
        super().__init__(n, arcLen, startPoint, endPoint)

    def calculate_startPoint(self):
        return (self.initialRadius, 0)

    def calculate_endPoint(self):
        return (0, self.finalRadius)

    def get_xy_n(self, point, dim):
        pass
    
    def get_neutralSurface(self, resl):
        pass

    def get_dxi_n(self, point):
        pass

    def get_rn(self, point):
        pass

    def get_drn(self, point):
        pass

    def get_alpha(self, point):
        pass

    def get_dalpha(self, point):
        pass

    def get_d2alpha(self, point):
        pass

    def measure_length(self):
        pass

class LinearRnSpiral(Path):
    def __init__(self, n: int, arcLen: float,
                 # Define both of these
                 startPoint: tuple = None, endPoint: tuple = None,
                 # OR Define both of these, OR Define all four
                 initialRadius: float = None, finalRadius: float = None):
        # Handle whenter rn were given 
        self.startPoint = startPoint
        self.endPoint   = endPoint
        self.calculate_geometryParams(arcLen, startPoint, endPoint, initialRadius, finalRadius)
        super().__init__(n, self.arcLen, startPoint, endPoint)
        self.innerRadius = lin.norm(self.startPoint)
        self.outerRadius = lin.norm(self.endPoint)
        self.startingAngle = np.atan2(self.startPoint[1], self.startPoint[0])

    def calculate_geometryParams(self, arcLen, startPoint, endPoint, initialRadius, finalRadius):
        self.centerOfCurvature=(0,0)
        self.arcLen=arcLen
        if (initialRadius is not None) and (finalRadius is not None):
            self.initialRadius = initialRadius
            self.finalRadius   = finalRadius
            if startPoint is not None:
                self.centerOfCurvature = circle_intersection(initialRadius, finalRadius, startPoint, endPoint)
                vec1 = np.subtract(self.endPoint,self.centerOfCurvature)
                vec2 = np.subtract(self.startPoint,self.centerOfCurvature)
                self.arcLen = np.arccos((vec1.dot(vec2))/(initialRadius*finalRadius))
        elif (initialRadius is not None) ^ (finalRadius is not None):
            raise ValueError("Gotta define both initial AND final radius")
        else:
            self.initialRadius=lin.norm(startPoint)
            self.finalRadius=lin.norm(endPoint)
            self.arcLen = np.atan2(self.endPoint[1], self.endPoint[0])-np.atan2(self.startPoint[1], self.startPoint[0])
         
    def calculate_startPoint(self):
        return (self.initialRadius, 0)

    def calculate_endPoint(self):
        return (0, self.finalRadius)

    def get_xy_n(self, point, dim):
        if dim=='x':
            return self.centerOfCurvature[0]+self.get_rn(point)*np.cos(point+self.startingAngle)
        if dim=='y':
            return self.centerOfCurvature[1]+self.get_rn(point)*np.sin(point+self.startingAngle)
    
    def get_dxdy_n(self, point, dim):
        if dim=='x':
            return -self.get_rn(point)*np.sin(point+self.startingAngle)
        if dim=='y':
            return self.get_rn(point)*np.cos(point+self.startingAngle)

    def get_neutralSurface(self, resl):
        return super().get_neutralSurface(resl)

    def get_dxi_n(self, point):
        dxdxi = self.get_dxdy_n(point, 'x')
        dydxi = self.get_dxdy_n(point, 'y')
        dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
        return dxids

    def get_rn(self, point):
        return self.initialRadius+(self.finalRadius-self.initialRadius)*(point)/self.arcLen

    def get_drn(self, point):
        if hasattr(point, "__len"):
            drndxi = np.ones_like(point)*(self.finalRadius-self.initialRadius)/self.arcLen
        else:
            drndxi = (self.finalRadius-self.initialRadius)/self.arcLen
        return drndxi

    def get_alpha(self, point):
        dydxi = self.get_dxdy_n(point, 'y')
        dxdxi = self.get_dxdy_n(point, 'x')
        return np.atan2(dydxi,dxdxi)

    def get_dalpha(self, point):
        return 1/self.get_rn(point)

    def get_d2alpha(self, point):
        return -self.get_drn(point)/(self.get_rn(point)**2)

    def measure_length(self):
        self.measureResolution = 200
        ximesh            = np.linspace(0,self.arcLen,self.measureResolution+1)

        # calcualte derivative of x and y with respect to xi (internal parameter)
        self.dxdxi = self.get_dxdy_n(ximesh, 'x')
        self.dydxi = self.get_dxdy_n(ximesh, 'y')
        integrand = np.sqrt(self.dxdxi**2+self.dydxi**2)
        # create a mesh in s according by integrating along the parameter to find
        # arc length
        self.smesh = np.zeros(len(ximesh))
        self.smesh[1:len(self.smesh)] = scipy.integrate.cumulative_trapezoid(
                                         integrand, ximesh, ximesh[1]-ximesh[0])
        # print(self.smesh)
        assert(self.smesh[0]==0)
        # return the overall length of the spring
        return self.smesh[-1]

class RadiallyEndedPolynomial(Path):
    def __init__(self, n: int, arcLen: float,
                 # whole bunch of default values:
                 radii: float        = np.array([1.1,2.025,2.025,2.95]),
                 ffradii: float      = np.array([1, 2.5]),
                 alphaAngles         = np.array([45,0])*deg2rad,
                 betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
                 XYFactors         = np.array([0.333,0.667])):
        self.alphaAngles = alphaAngles
        self.init(n, arcLen, radii,
                  betaAngles, XYFactors)
        # iterate spring initialization until xi mesh approximates s mesh as well as possible
        conv=1
        err = 1
        i = 0
        # iterate until you converge on a value
        while conv>10e-6:
            SSProfile("reinit", loud=False).tic()
            errPrev = err
            # measure real arc length each time
            correctedFullLength = self.measure_length()
            # initialize the spring with the remeasured legnth
            self.init(n, correctedFullLength, ### WE ONLY CHANGE THIS
                      radii, betaAngles, XYFactors)
            ximesh = np.linspace(0,correctedFullLength,self.measureResolution+1)
            dxids  = self.get_dxi_n(ximesh)
            err = lin.norm(dxids-np.ones(self.measureResolution+1))
            conv = abs(err-errPrev)
            i+=1
            SSProfile("reinit").toc()
        self.innerRadius=ffradii[0]
        self.outerRadius=ffradii[0]
        
    def init(self,n,arcLen,radii,betaAngles,XYFactors):

        self.pts = np.empty((len(radii),2))
        for i in range(len(radii)):
            self.pts[i,:] = [radii[i]*np.cos(betaAngles[i]),radii[i]*np.sin(betaAngles[i])]
        # print(self.pts)
        startPoint = tuple(self.pts[0,:])
        endPoint = tuple(self.pts[-1,:])

        # Cretae parameter control points for path based on passed factors
        self.XYParamLens = np.empty(len(radii))
        for i in range(len(radii)):
            if i==0:
                self.XYParamLens[i] = 0
            elif i==len(radii)-1:
                self.XYParamLens[i] = arcLen
            else:
                self.XYParamLens[i] = arcLen*XYFactors[i-1]

        super().__init__(n, arcLen, startPoint, endPoint)

        self.n = n
        self.radii = radii
        self.betaAngles = betaAngles
        self.XYFactors = XYFactors

        # generate the coefficients for polynomials
        self.XCoeffs, self.YCoeffs  = self.xy_poly()

    def measure_length(self):
        self.measureResolution = 200
        ximesh            = np.linspace(0,self.arcLen,self.measureResolution+1)

        # calcualte derivative of x and y with respect to xi (internal parameter)
        self.dxdxi = self.get_dxdy_n(ximesh, 'x')
        self.dydxi = self.get_dxdy_n(ximesh, 'y')
        integrand = np.sqrt(self.dxdxi**2+self.dydxi**2)
        # create a mesh in s according by integrating along the parameter to find
        # arc length
        self.smesh = np.zeros(len(ximesh))
        self.smesh[1:len(self.smesh)] = scipy.integrate.cumulative_trapezoid(
                                         integrand, ximesh, ximesh[1]-ximesh[0])
        # print(self.smesh)
        assert(self.smesh[0]==0)
        # return the overall length of the spring
        return self.smesh[-1]

    def prepare_rn_expr(self):
        XCoeffs = self.XCoeffs.flatten()
        YCoeffs = self.YCoeffs.flatten()

        s = sp.symbols('s')
        yPoly = sp.Poly(YCoeffs, s).as_expr()
        xPoly = sp.Poly(XCoeffs, s).as_expr()

        dyds   = sp.diff(yPoly, s)
        dxds   = sp.diff(xPoly, s)
        d2yds2 = sp.diff(yPoly, s, 2)
        d2xds2 = sp.diff(xPoly, s, 2)

        rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
        # print(rn)
        # print("---------------------------------")
        self.drnds_expr = sp.diff(rn, s)
        # print(self.drnds_expr)
        # print(self.fullParamLength)
        self.drnds = sp.lambdify(s, self.drnds_expr, "numpy")

    def xy_poly(self):

        """
        THIS IS THE ONLY FUNCTION THAT CHANGES:

        Instead of being locally radial at either end, we want to be locally
        tangent to the circle
        """

        # if we have n points we need:
        # n constraints on zeroth derivative
        # 2 constraints on first derivative (angle at start and end)
        # 2 constraints on second derivative (straightness at start and end)
        # n + 4 constraints
        # n + 3 degree polynomial

        # initialize matrix sizes
        nConstraints = len(self.pts)+4
        Mat = np.empty((nConstraints,nConstraints))
        YTarg = np.empty((nConstraints,1))
        XTarg = np.empty((nConstraints,1))
        # FIRST ASSIGN 0th DERIVATIVE CONSTRAINTS (c(s)=p)
        for i in range(len(self.XYParamLens)):
            # target correct x-y value
            XTarg[i]=self.pts[i,0]
            YTarg[i]=self.pts[i,1]
            # at associated s value
            for j in range(nConstraints):
                Mat[i,j]   = self.XYParamLens[i]**(nConstraints-1-j)
        # NOW ASSIGN FIRST AND SECOND DERIVATIVE CONSTRAINTS
        for i in range(-1,-5,-1):
            # target correct index (first or last)
            index = (i % 2)*(len(self.XYParamLens)-1)
            # print(index)
            if i < -2:
                # with correct first derivative (radial direction)
                XTarg1 = self.pts[index,0]/lin.norm(self.pts[index,:])
                YTarg1 = self.pts[index,1]/lin.norm(self.pts[index,:])
                TargVect = np.array([[XTarg1],[YTarg1]])
                if index == 0:
                    Mat2 = np.array([[np.cos(self.alphaAngles[index]),-np.sin(self.alphaAngles[index])],
                                    [np.sin(self.alphaAngles[index]), np.cos(self.alphaAngles[index])]])
                    TargVect = Mat2.dot(TargVect)
                else:
                    Mat2 = np.array([[np.cos(self.alphaAngles[-1]),-np.sin(self.alphaAngles[-1])],
                                    [np.sin(self.alphaAngles[-1]), np.cos(self.alphaAngles[-1])]])
                    TargVect = Mat2.dot(TargVect)
                # print(TargVect)
                XTarg[i] = TargVect[0]
                YTarg[i] = TargVect[1]
                # at correct s value
                for j in range(nConstraints):
                    Mat[i,j] = (nConstraints-1-j)* \
                                self.XYParamLens[index]**max(nConstraints-2-j,0) ## FIRST TWO ARE FIRST DERIVATIVES
            else:
                # and with zero second derivative at both points
                XTarg[i] = 0
                YTarg[i] = 0
                for j in range(nConstraints):
                    Mat[i,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)* \
                                self.XYParamLens[index]**max(nConstraints-3-j,0) ## LAST TWO ARE SECOND DERIVATIVES
        # print(XTarg)
        # print(YTarg)
        # print(Mat)
        XCoeffs = lin.solve(Mat, XTarg)
        YCoeffs = lin.solve(Mat, YTarg)
        # print(XCoeffs)
        # print(YCoeffs)
        # print(lin.cond(Mat))
        for i in range(len(self.XYParamLens)):
            diffX = PPoly_Eval(self.XYParamLens[i],XCoeffs)-self.pts[i,0]
            diffY = PPoly_Eval(self.XYParamLens[i],YCoeffs)-self.pts[i,1]
            diff = lin.norm([diffX,diffY])
            # print(diff)
            assert(np.isclose(diff,0))
        # TODO: assert(np.isclose(diff,0))
        # print(XCoeffs, YCoeffs)
        return XCoeffs, YCoeffs

    def get_xy_n(self, point, dim):
        if dim=='x':
            out = PPoly_Eval(point, self.XCoeffs)
        if dim=='y':
            out = PPoly_Eval(point, self.YCoeffs)
        return out

    def get_dxdy_n(self, point, dim):
        if dim=='x':
            out = PPoly_Eval(point, self.XCoeffs, deriv=1)
        if dim=='y':
            out = PPoly_Eval(point, self.YCoeffs, deriv=1)
        return out

    def get_neutralSurface(self, resl):
        return super().get_neutralSurface(resl)

    def get_dxi_n(self, point):
        dxdxi = self.get_dxdy_n(point, 'x')
        dydxi = self.get_dxdy_n(point, 'y')
        dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
        return dxids

    def get_rn(self, point):
        point=np.atleast_1d(point)
        d2yds2 = PPoly_Eval(point, self.YCoeffs, deriv=2)
        d2xds2 = PPoly_Eval(point, self.XCoeffs, deriv=2)
        dyds   = PPoly_Eval(point, self.YCoeffs, deriv=1)
        dxds   = PPoly_Eval(point, self.XCoeffs, deriv=1)
        # dAlphadS = (d2yds2/dxds-d2xds2*dyds/dxds)/(1+dyds**2/dxds)
        rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
        for i in range(len(rn)):
            if not np.isfinite(rn[i]):
                rn[i] = float('inf')
            if abs((d2yds2[i]/dxds[i]-d2xds2[i]*dyds[i]/dxds[i]**2)) <= 10**-13:
                rn[i] = float('inf')*np.sign(rn[i-1])
        if len(rn)==1:
            return rn[0]
        else:
            return rn

    def get_drn(self, point):
        point = np.atleast_1d(point)
        out = self.drnds(point)
        if len(out)==1:
            out = out[0]
        if np.isnan(out):
            out = 0
        return out
    
    def get_alpha(self, point):
        point = np.atleast_1d(point)
        # print(point)
        alphaList = np.arctan2(PPoly_Eval(point, self.YCoeffs, deriv=1),
                               PPoly_Eval(point, self.XCoeffs, deriv=1))
        
        for i in range(len(alphaList)):
            if alphaList[i]<0:
                alphaList[i]=alphaList[i]+2*np.pi
        if len(alphaList)==1:
            return alphaList[0]
        else:
            return alphaList

    def get_dalpha(self, point):
        if hasattr(point, "__len__"):
            alphaList = np.arctan2(PPoly_Eval(point, self.YCoeffs, deriv=1),
                                       PPoly_Eval(point, self.XCoeffs, deriv=1))
            for i in range(len(alphaList)):
                if alphaList[i]<0:
                    alphaList[i]=alphaList[i]+2*np.pi
            return alphaList
        else:
            alpha = np.arctan2(PPoly_Eval(point, self.YCoeffs, deriv=1),
                                       PPoly_Eval(point, self.XCoeffs, deriv=1))
            if alpha<0:
                alpha=alpha+2*np.pi
            return alpha

    def get_d2alpha(self, point):
        d2ads2 = -self.get_dalpha(point)**2*self.get_drn(point)
        if np.isnan(d2ads2):
            d2ads2 = 0
        return d2ads2

## LEGACY BELOW THIS POINT ##

"""
##############################################
# STANDARD INTERFACE FOR ALL PalongTHDEF OBJECTS #
##############################################
# Methods to include:
#
#   get_crscRef(crscRef):
#           crscRef as a CRSCDEF object
#   
#   get_xy_n(coord, dim):
#           coord as scalar or vector float
#           dim as str 'x' or 'y'
#
#   get_dxdy_n(coord, dim):
#
#   get_dxi_n(coord):
#
#   get_rn(coord):
#
#   get_drn(coord):
#
#   get_alpha(coord):
#
#   get_dalpha(coord):
#
#   get_d2alpha(coord):
#
#   get_neutralSurface(resolution):
#           resolution as int
#
#   get_centroidalSurface(resolution):
#
# Attributes to include:
#
#   self.x0
#        y0
#        xL
#        yL
#        momentArmX
#        momentArmY
#
#   self.n
#
#   self.fullParamLength
#
#   self.innerRadius
#        outerRadius
#
#
"""

# class Empty:
#     def __init__(self, n, x0, y0, arcLen):
#         self.returnValue = 0

#         self.n = n

#         self.x0 = x0
#         self.y0 = y0
#         self.xL = x0
#         self.yL = y0+arcLen
#         self.momentArmX = self.xL - self.x0
#         self.momentArmY = self.yL - self.y0

#         self.fullParamLength = arcLen
#         self.innerRadius     = lin.norm([x0,y0])
#         self.outerRadius     = lin.norm([self.xL,self.yL])

#     def measure_length(self):
#         return self.returnValue

#     ########################
#     ## Standard Interface ##
#     ########################

#     def get_crscRef(self, crscRef):
#         self.crsc = crscRef

#     def get_xy_n(self, xi, dim):
#         return self.returnValue

#     def get_dxdy_n(self, xi, dim):
#         return self.returnValue

#     def get_rn(self, xi):
#         return self.returnValue

#     def get_drn(self, xi):
#         return self.returnValue
    
#     def get_dxi_n(self, xi):
#         # what am I going to do here
#         return 1
    
#     def get_alpha(self, xi):
#         return self.returnValue

#     def get_dalpha(self, xi):
#         return self.returnValue
    
#     def get_d2alpha(self, xi):
#         return self.returnValue
    
#     def get_neutralSurface(self, resolution):
#         return self.returnValue
    
#     def get_centroidalSurface(self):
#         return self.returnValue

# class Linear_Rn_Spiral:
#     def __init__(self, start=1, end=2, xiRange = np.pi/2, n=2):

#         self.n = n

#         self.start = start
#         self.end   = end
#         self.xiRange = xiRange

#         self.drndxi = (self.end-self.start)/(xiRange)

#         self.x0 = self.get_xy_n(0, 'x')
#         self.y0 = self.get_xy_n(0, 'y')
#         self.xL = self.get_xy_n(xiRange, 'x')
#         self.yL = self.get_xy_n(xiRange, 'y')
#         self.momentArmX = self.xL - self.x0
#         self.momentArmY = self.yL - self.y0

#         self.fullParamLength = xiRange
#         self.innerRadius  = start
#         self.outerRadius  = end

#     def measure_length(self):
#         self.measureResolution = 200
#         ximesh            = np.linspace(0,self.xiRange,self.measureResolution+1)

#         # calcualte derivative of x and y with respect to xi (internal parameter)
#         self.dxdxi = self.get_dxi_n(ximesh, 'x')
#         self.dydxi = self.get_dxi_n(ximesh, 'y')
#         integrand = np.sqrt(self.dxdxi**2+self.dydxi**2)
#         # create a mesh in s according by integrating along the parameter to find
#         # arc length
#         self.smesh = np.zeros(len(ximesh))
#         self.smesh[1:len(self.smesh)] = scipy.integrate.cumulative_trapezoid(
#                                          integrand, ximesh, ximesh[1]-ximesh[0])
#         # print(self.smesh)
#         assert(self.smesh[0]==0)
#         # return the overall length of the spring
#         return self.smesh[-1]

#     ########################
#     ## Standard Interface ##
#     ########################

#     def get_crscRef(self, crscRef):
#         self.crsc = crscRef

#     def get_xy_n(self, xi, dim):
#         if dim=='x':
#             return self.get_rn(xi)*np.cos(xi)
#         if dim=='y':
#             return self.get_rn(xi)*np.sin(xi)

#     def get_dxdy_n(self, xi, dim):
#         if dim=='x':
#             return -self.get_rn(xi)*np.sin(xi)
#         if dim=='y':
#             return self.get_rn(xi)*np.cos(xi)

#     def get_rn(self, xi):
#         return self.start+self.drndxi*xi

#     def get_drn(self, xi):
#         return self.drndxi
    
#     def get_dxi_n(self, xi):
#         dxdxi = self.get_dxdy_n(xi, 'x')
#         dydxi = self.get_dxdy_n(xi, 'y')
#         dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
#         return dxids
    
#     def get_alpha(self, xi):
#         dydxi = self.get_dxdy_n(xi, 'y')
#         dxdxi = self.get_dxdy_n(xi, 'x')
#         return np.atan2(dydxi,dxdxi)

#     def get_dalpha(self, xi):
#         return 1
    
#     def get_d2alpha(self, xi):
#         return 0
    
#     def get_neutralSurface(self, resolution):
#         ximesh = np.linspace(0,self.xiRange,resolution+1)
#         self.undeformedNeutralSurface = np.hstack(
#                                   (np.atleast_2d(self.get_xy_n(ximesh, 'x')).T,
#                                    np.atleast_2d(self.get_xy_n(ximesh, 'y')).T))
#         return self.undeformedNeutralSurface
    
#     def get_centroidalSurface(self, resolution):
#         ximesh = np.linspace(0,self.fullParamLength,resolution+1)
#         Ic = self.crsc.get_Ic(ximesh)
#         t = self.crsc.t
#         h = self.crsc.get_Thk(ximesh)
#         rn = self.get_rn(ximesh)
#         ecc = Ic/(t*h*rn)
#         alpha = self.get_alpha(ximesh)

#         neutralSurface = self.get_neutralSurface(resolution)
#         self.undeformedCentroidalSurface = neutralSurface + \
#                        np.hstack((np.atleast_2d(ecc*np.sin(alpha)).T,
#                                   np.atleast_2d(ecc*np.cos(alpha)).T))

#         return self.undeformedCentroidalSurface

# class Minimal_Polynomial_Definition:
#     def __init__(self,
#                  # whole bunch of default values:
#                  n                   = 2,
#                  fullParamLength     = 6,
#                  radii               = np.array([1.1,2.025,2.025,2.95]),
#                  betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
#                  XYFactors         = np.array([0.333,0.667])
#                  ):

#         self.init(n, fullParamLength, radii,
#                   betaAngles, XYFactors)
#         # iterate spring initialization until xi mesh approximates s mesh as well as possible
#         conv=1
#         err = 1
#         i = 0
#         # iterate until you converge on a value
#         while conv>10e-6:
#             SSProfile("reinit").tic()
#             errPrev = err
#             # measure real arc length each time
#             correctedFullLength = self.measure_length()
#             # initialize the spring with the remeasured legnth
#             self.init(n, correctedFullLength, ### WE ONLY CHANGE THIS
#                       radii, betaAngles, XYFactors)
#             ximesh = np.linspace(0,correctedFullLength,self.measureResolution+1)
#             dxids  = self.get_dxi_n(ximesh)
#             err = lin.norm(dxids-np.ones(self.measureResolution+1))
#             # err = self.fullParamLength-self.measure_length() # These two error values seem equivalent
#             conv = abs(err-errPrev)
#             i+=1
#             SSProfile("reinit").toc()
#         # if the spring was given a name, save its input parameters in a filej

#         # this takes way too many lines but its better to see the math written
#         # out:

#         # Find coordinates of frames:
#         # at base of spring
#         self.x0 = PPoly_Eval(0,self.XCoeffs)
#         self.y0 = PPoly_Eval(0,self.YCoeffs)
#         # at end of spring
#         self.xL = PPoly_Eval(self.fullParamLength,self.XCoeffs)
#         self.yL = PPoly_Eval(self.fullParamLength,self.YCoeffs)
#         # at end of spring, interms of frame at base of spring
#         self.momentArmX = self.xL - self.x0
#         self.momentArmY = self.yL - self.y0

#         self.innerRadius= radii[0]
#         self.outerRadius= radii[-1]

#         self.saveVariables = {'n':     n,     'length':      correctedFullLength,
#                               'radii': radii, 'beta angles': betaAngles,
#                               'control factors': XYFactors}
#         # np.savetxt("springs/")

#     def init(self,
#              n                   = 2,
#              fullParamLength     = 6,
#              radii               = np.array([1.1,2.025,2.025,2.95]),
#              betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
#              XYParamLens           = np.array([0.333,0.667])):       # radii - 2

#         """
#         THIS FUNCTION:
#             Takes in default or custom parameters and
#             defines the quantities and function needed for analysis
#         """

#         ############################ PARAMETER KEY ############################
#         # E                 - youngs modulus                                  #
#         # designStress      - max allowable stress, currently 80% yield       #
#         # n                 - # of flexures                                   #
#         # fullParamLength   - maximum value of internal parameter             #
#         # outPlaneThickness - axial thickness of spring (in)                  #
#         # radii             - radii of control points of netural surface (in) #
#         # betaAngles        - angle of control points of neutral surface (rad)#
#         # IcPts             - curved second moment of area of beam at points  #
#         #                     along beam (in^4)                               #
#         # IcParamLens       - internal parameter values at which Ic control   #
#         #                     points apply                                    #
#         # IcParamLens       - internal parameter values at which x-y          #
#         #                     (radii-angle) control points apply              #
#         # resolution        - number of points in mesh                        #
#         ########################################################################

#         # stick all the arguments in the object
#         self.n = n
#         self.fullParamLength = fullParamLength
#         self.radii = radii
#         self.betaAngles = betaAngles
#         self.XYFactors = XYParamLens

#         # create control points for polynomials

#         # x-y control points
#         self.pts = np.empty((len(radii),2))
#         for i in range(len(radii)):
#             self.pts[i,:] = [self.radii[i]*np.cos(self.betaAngles[i]),self.radii[i]*np.sin(self.betaAngles[i])]

#         # Cretae parameter control points for path based on passed factors
#         self.XYParamLens = np.empty(len(radii))
#         for i in range(len(radii)):
#             if i==0:
#                 self.XYParamLens[i] = 0
#             elif i==len(radii)-1:
#                 self.XYParamLens[i] = self.fullParamLength
#             else:
#                 self.XYParamLens[i] = self.fullParamLength*XYParamLens[i-1]

#         # create a vector of all mutable parameters
#         self.parameterVector = np.concatenate((self.radii, self.betaAngles,
#                                               self.XYFactors))

#         # generate the coefficients for polynomials
#         self.XCoeffs, self.YCoeffs  = self.xy_poly()

#     def xy_poly(self):

#         # if we have n points we need:
#         # n constraints on zeroth derivative
#         # 2 constraints on first derivative (angle at start and end)
#         # 2 constraints on second derivative (straightness at start and end)
#         # n + 4 constraints
#         # n + 3 degree polynomial

#         # initialize matrix sizes
#         nConstraints = len(self.pts)+4
#         Mat = np.empty((nConstraints,nConstraints))
#         YTarg = np.empty((nConstraints,1))
#         XTarg = np.empty((nConstraints,1))
#         # FIRST ASSIGN 0th DERIVATIVE CONSTRAINTS (c(s)=p)
#         for i in range(len(self.XYParamLens)):
#             # target correct x-y value
#             XTarg[i]=self.pts[i,0]
#             YTarg[i]=self.pts[i,1]
#             # at associated s value

#             for j in range(nConstraints):
#                 Mat[i,j]   = self.XYParamLens[i]**(nConstraints-1-j)
#         # NOW ASSIGN FIRST AND SECOND DERIVATIVE CONSTRAINTS
#         for i in range(-1,-5,-1):
#             # target correct index (first or last)
#             index = (i % 2)*(len(self.XYParamLens)-1)
#             # print(index)
#             if i < -2:
#                 # with correct first derivative (radial direction)
#                 XTarg[i] = self.pts[index,0]/lin.norm(self.pts[index,:])
#                 YTarg[i] = self.pts[index,1]/lin.norm(self.pts[index,:])
#                 # at correct s value
#                 for j in range(nConstraints):
#                     Mat[i,j] = (nConstraints-1-j)*self.XYParamLens[index]**max(nConstraints-2-j,0) ## FIRST TWO ARE FIRST DERIVATIVES
#             else:
#                 # and with zero second derivative at both points
#                 XTarg[i] = 0
#                 YTarg[i] = 0
#                 for j in range(nConstraints):
#                     Mat[i,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)*self.XYParamLens[index]**max(nConstraints-3-j,0) ## LAST TWO ARE SECOND DERIVATIVES
#         # print(XTarg)
#         # print(YTarg)
#         # print(Mat)
#         XCoeffs = lin.solve(Mat, XTarg)
#         YCoeffs = lin.solve(Mat, YTarg)
#         # print(lin.cond(Mat))
#         for i in range(len(self.XYParamLens)):
#             diffX = PPoly_Eval(self.XYParamLens[i],XCoeffs)-self.pts[i,0]
#             diffY = PPoly_Eval(self.XYParamLens[i],YCoeffs)-self.pts[i,1]
#             diff = lin.norm([diffX,diffY])
#             assert(np.isclose(diff,0))
#         # TODO: assert(np.isclose(diff,0))
#         # print(XCoeffs, YCoeffs)
#         return XCoeffs, YCoeffs

#     def measure_length(self):

#         """
#         THIS FUNCTION:
#             measures the length of an initialized spring
#         """

#         self.measureResolution = 200
#         ximesh            = np.linspace(0,self.fullParamLength,self.measureResolution+1)

#         # calcualte derivative of x and y with respect to xi (internal parameter)
#         self.dxdxi = PPoly_Eval(ximesh, self.XCoeffs, deriv = 1)
#         self.dydxi = PPoly_Eval(ximesh, self.YCoeffs, deriv = 1)
#         integrand = np.sqrt(self.dxdxi**2+self.dydxi**2)
#         # create a mesh in s according by integrating along the parameter to find
#         # arc length
#         self.smesh = np.zeros(len(ximesh))
#         self.smesh[1:len(self.smesh)] = scipy.integrate.cumulative_trapezoid(
#                                          integrand, ximesh, ximesh[1]-ximesh[0])
#         # print(self.smesh)
#         assert(self.smesh[0]==0)
#         # return the overall length of the spring
#         return self.smesh[-1]

# #### STANDARD INTERFACE ####

#     def get_crscRef(self, crscRef):
#         self.crsc = crscRef

#     def get_xy_n(self, coord, dim):
#         if dim=='x':
#             out = PPoly_Eval(coord, self.XCoeffs)
#         if dim=='y':
#             out = PPoly_Eval(coord, self.YCoeffs)
#         return out

#     def get_dxdy_n(self, coord, dim):
#         if dim=='x':
#             out = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         if dim=='y':
#             out = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         return out

#     def get_dxi_n(self, coord):
#         dxdxi = self.get_dxdy_n(coord, 'x')
#         dydxi = self.get_dxdy_n(coord, 'y')
#         dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
#         return dxids

#     def get_rn(self, coord):
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         # dAlphadS = (d2yds2/dxds-d2xds2*dyds/dxds)/(1+dyds**2/dxds)
#         if hasattr(coord, "__len__"):
#             rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#             for i in range(len(rn)):
#                 if not np.isfinite(rn[i]):
#                     rn[i] = float('inf')
#                 if abs((d2yds2[i]/dxds[i]-d2xds2[i]*dyds[i]/dxds[i]**2)) <= 10**-13:
#                     rn[i] = float('inf')*np.sign(rn[i-1])
#         else:
#             if abs(d2yds2/dxds-d2xds2*dyds/dxds**2)<=10**-13:
#                 rn = float('inf')
#             else:
#                 rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#         return rn

#     def get_drn(self, coord):
#         d3yds3 = PPoly_Eval(coord, self.YCoeffs, deriv=3)
#         d3xds3 = PPoly_Eval(coord, self.XCoeffs, deriv=3)
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)

#         denominator = (d2yds2*dxds-d2xds2*dyds)**2
#         numerator   = ( dxds**3*d3yds3 - dyds**3*d3xds3 +
#                         2*dyds*d2xds2**2*dxds - 2*dyds*d2yds2**2*dxds - dxds**2*d3xds3*dyds -
#                         2*dxds**2*d2yds2*d2xds2 + 2*dyds**2*d2yds2*d2xds2 +
#                         dyds**2*d3yds3*dxds )

#         drnds = -numerator/denominator

#         return drnds

#     def get_alpha(self, coord):
#         if hasattr(coord, "__len__"):
#             alphaList = np.arctan2(PPoly_Eval(coord, self.YCoeffs, deriv=1),
#                                        PPoly_Eval(coord, self.XCoeffs, deriv=1))
#             for i in range(len(alphaList)):
#                 if alphaList[i]<0:
#                     alphaList[i]=alphaList[i]+2*np.pi
#             return alphaList
#         else:
#             alpha = np.arctan2(PPoly_Eval(coord, self.YCoeffs, deriv=1),
#                                        PPoly_Eval(coord, self.XCoeffs, deriv=1))
#             if alpha<0:
#                 alpha=alpha+2*np.pi
#             return alpha

#     def get_dalpha(self, coord):
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         dads = ((d2yds2/dxds-d2xds2*dyds/dxds**2)/(1+dyds**2/dxds**2))
#         return dads

#     def get_d2alpha(self, coord):
#         d2ads2 = -self.get_dalpha(coord)**2*self.get_drn(coord)
#         if np.isnan(d2ads2):
#             d2ads2 = 0
#         return d2ads2

#     def get_neutralSurface(self, resolution):

#         ximesh = np.linspace(0,self.xiRange,resolution+1)
#         self.undeformedNeutralSurface = np.hstack(
#                                   (np.atleast_2d(self.get_xy_n(ximesh, 'x')).T,
#                                    np.atleast_2d(self.get_xy_n(ximesh, 'y')).T))

#         return self.undeformedNeutralSurface

#     def get_centroidalSurface(self, resolution):
#         ximesh = np.linspace(0,self.fullParamLength,resolution+1)
#         Ic = self.crsc.get_Ic(ximesh)
#         t = self.crsc.t
#         h = self.crsc.get_Thk(ximesh)
#         rn = self.get_rn(ximesh)
#         ecc = Ic/(t*h*rn)
#         alpha = self.get_alpha(ximesh)

#         neutralSurface = self.get_neutralSurface(resolution)
#         self.undeformedCentroidalSurface = neutralSurface + \
#                        np.hstack((np.atleast_2d(ecc*np.sin(alpha)).T,
#                                   np.atleast_2d(ecc*np.cos(alpha)).T))

#         return self.undeformedCentroidalSurface

# class Minimal_Polynomial_Definition2:
#     def __init__(self,
#                  # whole bunch of default values:
#                  n                   = 2,
#                  fullParamLength     = 6,
#                  radii               = np.array([1.1,2.025,2.025,2.95]),
#                  ffradii             = np.array([1, 2.5]),
#                  betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
#                  XYFactors         = np.array([0.333,0.667])
#                  ):

#         self.init(n, fullParamLength, radii,
#                   betaAngles, XYFactors)
#         # iterate spring initialization until xi mesh approximates s mesh as well as possible
#         conv=1
#         err = 1
#         i = 0
#         # iterate until you converge on a value
#         while conv>10e-6:
#             SSProfile("reinit").tic()
#             errPrev = err
#             # measure real arc length each time
#             correctedFullLength = self.measure_length()
#             # initialize the spring with the remeasured legnth
#             self.init(n, correctedFullLength, ### WE ONLY CHANGE THIS
#                       radii, betaAngles, XYFactors)
#             ximesh = np.linspace(0,correctedFullLength,self.measureResolution+1)
#             dxids  = self.get_dxi_n(ximesh)
#             err = lin.norm(dxids-np.ones(self.measureResolution+1))
#             # err = self.fullParamLength-self.measure_length() # These two error values seem equivalent
#             conv = abs(err-errPrev)
#             i+=1
#             SSProfile("reinit").toc()
#         # if the spring was given a name, save its input parameters in a filej

#         # this takes way too many lines but its better to see the math written
#         # out:

#         # Find coordinates of frames:
#         # at base of spring
#         self.x0 = PPoly_Eval(0,self.XCoeffs)
#         self.y0 = PPoly_Eval(0,self.YCoeffs)
#         # at end of spring
#         self.xL = PPoly_Eval(self.fullParamLength,self.XCoeffs)
#         self.yL = PPoly_Eval(self.fullParamLength,self.YCoeffs)
#         # at end of spring, interms of frame at base of spring
#         self.momentArmX = self.xL - self.x0
#         self.momentArmY = self.yL - self.y0

#         self.innerRadius= ffradii[0]
#         self.outerRadius= ffradii[-1]

#         self.saveVariables = {'n':     n,     'length':      correctedFullLength,
#                               'radii': radii, 'beta angles': betaAngles,
#                               'control factors': XYFactors}

#         # print(self.XCoeffs)
#         # print(self.YCoeffs)

#         # self.get_neutralSurface(200)
#         # plt.plot(self.undeformedNeutralSurface[:,0], self.undeformedNeutralSurface[:,1])
#         # plt.show()
#         # np.savetxt("springs/")

#     def init(self,
#              n                   = 2,
#              fullParamLength     = 6,
#              radii               = np.array([1.1,2.025,2.025,2.95]),
#              betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
#              XYParamLens           = np.array([0.333,0.667])):       # radii - 2

#         """
#         THIS FUNCTION:
#             Takes in default or custom parameters and
#             defines the quantities and function needed for analysis
#         """

#         ############################ PARAMETER KEY ############################
#         # E                 - youngs modulus                                  #
#         # designStress      - max allowable stress, currently 80% yield       #
#         # n                 - # of flexures                                   #
#         # fullParamLength   - maximum value of internal parameter             #
#         # outPlaneThickness - axial thickness of spring (in)                  #
#         # radii             - radii of control points of netural surface (in) #
#         # betaAngles        - angle of control points of neutral surface (rad)#
#         # IcPts             - curved second moment of area of beam at points  #
#         #                     along beam (in^4)                               #
#         # IcParamLens       - internal parameter values at which Ic control   #
#         #                     points apply                                    #
#         # IcParamLens       - internal parameter values at which x-y          #
#         #                     (radii-angle) control points apply              #
#         # resolution        - number of points in mesh                        #
#         ########################################################################

#         # stick all the arguments in the object
#         self.n = n
#         self.fullParamLength = fullParamLength
#         self.radii = radii
#         self.betaAngles = betaAngles
#         self.XYFactors = XYParamLens

#         # create control points for polynomials

#         # x-y control points
#         self.pts = np.empty((len(radii),2))
#         for i in range(len(radii)):
#             self.pts[i,:] = [self.radii[i]*np.cos(self.betaAngles[i]),self.radii[i]*np.sin(self.betaAngles[i])]

#         # Cretae parameter control points for path based on passed factors
#         self.XYParamLens = np.empty(len(radii))
#         for i in range(len(radii)):
#             if i==0:
#                 self.XYParamLens[i] = 0
#             elif i==len(radii)-1:
#                 self.XYParamLens[i] = self.fullParamLength
#             else:
#                 self.XYParamLens[i] = self.fullParamLength*XYParamLens[i-1]

#         # create a vector of all mutable parameters
#         self.parameterVector = np.concatenate((self.radii, self.betaAngles,
#                                               self.XYFactors))

#         # generate the coefficients for polynomials
#         self.XCoeffs, self.YCoeffs  = self.xy_poly()

#     def xy_poly(self):

#         """
#         THIS IS THE ONLY FUNCTION THAT CHANGES:

#         Instead of being locally radial at either end, we want to be locally
#         tangent to the circle
#         """

#         # if we have n points we need:
#         # n constraints on zeroth derivative
#         # 2 constraints on first derivative (angle at start and end)
#         # 2 constraints on second derivative (straightness at start and end)
#         # n + 4 constraints
#         # n + 3 degree polynomial

#         # initialize matrix sizes
#         nConstraints = len(self.pts)+4
#         Mat = np.empty((nConstraints,nConstraints))
#         YTarg = np.empty((nConstraints,1))
#         XTarg = np.empty((nConstraints,1))
#         # FIRST ASSIGN 0th DERIVATIVE CONSTRAINTS (c(s)=p)
#         for i in range(len(self.XYParamLens)):
#             # target correct x-y value
#             XTarg[i]=self.pts[i,0]
#             YTarg[i]=self.pts[i,1]
#             # at associated s value
#             for j in range(nConstraints):
#                 Mat[i,j]   = self.XYParamLens[i]**(nConstraints-1-j)
#         # NOW ASSIGN FIRST AND SECOND DERIVATIVE CONSTRAINTS
#         for i in range(-1,-5,-1):
#             # target correct index (first or last)
#             index = (i % 2)*(len(self.XYParamLens)-1)
#             # print(index)
#             if i < -2:
#                 # with correct first derivative (radial direction)
#                 XTarg1 = self.pts[index,0]/lin.norm(self.pts[index,:])
#                 YTarg1 = self.pts[index,1]/lin.norm(self.pts[index,:])
#                 # rotate 90 degrees ccw
#                 XTarg[i] = -1*YTarg1
#                 YTarg[i] =  1*XTarg1
#                 # at correct s value
#                 for j in range(nConstraints):
#                     Mat[i,j] = (nConstraints-1-j)* \
#                                 self.XYParamLens[index]**max(nConstraints-2-j,0) ## FIRST TWO ARE FIRST DERIVATIVES
#             else:
#                 # and with zero second derivative at both points
#                 XTarg[i] = 0
#                 YTarg[i] = 0
#                 for j in range(nConstraints):
#                     Mat[i,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)* \
#                                 self.XYParamLens[index]**max(nConstraints-3-j,0) ## LAST TWO ARE SECOND DERIVATIVES
#         # print(XTarg)
#         # print(YTarg)
#         # print(Mat)
#         XCoeffs = lin.solve(Mat, XTarg)
#         YCoeffs = lin.solve(Mat, YTarg)
#         # print(lin.cond(Mat))
#         for i in range(len(self.XYParamLens)):
#             diffX = PPoly_Eval(self.XYParamLens[i],XCoeffs)-self.pts[i,0]
#             diffY = PPoly_Eval(self.XYParamLens[i],YCoeffs)-self.pts[i,1]
#             diff = lin.norm([diffX,diffY])
#             assert(np.isclose(diff,0))
#         # TODO: assert(np.isclose(diff,0))
#         # print(XCoeffs, YCoeffs)
#         return XCoeffs, YCoeffs

#     def measure_length(self):

#         """
#         THIS FUNCTION:
#             measures the length of an initialized spring
#         """

#         self.measureResolution = 200
#         ximesh            = np.linspace(0,self.fullParamLength,self.measureResolution+1)

#         # calcualte derivative of x and y with respect to xi (internal parameter)
#         self.dxdxi = PPoly_Eval(ximesh, self.XCoeffs, deriv = 1)
#         self.dydxi = PPoly_Eval(ximesh, self.YCoeffs, deriv = 1)
#         integrand = np.sqrt(self.dxdxi**2+self.dydxi**2)
#         # create a mesh in s according by integrating along the parameter to find
#         # arc length
#         self.smesh = np.zeros(len(ximesh))
#         self.smesh[1:len(self.smesh)] = scipy.integrate.cumulative_trapezoid(
#                                          integrand, ximesh, ximesh[1]-ximesh[0])
#         # print(self.smesh)
#         assert(self.smesh[0]==0)
#         # return the overall length of the spring
#         return self.smesh[-1]

# #### STANDARD INTERFACE ####

#     def get_crscRef(self, crscRef):
#         self.crsc = crscRef

#     def get_xy_n(self, coord, dim):
#         if dim=='x':
#             out = PPoly_Eval(coord, self.XCoeffs)
#         if dim=='y':
#             out = PPoly_Eval(coord, self.YCoeffs)
#         return out

#     def get_dxdy_n(self, coord, dim):
#         if dim=='x':
#             out = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         if dim=='y':
#             out = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         return out

#     def get_dxi_n(self, coord):
#         dxdxi = self.get_dxdy_n(coord, 'x')
#         dydxi = self.get_dxdy_n(coord, 'y')
#         dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
#         return dxids

#     def get_rn(self, coord):
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         # dAlphadS = (d2yds2/dxds-d2xds2*dyds/dxds)/(1+dyds**2/dxds)
#         if hasattr(coord, "__len__"):
#             rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#             for i in range(len(rn)):
#                 if not np.isfinite(rn[i]):
#                     rn[i] = float('inf')
#                 if abs((d2yds2[i]/dxds[i]-d2xds2[i]*dyds[i]/dxds[i]**2)) <= 10**-13:
#                     rn[i] = float('inf')*np.sign(rn[i-1])
#         else:
#             if abs(d2yds2/dxds-d2xds2*dyds/dxds**2)<=10**-13 or not np.isfinite(abs(d2yds2/dxds-d2xds2*dyds/dxds**2)):
#                 rn = float('inf')
#             else:
#                 rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#         return rn

#     def get_drn(self, coord):
#         d3yds3 = PPoly_Eval(coord, self.YCoeffs, deriv=3)
#         d3xds3 = PPoly_Eval(coord, self.XCoeffs, deriv=3)
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)

#         denominator = (d2yds2*dxds-d2xds2*dyds)**2
#         numerator   = ( dxds**3*d3yds3 - dyds**3*d3xds3 +
#                         2*dyds*d2xds2**2*dxds - 2*dyds*d2yds2**2*dxds - dxds**2*d3xds3*dyds -
#                         2*dxds**2*d2yds2*d2xds2 + 2*dyds**2*d2yds2*d2xds2 +
#                         dyds**2*d3yds3*dxds )

#         drnds = -numerator/denominator

#         return drnds

#     def get_alpha(self, coord):
#         if hasattr(coord, "__len__"):
#             alphaList = np.arctan2(PPoly_Eval(coord, self.YCoeffs, deriv=1),
#                                        PPoly_Eval(coord, self.XCoeffs, deriv=1))
#             for i in range(len(alphaList)):
#                 if alphaList[i]<0:
#                     alphaList[i]=alphaList[i]+2*np.pi
#             return alphaList
#         else:
#             alpha = np.arctan2(PPoly_Eval(coord, self.YCoeffs, deriv=1),
#                                        PPoly_Eval(coord, self.XCoeffs, deriv=1))
#             if alpha<0:
#                 alpha=alpha+2*np.pi
#             return alpha

#     def get_dalpha(self, coord):
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         dads = ((d2yds2/dxds-d2xds2*dyds/dxds**2)/(1+dyds**2/dxds**2))
#         return dads

#     def get_d2alpha(self, coord):
#         d2ads2 = -self.get_dalpha(coord)**2*self.get_drn(coord)
#         if np.isnan(d2ads2):
#             d2ads2 = 0
#         return d2ads2

#     def get_neutralSurface(self, resolution):

#         ximesh = np.linspace(0,self.fullParamLength,resolution+1)
#         self.undeformedNeutralSurface = np.hstack(
#                                   (np.atleast_2d(self.get_xy_n(ximesh, 'x')).T,
#                                    np.atleast_2d(self.get_xy_n(ximesh, 'y')).T))

#         return self.undeformedNeutralSurface

#     def get_centroidalSurface(self, resolution):
#         ximesh = np.linspace(0,self.fullParamLength,resolution+1)
#         Ic = self.crsc.get_Ic(ximesh)
#         t = self.crsc.t
#         h = self.crsc.get_Thk(ximesh, hasPrev=False)
#         rn = self.get_rn(ximesh)
#         ecc = Ic/(t*h*rn)
#         alpha = self.get_alpha(ximesh)

#         neutralSurface = self.get_neutralSurface(resolution)
#         self.undeformedCentroidalSurface = neutralSurface + \
#                        np.hstack((np.atleast_2d(ecc*np.sin(alpha)).T,
#                                   np.atleast_2d(ecc*np.cos(alpha)).T))

#         return self.undeformedCentroidalSurface

# class Minimal_Polynomial_Definition3:
#     def __init__(self,
#                  # whole bunch of default values:
#                  n                   = 2,
#                  fullParamLength     = 6,
#                  radii               = np.array([1.1,2.025,2.025,2.95]),
#                  ffradii             = np.array([1, 2.5]),
#                  betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
#                  XYFactors         = np.array([0.333,0.667])
#                  ):

#         self.init(n, fullParamLength, radii,
#                   betaAngles, XYFactors)
#         # iterate spring initialization until xi mesh approximates s mesh as well as possible
#         conv=1
#         err = 1
#         i = 0
#         # iterate until you converge on a value
#         while conv>10e-6:
#             SSProfile("reinit").tic()
#             errPrev = err
#             # measure real arc length each time
#             correctedFullLength = self.measure_length()
#             # initialize the spring with the remeasured legnth
#             self.init(n, correctedFullLength, ### WE ONLY CHANGE THIS
#                       radii, betaAngles, XYFactors)
#             ximesh = np.linspace(0,correctedFullLength,self.measureResolution+1)
#             dxids  = self.get_dxi_n(ximesh)
#             err = lin.norm(dxids-np.ones(self.measureResolution+1))
#             # err = self.fullParamLength-self.measure_length() # These two error values seem equivalent
#             conv = abs(err-errPrev)
#             i+=1
#             SSProfile("reinit").toc()
#         # if the spring was given a name, save its input parameters in a filej

#         # this takes way too many lines but its better to see the math written
#         # out:

#         # Find coordinates of frames:
#         # at base of spring
#         self.x0 = PPoly_Eval(0,self.XCoeffs)
#         self.y0 = PPoly_Eval(0,self.YCoeffs)
#         # at end of spring
#         self.xL = PPoly_Eval(self.fullParamLength,self.XCoeffs)
#         self.yL = PPoly_Eval(self.fullParamLength,self.YCoeffs)
#         # at end of spring, interms of frame at base of spring
#         self.momentArmX = self.xL - self.x0
#         self.momentArmY = self.yL - self.y0

#         self.innerRadius= ffradii[0]
#         self.outerRadius= ffradii[-1]

#         self.saveVariables = {'n':     n,     'length':      correctedFullLength,
#                               'radii': radii, 'beta angles': betaAngles,
#                               'control factors': XYFactors}

#         # print(self.XCoeffs)
#         # print(self.YCoeffs)

#         # self.get_neutralSurface(200)
#         # plt.plot(self.undeformedNeutralSurface[:,0], self.undeformedNeutralSurface[:,1])
#         # plt.show()
#         # np.savetxt("springs/")

#     def init(self,
#              n                   = 2,
#              fullParamLength     = 6,
#              radii               = np.array([1.1,2.025,2.025,2.95]),
#              betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
#              XYParamLens           = np.array([0.333,0.667])):       # radii - 2

#         """
#         THIS FUNCTION:
#             Takes in default or custom parameters and
#             defines the quantities and function needed for analysis
#         """

#         ############################ PARAMETER KEY ############################
#         # E                 - youngs modulus                                  #
#         # designStress      - max allowable stress, currently 80% yield       #
#         # n                 - # of flexures                                   #
#         # fullParamLength   - maximum value of internal parameter             #
#         # outPlaneThickness - axial thickness of spring (in)                  #
#         # radii             - radii of control points of netural surface (in) #
#         # betaAngles        - angle of control points of neutral surface (rad)#
#         # IcPts             - curved second moment of area of beam at points  #
#         #                     along beam (in^4)                               #
#         # IcParamLens       - internal parameter values at which Ic control   #
#         #                     points apply                                    #
#         # IcParamLens       - internal parameter values at which x-y          #
#         #                     (radii-angle) control points apply              #
#         # resolution        - number of points in mesh                        #
#         ########################################################################

#         # stick all the arguments in the object
#         self.n = n
#         self.fullParamLength = fullParamLength
#         self.radii = radii
#         self.betaAngles = betaAngles
#         self.XYFactors = XYParamLens

#         # create control points for polynomials

#         # x-y control points
#         self.pts = np.empty((len(radii),2))
#         for i in range(len(radii)):
#             self.pts[i,:] = [self.radii[i]*np.cos(self.betaAngles[i]),self.radii[i]*np.sin(self.betaAngles[i])]

#         # Cretae parameter control points for path based on passed factors
#         self.XYParamLens = np.empty(len(radii))
#         for i in range(len(radii)):
#             if i==0:
#                 self.XYParamLens[i] = 0
#             elif i==len(radii)-1:
#                 self.XYParamLens[i] = self.fullParamLength
#             else:
#                 self.XYParamLens[i] = self.fullParamLength*XYParamLens[i-1]

#         # create a vector of all mutable parameters
#         self.parameterVector = np.concatenate((self.radii, self.betaAngles,
#                                               self.XYFactors))

#         # generate the coefficients for polynomials
#         self.XCoeffs, self.YCoeffs  = self.xy_poly()

#     def xy_poly(self):

#         """
#         THIS IS THE ONLY FUNCTION THAT CHANGES:

#         Instead of being locally radial at either end, we want to be locally
#         tangent to the circle
#         """

#         # if we have n points we need:
#         # n constraints on zeroth derivative
#         # 2 constraints on first derivative (angle at start and end)
#         # 2 constraints on second derivative (straightness at start and end)
#         # n + 4 constraints
#         # n + 3 degree polynomial

#         # initialize matrix sizes
#         nConstraints = len(self.pts)+4
#         Mat = np.empty((nConstraints,nConstraints))
#         YTarg = np.empty((nConstraints,1))
#         XTarg = np.empty((nConstraints,1))
#         # FIRST ASSIGN 0th DERIVATIVE CONSTRAINTS (c(s)=p)
#         for i in range(len(self.XYParamLens)):
#             # target correct x-y value
#             XTarg[i]=self.pts[i,0]
#             YTarg[i]=self.pts[i,1]
#             # at associated s value
#             for j in range(nConstraints):
#                 Mat[i,j]   = self.XYParamLens[i]**(nConstraints-1-j)
#         # NOW ASSIGN FIRST AND SECOND DERIVATIVE CONSTRAINTS
#         for i in range(-1,-5,-1):
#             # target correct index (first or last)
#             index = (i % 2)*(len(self.XYParamLens)-1)
#             print(index)
#             if i < -2:
#                 # with correct first derivative (radial direction)
#                 XTarg1 = self.pts[index,0]/lin.norm(self.pts[index,:])
#                 YTarg1 = self.pts[index,1]/lin.norm(self.pts[index,:])
#                 if not index==0:# rotate 90 degrees ccw
#                     XTarg[i] = -1*YTarg1
#                     YTarg[i] =  1*XTarg1
#                 else:
#                     XTarg[i] = XTarg1
#                     YTarg[i] = YTarg1
#                 # at correct s value
#                 for j in range(nConstraints):
#                     Mat[i,j] = (nConstraints-1-j)* \
#                                 self.XYParamLens[index]**max(nConstraints-2-j,0) ## FIRST TWO ARE FIRST DERIVATIVES
#             else:
#                 # and with zero second derivative at both points
#                 XTarg[i] = 0
#                 YTarg[i] = 0
#                 for j in range(nConstraints):
#                     Mat[i,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)* \
#                                 self.XYParamLens[index]**max(nConstraints-3-j,0) ## LAST TWO ARE SECOND DERIVATIVES
#         # print(XTarg)
#         # print(YTarg)
#         # print(Mat)
#         XCoeffs = lin.solve(Mat, XTarg)
#         YCoeffs = lin.solve(Mat, YTarg)
#         print(XCoeffs)
#         print(YCoeffs)
#         # print(lin.cond(Mat))
#         for i in range(len(self.XYParamLens)):
#             diffX = PPoly_Eval(self.XYParamLens[i],XCoeffs)-self.pts[i,0]
#             diffY = PPoly_Eval(self.XYParamLens[i],YCoeffs)-self.pts[i,1]
#             diff = lin.norm([diffX,diffY])
#             print(diff)
#             assert(np.isclose(diff,0))
#         # TODO: assert(np.isclose(diff,0))
#         # print(XCoeffs, YCoeffs)
#         return XCoeffs, YCoeffs

#     def measure_length(self):

#         """
#         THIS FUNCTION:
#             measures the length of an initialized spring
#         """

#         self.measureResolution = 200
#         ximesh            = np.linspace(0,self.fullParamLength,self.measureResolution+1)

#         # calcualte derivative of x and y with respect to xi (internal parameter)
#         self.dxdxi = PPoly_Eval(ximesh, self.XCoeffs, deriv = 1)
#         self.dydxi = PPoly_Eval(ximesh, self.YCoeffs, deriv = 1)
#         integrand = np.sqrt(self.dxdxi**2+self.dydxi**2)
#         # create a mesh in s according by integrating along the parameter to find
#         # arc length
#         self.smesh = np.zeros(len(ximesh))
#         self.smesh[1:len(self.smesh)] = scipy.integrate.cumulative_trapezoid(
#                                          integrand, ximesh, ximesh[1]-ximesh[0])
#         # print(self.smesh)
#         assert(self.smesh[0]==0)
#         # return the overall length of the spring
#         return self.smesh[-1]

# #### STANDARD INTERFACE ####

#     def get_crscRef(self, crscRef):
#         self.crsc = crscRef

#     def get_xy_n(self, coord, dim):
#         if dim=='x':
#             out = PPoly_Eval(coord, self.XCoeffs)
#         if dim=='y':
#             out = PPoly_Eval(coord, self.YCoeffs)
#         return out

#     def get_dxdy_n(self, coord, dim):
#         if dim=='x':
#             out = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         if dim=='y':
#             out = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         return out

#     def get_dxi_n(self, coord):
#         dxdxi = self.get_dxdy_n(coord, 'x')
#         dydxi = self.get_dxdy_n(coord, 'y')
#         dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
#         return dxids

#     def get_rn(self, coord):
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         # dAlphadS = (d2yds2/dxds-d2xds2*dyds/dxds)/(1+dyds**2/dxds)
#         if hasattr(coord, "__len__"):
#             rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#             for i in range(len(rn)):
#                 if not np.isfinite(rn[i]):
#                     rn[i] = float('inf')
#                 if abs((d2yds2[i]/dxds[i]-d2xds2[i]*dyds[i]/dxds[i]**2)) <= 10**-13:
#                     rn[i] = float('inf')*np.sign(rn[i-1])
#         else:
#             if abs(d2yds2/dxds-d2xds2*dyds/dxds**2)<=10**-13 or not np.isfinite(abs(d2yds2/dxds-d2xds2*dyds/dxds**2)):
#                 rn = float('inf')
#             else:
#                 rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#         return rn

#     def get_drn(self, coord):
#         d3yds3 = PPoly_Eval(coord, self.YCoeffs, deriv=3)
#         d3xds3 = PPoly_Eval(coord, self.XCoeffs, deriv=3)
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)

#         # denominator = (d2yds2*dxds-d2xds2*dyds)**2
#         # numerator   = ( dxds**3*d3yds3 - dyds**3*d3xds3 +
#         #                 2*dyds*d2xds2**2*dxds - 2*dyds*d2yds2**2*dxds - dxds**2*d3xds3*dyds -
#         #                 2*dxds**2*d2yds2*d2xds2 + 2*dyds**2*d2yds2*d2xds2 +
#         #                 dyds**2*d3yds3*dxds )

#         # drnds = -numerator/denominator

#         num = (1+dyds**2/dxds**2)
#         den = (d2yds2/dxds-d2xds2*dyds/dxds**2)

#         dnumds = 2*dyds*(dxds*d2yds2-d2xds2*dyds)/dxds**3
#         ddends = (dxds**2*d3yds3-2*dxds*d2xds2*d2yds2+dyds*(2*d2xds2**2-dxds*d3xds3))/dxds**3

#         drnds = (dnumds*den-num*ddends)/den**2

#         return drnds

#     def get_alpha(self, coord):
#         if hasattr(coord, "__len__"):
#             alphaList = np.arctan2(PPoly_Eval(coord, self.YCoeffs, deriv=1),
#                                        PPoly_Eval(coord, self.XCoeffs, deriv=1))
#             for i in range(len(alphaList)):
#                 if alphaList[i]<0:
#                     alphaList[i]=alphaList[i]+2*np.pi
#             return alphaList
#         else:
#             alpha = np.arctan2(PPoly_Eval(coord, self.YCoeffs, deriv=1),
#                                        PPoly_Eval(coord, self.XCoeffs, deriv=1))
#             if alpha<0:
#                 alpha=alpha+2*np.pi
#             return alpha

#     def get_dalpha(self, coord):
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         dads = ((d2yds2/dxds-d2xds2*dyds/dxds**2)/(1+dyds**2/dxds**2))
#         return dads

#     def get_d2alpha(self, coord):
#         d2ads2 = -self.get_dalpha(coord)**2*self.get_drn(coord)
#         if np.isnan(d2ads2):
#             d2ads2 = 0
#         return d2ads2

#     def get_neutralSurface(self, resolution):

#         ximesh = np.linspace(0,self.fullParamLength,resolution+1)
#         self.undeformedNeutralSurface = np.hstack(
#                                   (np.atleast_2d(self.get_xy_n(ximesh, 'x')).T,
#                                    np.atleast_2d(self.get_xy_n(ximesh, 'y')).T))

#         return self.undeformedNeutralSurface

#     def get_centroidalSurface(self, resolution):
#         ximesh = np.linspace(0,self.fullParamLength,resolution+1)
#         Ic = self.crsc.get_Ic(ximesh)
#         t = self.crsc.t
#         h = self.crsc.get_Thk(ximesh, hasPrev=False)
#         rn = self.get_rn(ximesh)
#         ecc = Ic/(t*h*rn)
#         alpha = self.get_alpha(ximesh)

#         neutralSurface = self.get_neutralSurface(resolution)
#         self.undeformedCentroidalSurface = neutralSurface + \
#                        np.hstack((np.atleast_2d(ecc*np.sin(alpha)).T,
#                                   np.atleast_2d(ecc*np.cos(alpha)).T))

#         return self.undeformedCentroidalSurface

# class Minimal_Polynomial_Definition4:
#     def __init__(self,
#                  # whole bunch of default values:
#                  n                   = 2,
#                  fullParamLength     = 6,
#                  radii               = np.array([1.1,2.025,2.025,2.95]),
#                  ffradii             = np.array([1, 2.5]),
#                  alphaAngles         = np.array([45,90])*deg2rad,
#                  betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
#                  XYFactors         = np.array([0.333,0.667])
#                  ):

#         self.alphaAngles = alphaAngles

#         self.init(n, fullParamLength, radii,
#                   betaAngles, XYFactors)
#         # iterate spring initialization until xi mesh approximates s mesh as well as possible
#         conv=1
#         err = 1
#         i = 0
#         # iterate until you converge on a value
#         while conv>10e-6:
#             SSProfile("reinit").tic()
#             errPrev = err
#             # measure real arc length each time
#             correctedFullLength = self.measure_length()
#             # initialize the spring with the remeasured legnth
#             self.init(n, correctedFullLength, ### WE ONLY CHANGE THIS
#                       radii, betaAngles, XYFactors)
#             ximesh = np.linspace(0,correctedFullLength,self.measureResolution+1)
#             dxids  = self.get_dxi_n(ximesh)
#             err = lin.norm(dxids-np.ones(self.measureResolution+1))
#             # err = self.fullParamLength-self.measure_length() # These two error values seem equivalent
#             conv = abs(err-errPrev)
#             i+=1
#             SSProfile("reinit").toc()
#         # if the spring was given a name, save its input parameters in a filej

#         # do this so that symbolic algebra only happens once at beginning
#         self.prepare_rn_expr()

#         # this takes way too many lines but its better to see the math written
#         # out:

#         # Find coordinates of frames:
#         # at base of spring
#         self.x0 = PPoly_Eval(0,self.XCoeffs)
#         self.y0 = PPoly_Eval(0,self.YCoeffs)
#         # at end of spring
#         self.xL = PPoly_Eval(self.fullParamLength,self.XCoeffs)
#         self.yL = PPoly_Eval(self.fullParamLength,self.YCoeffs)
#         # at end of spring, interms of frame at base of spring
#         self.momentArmX = self.xL - self.x0
#         self.momentArmY = self.yL - self.y0

#         self.innerRadius= ffradii[0]
#         self.outerRadius= ffradii[-1]

#         self.saveVariables = {'n':     n,     'length':      correctedFullLength,
#                               'radii': radii, 'beta angles': betaAngles,
#                               'control factors': XYFactors}

#         # print(self.XCoeffs)
#         # print(self.YCoeffs)

#         # self.get_neutralSurface(200)
#         # plt.plot(self.undeformedNeutralSurface[:,0], self.undeformedNeutralSurface[:,1])
#         # plt.show()
#         # np.savetxt("springs/")

#     def init(self,
#              n                   = 2,
#              fullParamLength     = 6,
#              radii               = np.array([1.1,2.025,2.025,2.95]),
#              betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
#              XYParamLens           = np.array([0.333,0.667])):       # radii - 2

#         """
#         THIS FUNCTION:
#             Takes in default or custom parameters and
#             defines the quantities and function needed for analysis
#         """

#         ############################ PARAMETER KEY ############################
#         # E                 - youngs modulus                                  #
#         # designStress      - max allowable stress, currently 80% yield       #
#         # n                 - # of flexures                                   #
#         # fullParamLength   - maximum value of internal parameter             #
#         # outPlaneThickness - axial thickness of spring (in)                  #
#         # radii             - radii of control points of netural surface (in) #
#         # betaAngles        - angle of control points of neutral surface (rad)#
#         # IcPts             - curved second moment of area of beam at points  #
#         #                     along beam (in^4)                               #
#         # IcParamLens       - internal parameter values at which Ic control   #
#         #                     points apply                                    #
#         # IcParamLens       - internal parameter values at which x-y          #
#         #                     (radii-angle) control points apply              #
#         # resolution        - number of points in mesh                        #
#         ########################################################################

#         # stick all the arguments in the object
#         self.n = n
#         self.fullParamLength = fullParamLength
#         self.radii = radii
#         self.betaAngles = betaAngles
#         self.XYFactors = XYParamLens

#         # create control points for polynomials

#         # x-y control points
#         self.pts = np.empty((len(radii),2))
#         for i in range(len(radii)):
#             self.pts[i,:] = [self.radii[i]*np.cos(self.betaAngles[i]),self.radii[i]*np.sin(self.betaAngles[i])]

#         # Cretae parameter control points for path based on passed factors
#         self.XYParamLens = np.empty(len(radii))
#         for i in range(len(radii)):
#             if i==0:
#                 self.XYParamLens[i] = 0
#             elif i==len(radii)-1:
#                 self.XYParamLens[i] = self.fullParamLength
#             else:
#                 self.XYParamLens[i] = self.fullParamLength*XYParamLens[i-1]

#         # generate the coefficients for polynomials
#         self.XCoeffs, self.YCoeffs  = self.xy_poly()

#     def xy_poly(self):

#         """
#         THIS IS THE ONLY FUNCTION THAT CHANGES:

#         Instead of being locally radial at either end, we want to be locally
#         tangent to the circle
#         """

#         # if we have n points we need:
#         # n constraints on zeroth derivative
#         # 2 constraints on first derivative (angle at start and end)
#         # 2 constraints on second derivative (straightness at start and end)
#         # n + 4 constraints
#         # n + 3 degree polynomial

#         # initialize matrix sizes
#         nConstraints = len(self.pts)+4
#         Mat = np.empty((nConstraints,nConstraints))
#         YTarg = np.empty((nConstraints,1))
#         XTarg = np.empty((nConstraints,1))
#         # FIRST ASSIGN 0th DERIVATIVE CONSTRAINTS (c(s)=p)
#         for i in range(len(self.XYParamLens)):
#             # target correct x-y value
#             XTarg[i]=self.pts[i,0]
#             YTarg[i]=self.pts[i,1]
#             # at associated s value
#             for j in range(nConstraints):
#                 Mat[i,j]   = self.XYParamLens[i]**(nConstraints-1-j)
#         # NOW ASSIGN FIRST AND SECOND DERIVATIVE CONSTRAINTS
#         for i in range(-1,-5,-1):
#             # target correct index (first or last)
#             index = (i % 2)*(len(self.XYParamLens)-1)
#             # print(index)
#             if i < -2:
#                 # with correct first derivative (radial direction)
#                 XTarg1 = self.pts[index,0]/lin.norm(self.pts[index,:])
#                 YTarg1 = self.pts[index,1]/lin.norm(self.pts[index,:])
#                 TargVect = np.array([[XTarg1],[YTarg1]])
#                 if index == 0:
#                     Mat2 = np.array([[np.cos(self.alphaAngles[index]),-np.sin(self.alphaAngles[index])],
#                                     [np.sin(self.alphaAngles[index]), np.cos(self.alphaAngles[index])]])
#                     TargVect = Mat2.dot(TargVect)
#                 else:
#                     Mat2 = np.array([[np.cos(self.alphaAngles[-1]),-np.sin(self.alphaAngles[-1])],
#                                     [np.sin(self.alphaAngles[-1]), np.cos(self.alphaAngles[-1])]])
#                     TargVect = Mat2.dot(TargVect)
#                 # print(TargVect)
#                 XTarg[i] = TargVect[0]
#                 YTarg[i] = TargVect[1]
#                 # at correct s value
#                 for j in range(nConstraints):
#                     Mat[i,j] = (nConstraints-1-j)* \
#                                 self.XYParamLens[index]**max(nConstraints-2-j,0) ## FIRST TWO ARE FIRST DERIVATIVES
#             else:
#                 # and with zero second derivative at both points
#                 XTarg[i] = 0
#                 YTarg[i] = 0
#                 for j in range(nConstraints):
#                     Mat[i,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)* \
#                                 self.XYParamLens[index]**max(nConstraints-3-j,0) ## LAST TWO ARE SECOND DERIVATIVES
#         # print(XTarg)
#         # print(YTarg)
#         # print(Mat)
#         XCoeffs = lin.solve(Mat, XTarg)
#         YCoeffs = lin.solve(Mat, YTarg)
#         # print(XCoeffs)
#         # print(YCoeffs)
#         # print(lin.cond(Mat))
#         for i in range(len(self.XYParamLens)):
#             diffX = PPoly_Eval(self.XYParamLens[i],XCoeffs)-self.pts[i,0]
#             diffY = PPoly_Eval(self.XYParamLens[i],YCoeffs)-self.pts[i,1]
#             diff = lin.norm([diffX,diffY])
#             # print(diff)
#             assert(np.isclose(diff,0))
#         # TODO: assert(np.isclose(diff,0))
#         # print(XCoeffs, YCoeffs)
#         return XCoeffs, YCoeffs

#     def measure_length(self):

#         """
#         THIS FUNCTION:
#             measures the length of an initialized spring
#         """

#         self.measureResolution = 200
#         ximesh            = np.linspace(0,self.fullParamLength,self.measureResolution+1)

#         # calcualte derivative of x and y with respect to xi (internal parameter)
#         self.dxdxi = PPoly_Eval(ximesh, self.XCoeffs, deriv = 1)
#         self.dydxi = PPoly_Eval(ximesh, self.YCoeffs, deriv = 1)
#         integrand = np.sqrt(self.dxdxi**2+self.dydxi**2)
#         # create a mesh in s according by integrating along the parameter to find
#         # arc length
#         self.smesh = np.zeros(len(ximesh))
#         self.smesh[1:len(self.smesh)] = scipy.integrate.cumulative_trapezoid(
#                                          integrand, ximesh, ximesh[1]-ximesh[0])
#         # print(self.smesh)
#         assert(self.smesh[0]==0)
#         # return the overall length of the spring
#         return self.smesh[-1]

#     def prepare_rn_expr(self):
#         XCoeffs = self.XCoeffs.flatten()
#         YCoeffs = self.YCoeffs.flatten()

#         s = sp.symbols('s')
#         yPoly = sp.Poly(YCoeffs, s).as_expr()
#         xPoly = sp.Poly(XCoeffs, s).as_expr()

#         dyds   = sp.diff(yPoly, s)
#         dxds   = sp.diff(xPoly, s)
#         d2yds2 = sp.diff(yPoly, s, 2)
#         d2xds2 = sp.diff(xPoly, s, 2)

#         rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#         # print(rn)
#         # print("---------------------------------")
#         self.drnds_expr = sp.diff(rn, s)
#         # print(self.drnds_expr)
#         # print(self.fullParamLength)
#         self.drnds = sp.lambdify(s, self.drnds_expr, "numpy")

# #### STANDARD INTERFACE ####

#     def get_crscRef(self, crscRef):
#         self.crsc = crscRef

#     def get_xy_n(self, coord, dim):
#         if dim=='x':
#             out = PPoly_Eval(coord, self.XCoeffs)
#         if dim=='y':
#             out = PPoly_Eval(coord, self.YCoeffs)
#         return out

#     def get_dxdy_n(self, coord, dim):
#         if dim=='x':
#             out = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         if dim=='y':
#             out = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         return out

#     def get_dxi_n(self, coord):
#         dxdxi = self.get_dxdy_n(coord, 'x')
#         dydxi = self.get_dxdy_n(coord, 'y')
#         dxids = 1/np.sqrt(dxdxi**2+dydxi**2)
#         return dxids

#     def get_rn(self, coord):
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         # dAlphadS = (d2yds2/dxds-d2xds2*dyds/dxds)/(1+dyds**2/dxds)
#         if hasattr(coord, "__len__"):
#             rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#             for i in range(len(rn)):
#                 if not np.isfinite(rn[i]):
#                     rn[i] = float('inf')
#                 if abs((d2yds2[i]/dxds[i]-d2xds2[i]*dyds[i]/dxds[i]**2)) <= 10**-13:
#                     rn[i] = float('inf')*np.sign(rn[i-1])
#         else:
#             if abs(d2yds2/dxds-d2xds2*dyds/dxds**2)<=10**-13 or not np.isfinite(abs(d2yds2/dxds-d2xds2*dyds/dxds**2)):
#                 rn = float('inf')
#             else:
#                 rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
#         return rn

#     def get_drn(self, coord):
#         coord = np.atleast_1d(coord)
#         out = self.drnds(coord)
#         if len(out)==1:
#             out = out[0]
#         if np.isnan(out):
#             out = 0
#         return out

#     def get_alpha(self, coord):
#         if hasattr(coord, "__len__"):
#             alphaList = np.arctan2(PPoly_Eval(coord, self.YCoeffs, deriv=1),
#                                        PPoly_Eval(coord, self.XCoeffs, deriv=1))
#             for i in range(len(alphaList)):
#                 if alphaList[i]<0:
#                     alphaList[i]=alphaList[i]+2*np.pi
#             return alphaList
#         else:
#             alpha = np.arctan2(PPoly_Eval(coord, self.YCoeffs, deriv=1),
#                                        PPoly_Eval(coord, self.XCoeffs, deriv=1))
#             if alpha<0:
#                 alpha=alpha+2*np.pi
#             return alpha

#     def get_dalpha(self, coord):
#         d2yds2 = PPoly_Eval(coord, self.YCoeffs, deriv=2)
#         d2xds2 = PPoly_Eval(coord, self.XCoeffs, deriv=2)
#         dyds   = PPoly_Eval(coord, self.YCoeffs, deriv=1)
#         dxds   = PPoly_Eval(coord, self.XCoeffs, deriv=1)
#         dads = ((d2yds2/dxds-d2xds2*dyds/dxds**2)/(1+dyds**2/dxds**2))
#         return dads

#     def get_d2alpha(self, coord):
#         d2ads2 = -self.get_dalpha(coord)**2*self.get_drn(coord)
#         if np.isnan(d2ads2):
#             d2ads2 = 0
#         return d2ads2

#     def get_neutralSurface(self, resolution):

#         ximesh = np.linspace(0,self.fullParamLength,resolution+1)
#         self.undeformedNeutralSurface = np.hstack(
#                                   (np.atleast_2d(self.get_xy_n(ximesh, 'x')).T,
#                                    np.atleast_2d(self.get_xy_n(ximesh, 'y')).T))

#         return self.undeformedNeutralSurface

#     def get_centroidalSurface(self, resolution):
#         ximesh = np.linspace(0,self.fullParamLength,resolution+1)
#         Ic = self.crsc.get_Ic(ximesh)
#         t = self.crsc.t
#         h = self.crsc.get_Thk(ximesh, hasPrev=False)
#         rn = self.get_rn(ximesh)
#         ecc = Ic/(t*h*rn)
#         alpha = self.get_alpha(ximesh)

#         neutralSurface = self.get_neutralSurface(resolution)
#         self.undeformedCentroidalSurface = neutralSurface + \
#                        np.hstack((np.atleast_2d(ecc*np.sin(alpha)).T,
#                                   np.atleast_2d(ecc*np.cos(alpha)).T))

#         return self.undeformedCentroidalSurface





