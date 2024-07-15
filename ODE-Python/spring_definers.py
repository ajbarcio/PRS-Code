#! Built in python dependencies
import numpy as np
import numpy.linalg as lin
import sympy as sp
#! Borrowed dependencies
from StatProfiler import SSProfile # function timing library from dr thomas
#! developed (included) dependencies
import materials

"""
Create PRS defined by a variety of parameterization strategies.
ALL UNITS ARE IN CUSTOMARY, IPS

See [fuck you] for [random bullshit]

Classes:

    Polynomial_Spring

Functions:

    none yet because fuck you

Variables:

    deg2rad
"""

deg2rad = np.pi/180

spring1 = Spring(al, 2, 200)
spring2 = Spring(steel, 3, 300)
# al, 2, 200
currentSpring = Silly_Spring()
# steel, 3, 300
otherSpring  =  Silly_Spring2()
....
Ic = currentSpring.get_Ic(xi)


class Spring:
    def __init__(self, material, outPlaneThickness, resolution):

        material = 
        outPlaneThickness = 
        resolution = 

class Polynomial_Spring:

    """
    Create PRS defined by arbitrary polynomials for path and geometry

    An arbitrary number of paths can be set
    Thickness (geometry) is indirectly defined through curved second moment area
        (Ic) and consists of an arbitraryt number of quadratic functions with
        first-order continutity

        Parameters:
            Too damn many to list here. Like just a whole bunch
    """

    def get_Ic(self, xi):
        Ic = self.PPoly_Eval(xi, self.IcCoeffs)
        return Ic

    def __init__(self, material, name=None,
                 # whole bunch of default values:
                 n                   = 2,
                 fullParamLength     = 6,
                 outPlaneThickness   = 0.375,
                 radii               = np.array([1.1,2.025,2.025,2.95]),
                 betaAngles          = np.array([0,50,100,150])*deg2rad, # FIRST VALUE IS ALWAYS 0 !!!!, same length as radii
                 IcPts               = np.array([0.008, 0.001, 0.008]),
                 IcParamLens         = np.array([0.5]),              # IcPts - 2
                 XYParamLens         = np.array([0.333,0.667]),      # radii - 2
                 resolution          = 200):
        """
        Iteratively initializes a spring such that the difference between the
        true arc length and the internal parameter is minimized.

            Parameters:
                material:          a material object from the materials module,
                                   containing relevant material properties
                name:              a name, used for saving relevant spring
                                   properties in an output file
                n:                 numer of paths (legs) for this spring
                fullParamLength:   overall length of spring in units of the
                                   internal parameter. This is somewhat
                                   arbitrary to specify due to the interation,
                                   but results improve if it is close to the
                                   overall true arc length
                outPlaneThickness: thickness of spring out of plane
                radii:             list of distance of each control point from
                                   center of spring
                betaAngles:        list of angles form horizontal of each
                                   control point
                IcPts:             list of Ic (curved I) to be controlled
                IcParamLengths:    internal parameter values at which IcPts are
                                   enforced
                XYParamLengths:    internal parameter values at which radii/beta
                                   angles are enforced
                resolution:        size of fixed mesh
        """

        # improperly initialize the spring (unable to start with good guess for arc length)
        E = material.E
        designStress = material.yieldStress*0.8
        self.reinit(E, designStress, n, fullParamLength, outPlaneThickness, radii,
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
            self.reinit(E, designStress, n, correctedFullLength, ### WE ONLY CHANGE THIS
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
            saveVariables = np.hstack([n,correctedFullLength,outPlaneThickness,
                                       radii,betaAngles,IcPts,IcParamLens,
                                       XYParamLens])
            filepath = "springs\\"+name
            np.savetxt(filepath, saveVariables, "%f", ",")

        # Save the names of the functions used for interface with deformer
        
        # np.savetxt("springs/")

    def reinit(self,
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
        One step in the iteration to bring fullParamLength nearer to true arc
        length

            Parameters:
                Same as __init__
        """

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
        self.x0 = self.PPoly_Eval(0,self.XCoeffs)
        self.y0 = self.PPoly_Eval(0,self.YCoeffs)
        # at end of spring
        self.xL = self.PPoly_Eval(self.fullParamLength,self.XCoeffs)
        self.yL = self.PPoly_Eval(self.fullParamLength,self.YCoeffs)
        # at end of spring, interms of frame at base of spring
        self.momentArmX = self.xL - self.x0
        self.momentArmY = self.yL - self.y0

    def generate_Ic_poly(self):

        """
        Generates a piecwise function of second-order polynomials which defines
        the Ic (curved I) along one path of the spring

            Parameters:
                none, uses IcPts and IcParamLens from initialization

            Returns:
                coeffs (double array, 2d): all coefficients of all polynomials
                                           which satisfy the control points
                ctrlX: (double array, 1d): starting points of each sub-domain
                                           where each polynomial is in effect
        """

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
        # for i in range(numPolys):
        #     Targ = np.array([[IcPts[i]],[],[]])

    def generate_xy_poly(self):

        """
        Generates a polynomial of arbitrary degree to define the path of one
        leg of the spring. Path is locally straight and radially pointed at
        either end.

            Parameters:
                none, uses pts and XYParamLens from initialization

            Returns:
                XCoeffs (double array): List of coefficients (in descending
                                        order) for polynomial defining x-coord-
                                        -inates for netural surface of spring
                YCoeffs (double array): List of coefficients (in descending
                                        order) for polynomial defining y-coord-
                                        -inates for netural surface of spring
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
                XTarg[i] = self.pts[index,0]/lin.norm(self.pts[index,:])
                YTarg[i] = self.pts[index,1]/lin.norm(self.pts[index,:])
                # at correct s value
                for j in range(nConstraints):
                    Mat[i,j] = (nConstraints-1-j)*self.XYParamLens[index]**max(nConstraints-2-j,0) ## FIRST TWO ARE FIRST DERIVATIVES
            else:
                # and with zero second derivative at both points
                XTarg[i] = 0
                YTarg[i] = 0
                for j in range(nConstraints):
                    Mat[i,j] = (nConstraints-1-j)*max(nConstraints-2-j,0)*self.XYParamLens[index]**max(nConstraints-3-j,0) ## LAST TWO ARE SECOND DERIVATIVES
        # print(XTarg)
        # print(YTarg)
        # print(Mat)
        XCoeffs = lin.solve(Mat, XTarg)
        YCoeffs = lin.solve(Mat, YTarg)
        # print(lin.cond(Mat))
        for i in range(len(self.XYParamLens)):
            diffX = self.PPoly_Eval(self.XYParamLens[i],XCoeffs)-self.pts[i,0]
            diffY = self.PPoly_Eval(self.XYParamLens[i],YCoeffs)-self.pts[i,1]
            diff = lin.norm([diffX,diffY])
            assert(np.isclose(diff,0))
        # TODO: assert(np.isclose(diff,0))
        # print(XCoeffs, YCoeffs)
        return XCoeffs, YCoeffs

    def PPoly_Eval(self, x, coeffs, deriv=0, ranges=0):

        """
        Evaluates piecewise polynomial of arbitrary order and number of domains,
        and any derivatives thereof

            Parameters:
                x (double, scalar or array-like):      Independent variable
                coeffs (double, up to 2d array-like):  Coefficients, in descend-
                                                       -ing order, of polynomial.
                                                       Each row applies over
                                                       each sub-domain
                deriv (int):                           Order of derivative to be
                                                       taken
                ranges (double, up to 1d array-like)): List of the start of each
                                                       subdomain
        """

        # Make ranges into a vector:
        ranges = np.atleast_1d(ranges)
        # Make x into a vector:
        x      = np.atleast_1d(x)
        # initialize output of polynomial function (y)
        y      = np.empty(len(x))
        i = 0
        # for each independent variable:
        for value in x:
            # initialize coefficient matrix
            U = np.empty(coeffs.shape[1])
            # for each coefficient
            for j in range(coeffs.shape[1]):
                # determine coefficient based on order of derivative being taken
                preCoeff = 1
                for k in range(deriv):
                    preCoeff = preCoeff*max(coeffs.shape[1]-j-(k+1),0)
                U[j] = preCoeff*value**max((coeffs.shape[1]-j-1-deriv),0)
            # for each subdomain
            for l in range(len(ranges)-1):
                # determine if x value lies in this subdomain
                if value >= ranges[l] and value < ranges[l+1]:
                    index = l
                    break
                index = len(ranges)-2
            # evaluate dependent variable
            y[i] = U.dot(coeffs[index,:])
            i+=1
        return y

    def alpha_xy(self, xi):

        """
        Gets alpha value at a given internal parameter (xi) value

            Parameters:
                xi (double): internal parameter (can be array)

            Returns:
                alphaList (double): alpha angle in radians (can be array)

        """

        # make xi into vector
        xi = np.atleast_1d(xi)
        # get alpha value
        alphaList = np.arctan2(self.PPoly_Eval(xi, self.YCoeffs, deriv=1),
                               self.PPoly_Eval(xi, self.XCoeffs, deriv=1))
        # no negative values
        for i in range(len(alphaList)):
            if alphaList[i]<0:
                alphaList[i]=alphaList[i]+2*np.pi
        return alphaList

    def d_alpha(self, xi):

        """
        Get derivative of alpha with respect to independent variable (nominally
        the internal parameter xi)

            Parameters:
                xi (double): independent variable (can be array)

            Returns:
                dadxi (double): derivative of alpha with respect to IV (can be
                                array)
        """

        #get derivatives of neutral surface
        d2ydxi2 = self.PPoly_Eval(xi, self.YCoeffs, deriv=2)
        d2xdxi2 = self.PPoly_Eval(xi, self.XCoeffs, deriv=2)
        dydxi   = self.PPoly_Eval(xi, self.YCoeffs, deriv=1)
        dxdxi   = self.PPoly_Eval(xi, self.XCoeffs, deriv=1)
        # analytical solution for alpha derivative
        dadxi = ((d2ydxi2/dxdxi-d2xdxi2*dydxi/dxdxi**2)/(1+dydxi**2/dxdxi**2))
        return dadxi

    def d_2_alpha(self, xi):

        """
        Get second derivative of alpha with respect to independent variable
        (nominally the internal parameter xi)

            Parameters:
                xi (double): independent variable (can be array)

            Returns:
                d2adxi2 (double): derivative of alpha with respect to IV (can be
                                  array)
        """

        d2adxi2 = -self.d_alpha(xi, self.XCoeffs, self.YCoeffs)**2* \
                   self.d_rn(xi, self.XCoeffs, self.YCoeffs)
        if np.isnan(d2adxi2):
            d2adxi2 = 0
        return d2adxi2

    def r_n(self, xi):

        """
        Get the value of the netural radius of curvature at some (perhaps list
        of) independent variable(s) (nominally the internal parameter xi)

            Parameters:
                xi (double): independent variable (can be array)

            Returns:
                rn (double): radius of curvature of neutral surface (can be
                             array)
        """

        d2ydxi2 = self.PPoly_Eval(xi, self.YCoeffs, deriv=2)
        d2xdxi2 = self.PPoly_Eval(xi, self.XCoeffs, deriv=2)
        dydxi   = self.PPoly_Eval(xi, self.YCoeffs, deriv=1)
        dxdxi   = self.PPoly_Eval(xi, self.XCoeffs, deriv=1)
        # make xi into an array
        xi = np.atleast_1d(xi)
        rn = ((1+dydxi**2/dxdxi**2)/(d2ydxi2/dxdxi-d2xdxi2*dydxi/dxdxi**2))
        for i in range(len(rn)):
            if not np.isfinite(rn[i]):
                rn[i] = float('inf')
            if abs((d2ydxi2[i]/dxdxi[i]-d2xdxi2[i]*dydxi[i]/dxdxi[i]**2)) <= 10**-13:
                rn[i] = float('inf')*np.sign(rn[i-1])

        return rn

    def d_rn(self, xi):

        """
        Get the derivative of the netural radius of curvature wrt some (perhaps
        list of) independent variable(s) (nominally the internal parameter xi)

            Parameters:
                xi (double): independent variable (can be array)

            Returns:
                drndxi (double): radius of curvature of neutral surface (can be
                             array)
        """

        d3ydxi3 = self.PPoly_Eval(xi, self.YCoeffs, deriv=3)
        d3xdxi3 = self.PPoly_Eval(xi, self.XCoeffs, deriv=3)
        d2ydxi2 = self.PPoly_Eval(xi, self.YCoeffs, deriv=2)
        d2xdxi2 = self.PPoly_Eval(xi, self.XCoeffs, deriv=2)
        dydxi   = self.PPoly_Eval(xi, self.YCoeffs, deriv=1)
        dxdxi   = self.PPoly_Eval(xi, self.YCoeffs, deriv=1)

        denominator = (d2ydxi2*dxdxi-d2xdxi2*dydxi)**2
        numerator   = ( dxdxi**3*d3ydxi3 - dydxi**3*d3xdxi3 +
                        2*dydxi*d2xdxi2**2*dxdxi - 2*dydxi*d2ydxi2**2*dxdxi - dxdxi**2*d3xdxi3*dydxi -
                        2*dxdxi**2*d2ydxi2*d2xdxi2 + 2*dydxi**2*d2ydxi2*d2xdxi2 +
                        dydxi**2*d3ydxi3*dxdxi )

        drndxi = -numerator/denominator

        return drndxi

    def geometry_coeffs(self):

        """
        Produces coefficients for the polynomials that define the
        geometry of the spring.

            Parameters:
                none

            Returns:
                geometeryDef (list): coefficients of polynomials of:
                                        neutral path (X, YCoeffs)
                                        Ic profile   (IcCoeffs)

        """
        # sticks the coefficients in the object
        self.IcCoeffs, self.domains = self.generate_Ic_poly()
        self.XCoeffs, self.YCoeffs  = self.generate_xy_poly()
        # compiles all the coeffs into one variable
        self.geometryDef = [self.XCoeffs, self.YCoeffs, self.IcCoeffs]
        return self.geometryDef

class Silly_Spring(Spring):

    def __init__(self, material, outPlaneThickness, resolution, *geometryArgs):
        super().__init__(material, outPlaneThickness, resolution)

    def get_Ic(self, xi):
        # do some stuff
        #call a whole bunch of functions
        return Ic

    def a_whole_bunch_of_Ic_Methods():
        pass

class Silly_Spring2(Spring):

    def get_Ic(self, xi):
        # do some other stuff
        return Ic

class Silly_Spring3(Spring):

    def get_Ic(self, xi):
        # do some other other stuff
        return Ic