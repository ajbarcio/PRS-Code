
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt

import math
from StatProfiler import SSProfile
import matplotlib.collections as mcoll
import sympy as sp

import subprocess as sb
import win32com.client
import pythoncom

deg2rad = np.pi/180

## Automate SW Verification?? (Currently not working)

def make_SW_part(Spring):
    startSW()
    sw = connectToSW()
    newFile(sw, "newPart.SLDPRT")

## Turn design variables into polynomials/meshes/surfaces

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

## Calculate spring characteristics based off of parameter-driven polynomials

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

## General function to differentiate an array of values in a fixed mesh

def numerical_fixed_mesh_diff(ymesh, xmesh):
    dydx = np.zeros(len(xmesh))
    step = xmesh[1]-xmesh[0]
    for i in range(len(xmesh)):
        if i==0:
            dydx[i] = (ymesh[i+1]-ymesh[i])/step
        elif i==len(xmesh)-1:
            dydx[i] = (ymesh[i]-ymesh[i-1])/step
        else:
            dydx[i] = (ymesh[i+1]-ymesh[i-1])/(2*step)
    return dydx

## Low-Level functions used to evaluate ODE's: Fixed mesh RK4

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

## Code from the internet used to graph stress

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

## Low level solidworks automation codes

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
