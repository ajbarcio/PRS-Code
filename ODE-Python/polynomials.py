
import numpy as np
import numpy.linalg as lin
import sympy as sp

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
