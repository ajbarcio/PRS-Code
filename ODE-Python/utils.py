import numpy as np
import matplotlib.collections as mcoll
from numpy import linalg as lin
import matplotlib.pyplot as plt

## Function to find intersection between two circles closer to origin

def circle_intersection(radius1, radius2, point1, point2):
    distance = lin.norm(np.subtract(point1, point2))
    if distance > radius1+radius2:
        raise ValueError("Impossible Circles")
    else:    
        # a is the distance from center 1 to the line connecting the intersection points
        a = (radius1**2 - radius2**2 + distance**2) / (2 * distance)
        # h is the height from the line to the intersection points
        h = np.sqrt(radius1**2 - a**2)
        # Find the midpoint between the two centers
        x0 = point1[0] + a * (point2[0] - point1[0]) / distance
        y0 = point1[1] + a * (point2[1] - point1[1]) / distance
        intersection1 = (x0 + h * (point2[1] - point1[1]) / distance, y0 - h * (point2[0] - point1[0]) / distance)
        intersection2 = (x0 - h * (point2[1] - point1[1]) / distance, y0 + h * (point2[0] - point1[0]) / distance)
        if lin.norm(intersection1)<lin.norm(intersection2):
            center=intersection1
        else:
            center=intersection2
    return center

## Function to evaluate piecewise polynomials

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
                stepRes = rk4_step(fun, xmesh[i-1], y0, step, *args)
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
                stepRes = rk4_step(fun, xmesh[i-1], y0, step, *args)
                res[i] = stepRes
                y0 = stepRes
    # print("gonna return")
    return res

def rk4_step(fun, x0, y0, dx, *args): # (fun, alphaCoeffs, cICoeffs, x0, y0, du, Fx, Fy)
    #
    #  Get four sample values of the derivative.
    #
    f1 = fun ( x0,            y0, *args)
    f2 = fun ( x0 + dx / 2.0, y0 + dx * f1 / 2.0, *args)
    f3 = fun ( x0 + dx / 2.0, y0 + dx * f2 / 2.0, *args)
    f4 = fun ( x0 + dx,       y0 + dx * f3, *args)
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