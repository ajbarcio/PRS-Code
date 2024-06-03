import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from polynomials import Ic_multiPoly, Ic_poly, Ic_spline, PPoly_Eval

################################################################################
# This script tests the functions that generate Ic
# (curved second moment of area) profiles. Change the variables
# IcPts and IcArcLens to test the limits of the functions.
################################################################################

# Input variables:
# IcPts:     values (in in^4) to be enforced for Ic at different points
#            along the beam
# IcArcLens: "arc lengths" (mathematically parameter lengths) at which
#            associated values in IcPts are to be enforced
#
# Lists should be able to be any (matching) length and
# include any (real positive) values and produce good results

IcPts = np.array([.008,.001,.006])
IcParamLens = np.array([0,3,5])

method = "multiPoly"
# the "multiPoly" method is the only current "supported" method for Ic generation
# the other two, poly and spline, are left here for posterity and record-keeping

## POLYNOMIAL METHOD
if method=="poly":
    # This method is deprecated, and so is not fully documented, but in essence
    # it uses the same strategy as the current method of generating polynomials
    # for X and Y coordinates of the netural surface
    IcCoeffs = Ic_poly(IcPts,IcParamLens)
    # print(IcCoeffs)
    for i in range(len(IcCoeffs)):
        print(IcCoeffs[i][0])

    mesh = np.linspace(IcParamLens[0],IcParamLens[-1],101)
    Ic = PPoly_Eval(mesh,IcCoeffs)
    dIcds = PPoly_Eval(mesh,IcCoeffs,deriv=1)
    d2Icds2 = PPoly_Eval(mesh,IcCoeffs,deriv=2)
    # d3Icds3 = Ic_s(smesh,IcCoeffs,deriv=3)
    plt.plot(mesh,Ic)
    plt.plot(mesh,dIcds)
    plt.plot(mesh,d2Icds2)
    # plt.plot(smesh,d3Icds3)
## SPLINE METHOD
elif method=="spline":
    # this method uses symbolic algebra to turn second order bezier splines into
    # quadratic polynomials. This is a bad idea, as the splines do not always
    # form quadratic polynomials with axes aligned with x and y. Therefore, the
    # "multipoly" method is used, which simply generates these quadratic
    # polynomials directly
    soln, diff, splPts = Ic_spline(IcPts,IcParamLens)
    print(soln)
    print(diff)
    x = sp.Symbol('x')
    sp.plotting.plot((soln[0],(x,splPts[0],splPts[1])),(soln[1],(x,splPts[1],splPts[2])),(soln[2],(x,splPts[2],splPts[3])),(soln[3],(x,splPts[3],splPts[4])))
    sp.plotting.plot((diff[0],(x,splPts[0],splPts[1])),(diff[1],(x,splPts[1],splPts[2])),(diff[2],(x,splPts[2],splPts[3])),(diff[3],(x,splPts[3],splPts[4])))
    # plt.plot(theList[:,0],theList[:,1])
## PIECEWISE POLYNOMIAL METHOD (This is the currently used version)
elif method=="multiPoly":
    # Generate the coefficients of multiple quadratic polynomials and the
    # domains over which they apply so that the input variables are satisfied
    IcCoeffs, domains = Ic_multiPoly(IcPts, IcParamLens)
    # evaluate the value of the profile at 0 (light check to indicate whether
    # domains are likely being enforced properly)
    Ic_test = PPoly_Eval(0,IcCoeffs,ranges=domains)
    print(IcCoeffs)

    print(Ic_test)
    # generate mesh over which to evaluate profile
    mesh = np.linspace(IcParamLens[0],IcParamLens[-1],101)
    # generate profile along mesh
    Ic = PPoly_Eval(mesh,IcCoeffs,ranges=domains)
    # plot the profile and its derivative
    plt.plot(mesh,Ic)
    dIcds = PPoly_Eval(mesh, IcCoeffs, ranges=domains, deriv=1)
    plt.plot(mesh,dIcds)

plt.show()