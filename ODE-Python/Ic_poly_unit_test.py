from stiffness_library import *
import numpy as np
import sympy as sp

IcPts = np.array([.008,.008,.008])
IcArcLens = np.array([0,2.5,5])

method = "multiPoly"

## POLYNOMIAL METHOD
if method=="poly":
    IcCoeffs = Ic_poly(IcPts,IcArcLens)
    # print(IcCoeffs)
    for i in range(len(IcCoeffs)):
        print(IcCoeffs[i][0])

    smesh = np.linspace(IcArcLens[0],IcArcLens[-1],101)
    Ic = PPoly_Eval(smesh,IcCoeffs)
    dIcds = PPoly_Eval(smesh,IcCoeffs,deriv=1)
    d2Icds2 = PPoly_Eval(smesh,IcCoeffs,deriv=2)
    # d3Icds3 = Ic_s(smesh,IcCoeffs,deriv=3)
    plt.plot(smesh,Ic)
    plt.plot(smesh,dIcds)
    plt.plot(smesh,d2Icds2)
    # plt.plot(smesh,d3Icds3)
## SPLINE METHOD
elif method=="spline":
    soln, diff, splPts = Ic_spline(IcPts,IcArcLens)
    print(soln)
    print(diff)
    x = sp.Symbol('x')
    sp.plotting.plot((soln[0],(x,splPts[0],splPts[1])),(soln[1],(x,splPts[1],splPts[2])),(soln[2],(x,splPts[2],splPts[3])),(soln[3],(x,splPts[3],splPts[4])))
    sp.plotting.plot((diff[0],(x,splPts[0],splPts[1])),(diff[1],(x,splPts[1],splPts[2])),(diff[2],(x,splPts[2],splPts[3])),(diff[3],(x,splPts[3],splPts[4])))
    # plt.plot(theList[:,0],theList[:,1])
elif method=="multiPoly":
    IcCoeffs, domains = Ic_multiPoly(IcPts, IcArcLens)
    # print(IcCoeffs, domains)
    plt.plot()

    Ic_test = PPoly_Eval(0,IcCoeffs,ranges=domains)
    print(IcCoeffs)
    print(Ic_test)

    smesh = np.linspace(IcArcLens[0],IcArcLens[-1],101)
    Ic = PPoly_Eval(smesh,IcCoeffs,ranges=domains)
    plt.plot(smesh,Ic)
    dIcds = PPoly_Eval(smesh, IcCoeffs, ranges=domains, deriv=1)
    plt.plot(smesh,dIcds)
    ic = PPoly_Eval(2.5,IcCoeffs, ranges=domains)
    print(ic)


plt.show()