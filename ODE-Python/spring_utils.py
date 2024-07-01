# TODO: Put spring utils in here (function is spring file which aren't in spring class)
import numpy as np
import math
from polynomials import PPoly_Eval

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

def d_alpha(s, xCoeffs, yCoeffs):
    d2yds2 = PPoly_Eval(s, yCoeffs, deriv=2)
    d2xds2 = PPoly_Eval(s, xCoeffs, deriv=2)
    dyds   = PPoly_Eval(s, yCoeffs, deriv=1)
    dxds   = PPoly_Eval(s, xCoeffs, deriv=1)
    dads = ((d2yds2/dxds-d2xds2*dyds/dxds**2)/(1+dyds**2/dxds**2))
    return dads

def d_2_alpha(s, xCoeffs, yCoeffs):
    d2ads2 = -d_alpha(s, xCoeffs, yCoeffs)**2*d_rn(s, xCoeffs, yCoeffs)
    if np.isnan(d2ads2):
        d2ads2 = 0
    return d2ads2

def r_n(s, xCoeffs, yCoeffs):
    d2yds2 = PPoly_Eval(s, yCoeffs, deriv=2)
    d2xds2 = PPoly_Eval(s, xCoeffs, deriv=2)
    dyds   = PPoly_Eval(s, yCoeffs, deriv=1)
    dxds   = PPoly_Eval(s, xCoeffs, deriv=1)
    # dAlphadS = (d2yds2/dxds-d2xds2*dyds/dxds)/(1+dyds**2/dxds)
    if hasattr(s, "__len__"):
        rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
        for i in range(len(rn)):
            if not np.isfinite(rn[i]):
                rn[i] = float('inf')
            if abs((d2yds2[i]/dxds[i]-d2xds2[i]*dyds[i]/dxds[i]**2)) <= 10**-13:
                rn[i] = float('inf')*np.sign(rn[i-1])
    else:
        if abs(d2yds2/dxds-d2xds2*dyds/dxds**2)<=10**-13:
            rn = float('inf')
        else:
            rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
    return rn

def d_rn(s, xCoeffs, yCoeffs):
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