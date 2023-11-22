import numpy as np
import scipy as sp
from scipy import optimize as op
import numpy.linalg as lin
import scipy.linalg as slin
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from materials import *
from openpyxl import *

deg2rad = np.pi/180

def deformSpring(material **kwargs):
        ## PARAMETERS TO ITERATE:
            # these are default values, which can be changed by passing named arguments to the class

        Vparams = {'rotorThickness': 0.88,              \
                        'ctrlThickness':  0.27,             \
                        'minThickness':   0.23,             \
                        'p3radius':       2.13,             \
                        'p6radius':       2,              \
                        'p3angle':        deg2rad*45,       \
                        'p6angle':        deg2rad*110,      \
                        'p9angle':        (180)*deg2rad, \
                        'tangentProp':    1                 }

        ## INVARIANT PARAMETERS

        innerRadius = 1.62/2
        outerRadius = 5/2
        
        # Set variable parameters to specified values

        for key, value in kwargs.items():
            Vparams[key] = value

        # Parameters set by req. arg or default

        material = material
        nanFlag = 1

        # This is really stupid but I started using a Dict AFTER I 
        # wrote all my code for generating geometry so this is what you get

        rotorThickness = Vparams["rotorThickness"]
        ctrlThickness  = Vparams["ctrlThickness"]
        minThickness   = Vparams["minThickness"]
        p3radius       = Vparams["p3radius"]
        p6radius       = Vparams["p6radius"]
        p3angle        = Vparams["p3angle"]
        p6angle        = Vparams["p6angle"]
        p9angle        = Vparams["p9angle"]
        tangentProp    = Vparams["tangentProp"]

    def generate_ctrl_points(self):
        print("control points")
        # Points predefined by parameters and arbitrary constraints:

        # p9 on outer radius at angle
        # p8 radial with p9, 1 'rotorThickness' away 
        #       (scaled by tangent prop parameter to keep neutral radius from blowing up)
        # p6 at angle and distance parameters
        # p3 at angle and distance parameters
        # p0 on inner radius at horizontal 
        # p1 radial with p0, 1 'rotorThickness' away
        #       (scaled by tangent prop parameter to keep neutral radius from blowing up)

        p9 = np.array([outerRadius*np.cos(p9angle), outerRadius*np.sin(p9angle)])
        p8 = np.array([p9[0]-tangentProp*rotorThickness*np.cos(p9angle), \
                       p9[1]-tangentProp*rotorThickness*np.sin(p9angle) ])
        p6 = np.array([p6radius*np.cos(p6angle), p6radius*np.sin(p6angle)])
        p3 = np.array([p3radius*np.cos(p3angle), p3radius*np.sin(p3angle)])
        p0 = np.array([innerRadius, 0])
        p1 = np.array(p0+[tangentProp*rotorThickness,0])
        
        # Rules for the spline:

        # C1:

        # P4 = P3 + (P3-P2)
        # P7 = P6 + (P6-P5)

        # C2:

        # p5 = p1 + 4*(p3-p2)
        # p8 = p4 + 4*(p6-p5)

        # Boundary Conditions

        # knowns: 0, 1, 3, 6, 8, 9
        # unknowns: 2, 4, 5, 7

        A = np.array([[4,0,1,0],[1,1,0,0],[0,-1,4,0],[0,0,1,1]])
        B = np.array([4*p3+p1,2*p3,4*p6-p8,2*p6])

        result = lin.solve(A,B)
        pts = np.array([p0,p1,result[0],p3,result[1],result[2],p6,result[3],p8,p9])
        generate_centroidal_profile()
    
    def generate_centroidal_profile(self):
            print("profile")
            u = np.linspace(0,3,3001)
            offset = 0
            # print(u)
            u2 = u - np.floor(u)
            xyc = np.empty([len(u), 2])
            ang   = np.empty(len(u))
            norm = np.empty([len(u), 2])
            ck = np.empty([len(u), 2])
            rc   = np.empty(len(u))
            d    = np.empty(len(u))
            dd   = np.empty(len(u))

            print(pts)
    
            for i in range(len(u)):
                ptndx = int(np.floor(u[i]))
                t = u2[i]
                #cart[i] = pts[3*ptndx]*(1-t)**3+3*t*pts[3*ptndx+1]*(1-t)**2+3*pts[3*ptndx+2]*(1-t)**2+pts[3*ptndx+3]*t**3
                if (ptndx == 3):
                    t = 1
                    ptndx = 2
                    U = np.array([t**3, t**2, t, 1])
                    D = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
                    P = np.array([pts[3*ptndx],pts[3*ptndx+1],pts[3*ptndx+2],pts[3*ptndx+3]])
                    xyc[i] = U.dot(D).dot(P)-[0,offset]
                    ang[i] = np.arctan2(xyc[i,1],xyc[i,0])
                else:
                    U = np.array([t**3, t**2, t, 1])
                    D = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
                    P = np.array([pts[3*ptndx],pts[3*ptndx+1],pts[3*ptndx+2],pts[3*ptndx+3]])
                    xyc[i] = U.dot(D).dot(P)-[0,offset]
                    ang[i] = np.arctan2(xyc[i,1],xyc[i,0])
                dU = np.array([3*t**2, 2*t, 1, 0])
                ddU = np.array([6*t, 2, 0, 0])
                coeffs = D.dot(P)
                if i == 0: 
                    # tan  = [cart[i+1,0]-cart[i,0],cart[i+1,1]-cart[i,1]]/np.sqrt((cart[i+1,0]-cart[i,0])**2+(cart[i+1,1]-cart[i,1])**2)
                    dxydt = dU.dot(coeffs)
                    dydx = dxydt[1]/dxydt[0]
                    tan = [dxydt[0],dxydt[1]]
                    #dydx = (cart[i+1,1]-cart[i,1])/(cart[i+1,0]-cart[i,0])
                    # dydxprev = 0
                elif i == 3000:
                    # dydxprev = dydx
                    # tan  = [cart[i,0]-cart[i-1,0],cart[i,1]-cart[i-1,1]]/np.sqrt((cart[i,0]-cart[i-1,0])**2+(cart[i,1]-cart[i-1,1])**2)
                    dxydt = dU.dot(coeffs)
                    dydx = dxydt[1]/dxydt[0]
                    tan = [dxydt[0],dxydt[1]]
                    #dydx = (cart[i,1]-cart[i-1,1])/(cart[i,0]-cart[i-1,0])
                else:
                    # dydxprev = dydx
                    # tan  = [cart[i+1,0]-cart[i-1,0],cart[i+1,1]-cart[i-1,1]]/np.sqrt((cart[i+1,0]-cart[i-1,0])**2+(cart[i+1,1]-cart[i-1,1])**2)
                    dxydt = dU.dot(coeffs)
                    dydx = dxydt[1]/dxydt[0]
                    tan = [dxydt[0],dxydt[1]]
                    # dydx = (cart[i+1,1]-cart[i-1,1])/(cart[i+1,0]-cart[i-1,0])
                norm[i] = [-tan[1],tan[0]]/lin.norm([-tan[1],tan[0]])
                d2xydt2 = ddU.dot(coeffs)
                d2ydx2 = (d2xydt2[1] - dydx*d2xydt2[0])/(dxydt[0])
                # d2ydx2 = (dydx-dydxprev)/(cart[i,0]-cart[i-1,0])
                rc[i] = ((1+dydx**2)**(1.5)/(d2ydx2))
                ck[i] = xyc[i]-norm[i]*rc[i]
                d[i] = dydx
                dd[i] = d2ydx2
            #print("about to call next")
            generate_thickness_profile()
            return np.min(rc)
        #print(cart)

    # def generate_centroidal_radius(self, cartc):
        
    #     rc = np.empty([len(cartc)])
    #     for i in range(len(rc)):
    #         rc[i] = np.sqrt(cartc[i,0]**2+cartc[i,1]**2)
    #     return rc
    
    def generate_thickness_profile(self):
        print("thickness")
        u = np.linspace(0,3,3001)
        # a b c d e f g
        # rotorThickness = g
        # ctrlThickness  = a + b + c + d + e + f + g
        # minThickness   = a*1.5**6+b*1.5**5+c*1.5**4+d*1.5**3+e*1.5**2+f*1.5+g
        # ctrlThickness  = a*2**6+b*2**5+c*2**4+d*2**3+e*2**2+f*2+g
        # rotorThickness = a*3**6+b*3**5+c*3**4+d*3**3+e*3**2+f*3+g
        # 0 = f
        # 0 = a*3**5+b*3**4+c*3**3+d*3**2+e*3+f
        REG = np.array([[0,0,0,0,0,0,1],[1,1,1,1,1,1,1], \
                        [1.5**6,1.5**5,1.5**4,1.5**3,1.5**2,1.5,1], \
                        [2**6,2**5,2**4,2**3,2**2,2,1], \
                        [3**6,3**5,3**4,3**3,3**2,3,1], \
                        [0,0,0,0,0,1,0],[6*3**5,5*3**4,4*3**3,3*3**2,2*3,1,0]])
        Targ = np.array([[rotorThickness],[ctrlThickness], \
                         [minThickness],[ctrlThickness], \
                         [rotorThickness],[0],[0]])
        coeffs = lin.solve(REG,Targ)
        u = np.linspace(0,3,3001)
        thks = coeffs[0]*u**6 + coeffs[1]*u**5 + coeffs[2]*u**4 + coeffs[3]*u**3 + coeffs[4]*u**2 + coeffs[5]*u + coeffs[6]
        # coeffs = lin.inv(REG.T.dot(REG)).dot(REG.T).dot(Targ)
        # generate_neutral_profile()
    
    # def generate_neutral_radius(self, rc, t):
    #     rn = t/np.log((rc+.5*t)/(rc-.5*t))
    #     return rn
    
    def generate_neutral_profile(self):
        print("neutral")
        u = np.linspace(0,3,3001)
        xyn = np.empty([len(u), 2])
        rn  = np.empty(len(u))
        for i in range(len(u)):
            if np.isinf(rc[i]):
                rn[i] = rc[i]
            else:
                rn[i] = np.sign(rc[i])*thks[i]/np.log(((abs(rc[i])+.5*thks[i])/(abs(rc[i])-.5*thks[i])))
            diff = rc[i] - rn[i]
            if np.isnan(diff):
                b = rc[i]+.5*thks[i]
                a = rc[i]-.5*thks[i]
                rn[i] = (b-a)/np.log(b/a)
            if np.isnan(rn[i]) and not(tangentProp>2):
                 print(tangentProp)
                 tangentProp = tangentProp*1.01
                 print(tangentProp)
                 return
            # if i%1 == 0:
                # print(diff, rc[i], rn)
            xyn[i] = xyc[i]+norm[i]*diff
            norm[i] = norm[i]*diff
            # print(lin.norm(norm[i]))
        # print(tangentProp)
        if (not(np.isnan(rn).any()) or tangentProp > 2):
            plotResults()
            nanFlag = 0
        
    
    # def generate_neutral_profile2(self, rn, th):
    #     u = np.linspace(0,3,3001)
    #     cart = np.empty([len(u), 2])
    #     for i in range(len(u)):
    #         cart[i,0] = rn[i]*np.cos(th[i])
    #         cart[i,1] = rn[i]*np.sin(th[i])
    #     return cart 
    
    def plotResults(self):
        print("plot")
        plt.figure(1)
        plt.clf()
        plt.plot(xyc[:,0], xyc[:,1])
        plt.plot(xyn[:,0], xyn[:,1])
        plt.plot(pts[:,0], pts[:,1])
        for i in range(len(np.linspace(0,3,3001))):
            if i%25 == 0:
                #plt.arrow(xyc[i,0], xyc[i,1], norm2[i,0], norm2[i,1])
                plt.arrow(xyc[i,0], xyc[i,1], norm[i,0], norm[i,1])
        plt.figure(2)
        plt.clf()
        # plt.plot(norm[:,0],norm[:,1])
        normnorm = np.empty(len(norm))
        for i in range(len(np.linspace(0,3,3001))):
            normnorm[i] = lin.norm(norm[i])
        plt.plot(np.linspace(0,3,3001), rc)
        plt.plot(np.linspace(0,3,3001), normnorm)
        plt.plot(np.linspace(0,3,3001), d)
        plt.plot(np.linspace(0,3,3001), thks)

def minimize_curvature(factors):
    #print("profile")
    u = np.linspace(2.25,3,751)
    u2 = u - np.floor(u)
    rc   = np.empty(len(u))
    pts = np.array([[ 8.10000000e-01,  0.00000000e+00], \
                    [ 2.01300190e+00,  0.00000000e+00], \
                    [ 2.03847299e+00 , 9.04558883e-01], \
                    [ 1.50613744e+00,  1.50613744e+00], \
                    [ 9.73801898e-01,  2.10771600e+00], \
                    [-1.16340286e-01 , 2.40631424e+00], \
                    [-6.84040287e-01 , 1.87938524e+00], \
                    [-1.25174029e+00 , 1.35245624e+00], \
                    [-1.29699810e+00 , 1.58836458e-16], \
                    [-2.50000000e+00 , 3.06161700e-16]])
    tempPts = pts
    tempPts[1,:] = pts[1,:]*factors[0]
    tempPts[8,:] = pts[8,:]*factors[1]

    for i in range(len(u)):
        ptndx = int(np.floor(u[i]))
        t = u2[i]
        D = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
        if (ptndx == 3):
            t = 1
            ptndx = 2
            P = np.array([tempPts[3*ptndx],tempPts[3*ptndx+1],tempPts[3*ptndx+2],tempPts[3*ptndx+3]])
        else:
            P = np.array([tempPts[3*ptndx],tempPts[3*ptndx+1],tempPts[3*ptndx+2],tempPts[3*ptndx+3]])
        
        dU = np.array([3*t**2, 2*t, 1, 0])
        ddU = np.array([6*t, 2, 0, 0])
        
        coeffs = D.dot(P)
        dxydt = dU.dot(coeffs)

        dydx = dxydt[1]/dxydt[0]
        d2xydt2 = ddU.dot(coeffs)

        d2ydx2 = (d2xydt2[1] - dydx*d2xydt2[0])/(dxydt[0])
        rc[i] = ((1+dydx**2)**(1.5)/(d2ydx2))
    return 1/abs(np.min(rc))

def main():

    material = MaraginSteelC300()
    startingParameters = spring(material, rotorThickness=.8)
    while startingParameters.nanFlag == 1:
        startingParameters.generate_ctrl_points()
        startingParameters.generate_neutral_profile()
    factors = np.array([1,1])
    # plt.show()
    pts = startingParameters.pts
    print(pts)
    nPts = startingParameters.pts
    result = op.minimize(minimize_curvature, factors, method='COBYLA')
    print(result)
    # print(nPts)
    # print(nPts[1,:])
    newFactors = result.x
    nPts[1] = newFactors[0]*pts[1]
    nPts[8] = newFactors[1]*pts[8]
    print(nPts)
    startingParameters.pts = nPts
    print(startingParameters.pts)
    startingParameters.generate_centroidal_profile()
    startingParameters.generate_neutral_profile()
    # plt.clf()
    # plt.plot(startingParameters.xyc[:,0], startingParameters.xyc[:,1])
    # plt.plot(startingParameters.xyn[:,0], startingParameters.xyn[:,1])
    # plt.plot(startingParameters.pts[:,0], startingParameters.pts[:,1])
    plt.show()
if __name__ == '__main__':
    main()

