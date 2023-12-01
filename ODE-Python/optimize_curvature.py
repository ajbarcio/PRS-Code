import time
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

def deformSpring(springObject, angle):
    stuff = things

class spring(object):

    def __init__(self, material, **kwargs):
        ## PARAMETERS TO ITERATE:
            # these are default values, which can be changed by passing named arguments to the class

        self.Vparams = {'rotorThickness': 0.88,              \
                        'ctrlThickness':  0.27,             \
                        'minThickness':   0.23,             \
                        'p3radius':       2.13,             \
                        'p6radius':       2,              \
                        'p3angle':        deg2rad*45,       \
                        'p6angle':        deg2rad*110,      \
                        'p9angle':        (180)*deg2rad, \
                        'tangentProp':    1                 }

        ## INVARIANT PARAMETERS

        self.innerRadius = 1.62/2
        self.outerRadius = 5/2
        
        # Set variable parameters to specified values

        for key, value in kwargs.items():
            self.Vparams[key] = value

        # Parameters set by req. arg or default

        self.material = material
        self.nanFlag = 1

        # This is really stupid but I started using a Dict AFTER I 
        # wrote all my code for generating geometry so this is what you get

        self.rotorThickness = self.Vparams["rotorThickness"]
        self.ctrlThickness  = self.Vparams["ctrlThickness"]
        self.minThickness   = self.Vparams["minThickness"]
        self.p3radius       = self.Vparams["p3radius"]
        self.p6radius       = self.Vparams["p6radius"]
        self.p3angle        = self.Vparams["p3angle"]
        self.p6angle        = self.Vparams["p6angle"]
        self.p9angle        = self.Vparams["p9angle"]
        self.tangentProp    = self.Vparams["tangentProp"]

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

        p9 = np.array([self.outerRadius*np.cos(self.p9angle), self.outerRadius*np.sin(self.p9angle)])
        p8 = np.array([p9[0]-self.tangentProp*self.rotorThickness*np.cos(self.p9angle), \
                       p9[1]-self.tangentProp*self.rotorThickness*np.sin(self.p9angle) ])
        p6 = np.array([self.p6radius*np.cos(self.p6angle), self.p6radius*np.sin(self.p6angle)])
        p3 = np.array([self.p3radius*np.cos(self.p3angle), self.p3radius*np.sin(self.p3angle)])
        p0 = np.array([self.innerRadius, 0])
        p1 = np.array(p0+[self.rotorThickness,0])
        
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
        self.pts = np.array([p0,p1,result[0],p3,result[1],result[2],p6,result[3],p8,p9])
        return self.pts
    
    def reassign_points(self, freePts):
        
        #input points

        p1 = freePts[0]
        p8 = freePts[-1]

        p3 = freePts[1]
        p6 = freePts[2]

        # fixed points

        p0 = self.pts[0]
        p9 = self.pts[9]

        A = np.array([[4,0,1,0],[1,1,0,0],[0,-1,4,0],[0,0,1,1]])
        B = np.array([4*p3+p1,2*p3,4*p6-p8,2*p6])

        result = lin.solve(A,B)
        self.pts = np.array([p0,p1,result[0],p3,result[1],result[2],p6,result[3],p8,p9])

    def generate_centroidal_profile(self):
            # print("profile")
            u = np.linspace(0,3,3001)
            offset = 0
            # print(u)
            u2 = u - np.floor(u)
            self.xyc = np.empty([len(u), 2])
            self.ang   = np.empty(len(u))
            self.norm = np.empty([len(u), 2])
            self.ck = np.empty([len(u), 2])
            self.rc   = np.empty(len(u))
            self.d    = np.empty(len(u))
            self.dd   = np.empty(len(u))

            #  print(self.pts)
    
            for i in range(len(u)):
                ptndx = int(np.floor(u[i]))
                t = u2[i]
                #cart[i] = pts[3*ptndx]*(1-t)**3+3*t*pts[3*ptndx+1]*(1-t)**2+3*pts[3*ptndx+2]*(1-t)**2+pts[3*ptndx+3]*t**3
                if (ptndx == 3):
                    t = 1
                    ptndx = 2
                    U = np.array([t**3, t**2, t, 1])
                    D = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
                    P = np.array([self.pts[3*ptndx],self.pts[3*ptndx+1],self.pts[3*ptndx+2],self.pts[3*ptndx+3]])
                    self.xyc[i] = U.dot(D).dot(P)-[0,offset]
                    self.ang[i] = np.arctan2(self.xyc[i,1],self.xyc[i,0])
                else:
                    U = np.array([t**3, t**2, t, 1])
                    D = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
                    P = np.array([self.pts[3*ptndx],self.pts[3*ptndx+1],self.pts[3*ptndx+2],self.pts[3*ptndx+3]])
                    self.xyc[i] = U.dot(D).dot(P)-[0,offset]
                    self.ang[i] = np.arctan2(self.xyc[i,1],self.xyc[i,0])
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
                self.norm[i] = [-tan[1],tan[0]]/lin.norm([-tan[1],tan[0]])
                d2xydt2 = ddU.dot(coeffs)
                d2ydx2 = (d2xydt2[1] - dydx*d2xydt2[0])/(dxydt[0])
                # d2ydx2 = (dydx-dydxprev)/(cart[i,0]-cart[i-1,0])
                self.rc[i] = ((1+(dydx)**2)**(1.5)/((d2ydx2)))
                self.ck[i] = self.xyc[i]-self.norm[i]*self.rc[i]
                self.d[i] = (dydx)
                self.dd[i] = (d2ydx2)
            #print("about to call next")
            #self.generate_thickness_profile()
            return self.xyc, self.norm, self.rc
        #print(cart)

    # def generate_centroidal_radius(self, cartc):
        
    #     rc = np.empty([len(cartc)])
    #     for i in range(len(rc)):
    #         rc[i] = np.sqrt(cartc[i,0]**2+cartc[i,1]**2)
    #     return rc
    
    def generate_thickness_profile(self):
        # u = np.linspace(0,3,3001)
        # self.ttop = np.empty([len(u), 2])
        # self.tbottom = np.empty([len(u), 2])
        # print("thickness")
        u = np.linspace(0,3,3001)
        self.ttop = np.empty([len(u), 2])
        self.tbottom = np.empty([len(u), 2])
        
        # a b c d e f g
        # self.rotorThickness = g
        # self.ctrlThickness  = a + b + c + d + e + f + g
        # self.minThickness   = a*1.5**6+b*1.5**5+c*1.5**4+d*1.5**3+e*1.5**2+f*1.5+g
        # self.ctrlThickness  = a*2**6+b*2**5+c*2**4+d*2**3+e*2**2+f*2+g
        # self.rotorThickness = a*3**6+b*3**5+c*3**4+d*3**3+e*3**2+f*3+g
        # 0 = f
        # 0 = a*3**5+b*3**4+c*3**3+d*3**2+e*3+f
        REG = np.array([[0,0,0,0,0,0,1],[1,1,1,1,1,1,1], \
                        [1.5**6,1.5**5,1.5**4,1.5**3,1.5**2,1.5,1], \
                        [2**6,2**5,2**4,2**3,2**2,2,1], \
                        [3**6,3**5,3**4,3**3,3**2,3,1], \
                        [0,0,0,0,0,1,0],[6*3**5,5*3**4,4*3**3,3*3**2,2*3,1,0]])
        Targ = np.array([[self.rotorThickness],[self.ctrlThickness], \
                         [self.minThickness],[self.ctrlThickness], \
                         [self.rotorThickness],[0],[0]])
        coeffs = lin.solve(REG,Targ)
        u = np.linspace(0,3,3001)
        self.thks = coeffs[0]*u**6 + coeffs[1]*u**5 + coeffs[2]*u**4 + coeffs[3]*u**3 + coeffs[4]*u**2 + coeffs[5]*u + coeffs[6]
        # self.thks = u/u*self.ctrlThickness
        # coeffs = lin.inv(REG.T.dot(REG)).dot(REG.T).dot(Targ)
        # self.generate_neutral_profile()
        for i in range(len(u)):
            # if abs(self.rc[i]) < .5*self.thks[i]:
            #     self.rotorThickness = self.rotorThickness*0.99
            #     self.generate_thickness_profile()
            self.ttop[i] = self.xyc[i]+self.norm[i]*0.5*self.thks[i]
            self.tbottom[i] = self.xyc[i]-self.norm[i]*0.5*self.thks[i]
        return self.thks
    
    # def generate_neutral_radius(self, rc, t):
    #     rn = t/np.log((rc+.5*t)/(rc-.5*t))
    #     return rn
    
    def generate_neutral_profile(self):
        # print("neutral")
        u = np.linspace(0,3,3001)
        self.xyn = np.empty([len(u), 2])
        self.rn  = np.empty(len(u))
        for i in range(len(u)):
            if np.isinf(self.rc[i]):
                self.rn[i] = self.rc[i]
            else:
                self.rn[i] = np.sign(self.rc[i])*self.thks[i]/abs(np.log(abs(self.rc[i])+.5*self.thks[i])-np.log(abs(self.rc[i])-.5*self.thks[i]))
                # if np.isnan(self.rn[i]):
                    # print(i, self.rn[i])
                if not np.sign(self.rc[i])==np.sign(self.rn[i]):
                    self.rn[i] = -1*self.rn[i]
                    # print("tripped", i)
            diff = self.rc[i] - self.rn[i]
            if np.isnan(diff):
                diff = 1
            # if np.isnan(self.rn[i]) and not(self.tangentProp>2):
            #      print(self.tangentProp)
            #      self.tangentProp = self.tangentProp*1.01
            #      print(self.tangentProp)
            #      return
            # if i%1 == 0:
                # print(diff, self.rc[i], rn)
            self.xyn[i] = self.xyc[i]+self.norm[i]*diff
            self.norm[i] = self.norm[i]*diff
            # print(lin.norm(self.norm[i]))
        # print(self.tangentProp)
        # if (not(np.isnan(self.rn).any()) or self.tangentProp > 2):
        #    self.plotResults()
        #   self.nanFlag = 0
        # self.plotResults()
        return self.norm, self.rn
        
    
    # def generate_neutral_profile2(self, rn, th):
    #     u = np.linspace(0,3,3001)
    #     cart = np.empty([len(u), 2])
    #     for i in range(len(u)):
    #         cart[i,0] = rn[i]*np.cos(th[i])
    #         cart[i,1] = rn[i]*np.sin(th[i])
    #     return cart 
    
    def plotResults(self, oldPts):
        print("plot")
        plt.figure(1)
        plt.clf()
        plt.plot(self.xyc[:,0], self.xyc[:,1])
        plt.plot(self.xyn[:,0], self.xyn[:,1])
        plt.plot(self.pts[:,0], self.pts[:,1])
        plt.plot(oldPts[:,0], oldPts[:,1])
        plt.plot(self.ttop[:,0], self.ttop[:,1])
        plt.plot(self.tbottom[:,0], self.tbottom[:,1])
        for i in range(len(np.linspace(0,3,3001))):
            if i%25 == 0:
                #plt.arrow(xyc[i,0], xyc[i,1], norm2[i,0], norm2[i,1])
                plt.arrow(self.xyc[i,0], self.xyc[i,1], self.norm[i,0], self.norm[i,1])
        plt.figure(2)
        plt.clf()
        # plt.plot(norm[:,0],norm[:,1])
        normnorm = np.empty(len(self.norm))
        for i in range(len(np.linspace(0,3,3001))):
            normnorm[i] = lin.norm(self.norm[i])
        plt.plot(np.linspace(0,3,3001), (self.rn))
        plt.plot(np.linspace(0,3,3001), (self.rc))
        plt.plot(np.linspace(0,3,3001), -.5*(self.thks))
        plt.plot(np.linspace(0,3,3001), .5*(self.thks))
        # plt.plot(np.linspace(0,3,3001), normnorm)
        # plt.plot(np.linspace(0,3,3001), self.d)
        # plt.plot(np.linspace(0,3,3001), self.thks)
        #  plt.plot([-1,2])
        plt.figure(3)
        plt.clf()
        # self.rc[i] = ((1+(dydx)**2)**(1.5)/((d2ydx2)))
        plt.plot(np.linspace(0,3,3001), (1+(self.d)**2)**1.5)
        plt.plot(np.linspace(0,3,3001), self.dd)
        plt.ylim([-15,50])

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
    def full_spline(freePts):
        # print("in full spline")
        fixedPts = startingParameters.pts
        # free points are 1, 3, 6, 8
        
        stackFreePts = np.array([[freePts[0],0],[freePts[1],freePts[2]],[freePts[3],freePts[4]],[freePts[-1]*np.cos(startingParameters.p9angle),freePts[-1]*np.sin(startingParameters.p9angle)]])
        # p2 = 2*fixedPts[3]-stackFreePts[1]
        # p7 = 2*fixedPts[6]-stackFreePts[2]
        # freePts[7] = (fixedPts[9,1]/fixedPts[9,0])*freePts[6]
        # startingParameters.pts = np.array([fixedPts[0],stackFreePts[0],p2,fixedPts[3],stackFreePts[1],stackFreePts[2],fixedPts[6],p7,stackFreePts[3],fixedPts[9]])
        startingParameters.reassign_points(stackFreePts)

        [xyc, rc, norm] = startingParameters.generate_centroidal_profile()
        thks = startingParameters.generate_thickness_profile()
        [norm, rn] = startingParameters.generate_neutral_profile()
        # print(norm)
        return 1/lin.norm(abs(rn), -1)

    material = MaraginSteelC300()
    startingParameters = spring(material, rotorThickness=.88, p9angle=(180)*deg2rad)
    exitflag = 0
    j = 0
    while not exitflag == 1:
        print("entered")
        pts = startingParameters.generate_ctrl_points()
        [xycold, rcold, normold] = startingParameters.generate_centroidal_profile()
        thksold = startingParameters.generate_thickness_profile()
        for i in range(len(np.linspace(0,3,3001))):
                if abs(startingParameters.rc[i]) < .5*startingParameters.thks[i]:
                    exitflag = 2
        if exitflag == 0:
            exitflag = 1
        if exitflag == 2:
            startingParameters.tangentProp = startingParameters.tangentProp+0.001
            exitflag = 0
        j+=1
        print(j)
    [normold, rnold] = startingParameters.generate_neutral_profile()
    #startingParameters.plotResults()
    #plt.show(block=False)
    #time.sleep(15)

    print(startingParameters.pts)
    oldPts = pts

    # Bnds = ((np.sqrt(pts[1,0]**2+pts[1,1]**2),None), (None,None), (None,None), (None,None), (None,None), (None,np.sqrt(pts[8,0]**2+pts[8,1]**2)))
    Bnds = ((None,None), (None,None), (None,None), (None,None), (None,None), (None,None))

    freePts = np.array([np.sqrt(pts[1,0]**2+pts[1,1]**2), pts[3,0], pts[3,1], pts[6,0], pts[6,1], np.sqrt(pts[8,0]**2+pts[8,1]**2)])
    # full_spline(freePts)
    res = op.minimize(full_spline, freePts, method='Nelder-Mead', bounds=Bnds)
    print(res)
    stackResults = np.array([[res.x[0],0],[res.x[1],res.x[2]],[res.x[3],res.x[4]],[res.x[-1]*np.cos(startingParameters.p9angle),res.x[-1]*np.sin(startingParameters.p9angle)]])
    startingParameters.reassign_points(stackResults)
    [xyc, rc, norm] = startingParameters.generate_centroidal_profile()
    thks = startingParameters.generate_thickness_profile()
    [norm, rn] = startingParameters.generate_neutral_profile()
    print(startingParameters.pts)
    plt.close()
    startingParameters.plotResults(oldPts)
    print(startingParameters.rotorThickness)
    plt.show()

if __name__ == '__main__':
    main()

