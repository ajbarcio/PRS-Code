import numpy as np
import scipy as sp
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

class parameters(object):

    def __init__(self, material, **kwargs):
        ## PARAMETERS TO ITERATE:
            # these are default values, which can be changed by passing named arguments to the class
        # Thickness parameters

        self.rotorThickness = .8
        self.ctrlThickness = 0.27
        self.minThickness = .23

        # Spline geometry parameters

        self.p3radius = 1.9
        self.p6radius = 1.8
        self.p3angle = deg2rad*45
        self.p6angle = deg2rad*115
        self.p9angle = (180-10)*deg2rad

        # This is a MINIMUM STARTING VALUE that MAY be updated as time goes on

        self.tangentProp = 2

        ## INVARIANT PARAMETERS

        self.innerRadius = 1.62/2
        self.outerRadius = 5/2
        print(self.rotorThickness)

        for key, value in kwargs.items():
            keystr = str(key)
            self.keystr = value
            print("%s == %s" % (key, value))
            print(self.keystr)
            print(self.rotorThickness)

        # Parameters set by req. arg or default

        self.material = material
        self.nanFlag = 1

# def ODE():

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
        p1 = np.array(p0+[self.tangentProp*self.rotorThickness,0])
        
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
        self.generate_centroidal_profile()
    
    def generate_centroidal_profile(self):
            print("profile")
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
                self.rc[i] = ((1+dydx**2)**(1.5)/(d2ydx2))
                self.ck[i] = self.xyc[i]-self.norm[i]*self.rc[i]
                self.d[i] = dydx
                self.dd[i] = d2ydx2
            print("about to call next")
            self.generate_thickness_profile()
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
        # coeffs = lin.inv(REG.T.dot(REG)).dot(REG.T).dot(Targ)
        # self.generate_neutral_profile()
    
    # def generate_neutral_radius(self, rc, t):
    #     rn = t/np.log((rc+.5*t)/(rc-.5*t))
    #     return rn
    
    def generate_neutral_profile(self):
        print("neutral")
        u = np.linspace(0,3,3001)
        self.xyn = np.empty([len(u), 2])
        self.rn  = np.empty(len(u))
        for i in range(len(u)):
            if np.isinf(self.rc[i]):
                self.rn[i] = self.rc[i]
            else:
                self.rn[i] = np.sign(self.rc[i])*self.thks[i]/np.log(((abs(self.rc[i])+.5*self.thks[i])/(abs(self.rc[i])-.5*self.thks[i])))
            diff = self.rc[i] - self.rn[i]
            if np.isnan(diff):
                b = self.rc[i]+.5*self.thks[i]
                a = self.rc[i]-.5*self.thks[i]
                self.rn[i] = (b-a)/np.log(b/a)
            if np.isnan(self.rn[i]):
                self.tangentProp = self.tangentProp*1.05
                return
            # if i%1 == 0:
                # print(diff, self.rc[i], rn)
            self.xyn[i] = self.xyc[i]+self.norm[i]*diff
            self.norm[i] = self.norm[i]*diff
            # print(lin.norm(self.norm[i]))
        # print(self.tangentProp)
        if(not np.isnan(self.rn).any()):
            self.plotResults()
            self.nanFlag = 0
        
    
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
        plt.plot(self.xyc[:,0], self.xyc[:,1])
        plt.plot(self.xyn[:,0], self.xyn[:,1])
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
        plt.plot(np.linspace(0,3,3001), self.rc)
        plt.plot(np.linspace(0,3,3001), normnorm)
        plt.plot(np.linspace(0,3,3001), self.d)
        plt.plot(np.linspace(0,3,3001), self.thks)



def main():

    material = MaraginSteelC300()
    startingParameters = parameters(material, rotorThickness=.9)
    while startingParameters.nanFlag == 1:
        startingParameters.generate_ctrl_points()
        startingParameters.generate_neutral_profile()
    plt.show()
    print(startingParameters.tangentProp)
    print(startingParameters.rotorThickness)



if __name__ == '__main__':
    main()

