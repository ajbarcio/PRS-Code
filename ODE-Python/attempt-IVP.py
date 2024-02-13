import time
import math
import numpy as np
import scipy as sp
from scipy import integrate
from scipy import optimize as op
import numpy.linalg as lin
import scipy.linalg as slin
import matplotlib.pyplot as plt
import matplotlib.transforms as tfm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from materials import *
from openpyxl import *

deg2rad = np.pi/180

EPS = np.finfo(float).eps
funnyNumber = 1

globalRes = 500
globalLen = 3*globalRes+1
globalMaxIndex = 3*globalRes
globalBCresLimit = 2
globalRelresLimit = 10**6

def get_derivative(vec1, vec2):
    dv1dv2 = np.empty(len(vec1))
    for i in range(len(dv1dv2)):
        if i == 0:
            dv1dv2[i] = 1
        elif i == 1:
            dv1dv2[i-1] = (vec1[i]-vec1[i-1])/(vec2[i]-vec2[i-1])
        elif i == len(vec1)-1:
            dv1dv2[i-1] = (vec1[i]-vec1[i-1])/(vec2[i]-vec2[i-1])
            dv1dv2[i] = (vec1[i]-vec1[i-1])/(vec2[i]-vec2[i-1])
        else:
            dv1dv2[i-1] = (vec1[i]-vec1[i-2])/(vec2[i]-vec2[i-2])
    return dv1dv2

class spring(object):

    def __init__(self, material, **kwargs):
        ## PARAMETERS TO ITERATE:
            # these are default values, which can be changed by passing named arguments to the class

        self.Vparams = {'rotorThickness':      0.5,             \
                        'ctrlThickness':       0.5,             \
                        'minThickness':        0.5,             \
                        'p3radius':            .5+2.5/3,             \
                        'p6radius':            .5+5/3,                \
                        'p3angle':             5*deg2rad,       \
                        'p6angle':             10*deg2rad,      \
                        'p9angle':             15*deg2rad,    \
                        'tangentPropStart':    1,           \
                        'tangentPropEnd':      1 }

        ## INVARIANT PARAMETERS

        self.innerRadius = 1/2
        self.outerRadius = 6/2
        
        # Set variable parameters to specified values

        for key, value in kwargs.items():
            self.Vparams[key] = value

        # Parameters set by req. arg or default

        self.material = material
        self.nanFlag = 1
        self.outPlaneThickness = 0.25

        # This is really stupid but I started using a Dict AFTER I 
        # wrote all my code for generating geometry so this is what you get

        self.rotorThickness   = self.Vparams["rotorThickness"]
        self.ctrlThickness    = self.Vparams["ctrlThickness"]
        self.minThickness     = self.Vparams["minThickness"]
        self.p3radius         = self.Vparams["p3radius"]
        self.p6radius         = self.Vparams["p6radius"]
        self.p3angle          = self.Vparams["p3angle"]
        self.p6angle          = self.Vparams["p6angle"]
        self.p9angle          = self.Vparams["p9angle"]
        self.tangentPropStart = self.Vparams["tangentPropStart"]
        self.tangentPropEnd   = self.Vparams["tangentPropEnd"]

    def print_parameters(self):

        print(self.rotorThickness  )
        print(self.ctrlThickness   )
        print(self.minThickness    )
        print(self.p3radius        )
        print(self.p6radius        )
        print(self.p3angle         )
        print(self.p6angle         )
        print(self.p9angle         )
        print(self.tangentPropStart)
        print(self.tangentPropEnd  )

    def deform_force(self, force, test_init):

        u = np.linspace(0,3,globalLen)

        def rk4_deform (u0, gamma0, du):

            #
            #  Get four sample values of the derivative.
            #
            f1 = deformation_function ( u0,            gamma0 )
            f2 = deformation_function ( u0 + du / 2.0, gamma0 + du * f1 / 2.0 )
            f3 = deformation_function ( u0 + du / 2.0, gamma0 + du * f2 / 2.0 )
            f4 = deformation_function ( u0 + du,       gamma0 + du * f3 )
            #
            #  Combine them to estimate the solution U1 at time T1 = T0 + DT.
            #
            gamma1 = gamma0 + du * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0

            return gamma1
        
        largeCouple = np.empty(len(u))

        for i in range(len(u)):
            if np.isinf(self.rn[i]) and self.ecc[i]==0:
                largeCouple[i] = 0
            else:
                largeCouple[i] = self.outPlaneThickness*self.thks[i]*self.material.E*(self.ecc[i])*self.rn[i]

        # largeCouple = np.transpose(largeCouple)
        
        dlCdu       = np.transpose(get_derivative(largeCouple, np.linspace(0,3,globalLen)))
        dxdu        = np.transpose(self.dxyndu[:,0])
        dydu        = np.transpose(self.dxyndu[:,1])
        d2xdu2      = np.transpose(get_derivative(dxdu, np.linspace(0,3,globalLen)))
        d2ydu2      = np.transpose(get_derivative(dydu, np.linspace(0,3,globalLen)))

        sqrt_replace = np.sqrt(np.square(dxdu)+np.square(dydu))
        Fx = force*np.sin(self.p9angle)
        Fy = force*np.cos(self.p9angle)

        def deformation_function(u, gamma0):
            du = 1/globalRes
            i = int(u/du)
            return np.array([gamma0[1], (Fx*np.sin(self.tann[i]+gamma0[0])-Fy*np.cos(self.tann[i]+gamma0[0]))*sqrt_replace[i]/largeCouple[i] \
                              +gamma0[1]*((dxdu[i]*d2xdu2[i]+dydu[i]*d2ydu2[i])-dlCdu[i]/sqrt_replace[i])])

        u = np.linspace(0,3,int((globalLen-1)/2+1))
        u0 = 0
        gamma0 = np.array([0, test_init])
        self.gammaSol = np.empty((u.shape[0], 2))
        du = 1/globalRes*2
        for i in range(len(u)-1):
            if i == 0:
                self.gammaSol[i] = gamma0
                self.gammaSol[i+1] = rk4_deform(u[i], gamma0, du)
                gamma0 = self.gammaSol[i+1]
            else:
                self.gammaSol[i+1] = rk4_deform(u[i], gamma0, du)
                gamma0 = self.gammaSol[i+1]
            print(i,"success", gamma0)

        
        return self.gammaSol[-1,1], d2xdu2, d2ydu2

    def deform_torque(self, torqueIn):
        pass

    def get_max_stress(self):
        pass

    def get_stiffness(self):
        pass
    
    def generate_ctrl_points(self):

        p9 = np.array([self.outerRadius*np.cos(self.p9angle), self.outerRadius*np.sin(self.p9angle)])
        p6 = np.array([self.p6radius*np.cos(self.p6angle), self.p6radius*np.sin(self.p6angle)])
        p3 = np.array([self.p3radius*np.cos(self.p3angle), self.p3radius*np.sin(self.p3angle)])
        p0 = np.array([self.innerRadius, 0])

        p1 = np.array([self.innerRadius+1*(self.outerRadius-self.innerRadius)/9,0])
        p8 = np.array([p9[0]-(self.outerRadius-self.innerRadius)/9*np.cos(self.p9angle), \
                       p9[1]-(self.outerRadius-self.innerRadius)/9*np.sin(self.p9angle) ])
        
        A = np.array([[4,0,1,0],[1,1,0,0],[0,-1,4,0],[0,0,1,1]])
        B = np.array([4*p3+p1,2*p3,4*p6-p8,2*p6])

        result = lin.solve(A,B)
        self.pts = np.array([p0,p1,result[0],p3,result[1],result[2],p6,result[3],p8,p9])
        return self.pts
    
    def generate_centroidal_profile(self, leng, **kwargs):
            # print("profile")
            u = kwargs.get('mesh', [None, None])
            if np.any(u) == None:
                u = np.linspace(0,3,leng)
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
            self.tan  = np.empty(len(u))
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
                dxydt = dU.dot(coeffs)
                dydx = dxydt[1]/dxydt[0]
                tan = np.array([dxydt[0],dxydt[1]])
                # dydx = (cart[i+1,1]-cart[i-1,1])/(cart[i+1,0]-cart[i-1,0])
                self.norm[i] = [-tan[1],tan[0]]/lin.norm([-tan[1],tan[0]])
                d2xydt2 = ddU.dot(coeffs)
                d2ydx2 = (d2xydt2[1] - dydx*d2xydt2[0])/(dxydt[0])
                # d2ydx2 = (dydx-dydxprev)/(cart[i,0]-cart[i-1,0])
                self.rc[i] = ((1+(dydx)**2)**(1.5)/((d2ydx2)))
                self.ck[i] = self.xyc[i]-self.norm[i]*self.rc[i]
                self.d[i] = (dydx)
                self.dd[i] = (d2ydx2)
                self.tan[i] = np.arctan(tan[1]/tan[0])
            #print("about to call next")
            #self.generate_thickness_profile()
            return self.xyc, self.norm, self.rc
        #print(cart)
    
    def generate_thickness_profile(self, leng, **kwargs):
        
        u = kwargs.get('mesh', [None, None])
        if np.any(u) == None:
            u = np.linspace(0,3,leng)
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
        u = np.linspace(0,3,leng)
        self.thks = coeffs[0]*u**6 + coeffs[1]*u**5 + coeffs[2]*u**4 + coeffs[3]*u**3 + coeffs[4]*u**2 + coeffs[5]*u + coeffs[6]
        # self.thks = u/u*self.ctrlThickness
        # coeffs = lin.inv(REG.T.dot(REG)).dot(REG.T).dot(Targ)
        # self.generate_neutral_profile()
        nanerrorflag = 0
        for i in range(len(u)):
            if abs(self.rc[i]) < .5*self.thks[i]:
            #     self.rotorThickness = self.rotorThickness*0.99
            #     self.generate_thickness_profile()
               nanerrorflag = 1
               print("nanerror")
            self.ttop[i] = self.xyc[i]+self.norm[i]*0.5*self.thks[i]
            self.tbottom[i] = self.xyc[i]-self.norm[i]*0.5*self.thks[i]
        # if nanerrorflag == 0:
        #     print("nanerror")
        return self.thks
        
    def generate_neutral_profile(self, leng, **kwargs):
        # print("neutral")
        self.curveError = False
        # write mesh
        u = kwargs.get('mesh', [None, None])
        if np.any(u) == None:
            u = np.linspace(0,3,leng)
        # allocate arrays for outputs: x,y coords for neutral axis, neutral radius afo time,
        # chord length afo time, derivative of x,y of neutral axis afo time, tangent angle of neutral axis afo time
        self.xyn = np.empty([len(u), 2])
        self.rn  = np.empty(len(u))
        self.s   = np.empty(len(u))
        self.dxyndu     = np.empty([len(u), 2])
        self.tann       = np.empty(len(u))
        self.derivCheck = np.empty(len(u))
        self.ecc        = np.empty(len(u))
        # across mesh:
        for i in range(len(u)):
            # if the centroidal axis is flat
            if np.isinf(self.rc[i]):
                # then the neutral axis is also flat
                self.rn[i] = self.rc[i]
            # otherwise
            else:
                # calculate neutral radius
                self.rn[i] = np.sign(self.rc[i])*self.thks[i]/abs(np.log(abs(self.rc[i])+.5*self.thks[i])-np.log(abs(self.rc[i])-.5*self.thks[i]))
                # funny stuff to make sure sign of stuff is right
                if not np.sign(self.rc[i])==np.sign(self.rn[i]):
                    self.rn[i] = -1*self.rn[i]
            # calculate the difference between neutral and centroidal axis
            diff = self.rc[i] - self.rn[i]
            # get the nan error out
            if np.isnan(diff):
                if np.isinf(self.rn[i]) and np.isinf(self.rc[i]):
                    diff = 0
                else:
                    diff = 1
            # record that a nan error occurred
            if np.isnan(self.rn[i]):
                self.curveError = True
            # calculate x/y coordinate of neutral axis
            self.xyn[i] = self.xyc[i]+self.norm[i]*diff

            # calculate deriatives of neutral axis with time
            if i == 0:
                self.dxyndu[i] = 1
                # this is a temporary value that will be overwritten with an accurate one
            elif i == 1:
                # calculate derivative at 0
                self.dxyndu[i-1] = [(self.xyn[i,0]-self.xyn[i-1,0])/(1/globalLen),(self.xyn[i,1]-self.xyn[i-1,1])/(1/globalLen)]
                self.derivCheck[i-1] = 1
            elif i == len(u)-1:
                self.dxyndu[i-1] = [(self.xyn[i,0]-self.xyn[i-2,0])/(2/globalLen),(self.xyn[i,1]-self.xyn[i-2,1])/(2/globalLen)]
                self.derivCheck[i-1] = 4
                self.dxyndu[i] = [(self.xyn[i,0]-self.xyn[i-1,0])/(1/globalLen),(self.xyn[i,1]-self.xyn[i-1,1])/(1/globalLen)]
                self.derivCheck[i] = 3
            else:
                self.dxyndu[i-1] = [(self.xyn[i,0]-self.xyn[i-2,0])/(2/globalLen),(self.xyn[i,1]-self.xyn[i-2,1])/(2/globalLen)]
                self.derivCheck[i-1] = 2
        
            if i == 0:
                self.s[i] = 0
            else:
                integrand = np.sqrt(self.dxyndu[0:i,0]**2+self.dxyndu[0:i,1]**2)*(1/globalLen)
                self.s[i] = np.sum(integrand)
            
            self.norm[i] = self.norm[i]*diff
            self.tann[i] = np.arctan2(self.dxyndu[i,1],self.dxyndu[i,0])

            if np.isinf(self.rn[i]) and np.isinf(self.rc[i]):
                self.ecc[i] = 0
            else:
                self.ecc[i] = self.rc[i]-self.rn[i]

        # if self.curveError == True:
        #     print(self.tangentPropEnd)
        #     print(self.tangentPropStart)
        return self.norm, self.rn, self.curveError
          
    def plotResults(self, oldPts):
        print("plot")
        # Figure 1 holds spring reepresentation
        plt.figure(1)

        # base = plt.gca().transData
        # rot = tfm.Affine2D().rotate_deg(90)

        plt.clf()
        plt.plot(self.xyc[:,0], self.xyc[:,1])
        plt.plot(self.xyn[:,0], self.xyn[:,1])
        plt.plot(self.pts[:,0], self.pts[:,1])
        plt.plot(oldPts[:,0], oldPts[:,1])
        plt.plot(self.ttop[:,0], self.ttop[:,1])
        plt.plot(self.tbottom[:,0], self.tbottom[:,1])

        plt.plot(-self.ttop[:,0], -self.ttop[:,1])
        plt.plot(-self.tbottom[:,0], -self.tbottom[:,1])

        theta = np.linspace(0, 2*np.pi, 100)
        outerCircleX = self.outerRadius*np.cos(theta)
        outerCircleY = self.outerRadius*np.sin(theta)
        innerCircleX = self.innerRadius*np.cos(theta)
        innerCircleY = self.innerRadius*np.sin(theta)
        plt.plot(outerCircleX,outerCircleY)
        plt.plot(innerCircleX,innerCircleY)

        for i in range(len(np.linspace(0,3,globalLen))):
            if i%25 == 0:
                #plt.arrow(xyc[i,0], xyc[i,1], norm2[i,0], norm2[i,1])
                plt.arrow(self.xyc[i,0], self.xyc[i,1], self.norm[i,0], self.norm[i,1])
                # plt.arrow(self.xyc[i,0], self.xyc[i,1], self.ck[i,0], self.ck[i,1])
        self.stressConcentration = np.argmin(abs(self.rn))
        plt.arrow(0,0,self.xyc[self.stressConcentration,0],self.xyc[self.stressConcentration,1])
        # plt.ylim([-1,3])
        plt.xlim([-3,3])
        plt.gca().set_aspect('equal')

        # Figure 2 holds plots of properties
        plt.figure(2)
        plt.clf()
        # plt.plot(norm[:,0],norm[:,1])
        normnorm = np.empty(len(self.norm))
        for i in range(len(np.linspace(0,3,globalLen))):
            normnorm[i] = lin.norm(self.norm[i])
        plt.plot(np.linspace(0,3,globalLen), abs(self.rn))
        plt.plot(np.linspace(0,3,globalLen), (self.rc))
        plt.plot(np.linspace(0,3,globalLen), -.5*(self.thks))
        plt.plot(np.linspace(0,3,globalLen), .5*(self.thks))
        plt.axhline(np.min(abs((self.rn))))
        # plt.plot(np.linspace(0,3,3001), normnorm)
        # plt.plot(np.linspace(0,3,3001), self.d)
        # plt.plot(np.linspace(0,3,3001), self.thks)
        #  plt.plot([-1,2])
        plt.ylim([-1.5,1.5])
        plt.plot(np.linspace(0,3,globalLen), self.s)
        
        plt.plot(np.linspace(0,3,globalLen), np.sqrt(self.xyc[:,0]**2+self.xyc[:,1]**2))

        # Figure 3 holds derivatives
        plt.figure(3)
        plt.clf()
        # self.rc[i] = ((1+(dydx)**2)**(1.5)/((d2ydx2)))
        plt.plot(np.linspace(0,3,globalLen), ((self.d)**2)**1.5)
        plt.plot(np.linspace(0,3,globalLen), self.dd)
        plt.ylim([-15,50])
        print(np.min(abs(self.rn)))
        print(self.stressConcentration)
        print(self.s[-1])

def main():
    def drag_end_points(dragVector):
        startingParameters.tangentPropStart = dragVector[0]
        startingParameters.tangentPropEnd = dragVector[1]   
        startingParameters.p3radius = dragVector[2]
        startingParameters.p6radius = dragVector[3]
        startingParameters.p3angle = dragVector[4]
        startingParameters.p6angle = dragVector[5]     
        
        startingParameters.generate_ctrl_points()

        [xyc, rc, norm] = startingParameters.generate_centroidal_profile(globalLen)
        # sfor i in range(len(np.linspace(0,3,3001))):
            # if np.sqrt(xyc[i,0]**2+xyc[i,1]**2) > startingParameters.outerRadius:
            #     print("too-big-error-tripped")
            #     dragVector = prevDragVector
            #     [xyc, rc, norm] = startingParameters.generate_centroidal_profile()
            #     break
        thks = startingParameters.generate_thickness_profile(globalLen)
        [norm, rn, curveError] = startingParameters.generate_neutral_profile(globalLen)
        prevDragVector = dragVector
        # print(norm)
        # numberMins = np.array([i for i, x in enumerate(abs(startingParameters.rn)) if x == np.min(abs(startingParameters.rn))])
        # print(numberMins)
        # numberMins = 0
        # for i in range(len(np.linspace(0,3,3001))):
        #     if (abs(startingParameters.rn[i]) < min(abs(startingParameters.rn))+.00001):
        #         numberMins+=1
        # localMinDiff = startingParameters.rn[1000] - np.min(abs(startingParameters.rn))
        gainArray = np.array([1.5,1.5,1,2])
        gainArray = gainArray/lin.norm(gainArray,2)
        costArray = np.array([gainArray[0]*1/lin.norm(abs(rn[0:globalRes-1]), ord=-np.inf),gainArray[1]*1/lin.norm(abs(rn[globalRes:2*globalRes]), ord=-np.inf),gainArray[2]*1/lin.norm(abs(rn[globalRes+1:3*globalRes]), ord=-np.inf),gainArray[3]*1/startingParameters.s[-1]])
        
        #
        # print(numberMins)
        print(lin.norm(costArray, ord=2))
        return lin.norm(costArray, ord=2)

        #            {'rotorThickness': 0.88,             \
        #             'ctrlThickness':  0.27,             \
        #             'minThickness':   0.23,             \
        #             'p3radius':       2.13,             \
        #             'p6radius':       2,                \
        #             'p3angle':        deg2rad*45,       \
        #             'p6angle':        deg2rad*110,      \
        #             'p9angle':        (180)*deg2rad,    \
        #             'tangentPropStart':    1,           \
        #             'tangentPropEnd':      1 }

    material = MaraginSteelC300()

    startingParameters = spring(material)
    startingParameters.generate_ctrl_points()
    startingParameters.generate_centroidal_profile(globalLen)
    startingParameters.generate_thickness_profile(globalLen)
    [norm, rn, curveError] = startingParameters.generate_neutral_profile(globalLen)

    compareParameters = spring(material)
    compareParameters.generate_ctrl_points()
    compareParameters.generate_centroidal_profile(globalLen)
    compareParameters.generate_thickness_profile(globalLen)
    [norm, rn, curveError] = compareParameters.generate_neutral_profile(globalLen)
    
    if np.array_equal(compareParameters.xyn,startingParameters.xyn):
        print("neutral surface equal")
    if np.array_equal(compareParameters.xyn,startingParameters.xyn):
        print("derivatives equal")
    if np.array_equal(compareParameters.xyn,startingParameters.xyn):
        print("arc length equal")

    BC1Fun = np.empty((100))
    ang = 5*deg2rad

    [BC1Fun[0], d2xdu2, d2ydu2] = startingParameters.deform_force(100, 0)# -ang
    print(startingParameters.gammaSol)
    print(startingParameters.dxyndu[25])
    print(startingParameters.dxyndu[26])
    print(startingParameters.dxyndu[27])
    
    print(d2xdu2[25])
    print(d2xdu2[26])
    print(d2xdu2[27])

    print(d2ydu2[25])
    print(d2ydu2[26])
    print(d2ydu2[27])
    # step = 0.5*deg2rad
    # err = 1
    # i = 0
    # val = step
    # while err > 0.001 and i < BC1Fun.shape[0]:
    #     prevVal = val
    #     deriv = ((BC1Fun[i+1]-BC1Fun[i])/step)
    #     val = val - BC1Fun[i]/deriv
    #     print(prevVal+val)
    #     BC1Fun[i+1] = startingParameters.deform_force(100, prevVal+val)-ang
    #     err = prevVal-val
    #     i += 1

    # print(val)


if __name__ == '__main__':
    main()
