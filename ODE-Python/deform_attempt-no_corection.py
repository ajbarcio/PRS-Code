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

globalRes = 10
globalLen = 3*globalRes+1
globalMaxIndex = 3*globalRes
globalBCresLimit = 0.01
globalRelresLimit = .1

def get_derivative(vec1, vec2):
    dv1dv2 = np.empty(len(vec1))
    for i in range(len(dv1dv2)-1):
        if i == 0:
            dv1dv2[i] = 1
        elif i == 1:
            dv1dv2[i-1] = (vec1[i]-vec1[i-1])/(vec2[i]-vec2[i-1])
        elif i == len(vec1)-1:
            dv1dv2[i] = (vec1[i]-vec1[i-1])/(vec2[i]-vec2[i-1])
        else:
            dv1dv2[i-1] = (vec1[i]-vec1[i-2])/(vec2[i]-vec2[i-2])
    return dv1dv2

class spring(object):

    def __init__(self, material, **kwargs):
        ## PARAMETERS TO ITERATE:
            # these are default values, which can be changed by passing named arguments to the class

        self.Vparams = {'rotorThickness': 0.88,             \
                        'ctrlThickness':  0.27,             \
                        'minThickness':   0.23,             \
                        'p3radius':       2.13,             \
                        'p6radius':       2,                \
                        'p3angle':        deg2rad*45,       \
                        'p6angle':        deg2rad*110,      \
                        'p9angle':        (180)*deg2rad,    \
                        'tangentPropStart':    1,           \
                        'tangentPropEnd':      1 }

        ## INVARIANT PARAMETERS

        self.innerRadius = 1.62/2
        self.outerRadius = 5/2
        
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

    def deform_force(self, force):
        
        # ODE System
        # Parametrized version (no arc length)
        # g'  = g'
        # g'' = f(u)g' + h(u, g)

        def defomration_function(u, gamma, p):
            Fx = self.Fx
            Fy = self.Fy
            BetaMax = p[0]

            leng = len(u)

            self.generate_ctrl_points()
            self.generate_centroidal_profile(leng, mesh=u)
            self.generate_thickness_profile(leng, mesh=u)
            self.generate_neutral_profile(leng, mesh=u)

            largeCouple = np.transpose(self.outPlaneThickness*self.thks*self.material.E*(self.rc-self.rn)*self.rn)
            dlCdu       = np.transpose(get_derivative(largeCouple, np.linspace(0,3,leng)))
            dxdu        = np.transpose(self.dxyndu[:,0])
            dydu        = np.transpose(self.dxyndu[:,1])
            d2xdu2      = np.transpose(get_derivative(dxdu, np.linspace(0,3,leng)))
            d2ydu2      = np.transpose(get_derivative(dydu, np.linspace(0,3,leng)))

            # print(len(dlCdu), len(dxdu**2))

            # if np.isnan(dxdu.any):
            #     print("AAAAAA")
            # if np.isnan(dxdu.any):
            #     print("AAAAAA")    

            # u2 = u
            # u2 = u2.astype(int)
            #  print(np.sqrt(np.square(dxdu)+np.square(dydu)))
            # print(dxdu[-2], dydu[-2])

            sqrt_replace = np.sqrt(np.square(dxdu)+np.square(dydu))

            if 0 in sqrt_replace:
                # print("replaced")
                sqrt_replace = np.where(sqrt_replace != 0, sqrt_replace, .001)
                # print(sqrt_replace)
            if float('inf') in sqrt_replace:
                # print(sqrt_replace)
                sqrt_replace = np.where(sqrt_replace != float('inf'), sqrt_replace, 1.7976931348623158*10**308)


            
            # print("denom size",sqrt_replace.size)
            # print("num size",dlCdu.size)
            # print("div size",(dlCdu/sqrt_replace).size)
            # print(u)
           # print(sqrt_replace)
            # print("STOP")

            return np.vstack((gamma[1], \
                             (Fx*np.sin(self.tann+gamma[0])-Fy*np.cos(self.tann+gamma[0]))*np.sqrt(np.square(dxdu)+np.square(dydu))/largeCouple \
                              +gamma[1]*((dxdu*d2xdu2+dydu*d2ydu2)-dlCdu/sqrt_replace), \
                             np.cos(self.tann+gamma[0]), \
                             np.sin(self.tann+gamma[0])))

        def BC_function(left, right, p):
            BetaMax = p[0]
            Rn = self.rn[-1]
            return np.array([left[0], right[0]-BetaMax, left[2], right[2]-Rn*np.cos(beta0+BetaMax), right[3]-Rn*np.sin(beta0+BetaMax)])

        def rk4_deform (u0, gamma0, du):

            #
            #  Get four sample values of the derivative.
            #
            f1 = defomration_function ( u0,            gamma0 )
            f2 = defomration_function ( u0 + du / 2.0, gamma0 + du * f1 / 2.0 )
            f3 = defomration_function ( u0 + du / 2.0, gamma0 + du * f2 / 2.0 )
            f4 = defomration_function ( u0 + du,       gamma0 + du * f3 )
            #
            #  Combine them to estimate the solution U1 at time T1 = T0 + DT.
            #
            u1 = gamma0 + du * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0

            return u1
        
        self.Fx = force*np.sin(self.p9angle)
        self.Fy = force*np.cos(self.p9angle)
        BetaMax = 5*deg2rad
        beta0 = self.p9angle
        rn = self.rn
        Phi = BetaMax
        tol = 0.0001
        u = np.linspace(0,3,globalLen)
        u0 = u
        u2 = np.linspace(0,globalMaxIndex,globalLen)
        
        gamma = np.ones((4, u.shape[0]))
        print("solving BVP")
        # Finitial = np.sqrt((1530/self.outerRadius*25.4)**2/2)
        Pinitial = [BetaMax]
        flag = 1
        iterator = 0
        self.deformation = integrate.solve_bvp(defomration_function, BC_function, u, gamma, p=Pinitial, \
                                               verbose=2, max_nodes=50000, tol=globalRelresLimit, bc_tol=globalBCresLimit)

        # while flag==1:
        #     if iterator ==0:
        #         self.deformation = integrate.solve_bvp(defomration_function, BC_function, u, gamma, p=Pinitial, \
        #                                             verbose=2, max_nodes=600, tol=globalRelresLimit, bc_tol=globalBCresLimit)
        #     else:
        #         self.deformation = integrate.solve_bvp(defomration_function, BC_function, u, gamma, p=Pinitial, \
        #                                             verbose=2, max_nodes=100000, tol=globalRelresLimit, bc_tol=globalBCresLimit)
        #     flag = self.deformation.status
        #     gamma = self.deformation.sol(u0)
        #     Pinitial = self.deformation.p
        #     iterator+=1
        print(flag)


        # error = 1
        # boundConds = np.array([0, Phi])
        # gamma0 = np.array([0,0])
        
        # gamma = np.empty(len(u),2)
        # gamma[0] = gamma0
        # while error>tol:
        #     for i in range(len(u))-1:
        #         gamma[i+1] = rk4_deform(u[i], gamma[i], .001)
        #     errorBC = gamma[-1,1]-boundConds[1]
  
    def deform_torque(self, torqueIn):
        pass

    def get_max_stress(self):
        pass

    def get_stiffness(self):
        pass
    
    def generate_ctrl_points(self):
        # print("control points")
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
        p8 = np.array([p9[0]-self.tangentPropEnd*self.rotorThickness*np.cos(self.p9angle), \
                       p9[1]-self.tangentPropEnd*self.rotorThickness*np.sin(self.p9angle) ])
        p6 = np.array([self.p6radius*np.cos(self.p6angle), self.p6radius*np.sin(self.p6angle)])
        p3 = np.array([self.p3radius*np.cos(self.p3angle), self.p3radius*np.sin(self.p3angle)])
        p0 = np.array([self.innerRadius, 0])
        p1 = np.array(p0+[self.tangentPropStart*self.rotorThickness,0])
        
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
        self.dxyndu   = np.empty([len(u), 2])
        self.tann     = np.empty(len(u))
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
                self.dxyndu[i-1] = [(self.xyn[i,0]-self.xyn[i-1,0])/.001,(self.xyn[i,1]-self.xyn[i-1,1])/.001]

                
            elif i == len(u)-1:
                self.dxyndu[i] = [(self.xyn[i,0]-self.xyn[i-1,0])/.001,(self.xyn[i,1]-self.xyn[i-1,1])/.001]
                
            else:
                self.dxyndu[i-1] = [(self.xyn[i,0]-self.xyn[i-2,0])/.002,(self.xyn[i,1]-self.xyn[i-2,1])/.002]
                
            if i == 0:
                self.s[i] = 0
            else:
                integrand = np.sqrt(self.dxyndu[0:i,0]**2+self.dxyndu[0:i,1]**2)*.001
                # print(integrand)
                self.s[i] = np.sum(integrand)
            self.norm[i] = self.norm[i]*diff
            self.tann[i] = np.arctan2(self.dxyndu[i,1],self.dxyndu[i,0])
        # if self.curveError == True:
        #     print(self.tangentPropEnd)
        #     print(self.tangentPropStart)
        return self.norm, self.rn, self.curveError
          
    def plotResults(self, oldPts):
        print("plot")
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
        plt.figure(3)
        plt.clf()
        # self.rc[i] = ((1+(dydx)**2)**(1.5)/((d2ydx2)))
        plt.plot(np.linspace(0,3,globalLen), (1+(self.d)**2)**1.5)
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

    initialRotorThickness   = .85
    initialCtrlThickness    = .3
    initialMinThickness     = .3
    initialP3radius         = 2
    initialP6radius         = 1.466089122536943
    initialP3angle          = 0.6109352064805111
    initialP6angle          = 1.745535503850574
    initialP9angle          = 2.705260340591211
    initialTangentPropStart = 1.1491492041951368
    initialTangentPropEnd   = 1.2718246231042096

    startingParameters = spring(material, rotorThickness=initialRotorThickness, ctrlThickness=initialCtrlThickness, minThickness=initialMinThickness, \
                                p3radius=initialP3radius, p6radius=initialP6radius, p3angle=initialP3angle, p6angle=initialP6angle, \
                                p9angle=initialP9angle, tangentPropStart=initialTangentPropStart, tangentPropEnd=initialTangentPropEnd)

    startingParameters.generate_ctrl_points()
    startingParameters.generate_centroidal_profile(globalLen)
    startingParameters.generate_thickness_profile(globalLen)
    [norm, rn, curveError] = startingParameters.generate_neutral_profile(globalLen)

    pts = startingParameters.pts
    print(startingParameters.pts)
    print(startingParameters.s[-1])
    oldPts = pts
    # time.sleep(5)

    # Bnds = ((np.sqrt(pts[1,0]**2+pts[1,1]**2),None), (None,None), (None,None), (None,None), (None,None), (None,np.sqrt(pts[8,0]**2+pts[8,1]**2)))
    startingP3angle = startingParameters.p3angle
    startingP6angle = startingParameters.p6angle
    Bnds = ((1,2),(1,2), \
            (startingParameters.innerRadius,startingParameters.outerRadius-0.5),(startingParameters.innerRadius,startingParameters.outerRadius-0.5), \
            (startingP3angle-deg2rad*10,startingP3angle+deg2rad*10),(startingP6angle-deg2rad*10,startingP6angle+deg2rad*10))
    dragVector = np.array([startingParameters.tangentPropStart,startingParameters.tangentPropEnd, \
                           startingParameters.p3radius,startingParameters.p6radius, \
                           startingParameters.p3angle,startingParameters.p6angle])
    # print("minimizing")
    # res = op.minimize(drag_end_points, dragVector, method='Nelder-Mead', bounds=Bnds, options={'maxiter': 2000, 'adaptive': True, 'fatol': 0.01})
    # print("done minimizing")
    # print(res)
    
    startingParameters.generate_ctrl_points()
    [xyc, rc, norm] = startingParameters.generate_centroidal_profile(globalLen)
    thks = startingParameters.generate_thickness_profile(globalLen)
    [norm, rn, curveError] = startingParameters.generate_neutral_profile(globalLen)
    print(startingParameters.pts)
    plt.close()
    startingParameters.plotResults(oldPts)
    plt.show()

    startingParameters.print_parameters()

    startingParameters.deform_force(10000)


if __name__ == '__main__':
    main()

