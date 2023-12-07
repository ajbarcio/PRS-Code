import time
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


def get_derivative(vec1, vec2, s):
    dv1dv2 = np.empty(len(vec1))
    if i == 0:
        dv1dv2[i] = 0
    elif i == 1:
        dv1dv2[i-1] = [(vec1(i)-vec1(i-1))/(vec2(i)-vec2(i-1))]
    elif i == len(vec1)-1:
        dv1dv2[i] = [(vec1(i)-vec1(i-1))/(vec2(i)-vec2(i-1))]
    else:
        dv1dv2[i-1] = [(vec1(i)-vec1(i-2))/(vec2(i)-vec2(i-2))]

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

    def deform_angle(self, angle):
        dBetaMax = angle
        rn = self.rn
        alpha = self.tan
        E = self.material.E
        rc = self.rc
        dPhi = dBetaMax
        initial_conds = [0,dPhi]
        def bcs(gamma0, gammaL, p):
            F = p[0]
            return np.array([gamma0, gammaL-dPhi]), 
        def func(s, gamma, p):
            F = p[0]
            largeCouple = (self.outPlaneThickness)*self.thks*E*(self.rc-self.rn)*rn
            dlCds = get_derivative(largeCouple, self.s, s)
            
            return np.vstack(gamma[1],)
        sol = integrate.solve_bvp(func, bcs, self.s)
        ODE = lambda s, GAMMA:  \
                np.dot(np.array[[0,1],[0,1]],GAMMA)
        

    
    def deform_torque(self, torqueIn):
        stuff = things

    def get_max_stress(self):
        stuff = things

    def get_stiffness(self):
        stuff = things
    
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
                dxydt = dU.dot(coeffs)
                dydx = dxydt[1]/dxydt[0]
                self.tan = [dxydt[0],dxydt[1]]
                # dydx = (cart[i+1,1]-cart[i-1,1])/(cart[i+1,0]-cart[i-1,0])
                self.norm[i] = [-self.tan[1],self.tan[0]]/lin.norm([-self.tan[1],self.tan[0]])
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
    
    def generate_thickness_profile(self):
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
        nanerrorflag = 0
        for i in range(len(u)):
            if abs(self.rc[i]) < .5*self.thks[i]:
            #     self.rotorThickness = self.rotorThickness*0.99
            #     self.generate_thickness_profile()
               nanerrorflag = 1
            self.ttop[i] = self.xyc[i]+self.norm[i]*0.5*self.thks[i]
            self.tbottom[i] = self.xyc[i]-self.norm[i]*0.5*self.thks[i]
        if nanerrorflag == 0:
            print("nanerror")
        return self.thks
        
    def generate_neutral_profile(self):
        # print("neutral")
        self.curveError = False
        u = np.linspace(0,3,3001)
        self.xyn = np.empty([len(u), 2])
        self.rn  = np.empty(len(u))
        self.s   = np.empty(len(u))
        dxyndu   = np.empty([len(u), 2])
        for i in range(len(u)):
            if np.isinf(self.rc[i]):
                self.rn[i] = self.rc[i]
            else:
                self.rn[i] = np.sign(self.rc[i])*self.thks[i]/abs(np.log(abs(self.rc[i])+.5*self.thks[i])-np.log(abs(self.rc[i])-.5*self.thks[i]))
                if not np.sign(self.rc[i])==np.sign(self.rn[i]):
                    self.rn[i] = -1*self.rn[i]
            diff = self.rc[i] - self.rn[i]
            if np.isnan(diff):
                diff = 1
            if np.isnan(self.rn[i]):
                self.curveError = True
            self.xyn[i] = self.xyc[i]+self.norm[i]*diff
            if i == 0:
                dxyndu[i] = 0
            elif i == 1:
                dxyndu[i-1] = [(self.xyn[i,0]-self.xyn[i-1,0])/.001,(self.xyn[i,1]-self.xyn[i-1,1])/.001]
            elif i == 3000:
                dxyndu[i] = [(self.xyn[i,0]-self.xyn[i-1,0])/.001,(self.xyn[i,1]-self.xyn[i-1,1])/.001]
            else:
                dxyndu[i-1] = [(self.xyn[i,0]-self.xyn[i-2,0])/.002,(self.xyn[i,1]-self.xyn[i-2,1])/.002]
            if i == 0:
                self.s[i] = 0
            else:
                integrand = np.sqrt(dxyndu[0:i,0]**2+dxyndu[0:i,1]**2)*.001
                # print(integrand)
                self.s[i] = np.sum(integrand)
            self.norm[i] = self.norm[i]*diff
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

        for i in range(len(np.linspace(0,3,3001))):
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
        for i in range(len(np.linspace(0,3,3001))):
            normnorm[i] = lin.norm(self.norm[i])
        plt.plot(np.linspace(0,3,3001), abs(self.rn))
        plt.plot(np.linspace(0,3,3001), (self.rc))
        plt.plot(np.linspace(0,3,3001), -.5*(self.thks))
        plt.plot(np.linspace(0,3,3001), .5*(self.thks))
        plt.axhline(np.min(abs((self.rn))))
        # plt.plot(np.linspace(0,3,3001), normnorm)
        # plt.plot(np.linspace(0,3,3001), self.d)
        # plt.plot(np.linspace(0,3,3001), self.thks)
        #  plt.plot([-1,2])
        plt.ylim([-1.5,1.5])
        plt.plot(np.linspace(0,3,3001), self.s)
        plt.figure(3)
        plt.clf()
        # self.rc[i] = ((1+(dydx)**2)**(1.5)/((d2ydx2)))
        plt.plot(np.linspace(0,3,3001), (1+(self.d)**2)**1.5)
        plt.plot(np.linspace(0,3,3001), self.dd)
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

        [xyc, rc, norm] = startingParameters.generate_centroidal_profile()
        # sfor i in range(len(np.linspace(0,3,3001))):
            # if np.sqrt(xyc[i,0]**2+xyc[i,1]**2) > startingParameters.outerRadius:
            #     print("too-big-error-tripped")
            #     dragVector = prevDragVector
            #     [xyc, rc, norm] = startingParameters.generate_centroidal_profile()
            #     break
        thks = startingParameters.generate_thickness_profile()
        [norm, rn, curveError] = startingParameters.generate_neutral_profile()
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
        costArray = np.array([gainArray[0]*1/lin.norm(abs(rn[0:999]), ord=-np.inf),gainArray[1]*1/lin.norm(abs(rn[1000:2000]), ord=-np.inf),gainArray[2]*1/lin.norm(abs(rn[2001:3000]), ord=-np.inf),gainArray[3]*1/startingParameters.s[-1]])
        
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
    startingParameters = spring(material, rotorThickness=.85, p9angle=(180-25)*deg2rad, tangentPropStart=1.5, tangentPropEnd=1.5, p6radius=1.8, ctrlThickness=.3, minThickness=.3)

    startingParameters.generate_ctrl_points()
    startingParameters.generate_centroidal_profile()
    startingParameters.generate_thickness_profile()
    [norm, rn, curveError] = startingParameters.generate_neutral_profile()

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
    print("minimizing")
    res = op.minimize(drag_end_points, dragVector, method='Nelder-Mead', bounds=Bnds, options={'maxiter': 2000, 'adaptive': True, 'fatol': 0.01})
    print("done minimizing")
    print(res)
    
    startingParameters.generate_ctrl_points()
    [xyc, rc, norm] = startingParameters.generate_centroidal_profile()
    thks = startingParameters.generate_thickness_profile()
    [norm, rn, curveError] = startingParameters.generate_neutral_profile()
    print(startingParameters.pts)
    plt.close()
    startingParameters.plotResults(oldPts)
    plt.show()

if __name__ == '__main__':
    main()

