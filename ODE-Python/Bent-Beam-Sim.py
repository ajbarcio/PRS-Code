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

class parameters(object):
    def __init__(self, material):

        ## PARAMETERS TO ITERATE

        # Thickness parameters

        self.rotorThickness = .8
        self.ctrlThickness = 0.27
        self.minThickness = .23

        # Spline thickness parameters

        self.lp3 = 1.9
        self.lp6 = 1.8
        self.ap3 = deg2rad*45
        self.ap6 = deg2rad*110

        ## INVARIANT PARAMETERS

        self.material = material
        self.innerRadius = 1.62/2
        self.outerRadius = 5/2
        self.ap9 = (180-10)*deg2rad

# def ODE():

    def generate_ctrl_points(self):
        # Points predefined by parameters and arbitrary constraints:

        # p9 on outer radius at angle
        # p8 radial with p9, 1 'rotorThickness' away (maybe a parameter gets introduced to futz with this)
        # p6 at angle and distance parameters
        # p3 at angle and distance parameters
        # p0 on inner radius at horizontal 
        # p1 radial with p0, 1 'rotorThickness' away

        p9 = np.array([self.outerRadius*np.cos(self.ap9), self.outerRadius*np.sin(self.ap9)])
        p8 = np.array([p9[0]-self.rotorThickness*np.cos(self.ap9),p9[1]-self.rotorThickness*np.sin(self.ap9) ])
        p6 = np.array([self.lp6*np.cos(self.ap6), self.lp6*np.sin(self.ap6)])
        p3 = np.array([self.lp3*np.cos(self.ap3), self.lp3*np.sin(self.ap3)])
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
        pts = np.array([p0,p1,result[0],p3,result[1],result[2],p6,result[3],p8,p9])
        return pts
    
    def generate_centroidal_profile(self, pts):

            u = np.linspace(0,3,3001)
            offset = 0
            # print(u)
            u2 = u - np.floor(u)
            cart = np.empty([len(u), 2])
            th   = np.empty(len(u))
            norm = np.empty([len(u), 2])
            ck = np.empty([len(u), 2])
            rc   = np.empty(len(u))
    
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
                    cart[i] = U.dot(D).dot(P)-[0,offset]
                    th[i] = np.arctan2(cart[i,1],cart[i,0])
                else:
                    U = np.array([t**3, t**2, t, 1])
                    D = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
                    P = np.array([pts[3*ptndx],pts[3*ptndx+1],pts[3*ptndx+2],pts[3*ptndx+3]])
                    cart[i] = U.dot(D).dot(P)-[0,offset]
                    th[i] = np.arctan2(cart[i,1],cart[i,0])
                dU = np.array([3*t**2, 2*t, 1, 0])
                ddU = np.array([6*t, 2, 0, 0])
                coeffs = D.dot(P)
                if i == 0: 
                    # tan  = [cart[i+1,0]-cart[i,0],cart[i+1,1]-cart[i,1]]/np.sqrt((cart[i+1,0]-cart[i,0])**2+(cart[i+1,1]-cart[i,1])**2)
                    dxydt = dU.dot(coeffs)
                    dydx = dxydt[1]/dxydt[0]
                    tan = [dxydt[0],dxydt[0]]
                    #dydx = (cart[i+1,1]-cart[i,1])/(cart[i+1,0]-cart[i,0])
                    # dydxprev = 0
                elif i == 3000:
                    # dydxprev = dydx
                    # tan  = [cart[i,0]-cart[i-1,0],cart[i,1]-cart[i-1,1]]/np.sqrt((cart[i,0]-cart[i-1,0])**2+(cart[i,1]-cart[i-1,1])**2)
                    dxydt = dU.dot(coeffs)
                    dydx = dxydt[1]/dxydt[0]
                    tan = [dxydt[0],dxydt[0]]
                    #dydx = (cart[i,1]-cart[i-1,1])/(cart[i,0]-cart[i-1,0])
                else:
                    # dydxprev = dydx
                    # tan  = [cart[i+1,0]-cart[i-1,0],cart[i+1,1]-cart[i-1,1]]/np.sqrt((cart[i+1,0]-cart[i-1,0])**2+(cart[i+1,1]-cart[i-1,1])**2)
                    dxydt = dU.dot(coeffs)
                    dydx = dxydt[1]/dxydt[0]
                    tan = [dxydt[0],dxydt[0]]
                    # dydx = (cart[i+1,1]-cart[i-1,1])/(cart[i+1,0]-cart[i-1,0])
                norm[i] = [-tan[1],tan[0]]
                d2xydt2 = ddU.dot(coeffs)
                d2ydx2 = (d2xydt2[1] - dydx*d2xydt2[0])/(dxydt[0])
                # d2ydx2 = (dydx-dydxprev)/(cart[i,0]-cart[i-1,0])
                rc[i] = abs((1+dydx**2)**(1.5)/(d2ydx2))
                ck[i] = cart[i]-norm[i]*rc[i]

            return cart, norm, rc
        #print(cart)

    def generate_centroidal_radius(self, cartc):
        
        rc = np.empty([len(cartc)])
        for i in range(len(rc)):
            rc[i] = np.sqrt(cartc[i,0]**2+cartc[i,1]**2)
        return rc
    
    def generate_thickness_profile(self):

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
        t = coeffs[0]*u**6 + coeffs[1]*u**5 + coeffs[2]*u**4 + coeffs[3]*u**3 + coeffs[4]*u**2 + coeffs[5]*u + coeffs[6]
        # coeffs = lin.inv(REG.T.dot(REG)).dot(REG.T).dot(Targ)
        return t
    
    def generate_neutral_radius(self, rc, t):
        rn = t/np.log((rc+.5*t)/(rc-.5*t))
        return rn
    
    def generate_neutral_profile(self, rc, thks, norm, xyc):
        u = np.linspace(0,3,3001)
        cart = np.empty([len(u), 2])
        for i in range(len(u)):
            if np.isinf(rc[i]):
                rn = rc[i]
            else:
                rn = thks[i]/np.log((rc[i]+.5*thks[i])/(rc[i]-.5*thks[i]))
            diff = rc[i] - rn
            if np.isnan(diff):
                b = rc[i]+.5*thks[i]
                a = rc[i]-.5*thks[i]
                rn = (b-a)/np.log(b/a)
            if i%50 == 0:
                print(diff, rc[i], rn, rc[i]-.5*thks[i], .5*thks[i])
            cart[i] = xyc[i]+norm[i]*diff
            norm[i] = norm[i]*diff
        return cart, norm
    
    def generate_neutral_profile2(self, rn, th):
        u = np.linspace(0,3,3001)
        cart = np.empty([len(u), 2])
        for i in range(len(u)):
            cart[i,0] = rn[i]*np.cos(th[i])
            cart[i,1] = rn[i]*np.sin(th[i])
        return cart 



def main():

    #PATH = 'Spline_thickening-v3.xlsm'
    
    '''
    wb = load_workbook(filename = PATH)
    ws = wb.active
    pts = np.array([[ws['C2'].value,  ws['D2'].value], \
                    [ws['C3'].value,  ws['D3'].value], \
                    [ws['C4'].value,  ws['D4'].value], \
                    [ws['C5'].value,  ws['D5'].value], \
                    [ws['C6'].value,  ws['D6'].value], \
                    [ws['C7'].value,  ws['D7'].value], \
                    [ws['C8'].value,  ws['D8'].value], \
                    [ws['C9'].value,  ws['D9'].value], \
                    [ws['C10'].value, ws['D10'].value], \
                    [ws['C11'].value, ws['D11'].value], \
                    ])
    # print(pts)
    '''        
    material = MaraginSteelC300()
    startingParameters = parameters(material)

    pts  = startingParameters.generate_ctrl_points()
    [xyc, norm, rc]  = startingParameters.generate_centroidal_profile(pts)
    thks = startingParameters.generate_thickness_profile()
    [xyn, norm2] = startingParameters.generate_neutral_profile(rc, thks, norm, xyc)
    print(rc)
    # """ rc   = startingParameters.generate_centroidal_radius(xyc)
    # thks = startingParameters.generate_thickness_profile()
    # rn = startingParameters.generate_neutral_radius(rc, thks)
    # xyn  = startingParameters.generate_neutral_profile(rc, thks, xyc)
    # xyn2 = startingParameters.generate_neutral_profile2(rn, th) """

    # u = np.linspace(0,3,3001)
    # test = np.empty([len(u), 2])
    # for i in range(len(u)):
    #     test[i,0] = rc[i]*np.cos(u[i]*startingParameters.ap9/u[-1])
    #     test[i,1] = rc[i]*np.sin(u[i]*startingParameters.ap9/u[-1])

    #plt.plot(np.linspace(0,3,3001), rc)
    #plt.plot(np.linspace(0,3,3001), rn)
    # plt.figure(1,label="xyc")
    # plt.plot(xyc[:,0], xyc[:,1])
    # plt.figure(2,label="xyn")
    # plt.plot(xyn[:,0], xyn[:,1])
    # # plt.figure(3,label="test")
    # plt.plot(xyn[:,0], xyn[:,1])
    
    # plt.plot(xyn2[:,0], xyn[:,1])
    # plt.plot(0,0)

    # plt.plot(xyc[:,0],xyc[:,1])
    # plt.plot(ck[0:1001,1],ck[0:1001,1])
    # plt.figure(1,label="rc")
    # plt.plot(rc)
    # plt.plot(xyc[:,1])
    # plt.figure(2,label="xyc")
    plt.plot(xyc[:,0], xyc[:,1])
    plt.plot(xyn[:,0], xyn[:,1])
    for i in range(len(np.linspace(0,3,3001))):
        if i%25 == 0:
            #plt.arrow(xyc[i,0], xyc[i,1], norm2[i,0], norm2[i,1])
            plt.arrow(xyc[i,0], xyc[i,1], norm2[i,0], norm2[i,1])

    # plt.plot(ck[:,0], ck[:,1])

    plt.show()

    # print(rn)
    # print(rc)
    print("xyc:",xyc)
    # print("angs",angs)
    print("normdir",norm)
    # print("center",ck)
    print("radius",rc)




if __name__ == '__main__':
    main()

