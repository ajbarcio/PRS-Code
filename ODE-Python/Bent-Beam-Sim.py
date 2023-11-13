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
            # print(u)
            u2 = u - np.floor(u)
            cart = np.empty([len(u), 2])
    
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
                    cart[i] = U.dot(D).dot(P)
                else:
                    U = np.array([t**3, t**2, t, 1])
                    D = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
                    P = np.array([pts[3*ptndx],pts[3*ptndx+1],pts[3*ptndx+2],pts[3*ptndx+3]])
                    cart[i] = U.dot(D).dot(P)

            return cart
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
                        [0,0,0,0,0,1,0],[3**5,3**4,3**3,3**2,3,1,0]])
        Targ = np.array([[self.rotorThickness],[self.ctrlThickness], \
                         [self.minThickness],[self.ctrlThickness], \
                         [self.rotorThickness],[0]])


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

    pts = startingParameters.generate_ctrl_points()
    cartc  = startingParameters.generate_centroidal_profile(pts)
    rc     = startingParameters.generate_centroidal_radius(cartc)

    print(rc)


if __name__ == '__main__':
    main()

