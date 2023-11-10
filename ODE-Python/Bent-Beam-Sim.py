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

        # angle(P8) = angle(P9)

        # knowns: 0, 1, 3, 6, 8, 9
        # unknowns: 2, 4, 5, 7

        A = np.array([[4,0,1,0],[1,1,0,0],[0,-1,4,0],[0,0,1,1]])
        B = np.array([4*p3+p1,2*p3,4*p6-p8,2*p6])

        result = lin.solve(A,B)
        pts = np.array([p0,p1,result[0],p3,result[1],result[2],p6,result[3],p8,p9])
        return pts

def main():

    PATH = 'Spline_thickening-v3.xlsm'

    wb = load_workbook(filename = PATH)
    ws = wb.active
    '''
    pts = np.array([[ws['C2'].value, ws['D2'].value], \
                    [ws['C3'].value, ws['D3'].value], \
                    [ws['C4'].value, ws['D4'].value], \
                    [ws['C5'].value, ws['D5'].value], \
                    [ws['C6'].value, ws['D6'].value], \
                    [ws['C7'].value, ws['D7'].value], \
                    [ws['C8'].value, ws['D8'].value], \
                    [ws['C9'].value, ws['D9'].value], \
                    [ws['C10'].value, ws['D10'].value], \
                    [ws['C11'].value, ws['D11'].value], \
                    ])
    # print(pts)
    '''        
    material = MaraginSteelC300()
    startingParameters = parameters(material)

    u = np.linspace(0,3,3001)
    print(startingParameters.generate_ctrl_points())


if __name__ == '__main__':
    main()

