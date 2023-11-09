import numpy as np
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

    def generate_points(self):
        p9 = [self.outerRadius*np.cos(self.ap9), self.outerRadius*np.sin(self.ap9)]
        p6 = [self.lp6*np.cos(self.ap6), self.lp6*np.sin(self.ap6)]
        p3 = [self.lp3*np.cos(self.ap3), self.lp3*np.sin(self.ap3)]
        p0 = [self.innerRadius, 0]

        # Rules for the spline:

        # p5 = p1 + 4*(p3-p2)
        # p8 = p4 + 4*(p6-p5)




def main():

    PATH = 'C:\Stuff\SMM\Spline_thickening-v3.xlsm'

    wb = load_workbook(filename = PATH)
    ws = wb.active

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
        
    material = MaraginSteelC300()
    startingParameters = parameters(material)

    u = np.linspace(0,3,3001)
    x = 


if __name__ == '__main__':
    main()

