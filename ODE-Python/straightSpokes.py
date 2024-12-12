import numpy as np
from scipy import optimize as opt
import re

import pandas as pd

import modules.materials as materials
import modules.PATHDEF as PATHDEF
import modules.CRSCDEF as CRSCDEF
from modules.spring import Spring, determineFastestSolver
from modules.utils import deg2rad

import argparse

parser = argparse.ArgumentParser(description='decide which form factor to use')
parser.add_argument('-s', type=str, help='What form factor do you want to use? \
                                          Make this match any row title in  \
                                          Spring_Constraints.ods')
args = parser.parse_args()
sizeName = args.s

print("trying to make a spring to fit:", sizeName)
# Standard procedure: define a path, then define a thickness profile using the
#                     length of that path, then pass both of those objects into
#                     the overall spring object that contains deformation
#                     methods

springData = pd.read_excel('Spring_Constraints.ods', engine='odf', index_col=0)

IR = springData.loc[sizeName,'IR lim (in)']
OR = springData.loc[sizeName,'OR lim (in)']
testTorque = springData.loc[sizeName,'Max Torque (in.lbs)']/2

# Define these parameters first, hopefully variable names are clear
numberOfArms        = 3
outOfPlaneThickness = 0.375
inPlaneThickness    = 0.125

def defineSpring():
    # Define a spring based on user-supplied parameters
    pathDef = PATHDEF.RadiallyEndedPolynomial(n=numberOfArms, arcLen=OR-IR,
                radii = np.array([IR, OR]),
                ffradii = np.array([IR, OR]),
                alphaAngles = np.array([0,0]),
                betaAngles = np.array([0,0]),
                XYFactors = np.array([]))
    crscDef = CRSCDEF.Constant_Ic(path=pathDef, t = outOfPlaneThickness, h0 = inPlaneThickness)
    # Give it a material
    materialDef = materials.Maraging300Steel
    # Format export name properly:
    exportName = sizeName+" "+materialDef.name+" straight"
    exportName = re.sub(' ', '_', exportName)
    # Initialize spring
    manualSpring = Spring(pathDef, crscDef, materialDef, resolution=200, torqueCapacity=testTorque,
                        name=exportName)
    return manualSpring

def deformSpring(spring: Spring, ODE = None):
    # Deform in a (hopefully) robust way
    if ODE is None:
        ODE = spring.deform_ODE
    res, SF, divergeFlag, i = spring.deformMode(testTorque,ODE)
    if divergeFlag:
        # testSpring.detailed_deform_regression(testTorque, testSpring.deform_ODE,resl=10,degree=1)
        pass

    # Print total deformation
    print("total deformation:", spring.dBeta/deg2rad, u'\N{DEGREE SIGN}')
    # Keep track of whether or not this is a reasonable deformation (if this number is small you are good)
    print(spring.solnerr)
    if spring.solnerr > 0.01:
        print("This solution may not be valid, attempting to update deformMode")
        determineFastestSolver(spring)

    # Return the requried boundary conditions
    print(SF)

    return res, SF, divergeFlag, i

def showResults(spring: Spring):
    spring.plot_deform()

def exportResults(spring: Spring):
    spring.export_surfaces()
    spring.export_parameters()

def main():
    thisSpring = defineSpring()
    deformSpring(thisSpring, ODE = thisSpring.deform_withTension_ODE)
    showResults(thisSpring)
    exportResults(thisSpring)

if __name__ == "__main__":
    main()