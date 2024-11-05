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

# import os
# import json

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
testTorque = springData.loc[sizeName,'Max Torque (in.lbs)']

# Define these parameters first, hopefully variable names are clear
numberOfArms                  = 2
totalSweptAngle               = 155
beginningAndEndingAlphaAngles = np.array([20,0])*deg2rad

# Define these parameters for the thickness profile
# Currently, a piecewise quadratic polynomial defines the second moment 
# of area (I_c) which passes through the valeus of IcSetpoitns at certain points.
# The first value in IcSetpoints is the root of the spring
# The last  value in IcSetpoints is the tip  of the spring
# Any middle values occur at the proportions of the spring's arc length outlined
# in IcArcLens
outOfPlaneThickness           = .375
IcSetpoints                   = np.array([.00186, .00045, .000096, .00096, .00156])
IcArcLens                     = np.array([.2,.5,.8])

# These values are then calculated to account for the beginning and ending 
# alpha angles, thicknesses, and enforce the form factor constraints outlined 
# in the .ods
rootThk  = np.cbrt(12*IcSetpoints[0]/outOfPlaneThickness)
tipThk   = np.cbrt(12*IcSetpoints[-1]/outOfPlaneThickness)
thks = np.array(rootThk,tipThk)
offsets = thks*np.sin(beginningAndEndingAlphaAngles)/2
offsets[-1] = -offsets[-1]

# Define these parameters for the path
# radiiValues is a list of radii checkpoints through which the neutral surface 
# must pass
# betaAngleValues are the angles from the positive x axis at which these radii 
# checkpoints will be enforced
# radiiArcLens are the proportions of the springs arc length at which each 
# intermediate radius/angle checkpoint will be enforced
radiiValues = np.array([IR+offsets[0],(IR+OR)/2*.8,(IR+OR)/2*1.1,OR+offsets[1]])
betaAngleValues = np.array([0,totalSweptAngle*.333,totalSweptAngle*.666,totalSweptAngle])*deg2rad
radiiArcLens = np.array([0.4,0.6])


def defineSpring():
    # Define a spring based on user-supplied parameters
    pathDef = PATHDEF.RadiallyEndedPolynomial(n=numberOfArms, arcLen=6,
                radii = radiiValues,
                ffradii = np.array([IR, OR]),
                alphaAngles = beginningAndEndingAlphaAngles,
                betaAngles = betaAngleValues,
                XYFactors = radiiArcLens)
    crscDef = CRSCDEF.Piecewise_Ic_Control(path=pathDef,
                                        t = outOfPlaneThickness,
                                        IcPts = IcSetpoints,
                                        IcParamLens = IcArcLens)
    # Give it a material
    materialDef = materials.Maraging300Steel
    # Format export name properly:
    exportName = sizeName+" "+materialDef.name
    exportName = re.sub(' ', '_', exportName)
    # Initialize spring
    manualSpring = Spring(pathDef, crscDef, materialDef, resolution=500, torqueCapacity=testTorque,
                        name=exportName)
    return manualSpring

def deformSpring(spring):
    # Deform in a (hopefully) robust way
    res, SF, divergeFlag, i = spring.deformMode(testTorque,spring.deform_ODE)
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

def showResults(spring):
    spring.plot_deform()

def exportResults(spring):
    spring.export_surfaces()
    spring.export_parameters()

def main():
    thisSpring = defineSpring()
    deformSpring(thisSpring)
    showResults(thisSpring)
    exportResults(thisSpring)

if __name__ == "__main__":
    main()