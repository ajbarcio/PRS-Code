import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

import modules.PATHDEF as PATHDEF
import modules.CRSCDEF as CRSCDEF
from modules.utils import deg2rad
import modules.materials as materials
from modules.interactive import Interactive_Spring


parser = argparse.ArgumentParser(description='decide which form factor to use')
parser.add_argument('-s', type=str, help='What form factor do you want to use? \
                                          Make this match any row title in  \
                                          Spring_Constraints.ods')
parser.add_argument('-n', type=int, help='If you want to use a nonstandard  \
                                          number of spring legs (n !=2),  \
                                          specify this')
# parser.add_argument()
args = parser.parse_args()
sizeName = args.s
numLegs  = args.n

if sizeName is not None:
    print("trying to make a spring to fit:", sizeName)
    springData = pd.read_excel('Spring_Constraints.ods', engine='odf', index_col=0)
    IR = springData.loc[sizeName,'IR lim (in)']
    OR = springData.loc[sizeName,'OR lim (in)']
    testTorque = springData.loc[sizeName,'Max Torque (in.lbs)']
else:
    IR = input("Input inner radius of spring: ")
    OR = input("Input outer radius of spring: ")
    testTorque = input("How hard to torque the spring? (Nm): ")

    IR = float(IR)
    OR = float(OR)
    testTorque = float(testTorque)

# Define these parameters first, hopefully variable names are clear
if numLegs is not None:
    numberOfArms = numLegs
else:
    numberOfArms                  = 2
totalSweptAngle               = 360/numberOfArms*11/12
beginningAndEndingAlphaAngles = np.array([0,0])*deg2rad

# Define these parameters for the thickness profile
# Currently, a piecewise quadratic polynomial defines the second moment 
# of area (I_c) which passes through the valeus of IcSetpoitns at certain points.
# The first value in IcSetpoints is the root of the spring
# The last  value in IcSetpoints is the tip  of the spring
# Any middle values occur at the proportions of the spring's arc length outlined
# in IcArcLens
outOfPlaneThickness           = .375
IcSetpoints                   = np.array([.00186, .00045, .000096, .00096, .00156])*2
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
radiiValues = np.array([IR+offsets[0],(IR+OR)/2*1.1,(IR+OR)/2,OR+offsets[1]])
betaAngleValues = np.array([0,totalSweptAngle*.333,totalSweptAngle*.666,totalSweptAngle])*deg2rad
radiiArcLens = np.array([0.333,0.6666])

path = PATHDEF.RadiallyEndedPolynomial(n = numberOfArms, arcLen=6,
                                       radii = radiiValues,
                                       ffradii = np.array([IR, OR]),
                                       alphaAngles = beginningAndEndingAlphaAngles,
                                       betaAngles = betaAngleValues,
                                       XYFactors = radiiArcLens)
crsc = CRSCDEF.Piecewise_Ic_Control(path = path,
                                    t = outOfPlaneThickness,
                                    IcPts = IcSetpoints,
                                    IcParamLens = IcArcLens)

radius = path.outerRadius

# plot = DraggableSpring(x, y, np.linspace(0,path.arcLen,500), [1, 2], [3], semicircle)
plot = Interactive_Spring(path, crsc, materials.Maraging300Steel, torqueCapacity=testTorque, name="Interactive")
plot.ax.set_xlim(-radius*1.1,radius*1.1)
plot.ax.set_ylim(-radius*1.1,radius*1.1)
plot.ax.set_aspect('equal')
plt.show()