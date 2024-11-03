import numpy as np

from modules.utils import deg2rad
import modules.PATHDEF as PATHDEF
import modules.CRSCDEF as CRSCDEF
import modules.materials as materials
from modules.spring import Spring, determineFastestSolver

IR = 1.3
OR = 3.44
testTorque = 4549

# Define these parameters first, hopefully variable names are clear
numberOfArms                  = 3
totalSweptAngle               = 110
beginningAndEndingAlphaAngles = np.array([35,0])*deg2rad

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
radiiValues = np.array([IR+offsets[0],(IR+OR)/2*1.3,(IR+OR)/2*.75,OR+offsets[1]])
betaAngleValues = np.array([0,totalSweptAngle*.5,totalSweptAngle*.9,totalSweptAngle])*deg2rad
radiiArcLens = np.array([0.5,0.75])

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
manualSpring = Spring(pathDef, crscDef, materialDef, resolution=500, torqueCapacity=testTorque)

res, SF, divergeFlag, i = manualSpring.deformMode(testTorque,manualSpring.deform_ODE)
if divergeFlag:
    # testSpring.detailed_deform_regression(testTorque, testSpring.deform_ODE,resl=10,degree=1)
    pass

# Print total deformation
print("total deformation:", manualSpring.dBeta/deg2rad, u'\N{DEGREE SIGN}')
# Keep track of whether or not this is a reasonable deformation (if this number is small you are good)
print(manualSpring.solnerr)
if manualSpring.solnerr > 0.01:
    print("This solution may not be valid, attempting to update deformMode")
    determineFastestSolver(manualSpring)

# Return the requried boundary conditions
print(SF)

manualSpring.plot_deform()