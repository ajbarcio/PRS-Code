# built-in dependencies
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from os import remove
import glob

# Included Dependencies
# from spring import Spring
# from materials import Maraging300Steel

# Abaqqus interfacing dependencies
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup

## This is an example script in which I am playing around with using abaqus in
#  conjunction with this simulator

# get spring surfaces from txts
Aname      = "A_surface.txt"
Bname      = "B_surface.txt"
springName = "manual_spring"

A = np.genfromtxt("surfaces\\"+Aname, delimiter=',')
B = np.genfromtxt("surfaces\\"+Bname, delimiter=',')

# get spring parameters from txts
params = np.genfromtxt("springs\\"+springName, delimiter=',')
n = int(params[0])

# prepare spring surfaces for use in abaqus
A = A[:,0:2]
A = list(map(tuple,A))

B = B[:,0:2]
B = list(map(tuple,B))
# print(A)

# wtf does this do?
executeOnCaeStartup()

# should probably get rid of this
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)

# new model database
Mdb()
# let me select stuff using coordinates
session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)

# make the sketch
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
    sheetSize=200.0)
s.setPrimaryObject(option=STANDALONE)
# sketch both outer surfaces
s.Spline(points=A)
s.Spline(points=B)
# close the profile
s.Line(point1=A[0],point2=B[0])
s.Line(point1=A[-1],point2=B[-1])

# make the part
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
# this is probably unnecesary
p = mdb.models['Model-1'].parts['Part-1']

# make the assembly (prepare to pattern the thing)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
# this is probably still unnecesary
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Base Flexure', part=p, dependent=ON)

# make a copy in a revolved state
# this line is probably still unnecessary
a1 = mdb.models['Model-1'].rootAssembly
# this is the rotation step (IDK WHY THERE IS A POINT HERE)
a1.RadialInstancePattern(instanceList=('Base Flexure', ), point=(0.0,
    0.0, 0.0), axis=(0.0, 0.0, 1.0), number=2, totalAngle=360.0)

# very last thing you do is show it
# session.viewports['Viewport: 1'].setValues(displayedObject=a)
if False:
    for f in glob.glob("*.rpy.*"):
        remove(f)
    for f in glob.glob("*.rec"):
        remove(f)
    for f in glob.glob("*.exception"):
        remove(f)