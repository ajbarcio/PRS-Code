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
Aname      = "default_A_surface.txt"
Bname      = "default_B_surface.txt"
springName = "default_spring"

A = np.genfromtxt("..\\surfaces\\"+Aname, delimiter=',')
B = np.genfromtxt("..\\surfaces\\"+Bname, delimiter=',')

# get spring parameters from txts
params = np.genfromtxt("..\\springs\\"+springName, delimiter=',')
n = int(params[0])
thk = float(params[2])
ang = float(params[10])

A3D = list(map(tuple,A))
B3D = list(map(tuple,B))
# prepare spring surfaces for use in abaqus
A = A[:,0:2]
A = list(map(tuple,A))

B = B[:,0:2]
B = list(map(tuple,B))
# print(A)

# wtf does this do?
executeOnCaeStartup()

# should probably get rid of this
# session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
#     referenceRepresentation=ON)

# new model database
Mdb()
# let me select stuff using coordinates
session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)

# Prepare the material
mdb.models['Model-1'].Material(name='Maraging300')
mdb.models['Model-1'].materials['Maraging300'].Elastic(table=((27500000.0, 
    0.33), ))
mdb.models['Model-1'].materials['Maraging300'].Plastic(
    scaleStress=None, table=((309700.0, 0.0), ))

# sketch the inner circle
inner_radius = max(lin.norm(A[0]),lin.norm(B[0]))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
     sheetSize=8.0)
s1.setPrimaryObject(option=STANDALONE)
s1.CircleByCenterPerimeter(center=(0.0,0.0), point1=(inner_radius,0.0))
# Make a rigid inner body
p1 = mdb.models['Model-1'].Part(name='Inner-Boolean', dimensionality=THREE_D,
     type=DEFORMABLE_BODY)
p1 = mdb.models['Model-1'].parts['Inner-Boolean']
p1.BaseSolidExtrude(sketch=s1, depth=thk)
# Make a rigid inner surface
p2 = mdb.models['Model-1'].Part(name='Inner-Surface', dimensionality=THREE_D,
     type=DISCRETE_RIGID_SURFACE)
p2 = mdb.models['Model-1'].parts['Inner-Surface']
p2.BaseShellExtrude(sketch=s1, depth=thk)
s1.unsetPrimaryObject()

# Reference node for rigid inner surface
p2.ReferencePoint(point=(0.0, 0.0, 0.0))
p2.features.changeKey(fromName='RP', toName='Inner-Surface-Pt')

# sketch the outer circle
point = ((A[-1][0]+B[-1][0])/2, (A[-1][1]+B[-1][1])/2)
outer_radius = lin.norm(point)
s2 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
     sheetSize=8.0)
s2.setPrimaryObject(option=STANDALONE)
s2.CircleByCenterPerimeter(center=(0.0,0.0), point1=(outer_radius,0.0))
# Make a rigid outer surface
p2 = mdb.models['Model-1'].Part(name='Outer-Surface', dimensionality=THREE_D,
     type=DISCRETE_RIGID_SURFACE)
p2 = mdb.models['Model-1'].parts['Outer-Surface']
p2.BaseShellExtrude(sketch=s2, depth=thk)
s2.unsetPrimaryObject()
# Sketch again
s3 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
     sheetSize=8.0)
s3.setPrimaryObject(option=STANDALONE)
s3.CircleByCenterPerimeter(center=(0.0,0.0), point1=(outer_radius,0.0))
s3.CircleByCenterPerimeter(center=(0.0,0.0), point1=(outer_radius+0.25,0.0))
# Make an outer body
p1 = mdb.models['Model-1'].Part(name='Outer-Boolean', dimensionality=THREE_D,
     type=DEFORMABLE_BODY)
p1 = mdb.models['Model-1'].parts['Outer-Boolean']
p1.BaseSolidExtrude(sketch=s3, depth=thk)
s3.unsetPrimaryObject()

# Reference node for rigid outer surface
p2.ReferencePoint(point=(0.0, 0.0, thk/2))
p2.features.changeKey(fromName='RP', toName='Moment-Pt')

# sketch the selected geometry
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
    sheetSize=8.0)
s.setPrimaryObject(option=STANDALONE)
# sketch both outer surfaces
s.Spline(points=A)
s.Spline(points=B)
# close the profile
s.Line(point1=A[0],point2=B[0])
s.Line(point1=A[-1],point2=B[-1])

# make the part that deforms
p = mdb.models['Model-1'].Part(name='Elastic_Portion', dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Elastic_Portion']
p.BaseSolidExtrude(sketch=s, depth=thk)
s.unsetPrimaryObject()

# make the assembly (prepare to pattern the thing)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
# this is probably still unnecesary
p = mdb.models['Model-1'].parts['Elastic_Portion']
a.Instance(name='Base Flexure', part=p, dependent=OFF)

# make a copy in a revolved state
# this is the rotation step (IDK WHY THERE IS A POINT HERE)
a.RadialInstancePattern(instanceList=('Base Flexure', ), point=(0.0,
    0.0, 0.0), axis=(0.0, 0.0, 1.0), number=n, totalAngle=360.0)

# Insert all your booleans and rigid surfaces into the assembly
p = mdb.models['Model-1'].parts['Inner-Boolean']
a.Instance(name='Inner Boolean', part=p, dependent=OFF)
p = mdb.models['Model-1'].parts['Inner-Surface']
a.Instance(name='Inner Surface', part=p, dependent=OFF)
p = mdb.models['Model-1'].parts['Outer-Boolean']
a.Instance(name='Outer Boolean', part=p, dependent=OFF)
p = mdb.models['Model-1'].parts['Outer-Surface']
a.Instance(name='Outer Surface', part=p, dependent=OFF)

# Use booleans to cut
a.InstanceFromBooleanMerge(name='Flexure-Combine', instances=(
    a.instances['Base Flexure'], a.instances['Base Flexure-rad-2'], ),
    originalInstances=SUPPRESS, domain=GEOMETRY)
a.InstanceFromBooleanCut(name='Flexure-Inner-Cut',
    instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['Flexure-Combine-1'],
    cuttingInstances=(a.instances['Inner Boolean'], ),
    originalInstances=SUPPRESS)
a.InstanceFromBooleanCut(name='Flexure-Outer-Cut',
    instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['Flexure-Inner-Cut-1'],
    cuttingInstances=(a.instances['Outer Boolean'], ),
    originalInstances=SUPPRESS)

# Make sure the final part is independent
a.makeIndependent(instances=(a.instances['Flexure-Outer-Cut-1'], ))

# Assign Section
mdb.models['Model-1'].HomogeneousSolidSection(name='Deformable-Section',
    material='Maraging300Steel', thickness=None)

p = mdb.models['Model-1'].parts['Flexure-Outer-Cut']
c = p.cells
cells = c.getByBoundingBox(-6,-6,-6,6,6,6)
region = p.Set(cells=cells, name='Section-Set')
p.SectionAssignment(region=region, sectionName='Deformable-Section', offset=0.0,
    offsetType=MIDDLE_SURFACE, offsetField='',
    thicknessAssignment=FROM_SECTION)

# Partition Faces
f1 = a.instances['Outer Surface'].faces
f2 = a.instances['Inner Surface'].faces
f3 = a.instances['Flexure-Outer-Cut-1'].faces

outerPartitionFaces = f1.findAt(((outer_radius,0,0)),)
innerPartitionFaces = f2.findAt(((0,inner_radius,0)),)

innerRightFace       = f3.findAt(((inner_radius,0,thk/2),))
innerLeftFace       = f3.findAt(((-inner_radius,0,thk/2),))

outerLeftFacePoint   = (outer_radius*np.cos(ang),outer_radius*np.sin(ang),thk/2)
outerLeftFace       = f3.findAt(((outerLeftFacePoint),))
outerRightFacePoint   = (-outer_radius*np.cos(ang),-outer_radius*np.sin(ang),thk/2)
outerRightFace       = f3.findAt(((outerRightFacePoint),))

a.PartitionFaceByIntersectFace(faces=innerPartitionFaces, cuttingFaces=innerRightFace) #, name="right_inner")
mdb.models['Model-1'].rootAssembly.features.changeKey(
    fromName='Partition face-1', toName='Inner-Right')
a.PartitionFaceByIntersectFace(faces=outerPartitionFaces, cuttingFaces=outerLeftFace) #, name="left_outer")
mdb.models['Model-1'].rootAssembly.features.changeKey(
    fromName='Partition face-1', toName='Outer-Left')

a.PartitionFaceByIntersectFace(faces=innerPartitionFaces, cuttingFaces=innerLeftFace) #, name="left_inner")
mdb.models['Model-1'].rootAssembly.features.changeKey(
    fromName='Partition face-1', toName='Inner-Left')
a.PartitionFaceByIntersectFace(faces=outerPartitionFaces, cuttingFaces=outerRightFace) #, name="right_outer")
mdb.models['Model-1'].rootAssembly.features.changeKey(
    fromName='Partition face-1', toName='Outer-Right')

f1 = a.instances['Outer Surface'].faces
f2 = a.instances['Inner Surface'].faces

innerRightFace2 = f2.findAt(((inner_radius,0,thk/2),))
innerLeftFace2  = f2.findAt(((-inner_radius,0,thk/2),))
outerLeftFace2  = f1.findAt(((outerLeftFacePoint),))
outerRightFace2 = f1.findAt(((outerRightFacePoint),))

# INFINITE BULLSHIT TO CREATE TIE CONSTRAINTS

a.Surface(side1Faces=innerRightFace, name='innerRightSurfaceD')
a.Surface(side1Faces=innerLeftFace, name='innerLeftSurfaceD')
a.Surface(side1Faces=outerRightFace, name='outerRightSurfaceD')
a.Surface(side1Faces=outerLeftFace, name='outerLeftSurfaceD')

a.Surface(side1Faces=innerRightFace2, name='innerRightSurface')
a.Surface(side1Faces=innerLeftFace2, name='innerLeftSurface')
a.Surface(side1Faces=outerRightFace2, name='outerRightSurface')
a.Surface(side1Faces=outerLeftFace2, name='outerLeftSurface')

mainRegion1 = a.surfaces['innerRightSurface']
mainRegion2 = a.surfaces['innerLeftSurface']
mainRegion3 = a.surfaces['outerRightSurface']
mainRegion4 = a.surfaces['outerLeftSurface']

secondaryRegion1 = a.surfaces['innerRightSurfaceD']
secondaryRegion2 = a.surfaces['innerLeftSurfaceD']
secondaryRegion3 = a.surfaces['outerRightSurfaceD']
secondaryRegion4 = a.surfaces['outerLeftSurfaceD']

mdb.models['Model-1'].Tie(name='innerRightTie', main=mainRegion1,
                          secondary=secondaryRegion1,
                          positionToleranceMethod=COMPUTED, adjust=ON,
                          tieRotations=ON, thickness=ON)
mdb.models['Model-1'].Tie(name='innerLeftTie', main=mainRegion2,
                          secondary=secondaryRegion2,
                          positionToleranceMethod=COMPUTED, adjust=ON,
                          tieRotations=ON, thickness=ON)
mdb.models['Model-1'].Tie(name='outerRightTie', main=mainRegion3,
                          secondary=secondaryRegion3,
                          positionToleranceMethod=COMPUTED, adjust=ON,
                          tieRotations=ON, thickness=ON)
mdb.models['Model-1'].Tie(name='outerLeftTie', main=mainRegion4,
                          secondary=secondaryRegion4,
                          positionToleranceMethod=COMPUTED, adjust=ON,
                          tieRotations=ON, thickness=ON)

# Encastre inner edge

f1 = a.instances['Inner Surface'].faces
faces1 = f1.findAt(((0, inner_radius, thk/2), (0, -inner_radius, thk/2)))
region = a.Set(faces=faces1, name='Boundary-Fix')
mdb.models['Model-1'].EncastreBC(name='Inner-Encastre', createStepName='Initial',
    region=region, localCsys=None)


# a.ReferencePoint(point=(0.0, 0.0, 0.0))
# mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='RP-1',
#     toName='Moment-Pt')

# f1 = a.instances['Outer Surface'].faces
# faces1 = f1.findAt(((0, outer_radius, thk/2), ), ((0, -outer_radius, thk/2), ))
# region2=a.Set(faces=faces1, name='b_Set-1')

# r1 = a.referencePoints
# refPoints1 = (a.referencePoints[32],)
# region1=regionToolset.Region(referencePoints=refPoints1)

# mdb.models['Model-1'].RigidBody(name='Constraint-5', refPointRegion=region1,
#     bodyRegion=region2)
# mdb.models['Model-1'].constraints['Constraint-5'].setValues(tieRegion=region2)
# mdb.models['Model-1'].constraints.changeKey(fromName='Constraint-5',
#     toName='Moment-Tie')


# region1 = regionToolset.Region(side1Faces=(innerRightFace,))
# region2 = regionToolset.Region(side1Faces=innerRightFace2)
# mdb.models['Model-1'].Tie(name='Constraint-1', main=region1, secondary=region2,
#     positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
# # very last thing you do is show it
session.viewports['Viewport: 1'].setValues(displayedObject=a)
if False:
    for f in glob.glob("*.rpy.*"):
        remove(f)
    for f in glob.glob("*.rec"):
        remove(f)
    for f in glob.glob("*.exception"):
        remove(f)