# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2023 replay file
# Internal Version: 2022_09_28-13.11.55 183150
# Run by ajbarcio on Tue Jul 30 21:21:39 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=673.625549316406, 
    height=216.300003051758)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
execfile('Abaqus_Interface.py', __main__.__dict__)
#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.26702, 
    farPlane=5.16107, width=8.23935, height=2.53671, cameraPosition=(0.454462, 
    0.577462, 4.71405), cameraTarget=(0.454462, 0.577462, 0))
mdb.models['Model-1'].sketches.changeKey(fromName='__profile__', 
    toName='__save__')
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.retrieveSketch(sketch=mdb.models['Model-1'].sketches['__save__'])
del mdb.models['Model-1'].sketches['__save__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.72897, 
    farPlane=10.082, width=21.6849, height=6.89092, cameraPosition=(1.74309, 
    1.43215, 8.90547), cameraTarget=(1.74309, 1.43215, 0))
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseSolidExtrude(sketch=s, depth=0.375)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.44928, 
    farPlane=11.4256, width=10.4442, height=3.21552, cameraPosition=(5.46957, 
    5.89988, 5.09453), cameraTarget=(0.309526, 0.73983, -0.0655142))
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.08965, 
    farPlane=13.1529, width=9.82212, height=3.02401, viewOffsetX=0.356566, 
    viewOffsetY=-0.0919337)
a = mdb.models['Model-1'].rootAssembly
a.RadialInstancePattern(instanceList=('Part-1-1', ), point=(0.0, 0.0, 0.0), 
    axis=(0.0, 0.0, 1.0), number=2, totalAngle=360.0)
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.15443, 
    farPlane=13.0881, width=8.75813, height=2.69643, viewOffsetX=0.605678, 
    viewOffsetY=-0.206244)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=10.4214, 
    farPlane=12.4014, width=13.1398, height=4.04545, viewOffsetX=0.0531051, 
    viewOffsetY=0.127501)
session.viewports['Viewport: 1'].view.setProjection(projection=PARALLEL)
session.viewports['Viewport: 1'].view.setValues(nearPlane=10.5313, 
    farPlane=12.2914, width=12.7653, height=3.93015, cameraPosition=(0.474984, 
    -0.0382662, 11.5989), cameraTarget=(0.474984, -0.0382662, 0.1875))
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
a = mdb.models['Model-1'].rootAssembly
d11 = a.datums
a.ReferencePoint(point=d11[1].origin)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
session.viewports['Viewport: 1'].view.setValues(nearPlane=8.20953, 
    farPlane=14.6132, cameraPosition=(-6.82569, -4.94283, 7.88151), 
    cameraUpVector=(-0.509699, 0.704344, -0.494071), cameraTarget=(
    -1.93715e-06, -1.68383e-06, 0.1875))
a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[6], )
region1=a.Set(referencePoints=refPoints1, name='m_Set-1')
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-1-1-rad-2'].faces
faces1 = f1.findAt(((1.982571, -0.14831, 0.25), ))
f2 = a.instances['Part-1-1'].faces
faces2 = f2.findAt(((-1.982571, 0.14831, 0.25), ))
region2=a.Set(faces=faces1+faces2, name='s_Set-1')
mdb.models['Model-1'].MultipointConstraint(name='Constraint-1', 
    controlPoint=region1, surface=region2, mpcType=TIE_MPC, 
    userMode=DOF_MODE_MPC, userType=0, csys=None)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, interactions=OFF, constraints=OFF, 
    engineeringFeatures=OFF)
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
session.viewports['Viewport: 1'].view.setValues(nearPlane=10.3549, 
    farPlane=12.4678, width=16.0159, height=4.93094, cameraPosition=(0.0799735, 
    0.130918, 11.5989), cameraTarget=(0.0799735, 0.130918, 0.1875))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.26748, 
    farPlane=13.6776, cameraPosition=(2.63492, 3.90176, 10.6504), 
    cameraUpVector=(0.215229, 0.900777, -0.377198))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.48065, 
    farPlane=13.2948, cameraPosition=(-3.12673, -1.2234, 11.0696), 
    cameraUpVector=(0.092544, 0.984779, 0.147126), cameraTarget=(0.0492542, 
    0.103592, 0.189735))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.92839, 
    farPlane=12.9164, cameraPosition=(-0.69551, 1.55833, 11.4822), 
    cameraUpVector=(0.232231, 0.966446, -0.109775), cameraTarget=(0.0442107, 
    0.0978214, 0.188879))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.68659, 
    farPlane=13.1966, cameraPosition=(0.793988, 3.10779, 11.1708), 
    cameraUpVector=(0.269691, 0.923889, -0.271469), cameraTarget=(0.0456483, 
    0.0993168, 0.188578))
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-1-1'].faces
faces1 = f1.findAt(((1.34714, -0.04714, 0.25), ))
f2 = a.instances['Part-1-1-rad-2'].faces
faces2 = f2.findAt(((-1.34714, 0.04714, 0.25), ))
region = a.Set(faces=faces1+faces2, name='Set-3')
mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial', 
    region=region, localCsys=None)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
mdb.models['Model-1'].Material(name='M300')
mdb.models['Model-1'].materials['M300'].Elastic(table=((27500000.0, 0.33), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', 
    material='M300', thickness=None)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
    engineeringFeatures=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
a1 = mdb.models['Model-1'].rootAssembly
a1.InstanceFromBooleanMerge(name='Part-2', instances=(a1.instances['Part-1-1'], 
    a1.instances['Part-1-1-rad-2'], ), originalInstances=SUPPRESS, 
    domain=GEOMETRY)
p1 = mdb.models['Model-1'].parts['Part-2']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p = mdb.models['Model-1'].parts['Part-2']
c = p.cells
cells = c.findAt(((-0.908573, 1.535566, 0.1875), ), ((-1.297377, -1.175368, 
    0.1875), ))
region = p.Set(cells=cells, name='Set-1')
p = mdb.models['Model-1'].parts['Part-2']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
session.viewports['Viewport: 1'].setValues(displayedObject=a)
a1 = mdb.models['Model-1'].rootAssembly
a1.makeIndependent(instances=(a1.instances['Part-2-1'], ))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
a = mdb.models['Model-1'].rootAssembly
c1 = a.instances['Part-2-1'].cells
pickedRegions = c1.findAt(((-0.908573, 1.535566, 0.1875), ), ((-1.297377, 
    -1.175368, 0.1875), ))
a.setMeshControls(regions=pickedRegions, algorithm=MEDIAL_AXIS)
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['Part-2-1'], )
a.seedPartInstance(regions=partInstances, size=0.125, deviationFactor=0.01, 
    minSizeFactor=0.01)
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['Part-2-1'], )
a.generateMesh(regions=partInstances)
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.68818, 
    farPlane=13.1345, width=16.4286, height=5.06515, cameraPosition=(0.754341, 
    2.99404, 11.1732), cameraTarget=(0.00600398, -0.0144234, 0.191042))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, 
    adaptiveMeshConstraints=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', 
    initialInc=0.1)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[6], )
region = a.Set(referencePoints=refPoints1, name='Set-5')
mdb.models['Model-1'].Moment(name='Load-1', createStepName='Step-1', 
    region=region, cm3=4550.0, distributionType=UNIFORM, field='', 
    localCsys=None)
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.50612, 
    farPlane=13.3166, width=19.8075, height=6.0983, cameraPosition=(0.837379, 
    3.05312, 11.1514), cameraTarget=(0.0890423, 0.0446615, 0.169198))
session.viewports['Viewport: 1'].view.setValues(nearPlane=8.53324, 
    farPlane=14.2077, cameraPosition=(-4.98639, 3.64027, 9.73637), 
    cameraUpVector=(0.127109, 0.948805, -0.289157), cameraTarget=(0.0890423, 
    0.0446615, 0.169198))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
    multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.36344, 
    farPlane=13.3775, width=5.7698, height=1.77639, cameraPosition=(-4.24974, 
    3.74001, 10.0897), cameraTarget=(0.82569, 0.144404, 0.522508))
session.viewports['Viewport: 1'].view.setValues(nearPlane=10.0346, 
    farPlane=13.3542, cameraPosition=(-0.216227, 5.86379, 10.3382), 
    cameraUpVector=(0.291951, 0.839106, -0.458984), cameraTarget=(0.81119, 
    0.136769, 0.521614))
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-2-1'].faces
faces1 = f1.findAt(((1.34714, -0.04714, 0.25), ), ((-1.34714, 0.04714, 0.25), 
    ))
region = a.Set(faces=faces1, name='Set-6')
mdb.models['Model-1'].boundaryConditions['BC-1'].setValues(region=region)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON)
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.62323, 
    farPlane=13.7656, width=13.7205, height=4.22424, cameraPosition=(2.36912, 
    5.88595, 10.5959), cameraTarget=(3.39654, 0.158924, 0.779275))
session.viewports['Viewport: 1'].view.setValues(nearPlane=7.27098, 
    farPlane=12.473, cameraPosition=(-3.72775, 2.40311, 9.47705), 
    cameraUpVector=(0.0630258, 0.97567, -0.209991), cameraTarget=(3.24899, 
    0.0746344, 0.752198))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.97395, 
    farPlane=13.9091, cameraPosition=(3.69652, 3.98952, 11.0155), 
    cameraUpVector=(0.492186, 0.780883, -0.384676), cameraTarget=(2.09128, 
    -0.172745, 0.51229))
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-2-1'].faces
faces1 = f1.findAt(((-1.982571, 0.14831, 0.25), ), ((1.982571, -0.14831, 0.25), 
    ))
region2=a.Set(faces=faces1, name='s_Set-7')
mdb.models['Model-1'].constraints['Constraint-1'].setValues(surface=region2)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
    constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: The job input file "Job-1.inp" has been submitted for analysis.
#: Error in job Job-1: 1 nodes have inactive dof on which boundary conditions are specified. The nodes have been identified in node set ErrNodeBCInactiveDof.
#: Job Job-1: Analysis Input File Processor aborted due to errors.
#: Error in job Job-1: Analysis Input File Processor exited with an error - Please see the  Job-1.dat file for possible error messages if the file exists.
#: Job Job-1 aborted due to errors.
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.43121, 
    farPlane=14.4517, width=21.7357, height=6.90706, cameraPosition=(6.25852, 
    2.99097, 11.0197), cameraTarget=(4.65328, -1.1713, 0.516443))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.81551, 
    farPlane=14.0674, width=17.5157, height=5.39268, cameraPosition=(5.35047, 
    3.18986, 11.0796), cameraTarget=(3.74523, -0.9724, 0.576403))
session.viewports['Viewport: 1'].view.setValues(nearPlane=8.39768, 
    farPlane=12.6351, cameraPosition=(-0.13674, 0.49233, 11.2846), 
    cameraUpVector=(0.217362, 0.973559, -0.0702633), cameraTarget=(3.50164, 
    -1.09215, 0.585502))
session.viewports['Viewport: 1'].view.setValues(nearPlane=8.95639, 
    farPlane=12.0764, width=7.36575, height=2.26775, cameraPosition=(-1.42223, 
    1.03947, 10.7664), cameraTarget=(2.21615, -0.545009, 0.0673245))
session.viewports['Viewport: 1'].view.setValues(nearPlane=10.4858, 
    farPlane=13.3761, cameraPosition=(6.07313, -0.254066, 10.5712), 
    cameraUpVector=(0.248543, 0.960767, -0.123096), cameraTarget=(1.5783, 
    -0.43493, 0.0839369))
session.viewports['Viewport: 1'].view.setValues(nearPlane=10.169, 
    farPlane=13.6929, width=13.6753, height=4.21033, cameraPosition=(7.59426, 
    -0.415285, 9.92203), cameraTarget=(3.09943, -0.596149, -0.565237))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.71494, 
    farPlane=13.4487, cameraPosition=(6.56723, 1.0027, 10.1946), 
    cameraUpVector=(0.238413, 0.947426, -0.213411), cameraTarget=(3.0547, 
    -0.534394, -0.553367))
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.84868, 
    farPlane=13.315, width=10.677, height=3.28721, cameraPosition=(5.91619, 
    1.1021, 10.3932), cameraTarget=(2.40366, -0.434991, -0.354817))
o3 = session.openOdb(name='C:/Stuff/SMM/PRS-Code-Git/ODE-Python/Job-1.odb')
#: Model: C:/Stuff/SMM/PRS-Code-Git/ODE-Python/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       6
#: Number of Node Sets:          9
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.odbs['C:/Stuff/SMM/PRS-Code-Git/ODE-Python/Job-1.odb'])
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=8.15714, 
    farPlane=14.6649, width=13.3955, height=4.12418, cameraPosition=(7.3586, 
    6.89399, 5.6994), cameraTarget=(0.770433, 0.305825, -0.888758))
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.99242, 
    farPlane=13.1713, width=9.43421, height=2.90458, cameraPosition=(5.77779, 
    1.14382, 10.4324), cameraTarget=(2.26526, -0.39327, -0.315553))
a = mdb.models['Model-1'].rootAssembly
region = a.sets['Set-6']
mdb.models['Model-1'].boundaryConditions['BC-1'].setValues(region=region)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: The job input file "Job-1.inp" has been submitted for analysis.
#: Error in job Job-1: 1 nodes have inactive dof on which boundary conditions are specified. The nodes have been identified in node set ErrNodeBCInactiveDof.
#: Job Job-1: Analysis Input File Processor aborted due to errors.
#: Error in job Job-1: Analysis Input File Processor exited with an error - Please see the  Job-1.dat file for possible error messages if the file exists.
#: Job Job-1 aborted due to errors.
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.66458, 
    farPlane=13.4991, width=15.4769, height=4.91815, cameraPosition=(7.56288, 
    1.15194, 9.84787), cameraTarget=(4.05035, -0.38515, -0.900098))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
del mdb.models['Model-1'].boundaryConditions['BC-1']
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, interactions=ON, constraints=ON, 
    engineeringFeatures=ON)
del mdb.models['Model-1'].constraints['Constraint-1']
session.viewports['Viewport: 1'].view.setValues(nearPlane=9.32335, 
    farPlane=13.8404, width=18.7914, height=5.92132, cameraPosition=(8.39517, 
    0.814126, 9.62418), cameraTarget=(4.88264, -0.722967, -1.12379))
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, interactions=OFF, constraints=OFF, 
    engineeringFeatures=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=8.41976, 
    farPlane=14.403, cameraPosition=(-7.99594, 2.46207, 7.94782), 
    cameraUpVector=(0.0903641, 0.763614, -0.639318), cameraTarget=(
    -2.11596e-06, -5.96046e-07, 0.1875))
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-2-1'].faces
faces1 = f1.findAt(((1.34714, -0.04714, 0.25), ), ((-1.34714, 0.04714, 0.25), 
    ))
region = a.Set(faces=faces1, name='Set-8')
mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Initial', 
    region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, interactions=ON, constraints=ON, 
    engineeringFeatures=ON)
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
session.viewports['Viewport: 1'].view.setValues(nearPlane=8.2644, 
    farPlane=14.5583, cameraPosition=(9.62659, 1.82087, 6.03836), 
    cameraUpVector=(-0.399684, 0.877345, -0.265554), cameraTarget=(
    -1.90735e-06, -5.28991e-07, 0.1875))
a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[6], )
region1=a.Set(referencePoints=refPoints1, name='m_Set-9')
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-2-1'].faces
faces1 = f1.findAt(((-1.982571, 0.14831, 0.25), ), ((1.982571, -0.14831, 0.25), 
    ))
region2=a.Set(faces=faces1, name='s_Set-9')
mdb.models['Model-1'].MultipointConstraint(name='Constraint-1', 
    controlPoint=region1, surface=region2, mpcType=TIE_MPC, 
    userMode=DOF_MODE_MPC, userType=0, csys=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
    constraints=OFF, connectors=OFF, engineeringFeatures=OFF)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: The job input file "Job-1.inp" has been submitted for analysis.
#: Error in job Job-1: 1 nodes have inactive dof on which boundary conditions are specified. The nodes have been identified in node set ErrNodeBCInactiveDof.
#: Job Job-1: Analysis Input File Processor aborted due to errors.
#: Error in job Job-1: Analysis Input File Processor exited with an error - Please see the  Job-1.dat file for possible error messages if the file exists.
#: Job Job-1 aborted due to errors.
session.viewports['Viewport: 1'].view.setValues(nearPlane=8.61252, 
    farPlane=14.2102, width=10.219, height=3.32359, cameraPosition=(9.54515, 
    2.08845, 6.08908), cameraTarget=(-0.0814434, 0.267579, 0.238224))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON)
o3 = session.openOdb(name='C:/Stuff/SMM/PRS-Code-Git/ODE-Python/Job-1.odb')
#: Model: C:/Stuff/SMM/PRS-Code-Git/ODE-Python/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       8
#: Number of Node Sets:          13
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.odbs['C:/Stuff/SMM/PRS-Code-Git/ODE-Python/Job-1.odb'])
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
    predefinedFields=OFF, connectors=OFF)
