# built-in dependencies
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from os import remove
import glob
from utils import *

# Included Dependencies
# from spring import Spring
# from materials import Maraging300Steel

import subprocess
import sys

# Abaqqus interfacing dependencies
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup


class Abaqus_Spring:
    def __init__(self, name):
        self.Aname = name+"_A_surface.txt"
        self.Bname = name+"_A_surface.txt"
        self.springName = name

        self.A = np.genfromtxt("..\\surfaces\\"+self.Aname, delimiter=',')
        self.B = np.genfromtxt("..\\surfaces\\"+self.Bname, delimiter=',')
        self.params = csv_to_dict("..\\springs\\"+self.springName)

        self.A3D = list(map(tuple,self.A))
        self.B3D = list(map(tuple,self.B))
        self.A = self.A[:,0:2]
        self.A = list(map(tuple,self.A))
        self.B = self.B[:,0:2]
        self.B = list(map(tuple,self.B))

    def deform_with_abaqus(self):
        try:
            output = subprocess.check_output("abaqus cae script=Abaqus-Spring-Deform.py")
            print(output.decode('utf-8'))
        except:
            pass
        executeOnCaeStartup()

        n = self.params['n']
        thk = self.params['Out_Plane_Thk']
        ang = self.params[final_of_value(self.params,"beta_ang_")]
        overallLen = self.params['full_length']

        A = self.A
        B = self.B
        A3D = self.A3D
        B3D = self.B3D
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

        a.DatumCsysByThreePoints(name='Cylindrical-CSYS', coordSysType=CYLINDRICAL,
            origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 1.0, 0.0),
            isDependent=False)

        # Place an instance of the spring leg
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
            material='Maraging300', thickness=None)

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
        faces1 = f1.findAt(((0, inner_radius, thk/2), ))
        region = a.Set(faces=faces1, name='Boundary-Fix-a')
        mdb.models['Model-1'].EncastreBC(name='Upper-Encastre', createStepName='Initial',
            region=region, localCsys=None)
        faces1 = f1.findAt(((0, -inner_radius, thk/2), ))
        region = a.Set(faces=faces1, name='Boundary-Fix-2')
        mdb.models['Model-1'].EncastreBC(name='Lower-Encastre', createStepName='Initial',
            region=region, localCsys=None)

        # Mesh

        rigidSize  = 0.125
        deformSize = overallLen/200
        dev        = 0.01
        min        = 0.01

        c1 = a.instances['Flexure-Outer-Cut-1'].cells
        pickedRegions = c1.findAt(((inner_radius, 0, thk/2), ),
                                ((-inner_radius, 0, thk/2), ))
        a.setMeshControls(regions=pickedRegions, algorithm=MEDIAL_AXIS)

        partInstances =(a.instances['Outer Surface'], )
        a.seedPartInstance(regions=partInstances, size=rigidSize, deviationFactor=dev,
            minSizeFactor=min)
        partInstances =(a.instances['Inner Surface'], )
        a.seedPartInstance(regions=partInstances, size=rigidSize, deviationFactor=dev,
            minSizeFactor=min)
        partInstances =(a.instances['Flexure-Outer-Cut-1'], )
        a.seedPartInstance(regions=partInstances, size=deformSize, deviationFactor=dev,
            minSizeFactor=min)
        partInstances =(a.instances['Flexure-Outer-Cut-1'],
            a.instances['Inner Surface'], a.instances['Outer Surface'], )
        a.generateMesh(regions=partInstances)

        # make the load
        mdb.models['Model-1'].StaticStep(name='Moment-Step', previous='Initial',
            initialInc=0.1)
        r1 = a.instances['Outer Surface'].referencePoints
        refPoints1=(r1[2], )
        region = a.Set(referencePoints=refPoints1, name='Set-3')
        datum = mdb.models['Model-1'].rootAssembly.datums[2]
        mdb.models['Model-1'].Moment(name='Moment-Load', createStepName='Moment-Step',
            region=region, cm3=4550.0, distributionType=UNIFORM, field='',
            localCsys=datum)


        # make the job
        mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS,
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
            scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1,
            multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)

        # do the job
        myJob = mdb.Job(name='Job-1', model='Model-1')
        # mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
        myJob.submit()
        # mdb.jobs.waitForCompletion(timeOut=30)
        myJob.waitForCompletion()

        # # very last thing you do is show it
        o3 = session.openOdb(
            name='C:/Stuff/SMM/PRS-Code-Git/ODE-Python/ModelDBs/Job-1.odb')
        # Go to results view
        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
        session.viewports['Viewport: 1'].setValues(
            displayedObject=session.odbs['C:/Stuff/SMM/PRS-Code-Git/ODE-Python/ModelDBs/Job-1.odb'])
        session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
        session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
            meshTechnique=OFF)
        session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
            CONTOURS_ON_UNDEF, CONTOURS_ON_DEF, ))

        dtm = session.odbs['C:/Stuff/SMM/PRS-Code-Git/ODE-Python/ModelDBs/Job-1.odb'].rootAssembly.datumCsyses['ASSEMBLY__T-CYLINDRICAL-CSYS']
        session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
            transformationType=USER_SPECIFIED, datumCsys=dtm)
        # # make it not look like [expletive]
        # View tangential deformation
        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
            variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'), )
        # View real deflection
        session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
            deformationScaling=UNIFORM, uniformScaleFactor=1)
        session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
        session.viewports['Viewport: 1'].view.setProjection(projection=PARALLEL)
        session.viewports['Viewport: 1'].makeCurrent()

