# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2023 replay file
# Internal Version: 2022_09_28-13.11.55 183150
# Run by ajbarcio on Thu Jul 18 12:04:03 2024
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
execfile('./Deform_default_spring.py', __main__.__dict__)
#* AttributeError: 'module' object has no attribute 'integrate'
#* File "./Deform_default_spring.py", line 18, in <module>
#*     defaultSpring = spring.generate_default_spring(name="default_spring")
#* File "C:\Stuff\SMM\PRS-Code-Git\ODE-Python\spring.py", line 1053, in 
#* generate_default_spring
#*     name=name)
#* File "C:\Stuff\SMM\PRS-Code-Git\ODE-Python\spring.py", line 65, in __init__
#*     correctedFullLength = self.measure_length()
#* File "C:\Stuff\SMM\PRS-Code-Git\ODE-Python\spring.py", line 1011, in 
#* measure_length
#*     self.smesh[1:len(self.smesh)] = 
#* scipy.integrate.cumulative_trapezoid(integrand, self.ximesh, 
#* self.ximesh[1]-self.ximesh[0])
