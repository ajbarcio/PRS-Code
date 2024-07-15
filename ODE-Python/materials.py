import numpy as np

psi2Pa = 6894.76

# CUSTOMARY UNITS!!!
# PSI

class Material:
    def __init__(self, E, G, yieldStress):
        self.E = E
        self.G = G
        self.yieldStress = yieldStress

## whole bunch of options

Aluminum = Material(5, 78, 2)

Steel = Material(4, 79, 43)

~~~~ NEW FILE ~~~~~

# geometries module

class Geometry:
    def __init__(self, thickness, resolution):
        self.thickness = thickness
        self.resolution = resolution

class Silly_Spring(Geometry):

    def __init__(self, thickness, resolution, *geometryArgs):
        super().__init__(thickness, resolution)

    def get_Ic(self, xi):
        # do some stuff
        #call a whole bunch of functions
        return Ic

    def a_whole_bunch_of_Ic_Methods():
        pass

spring1 = Silly_Spring(3, 5, 3, 4, 6, 7, 4)


class Maraging300Steel:
    def __init__(self):
        self.E = 27500000
        self.G = self.E/(2*(1.33))
        self.yieldStress=309700


