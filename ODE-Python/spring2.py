from materials import Maraging300Steel

class Spring:
    def __init__(self, geoDef, material,
                 # The following are default values:
                 t          = 0.375,
                 resolution = 200):

        # Parameters and data structures from passed objects
        # Path accessor from path definition
        self.path     = geoDef.pathDef
        # Geometry accessor
        self.geom     = geoDef
        # Material info
        self.material = material

        self.fullArcLength = self.geom.fullParamLength

        # thickness arg
        self.t = t

        # Mesh information
        self.resl = resolution
        self.len = self.resl+1
        self.step = self.fullArcLength/self.resl
        self.endIndex = self.len-1

        # fuck your memory, extract commonly used attributes from passed
        # definition objects:
        self.E = self.material.E

    def deform_ODE(self, xi, deforms, *args):
        # Deal with ajrgs in a way that is hopefully better
        Fx    = args[0]
        Fy    = args[0]
        dgdx0 = args[0]
        # Prepare moment dimensionalizer
        Mdim = self.E*self.crsc.get_Ic(xi)