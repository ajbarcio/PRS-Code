from materials import Maraging300Steel

class Spring:
    def __init__(self, path_definition, thickness_definition, material,
                 # The following are default values:
                 t          = 0.375,
                 resolution = 200):
        
        self.t = t

        self.path = path_definition
        self.geometry = thickness_definition
        self.material = material

        self.fullArcLength = self.path.fullParamLength

        self.resl = resolution
        self.len = self.resl+1
        self.step = self.fullArcLength/self.resl
        self.endIndex = self.len-1
