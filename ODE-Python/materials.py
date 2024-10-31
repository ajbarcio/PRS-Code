import numpy as np

psi2Pa = 6894.76

# CUSTOMARY UNITS!!!
# PSI

class Material:
    def __init__(self, E, yieldStress, poisson=0.33, fatigueStress=0):
        self.E = E
        self.G = self.E/(2*(1+poisson))
        self.yieldStress = yieldStress
        if fatigueStress:
            self.fatigueStress = fatigueStress
    
    def ultimateStress_set(self, ultimateStress):
        self.ultimateStress = ultimateStress
        self.fatigueStress  = 0.5*ultimateStress

    def reset_designStrength(self):
        if self.ultimateStress:
            self.designStress=self.fatigueStress
        else:
            self.designStress=0.8*self.yieldStress
        pass

# Initialize other materials here

Maraging300Steel = Material(27500000, 309700)
Maraging300Steel.ultimateStress_set(314600)
TestMaterial     = Material(20000000, 200000)

# Aluminum         = Material()
