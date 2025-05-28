import numpy as np

psi2Pa = 6894.76

# CUSTOMARY UNITS!!!
# PSI

class Material:
    def __init__(self, name, E, yieldStress, poisson=0.33, fatigueStress=0):
        self.name=name
        self.E = E
        self.G = self.E/(2*(1+poisson))
        self.yieldStress = yieldStress
        if fatigueStress:
            self.fatigueStress = fatigueStress
        self.reset_designStress()
    
    def ultimateStress_set(self, ultimateStress):
        self.ultimateStress = ultimateStress
        self.fatigueStress  = 0.5*ultimateStress

    def reset_designStress(self):
        if hasattr(self, "ultimateStress"):
            self.designStress=self.fatigueStress
        else:
            self.designStress=0.8*self.yieldStress
        pass

# Initialize other materials here

Maraging300Steel = Material("Maraging300Steel", 27500000, 309700)
Maraging300Steel.ultimateStress_set(314600)
TestMaterial     = Material("TestMaterial", 100000, 100000)
Titanium5        = Material("Titanium5",113800*1e6/psi2Pa, 880*1e6/psi2Pa)

AL7075 = Material("AL7075", 10442717, 73244.1)
Delrin = Material("Delrin", 449617, 11385.46) # Should be most compliant
AISI4640 = Material("AISI4640", 209000*1000000/psi2Pa,1103*1000000/psi2Pa)
AISI4340Tempered = Material("AISI4340Tempered", 200000*1000/psi2Pa, 862*1000/psi2Pa)
# Aluminum         = Material()
