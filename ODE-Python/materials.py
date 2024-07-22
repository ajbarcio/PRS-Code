import numpy as np

psi2Pa = 6894.76

# CUSTOMARY UNITS!!!
# PSI

class Material:
    def __init__(self, E, yieldStress, poisson=0.33):
        self.E = E
        self.G = self.E/(2*(1+poisson))
        self.yieldStress = yieldStress

# Initialize other materials here and import where used

Maraging300Steel = Material(27500000, 309700)


