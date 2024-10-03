import numpy as np

psi2Pa = 6894.76

# CUSTOMARY UNITS!!!
# PSI

class Material:
    def __init__(self, E, yieldStress, poisson=0.33):
        self.E = E
        self.G = self.E/(2*(1+poisson))
        self.yieldStress = yieldStress

# Initialize other materials here

Maraging300Steel = Material(27500000, 309700)
# Maraging300Steel = Material(190000*1e6/psi2Pa, 2140*1e6/psi2Pa)
Aluminum7075     = Material(72000*1e6/psi2Pa, 505*1e6/psi2Pa)
AISI4340         = Material(200000*1e6/psi2Pa, 862*1e6/psi2Pa)
Titanium5        = Material(113800*1e6/psi2Pa, 880*1e6/psi2Pa)