import numpy as np

psi2Pa = 6894.76

# CUSTOMARY UNITS!!!
# PSI

class Maraging300Steel:
    def __init__(self):
        self.E = 27500000
        self.G = self.E/(2*(1.33))
        self.yieldStress=309700
