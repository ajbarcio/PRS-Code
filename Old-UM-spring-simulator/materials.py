import numpy as np

psi2Pa = 6894.76

class MaragingSteelC300(object):
    def __init__(self):
        # Material properties defined
        self.E = 27500000*psi2Pa
        self.G = self.E/(2*(1.33))
        self.yield_stress=309700*psi2Pa

class Aluminum6061(object):
    def __init__(self):
        # Material properties defined
        self.E = 68.9e9 # Pa (GPa) Young's modulus
        # self.eta = .33 # poisson's ratio
        self.G = 26e9# Pa (GPa) shear modulus
        self.yield_stress=276e6 # Pa (MPa) yield stress
        # self.shear_yield = 207e6 # Pa (MPa) shear yield stress
        # self.fatigue_stress = 97e6 # Pa (MPa) yield stress at 500e6 complete cycles

class Aluminum7075(object):
    def __init__(self):
        # Material properties defined
        self.E = 71.7e9 # Pa (GPa) Young's modulus
        # self.eta = .33 # poisson's ratio
        self.G = 26.9e9# Pa (GPa) shear modulus
        self.yield_stress=572e6 # Pa (MPa) yield stress
        # self.shear_yield = 331e6 # Pa (MPa) shear yield stress
        self.fatigue_stress = 159e6 # Pa (MPa) 


class Aluminum7075SW(object):
    def __init__(self):
        # Material properties defined
        self.E = 7.2e10 # Pa (GPa) Young's modulus
        self.eta = .33 # poisson's ratio
        self.G = 2.69e10# Pa (GPa) shear modulus
        self.yield_stress=505e6 # Pa (MPa) yield stress
        self.rho = 2810 # kg/m^3
        # self.shear_yield = 331e6 # Pa (MPa) shear yield stress
        self.fatigue_stress = 159e6 # Pa (MPa) 

class SS410(object):
    # Hardened and annealed at 350 degrees 
    def __init__(self):
        # Material properties defined 0 
        self.E = 200e9# Pa (GPa) Young's modulus
        # self.eta = .26 # poisson's ratio
        self.G = .38*self.E# Pa (GPa) shear modulus
        self.yield_stress=185e3*psi2Pa # Pa (ksi) yield stress
        # heated to 980°C (light orange hot) 30 minutes, then oil quenched
        self.yield_stress=721e6 # MPa, tough, 650°C 2H temper 
        self.yield_stress=807e6 # MPa, tough, 605°C 2H temper 
        self.yield_stress=835e6 # MPa, tough?, 508°C 2H temper 
        self.yield_stress=920e6 # MPa, less brittle, 453°C 2H temper 
        self.yield_stress=961e6 # MPa, less brittle, 343°C 2H temper (grey green temper)
        # self.yield_stress=1005e6 # MPa, less brittle, 233°C 2H temper (yellow temper)
        self.yield_stress=1225e6 # MPa, super brittle, 49°C 2H temper 
        # self.yield_stress=1018209985/0.9 # yield stress that results in a serpentine factor of 1.000000000
        self.rho = 7.74e3 # kg/m^3


class AISI1095(object):
    # Carbon Steel
    def __init__(self):
        # Material properties defined
        # self.rho = 7.85 # g/cm^3
        self.E = 205e9# Pa (GPa) Young's modulus
        # self.eta = .29 # poisson's ratio
        self.G = 80e9# Pa (GPa) shear modulus
        self.yield_stress = 703e6 # Pa (MPa) yield stress quenched 800 C (oil), tempered 480 C
        self.endurance_yield_stress = 400e6 # MPa (https://www.aksteel.com/sites/default/files/2018-11/410-stainless.pdf, Table 2)
        # self.yield_stress = 703e6 # Pa (MPa) yield stress quenched 800 C (oil), tempered 540 C
        # self.yield_stress = 525e6 # Pa yield stress
        # self.yield_stress = 495e6 # Pa yield stress air-cooled (matweb)
        # self.yield_stress=71800* psi2Pa  # #

class Titanium5(object):
    # Hardened and annelled at 350 degrees
    def __init__(self):
        # Material properties defined


        self.E = 113.8e9 # Pa (GPa) Young's modulus
        # self.eta = .342 # poisson's ratio
        self.G = 44e9# Pa (GPa) shear modulus
        self.yield_stress=880e6  # Pa (MPa) yield stress

class Titanium5STA(object):
    # Titanium Ti-6Al-4V (Grade 5), STA Bar 
    # http://www.matweb.com/search/DataSheet.aspx?MatGUID=f87a4a1c92d34da2b1ecde4e4dec7a73
    def __init__(self):
        # Material properties defined

        self.rho = 4.43e3 # kg/m^3
        self.E = 114e9 # Pa (GPa) Young's modulus
        # self.eta = .342 # poisson's ratio
        self.G = 44e9# Pa (GPa) shear modulus

        # Important to only use fatigue strength estimates evaluated at the
        # correct heat treatment settings. Most references I could find give
        # the much lower fatigue behavior of the annealed condition. Further
        # information might be obtainable through the references listed in the
        # matweb site, but they are all over 300 dollars in 2021.

        # Nominal:
        self.yield_stress = 965e6 # Pa (MPa) http://www.matweb.com/search/datasheet.aspx?MatGUID=f87a4a1c92d34da2b1ecde4e4dec7a73
        # Fatigue toughness:
        self.yield_stress = 700e6 # Pa (MPa) http://www.matweb.com/search/datasheet.aspx?MatGUID=f87a4a1c92d34da2b1ecde4e4dec7a73


if __name__ == '__main__2':
    # mat = Aluminum6061()
    steel = AISI1095()
    stainless = SS410()
    print(stainless.yield_stress/stainless.E)
    print(steel.yield_stress/steel.E)
    print(Titanium5().yield_stress/Titanium5().E)
    print(Aluminum7075().yield_stress/Aluminum7075().E)

    def energy_eval(mat):
        return mat.yield_stress**2/(6.* mat.E) * 1e-6

    print("Energy/Volume (mJ/mm^3):")
    print("1095   ", energy_eval(AISI1095()))
    print("SS410  ", energy_eval(SS410()))
    print("Ti 5   ", energy_eval(Titanium5()))
    print("Al 7075", energy_eval(Aluminum7075()))


    total_target_rate = 600 # Nm/rad
    spring_target = total_target_rate/6
    tooth_target = spring_target/24
    max_deflection = 15 * np.pi/180
    nominal_tooth_radius = 6e-3 # m
    FOS = 1.3
    tooth_normal_force_polar = max_deflection*tooth_target/nominal_tooth_radius * FOS
    print("tooth_normal_force_polar", tooth_normal_force_polar)
    tooth_contact_angle = 25.45*np.pi/180
    tooth_normal_force_radial = tooth_normal_force_polar*np.tan(tooth_contact_angle)
    tooth_normal_force_total = tooth_normal_force_polar/np.cos(tooth_contact_angle)
    print("tooth_normal_force_radial", tooth_normal_force_radial)
    print("tooth_normal_force_total", tooth_normal_force_total)
    rs = np.linspace(6, 32, 9)
    taus = (rs-6)*1e-3*tooth_normal_force_total # Nm

    
    strain = stainless.yield_stress/stainless.E
    print()
    print("strain goal", strain*100, "%")

    # strain = width/2*curvature
    # E I curvature = tauS
    # I = depth * width**3 /12
    # curvature = 1/R
    # width = sqrt(6*tau / (strain * E * depth))
    widths = np.sqrt(6* taus / (strain*stainless.E*5e-3)) *1e3 # mm

    # G*A*shear_strain = f
    # G * x * y * strain = f
    # x = f/(G*y*strain)

    min_width = tooth_normal_force_total / (strain* stainless.G*5e-3) *1e3

    # print(taus)

    # print()
    print("design guide:", list(zip(rs , widths)))
    print("\t".join(["%.5f"%w for w in widths]))
    print()
    # print(widths)
    print("min_width", min_width)


    R = 200 # mm
    width = 1 # mm
    print("strain", (width/2)/R, strain)
    print("radius", (width/2)/strain, R)

    # (R + w/2)/R = steel.yield_stress/steel.E



def old_pulley_math():
    mat = Aluminum7075()
    FoS = 3 # yield/working
    tau = 50 # Nm
    rad = .08 # m
    print("Main linkage components")
    force = tau/rad*FoS
    print("yield area tension, mm^2:", (force/mat.yield_stress)*1000**2)
    print("yield area shear, mm^2:", (force/mat.yield_stress*np.sqrt(3))*1000**2)
    length = 44 # mm
    print("buckling limit required area at %d mm length: "%length, np.sqrt( force /(np.pi**2 * mat.E) *12 )*(length*1e-3)*1000**2)

    print("pulley support componenets")
    force = tau/rad*2/5*FoS
    print("yield area tension, mm^2:", (force/mat.yield_stress)*1000**2)
    length = 20 # mm
    print("buckling limit required area at %d mm length: "%length, np.sqrt( force /(np.pi**2 * mat.E) *12 )*(length*1e-3)*1000**2)


def main():
    # 6 January 2023, math for Chris/Elliott dorsiflexion spring
    torque = 2 # Nm
    angle = 45 * np.pi/180 # rad
    # local minima of (x+1)^2/x = 4 at x=1
    energy = torque*angle*2
    print(f"{energy=} Joules")
    FOS=1.3
    mat = SS410()
    volume = energy * 3 * mat.E / (mat.yield_stress**2) * FOS**2
    print(f"{volume=} m^3")
    mass = mat.rho*volume
    print(f"{mass=} kg")

    area = volume/4e-3
    density_factor=0.25
    disk_area = 1/density_factor *100*100* area
    print(f"{density_factor=}")
    print(f"{disk_area=} cm^2" )
    radius = np.sqrt(disk_area/np.pi)
    print(f"{radius=} cm")








if __name__ == '__main__':
    main()