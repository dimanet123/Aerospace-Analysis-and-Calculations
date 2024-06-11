import numpy as np

class Rocket_Block:
    def __init__(self, name, mass_constr, mass_fuel, I_DU, F_DU, List_PN):
        self.name = name
        self.mass_constr = mass_constr
        self.mass_fuel = mass_fuel
        self.I_DU = I_DU
        self.F_DU = F_DU
        self.List_PN = List_PN
    
    def calculate_fullMass(self):
        massPn = sum(item.mass for item in self.List_PN)
        all_mass = self.mass_constr + self.mass_fuel + massPn
        return all_mass
        
    def rocket_impulse(self, deltaV):
        M_1 = self.calculate_fullMass()
        M_2 = M_1 * np.exp(-deltaV / self.I_DU)
        deltaM = M_1 - M_2
        self.mass_fuel -= deltaM
        
    def relase_effective_load(self, name):
        self.List_PN = [item for item in self.List_PN if item.name != name]