import numpy as np

class BaseCalculs:
    def __init__(self, M_fuel, K, n_prop, n_oxid, proppelant, oxidiser):
        self.M_fuel = M_fuel
        self.K = K
        self.n_prop = n_prop
        self.n_oxid = n_oxid
        self.M_prop = self.prop_mass()
        self.M_oxid = self.oxid_mass()
        self.proppelant = proppelant
        self.oxidiser = oxidiser
        self.V_prop = self.prop_volume()
        self.V_oxid = self.oxid_volume()
        
    def prop_mass(self):
        M_prop = self.M_fuel / (self.K + 1)
        return M_prop
    
    def oxid_mass(self):
        M_oxid = self.M_fuel/(1 + 1/self.K)
        return M_oxid
    
    def oneTank_propMass(self):
        return self.M_prop / self.n_prop
    
    def oneTank_oxidMass(self):
        return self.M_oxid / self.n_oxid
    
    def prop_volume(self):
        return self.M_prop / self.proppelant.density
    
    def oxid_volume(self):
        return self.M_oxid / self.oxidiser.density
    
    def oneTank_propVolume(self):
        return self.V_prop / self.n_prop
    
    def oneTank_oxidVolume(self):
        return self.V_oxid / self.n_oxid
    
    
    
