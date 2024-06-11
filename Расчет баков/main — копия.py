import numpy as np
import sys
from include import BaseCalculs_class
from include import Propellant_class
from include import Fueltank_class
from include import Material_class

M_fuel = 100
K = 1.85
n_prop = 1
n_oxid = 4
NDMG = Propellant_class.Propellant(density=793)
AT = Propellant_class.Propellant(density=1443)


alluminum = Material_class.Material(density = 2640, sigma_v = 310e6)
steel = Material_class.Material(density = 7800, sigma_v = 440e6)


First = BaseCalculs_class.BaseCalculs(M_fuel, K, n_prop, n_oxid, proppelant = NDMG, oxidiser = AT)
M_prop = First.oneTank_propMass()
M_oxid = First.oneTank_oxidMass()

V_prop = First.oneTank_propVolume()
V_oxid = First.oneTank_oxidVolume()

P_out = 1.03e6

propTank = Fueltank_class.SphericalTank(material = steel, n_air = 0.06, n_agr = 0.01, fuel_volume = V_prop, P_out = P_out, n_zap = 1.5, c = 0.85)
oxidTank = Fueltank_class.SphericalTank(material = steel, n_air = 0.06, n_agr = 0.01, fuel_volume = V_oxid, P_out = P_out, n_zap = 1.5, c = 0.85)

propTank.radius
oxidTank.radius

propTank.thikness
oxidTank.thikness

Full_volume = n_prop * propTank.full_volume + n_oxid * oxidTank.full_volume


P_sp = 10e6
V = P_out* Full_volume/(P_sp - P_out)

nadduvTank = Fueltank_class.SphericalTank(material = steel, n_air = 0, n_agr = 0, fuel_volume = V, P_out = P_sp, n_zap = 1.5, c = 0.85)

TotalMass = propTank.tank_mass * n_prop + propTank.tank_mass * n_oxid + nadduvTank.tank_mass

def output_info():
    print(f'Радиус бака топлива - {propTank.radius * 1000} мм')
    print(f'Радиус бака окислителя - {oxidTank.radius * 1000} мм')
    print(f'Радиус бака наддува - {nadduvTank.radius * 1000} мм')
    
output_info()