import numpy as np
from include import Effective_load_class
from include import Planet_class
from include import Rocket_Block_class
from include import deltaV_class

Earth = Planet_class.Planet(mu = 398616, eps = 26327850000, R = 6378.137, We = 7.29212e-5, sol_activity = 0, sol_activity_bool=False,  meters = False)

Calculator = deltaV_class.Calculate_Impulse(Earth)

Condor1 = Effective_load_class.Effective_load('Condor1', 2050 + 370)
listPN = [Condor1]

Fregat = Rocket_Block_class.Rocket_Block('Fregat', mass_constr = 950, mass_fuel = 4600, I_DU = 320 * 9.8, F_DU = 2030, List_PN=listPN)


if Fregat.calculate_fullMass() > 8000:
    print(f"Не летишь сегодня - {Fregat.calculate_fullMass()}")

print(Fregat.mass_fuel)
Fregat.rocket_impulse(2600)
print(Fregat.mass_fuel)



