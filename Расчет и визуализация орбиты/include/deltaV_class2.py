import numpy as np

class Calculate_Impulse:
    def __init__(self, planet):
        self.planet = planet
        
    def change_orbits(self, h_start, h_end, i_start, i_end, meters = False):
        
        mu = self.planet.mu
        R = self.planet.R
        
        # Вычисление начальной скорости на круговой орбите
        V0 = np.sqrt(mu/(R+h_start))

        # Вычисление радусов апогея и перигея
        R_apog = R + h_end
        R_perig = R + h_start

        # Вычисление эксцентриситета и фокального параметра
        e = (R_apog - R_perig)/(R_apog + R_perig)
        p = (2*R_apog*R_perig)/(R_apog + R_perig)

        # Вычисление скорсоти в перигее
        V_perig = np.sqrt((mu/p)*(e**2+2*e+1))

        # Вычисление изменения скорости для первого маневра
        deltaV1 = V_perig - V0
        
        # Вычисление скорости в апогее
        V_apog = np.sqrt((mu/p)*(e**2-2*e+1))
        
        # Вычисление конечной скорости на второй круговой орбите
        V1 = np.sqrt(mu/(R+h_end))
        
        # Вычисление прироста скорости для второго маневра
        deltaV2 = V1 - V_apog
        
        a = [deltaV1, deltaV2]
        
    
        return a