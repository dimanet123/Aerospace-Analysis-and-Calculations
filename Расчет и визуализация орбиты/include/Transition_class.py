import numpy as np

class Orbit_Transitions:
    def __init__(self, spacecraft, planet, h_start, angle_transition, n_vit):
        self.h_start = h_start
        self.angle_transition = angle_transition
        self.planet = planet
        self.n_vit = n_vit
        self.spacecraft = spacecraft
    
    def calculate_period(self, a):
        T = 2 * np.pi * np.sqrt(a**3 / self.planet.mu)
        return T
        
    def calculate_new_height(self):
        T1 = self.calculate_period(self.h_start + self.planet.R)
        
        T2 = (self.n_vit*T1 + T1 * self.angle_transition/360) / self.n_vit
        
        a_new = (T2 ** 2 / 4 / np.pi**2 * self.planet.mu)**(1/3)
        
        h_apog = 2 * a_new - self.h_start - 2 * self.planet.R
        
        return h_apog
    
    def deltaV(self, h_start, h_apog):
        # Вычисление начальной скорости на круговой орбите
        V0 = np.sqrt(self.planet.mu/(self.planet.R+h_start))

        # Вычисление радусов апогея и перигея
        R_apog = self.planet.R + h_apog
        R_perig = self.planet.R + h_start

        # Вычисление эксцентриситета и фокального параметра
        e = (R_apog - R_perig)/(R_apog + R_perig)
        p = (2*R_apog*R_perig)/(R_apog + R_perig)

        # Вычисление скорсоти в перигее
        V1 = np.sqrt((self.planet.mu/p)*(e**2+2*e+1))

        # Вычисление изменения скорости
        deltaV = V1 - V0
        
        return deltaV
    
    def get_fuel(self, deltaV):
        M = self.spacecraft.M_constr * (np.exp(deltaV/self.spacecraft.I) - 1)
        return M
    
    def get_time(self, deltaV):
        accelaration = self.spacecraft.mainF / self.spacecraft.M_constr
        t = deltaV / accelaration
        return t
        

