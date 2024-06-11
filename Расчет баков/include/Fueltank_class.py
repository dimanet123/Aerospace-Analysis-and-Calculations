import numpy as np
import math


class FuelTank:
    def __init__(self, material, wall_thickness = 0):
        self.material = material
        self.wall_thickness = wall_thickness

    def surface_area(self):
        # Этот метод будет переопределен в подклассах
        pass

class SphericalTank(FuelTank):
    def __init__(self, material, c, fuel_volume, n_air, n_agr, P_out, n_zap):
        super().__init__(material)
        self.n_air = n_air
        self.n_agr = n_agr
        self.fuel_volume = fuel_volume
        self.P_out = P_out
        self.full_volume = fuel_volume * (1 + n_agr + n_air)
        self.radius = (3 * self.full_volume/4/np.pi)**(1/3)
        self.thikness = n_zap * P_out * n_zap * self.radius / self.material.sigma_v / 2 / c
        self.thikness = math.ceil(self.thikness * 1000)/1000
        self.tank_mass = self.material.density * 4/3 * np.pi * ((self.radius + self.thikness)**3 - self.radius**3)

    def surface_area(self):
        # Расчет площади поверхности сферы1
        return 4 * np.pi * self.radius ** 2
    
class CylinderTank(FuelTank):
    def __init__(self, material, c, fuel_volume, n_air, n_agr, P_out, n_zap, radius):
        super().__init__(material)
        self.n_air = n_air
        self.n_agr = n_agr
        self.fuel_volume = fuel_volume
        self.P_out = P_out
        self.full_volume = fuel_volume * (1 + n_agr + n_air)
        self.radius = radius
        self.length = (4 * self.full_volume)/(np.pi * 2 * self.radius)
        self.thikness = n_zap * P_out * self.radius / self.material.sigma_v / c
        self.thikness = math.ceil(self.thikness * 1000)/1000
        self.tank_mass = self.material.density * 4/3 * np.pi * ((self.radius + self.thikness)**3 - self.radius**3)
 
    def surface_area(self):
        # Расчет площади поверхности сферы1
        return 2 * np.pi * self.radius * self.length
