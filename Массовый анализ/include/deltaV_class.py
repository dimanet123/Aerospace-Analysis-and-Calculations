import numpy as np

class Calculate_Impulse:
    def __init__(self, planet):
        self.planet = planet
        
    def change_orbits(self, h_start, h_end, i_start, i_end, meters = False):
        R = self.planet.R
        mu = self.planet.mu
        if meters == True:
            pass
        else:
            h_start = h_start * 1000
            h_end = h_end * 1000
        
        r1 = R + h_start
        r2 = R + h_end
        i1 = np.deg2rad(i_start)
        i2 = np.deg2rad(i_end)
        alpha = i2 - i1
        
        sumDeltaV = np.sqrt(mu/r1) * (np.sqrt(2*r2/(r2 + r1)) - 1) + np.sqrt((np.sqrt(mu/r2)*(np.cos(alpha) - np.sqrt(2*r1/(r1+r2))))**2 + (np.sqrt(mu/r2)*np.sin(alpha))**2)
        return sumDeltaV