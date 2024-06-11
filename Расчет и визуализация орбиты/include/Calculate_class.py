import numpy as np
import scipy
import pandas as pnd
from deltaV_class2 import Calculate_Impulse
# from scipy.optimize import minimize

class Calculate:
    def __init__(self, planet, spacecraft):
        self.planet = planet
        self.spacecraft = spacecraft
        self.results = 0
        
    # def find_omega_for_target_position(planet, spacecraft, h, i, Omega, target_latitude, target_longitude, time_target):
    #     """
    #     Находит значение omega, для которого спутник достигает заданной широты и долготы через определённое время.
        
    #     Аргументы:
    #     planet - объект планеты
    #     spacecraft - объект космического аппарата
    #     h - высота орбиты в км
    #     i - наклонение орбиты в градусах
    #     Omega - долгота восходящего узла в градусах
    #     target_latitude - целевая широта в градусах
    #     target_longitude - целевая долгота в градусах
    #     time_target - время в секундах, через которое должна быть достигнута цель
        
    #     Возвращает:
    #     Результат оптимизации, включая оптимальное значение omega.
    #     """

    #     def objective(omega):
    #         # Используем класс Calculate для получения траектории
    #         orbits = Calculate(planet, spacecraft)
    #         trajectory = orbits.orbit(h, i, omega[0], Omega, time_vis=time_target, one_lap=False)
    #         # Находим позицию спутника в заданное время
    #         idx = min(range(len(trajectory['time'])), key=lambda x: abs(trajectory['time'][x] - time_target))
    #         # Вычисляем разницу с целевыми координатами
    #         latitude_error = abs(trajectory['latitude'][idx] - target_latitude)
    #         longitude_error = abs(trajectory['longitude'][idx] - target_longitude)
    #         return latitude_error + longitude_error

    #     # Начальное предположение для omega
    #     omega_initial = [0]  
    #     result = minimize(objective, omega_initial, method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-6, 'disp': True})

    #     return result
    
    
    
    def orbit(self, h, i, omega = 0, Omega = 0, time_vis = 1000000, one_lap = True):
        
        global manuvr_var
        global sum_imp
        global h_now
        
        manuvr_var = 0
        sum_imp = [0,0]
        h_now = 0
        
        def initial_conditions(self, h, omega, i, Omega):
            # Constants
            mu = self.planet.mu  # Gravitational parameter of Earth
            R = self.planet.R    # Radius of Earth

            # Convert degrees to radians for trigonometric functions
            omega_rad = np.deg2rad(omega)
            i_rad = np.deg2rad(i)
            Omega_rad = np.deg2rad(Omega)

            # Semi-major axis a from altitude
            a = R + h

            # Orbital velocity at periapsis
            v = np.sqrt(mu / a)

            # Position in the orbital plane at periapsis
            r_pqw = np.array([a, 0, 0])

            # Velocity in the orbital plane at periapsis
            v_pqw = np.array([0, v, 0])

            # Rotation matrices for coordinate transformation
            R3_Omega = np.array([
                [np.cos(Omega_rad), -np.sin(Omega_rad), 0],
                [np.sin(Omega_rad), np.cos(Omega_rad), 0],
                [0, 0, 1]
            ])

            R1_i = np.array([
                [1, 0, 0],
                [0, np.cos(i_rad), -np.sin(i_rad)],
                [0, np.sin(i_rad), np.cos(i_rad)]
            ])

            R3_omega = np.array([
                [np.cos(omega_rad), -np.sin(omega_rad), 0],
                [np.sin(omega_rad), np.cos(omega_rad), 0],
                [0, 0, 1]
            ])

            # Complete rotation matrix from PQW to ECI
            R_eci = R3_omega @ R1_i @ R3_Omega  # Corrected order of multiplication

            # Position and velocity in ECI frame
            r_eci = R_eci @ r_pqw
            v_eci = R_eci @ v_pqw
            initial = [r_eci[0], v_eci[0], r_eci[1], v_eci[1], r_eci[2], v_eci[2]]

            return initial
        
        def dydt(t, y, self):
            global manuvr_var, sum_imp, h_now
            x, vx, y, vy, z, vz = y
            eps = self.planet.eps
            mu = self.planet.mu
            r = np.sqrt(x**2 + y**2 + z**2)
            v = np.sqrt(vx**2 + vy**2 + vz**2)
            ax = -mu * x / r**3 - eps * (x / r**5 - 5 * x * z**2 / r**7)
            ay = -mu * y / r**3 - eps * (y / r**5 - 5 * y * z**2 / r**7)
            az = -mu * z / r**3 - eps * (3 * z / r**5 - 5 * z**3 / r**7)
            
            altitude = r - self.planet.R
            
            rho = self.planet.calc_atm(altitude)
            Cd = self.spacecraft.Cd  
            S = self.spacecraft.S
            
            drag_coefficient = 0.5 * rho * v**2 * Cd * S / self.spacecraft.M_constr
            
            ax -= drag_coefficient * vx / v
            ay -= drag_coefficient * vy / v
            az -= drag_coefficient * vz / v
            
            
            
            return [vx, ax, vy, ay, vz, az]
        
        def event_h_eq_0(t, y):
            
            if t > 100:
                return y[4]  # This is the z-coordinate
            else:
                return 1  # Return non-zero to avoid triggering the event
        
        
        R = self.planet.R
        mu = self.planet.mu
        
        conds = initial_conditions(self, h, omega, i, Omega)
        
        event_h_eq_0.direction = 0
        # функция-детектор активна
        event_h_eq_0.terminal = True
        v = np.sqrt(mu / (R + h))
        initial_conditions = conds
        if one_lap:
            sol = scipy.integrate.solve_ivp(lambda t, param: dydt(t, param, self), [0, time_vis], initial_conditions, method='LSODA', events=event_h_eq_0, rtol=1e-8, max_step=1000)
        else:
            sol = scipy.integrate.solve_ivp(lambda t, param: dydt(t, param, self), [0, time_vis], initial_conditions, method='LSODA', rtol=1e-8, max_step=1000)
        
        a = pnd.DataFrame(columns=['x', 'vx', 'y', 'vy', 'z', 'vz', 'v'])
        a['x'] = sol.y[0]
        a['vx'] = sol.y[1]
        a['y'] = sol.y[2]
        a['vy'] = sol.y[3]
        a['z'] = sol.y[4]
        a['vz'] = sol.y[5]
        r = np.sqrt(sol.y[0] **2 + sol.y[2] **2 + sol.y[4] **2)
        a['v'] = np.sqrt(sol.y[1] **2 + sol.y[3] **2 + sol.y[5] **2)
        a['alt'] = r - R
        a['time'] = sol.t 
        a['latitude'], a['longitude'] = zip(*[self.cartesian_to_geographic(x, y, z, time) for x, y, z, time in zip(a['x'], a['y'], a['z'], a['time'])])
        
        a['longitude'] = a['longitude'].mask(a['longitude'] > 180, a['longitude'] - 360)

        # Коррекция значений долготы, которые ниже -180 градусов
        a['longitude'] = a['longitude'].mask(a['longitude'] < -180, a['longitude'] + 360)
        return a
    
    def cartesian_to_geographic(self, x, y, z, time):
        latitude = np.arctan2(z , np.sqrt(x**2 + y**2))  # Широта
        delta = time * self.planet.We
        delta_cell = delta // (2 * np.pi)
        longitude = np.arctan2(y , x) - (delta - delta_cell*2 *np.pi) # Долгота
        return np.rad2deg(latitude), np.rad2deg(longitude)
    
    def get_delta_lam(self, h, i):
        a = self.orbit(h, i)
        y = a.iloc[-1]['y']
        x = a.iloc[-1]['x']
        t = a.iloc[-1]['time']
        
        return np.rad2deg(np.arctan(y/x)) - t * np.rad2deg(self.planet.We)
    
    def zamik(self, h, turns, i):
        a = self.orbit(h, i)
        T = a.iloc[-1]['time']*turns
        cac = round(T/24/3600)
        delta_h = 1000
        h3 = h
        while abs((turns*(-1*self.get_delta_lam(h3, i))/360.0)-cac) > 1e-12:
            a = self.orbit(h, i)
            T = a.iloc[-1]['time']*turns
            cac = round(T/24/3600)
            h2 = h+delta_h
            a = -self.get_delta_lam(h2, i)
            b = -self.get_delta_lam(h, i)
            h3 = h2 - (a - cac*360.0/turns)*(delta_h)/(a-b)
            delta_h = h3-h2
            h = h2
        return h3
    
    def find_T(self, h,i):
        a = self.orbit(h, i)
        T = a.iloc[-1]['time']
        return T
    
    def sun_sync(self, h, turn, i):
        h1 = self.zamik(h, turn, i)
        T1 = self.find_T(h1,i)*turn/(3600*24)
        N1 = T1 - round(T1)
        i2 = i + 0.5
        h1 = self.zamik(h, turn, i2)
        T2 = self.find_T(h1,i)*turn/(3600*24)*turn/(3600*24)
        N2 = T2 - round(T2)
        delta = (N2 - N1)/(i2-i)
        i3 = i2 - N2/delta
        h1 = self.zamik(h, turn, i3)
        T3 = self.find_T(h1,i)*turn/(3600*24)
        N3 = T3 - round(T3)
        while abs(N3) > 1e-12:
            delta = (N3 - N2)/(i3-i2)
            i2 = i3
            i3 = i3 - N3/delta
            T2 = T3
            h1 = self.zamik(h, turn, i3)
            T3 = self.find_T(h1,i)*turn/(3600*24)
            N2 = N3
            N3 = T3 - round(T3)
        print(self.zamik(h, turn, i3))
        print(i3)
        print(T3)
            
    

