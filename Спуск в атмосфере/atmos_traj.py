import numpy as np
import pandas as pnd
import scipy
import matplotlib.pyplot as plt


class Planet:
    def __init__(self, h_atm, mu, R, meters=True):
        if meters != True:
            self.R = R*1000
            self.mu = mu*1e9
            self.h_atm = h_atm*1000
        else:
            self.R = R
            self.mu = mu
            self.h_atm = h_atm

    def calc_atm(self, h):
        a = np.array([-6.3759, -7.3012, -1.1817])
        b = np.array([-0.4754, -0.0096, -0.0068, -0.0120, 0.0042])
        c = np.array([0.1803, 0.0872, -0.0153, 0.0145, 0])

        # Функция вычисления плотности воздуха
        # высота h задаётся в километрах
        # результат -- плотность кг/м3
        h = h / 1000
        x = h/50.0-1
        sa = a[0] + a[1]*x + a[2]*x*x
        sbc = np.sum((b[i]*np.cos((i+1)*np.pi*x) + c[i] *
                     np.sin((i+1)*np.pi*x) for i in range(5)))
        return np.exp(sa + sbc)

    def g_acc(self, h):
        # Вычисление ускорения свободного падения м/с^2 на высоте h (в метрах)
        return self.mu/(self.R+h)**2


class SpaceCraft:

    def __init__(self, planet, thermo, a, e, v, deltaV, gamma, S, Cx, K, initial_mass, meters=True, mass_loss = True):
        self.planet = planet
        self.thermo = thermo
        self.mass_loss = mass_loss
        if meters == True:
            scale = 1
        else:
            scale = 1000
        self.a = a * scale
        self.e = e
        self.v = np.deg2rad(v)
        self.deltaV = deltaV * scale
        self.gamma = np.deg2rad(gamma)
        self.S = S
        self.Cx = Cx
        self.K = K
        self.Cy = self.K * self.Cx
        self.initial_mass = initial_mass
        self.V_enter, self.theta_enter = self.calc_enter_params(self)
        self.Rn = 1.5 * np.sqrt(S/np.pi)

    def calc_enter_params(self, h_atm):
        planet, a, e, v, deltaV, gamma = self.planet, self.a, self.e, self.v, self.deltaV, self.gamma
        r2 = planet.R + planet.h_atm

        p = a*(1-e**2)
        r = p/(1+e*np.cos(v))

        # Скорость в точке схода
        V1 = np.sqrt(planet.mu*(2/r - 1/a))

        # Угол наклона скорости к местному горизонту
        theta1 = np.arctan((e*np.sin(v))/(1+e*np.cos(v)))

        # Вектор суммарной скорости при выдаче тормозного импульса
        V2 = np.sqrt(V1**2 + deltaV**2 + 2*V1*deltaV*np.cos(gamma - theta1))

        # Угол наклона вектора скорости схода к горизонту
        theta2 = np.arctan((deltaV*np.sin(gamma) + V1*np.sin(theta1)) /
                           (V1*np.cos(theta1) + deltaV*np.cos(gamma)))

        # Большая полуось орбиты спуска
        a2 = (planet.mu*r)/(2*planet.mu - V2**2*r)

        # Параметр орбиты
        p2 = (V2*r*np.cos(theta2))**2/planet.mu

        # Эксцентриситет орбиты
        e2 = np.sqrt(1 - p2/a2)

        # Истинная аномалия точки схода

        buffer = (p2 - r2) / (e2 * r2)
        if buffer > 1 or buffer < -1:
            raise ValueError('Нет входа в атмосферу')
        else:
            v2 = np.arccos((p2 - r2) / (e2 * r2))
            Vin = np.sqrt(planet.mu * (2 / r2 - 1 / a2))
            thetain = np.arccos((V2 * r * np.cos(theta2)) / (Vin * r2))

        return Vin, thetain

class calculations:
    def __init__(self, planet, end_height, list_spacecrafts):
        self.planet = planet
        self.list_spacecrafts = list_spacecrafts
        self.results = []
        self.end_height = end_height
    
    
    
    def calc_params(self):
        global time
        time = 0
        def dydt(t, y, self):
            global time
            # Функция правых частей дифференциальных уравнений
            # вычисление правых частей дифференциальных уравнений для
            # момента времени t и вектора состояния y = [r(t), v(t), theta(t)]
            
            # Радиус-вектор точки в момент времени t
            r = y[0]
            # Скорость  в момент времени t
            v = y[1]
            # Угол наклона траекториии  в момент времени t
            theta = y[2]
            R = self.planet.R
            # Высота (км)
            h = (r - R)
            if t > 100:
                True
            # Скоростной напор Н/м^2
            q = self.planet.calc_atm(h)*v*v/2
            # Ускорение свободного падения
            g = self.planet.g_acc(r - R)
            Cx = self.spacecraft.Cx
            Cy = self.spacecraft.Cy
            S = self.spacecraft.S
            mass = y[5]
            dm = 0
            Tsurf = y[6]
            if self.spacecraft.mass_loss == True:
                qkonv_res = self.spacecraft.thermo.qkonv(h, v, self.spacecraft.Rn)
                qrad_res = self.spacecraft.thermo.qrad(h, v, self.spacecraft.Rn)
                is_Active = self.spacecraft.thermo.is_Active(y[6])
                
                dm = is_Active * 1 *  (qkonv_res + qrad_res - self.spacecraft.thermo.eps * self.spacecraft.thermo.delta * self.spacecraft.thermo.Tp**4)/self.spacecraft.thermo.nu
                
            Tsurf = self.spacecraft.thermo.temp(h, v, self.spacecraft.Rn)
            mass = mass - dm*S
            Fx = q*Cx*S
            Fy = q*Cy*S
            # dv/dt =
            dv =  -Fx/mass  - g*np.sin(theta)
            # dtheta/dt =
            dtheta = Fy/(mass * v) - g*np.cos(theta)/v + v/r*np.cos(theta)
            # dr/dt =
            dr = v*np.sin(theta)

            dl = v*np.cos(theta) * R/(h+R)
            
            overloads_x = Fx / (mass * 9.81)
            overloads_y = Fy / (mass * 9.81)
            overloads = np.sqrt(overloads_x**2 + overloads_y**2)
            # time_overloads = (t)

            return (dr, dv, dtheta, dl, overloads - y[4], -dm * S, Tsurf - y[6])


        def event_h_eq_0(t, y):
            # Функция-"детектор", передаваемая в интегратор (параметр events),
            # для определения времени достижения нулевой высоты и
            # остановки процесса интегрирования
            # функция возвращает высоту
            return y[0]-self.planet.R - self.end_height
        
        for spacecraft_one in list_spacecrafts:
            self.spacecraft = spacecraft_one
            # функция определяется условие h = 0 при движении "вниз"
            event_h_eq_0.direction = -1
            # функция-детектор активна
            event_h_eq_0.terminal = True
            
            R = self.planet.R
            # Высота (км)
            # Скоростной напор Н/м^2
            # Ускорение свободного падения
            Cx = self.spacecraft.Cx
            Cy = self.spacecraft.Cy
            S = self.spacecraft.S
            mass = self.spacecraft.initial_mass
            Tsurf = self.spacecraft.thermo.temp(self.planet.h_atm, self.spacecraft.V_enter, self.spacecraft.Rn)
            
            time = 0
            sol = scipy.integrate.solve_ivp(lambda t, y: dydt(t, y, self), [0, 1000000], [
                                            R+self.planet.h_atm, self.spacecraft.V_enter, -self.spacecraft.theta_enter, 0, 0, mass,Tsurf], method='RK45', events=event_h_eq_0, rtol=1e-8)
            a = pnd.DataFrame(columns=['Altitude', 'Speed', 'Overloads', 'Theta', 'Longtitude', 'time', 'Mass', 'Tempreture'])
            a['Altitude'] = sol.y[0] - R
            a['Speed'] = sol.y[1]
            a['Theta'] = sol.y[2]
            a['Longtitude'] = sol.y[3]
            a['Overloads'] = sol.y[4]
            a['time'] = sol.t
            a['Mass'] = self.spacecraft.initial_mass - sol.y[5]
            a['Tempreture'] = sol.y[6]
            
            self.results.append(a)
            
class PlotResults:
    def __init__(self, list_of_dataframes, names=None):
        self.list_of_dataframes = list_of_dataframes
        self.names = names
    
    
        
    def plot_single(self, parameter_index):
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'black', 'magenta']
        # Creating a figure for plotting
        plt.figure(figsize=(10, 6))
        
        # Titles for different parameters based on the index provided
        titles = ['Altitude (m)', 'Speed (m/s)', 'Theta (rad)', 'Longitude (m)', 'Overloads (g)', 'Mass (kg)', 'Temperature (K)']
        # The actual data column names in your DataFrames
        data_columns = ['Altitude', 'Speed', 'Theta', 'Longtitude', 'Overloads', 'Mass', 'Tempreture']
        
        # Check if the parameter_index is within the valid range
        if parameter_index < 0 or parameter_index >= len(data_columns):
            raise ValueError("parameter_index out of range")
        
        # Parameter to plot
        parameter = data_columns[parameter_index]
        
        # Colors for the lines
        
        
        # Plot each spacecraft's specified parameter
        for idx, df in enumerate(self.list_of_dataframes):
            spacecraft_name = self.names[idx] if self.names else f'Spacecraft {idx + 1}'
            plt.plot(df['time'], df[parameter], label=spacecraft_name, color=colors[idx % len(colors)])
        
        # Setting title, labels, and grid
        plt.title(titles[parameter_index])  # Set the title using the titles list
        plt.xlabel('Time (s)')  # X-axis label
        plt.ylabel(titles[parameter_index])  # Y-axis label based on the parameter
        plt.legend()  # Display legend to show spacecraft names
        plt.grid(True)  # Enable grid for better readability

        # Display the plot
        plt.show()
        
    def plot(self):
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'black', 'magenta']
        # Setting up a figure with multiple subplots
        fig, axs = plt.subplots(3, 3, figsize=(14, 10))
        fig.subplots_adjust(hspace=0.5, wspace=0.4)
        titles = ['Altitude (m)', 'Speed (m/s)', 'Theta (rad)', 'Longitude (m)', 'Overloads (g)', 'Mass (kg)', 'Tempreture (K)']

        # Iterate through each DataFrame and plot each parameter in a subplot
        for idx, df in enumerate(self.list_of_dataframes):
            spacecraft_name = self.names[idx] if self.names else f'Spacecraft {idx + 1}'
            axs[0, 0].plot(df['time'], df['Altitude'], label=spacecraft_name, color=colors[idx % len(colors)])
            axs[0, 1].plot(df['time'], df['Speed'], label=spacecraft_name, color=colors[idx % len(colors)])
            axs[0, 2].plot(df['time'], df['Theta'], label=spacecraft_name, color=colors[idx % len(colors)])
            axs[1, 0].plot(df['time'], df['Longtitude'], label=spacecraft_name, color=colors[idx % len(colors)])
            axs[1, 1].plot(df['time'], df['Overloads'], label=spacecraft_name, color=colors[idx % len(colors)])
            axs[1, 2].plot(df['time'], df['Mass'], label=spacecraft_name, color=colors[idx % len(colors)])
            axs[2, 0].plot(df['time'], df['Tempreture'], label=spacecraft_name, color=colors[idx % len(colors)])
            
        axs[2, 0].axhline(y=1650, color='grey', linestyle='--', linewidth=1)
        # Set titles, labels, and legends for each subplot
        for ax, title in zip(axs.flat, titles):
            ax.set_title(title)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel(title)
            ax.legend()
            ax.grid(True)

        # Hide any unused subplots
        axs[2, 1].axis('off')
        axs[2, 2].axis('off')

        # Show the plot
        plt.show()
    
class Thermo:
    def __init__(self, planet, Tp, ro_tz, CP, lam, nu, Ckonv, Crad):
        self.planet = planet
        self.Tp = Tp
        self.ro_tz = ro_tz
        self.CP = CP
        self.lam = lam
        self.nu = nu
        self.Ckonv = Ckonv
        self.Crad = Crad
        self.eps = 0.95
        self.delta = 5.67e-8
        self.Vfk = 7910
    
    
    
    def qkonv(self, h, V, Rn):
        result = self.Ckonv/(Rn**(1/2)) * (self.planet.calc_atm(h)/self.planet.calc_atm(0))**(1/2) * (V/self.Vfk)**3.2 
        return result
    
    def qrad(self, h, V, Rn):
        result = self.Crad*Rn * (self.planet.calc_atm(h)/self.planet.calc_atm(0))**(1.168) * (V/self.Vfk)**7.4
        return result
        
    def temp(self, h, V, Rn):
        result = ((self.qkonv(h, V, Rn) + self.qrad(h, V, Rn))/(self.eps*self.delta))**0.25
        return result
    
    def is_Active(self, temp):
        if temp > self.Tp:
            return 1
        else:
            return 0
        
class Parachut:
    def __init__(self, planet, spacecraft, Vp, Cph, np, kh, hpack, V_start, f, h_start):
        self.spacecraft = spacecraft
        self.planet = planet
        self.Vp = Vp
        self.Cph = Cph
        self.np = np
        self.kh = kh
        self.hpack = hpack
        self.V_start = V_start
        self.f = f
        self.P_KA = spacecraft.initial_mass * 9.8
        self.h_start = h_start
        self.M_par = f*self.P_KA/9.8
        self.F_potr = 2 * (self.M_par*9.8 + self.P_KA)/(planet.calc_atm(0) * Cph * Vp**2)
        self.Rpmax = ((np + 1) * (self.M_par*9.8 + self.P_KA)) / (1 + (spacecraft.Cx * spacecraft.S)/(2 * Cph * self.F_potr))
        
    def allowed_speed(self, h):
        result = np.sqrt((self.Rpmax*(self.kh*self.Vp**2/(self.planet.calc_atm(h)/self.planet.calc_atm(0))+2*np.sqrt(self.F_potr)))/(self.planet.calc_atm(0)*self.kh*(self.Vp**2)*self.Cph*self.F_potr)-np.sqrt(self.F_potr)/self.kh)
        return result
    
Earth = Planet(h_atm=100, mu=398600, R=6371, meters=False)
PICA_X = Thermo(Earth, Tp = 1650, ro_tz = 1600, CP = 1.667e3, lam = 4.3937e-1, nu = 4168e3, Ckonv = 2.24e8, Crad = 0.435)
Condor1 = SpaceCraft(Earth, PICA_X, a=6896.137, e=0, v=0, deltaV=0.13,
                    gamma=180, S=8 , Cx=1.2, K=0, initial_mass=900, meters=False, mass_loss = False)
# Condor2 = SpaceCraft(Earth, PICA_X, a=7200, e=0.006, v=140, deltaV=0.36,
#                     gamma=120, S=5.4 , Cx=1.3, K=0.16, initial_mass=5100, meters=False, mass_loss = True)
# Condor2 = SpaceCraft(Earth, PICA_X, a=7200, e=0.006, v=140, deltaV=0.36,
#                     gamma=120, S=5.4 , Cx=1.3, K=0.16, initial_mass=4600, meters=False, mass_loss = False)
# Condor3 = SpaceCraft(Earth, PICA_X, a=7200, e=0.006, v=140, deltaV=0.36,
#                     gamma=120, S=5.4 , Cx=1.3, K=0.16, initial_mass=4100, meters=False, mass_loss = False)
# Condor4 = SpaceCraft(Earth, PICA_X, a=7200, e=0.006, v=140, deltaV=0.36,
#                     gamma=120, S=5.4 , Cx=1.3, K=0.16, initial_mass=3600, meters=False, mass_loss = False)
# Condor5 = SpaceCraft(Earth, PICA_X, a=7200, e=0.006, v=140, deltaV=0.36,
#                     gamma=120, S=5.4 , Cx=1.3, K=0.16, initial_mass=5600, meters=False, mass_loss = False)
# Condor6 = SpaceCraft(Earth, PICA_X, a=7200, e=0.006, v=140, deltaV=0.36,
#                     gamma=120, S=5.4 , Cx=1.3, K=0.16, initial_mass=6100, meters=False, mass_loss = False)
# Condor7 = SpaceCraft(Earth, PICA_X, a=7200, e=0.006, v=140, deltaV=0.36,
#                     gamma=120, S=5.4 , Cx=1.3, K=0.16, initial_mass=6700, meters=False, mass_loss = False)

list_spacecrafts = [Condor1]
Calculations_all = calculations(Earth, 0, list_spacecrafts)
Calculations_all.calc_params()
# speed = Calculations_all.results[0].iloc[-1]['Speed']

# osn_parashut = Parachut(planet = Earth, spacecraft = Condor1, Vp = 8, Cph = 0.9, np = 3, kh = 0.027, hpack = 10000, V_start = speed, f = 0.07, h_start = 10000)
# dop_speed = osn_parashut.allowed_speed(10000)
spacecraft_names = ['Спуск в атмосфере', 'С уносом массы']
list_spacecrafts_dataframes = [Calculations_all.results[i] for i in range(len(Calculations_all.results))]

plotting = PlotResults(list_spacecrafts_dataframes, names=spacecraft_names)
plotting.plot()