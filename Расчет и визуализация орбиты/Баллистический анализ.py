import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import math
from Planet_class import Planet
from Spacecraft_class import Spacecraft
from Calculate_class import Calculate
from Visulise_class import VisualizeTrajectory
from Transition_class import Orbit_Transitions

Earth = Planet(mu = 398616, eps = 26327850000, R = 6378.137, We = 7.29212e-5, sol_activity = 150, sol_activity_bool=True,  meters = False)
Condor = Spacecraft(M_constr = 11000, M_fuel=100, I = 2000, mainF = 17, S = 4, Cd = 1.4)

orbits = Calculate(Earth, Condor)
# # orbits.sun_sync(500000, 243, 97.26690239412387)
time_vis = 3600 * 24 * 365
a = orbits.orbit(h = 513123.08798615524, i = 97.40350030687841,  time_vis = time_vis, one_lap = False, Omega = 0, omega = 0)
b = orbits.orbit(h = 513123.08798615524, i = 97.40350030687841, time_vis = time_vis, one_lap = True, Omega = 160.1, omega = 8.88)
lables = ['Первый спутник','Второй спутник']
list_dataframes = [a,b]
vis = VisualizeTrajectory(data = a, labels = lables, data_list = list_dataframes, planet = Earth)
vis.plot_trajectory(scale = 1000)

# b = orbits.orbit(h = 513123.08798615524, i = 97.40350030687841,  time_vis = time_vis, one_lap = False, Omega = 160.1, omega = 8.88)
# list_dataframes = [a,b]
# lables = ['Первый спутник','Второй спутник']
# vis = VisualizeTrajectory(data = a, labels = lables, data_list = list_dataframes, planet = Earth)
# vis.plot_trajectory_over_planet()
# result = Calculate.find_omega_for_target_position(Earth, Condor, 518000, 97.41405945929156, 8.88, 0, 0, 10000)
# print(result)
# print(orbits.zamik(500000, 243, 97.436)/1000)
# a = orbits.orbit(h = 500000, i = 97.26690239412387)
# orbits.sun_sync(512000, 243, 97)
# vis = VisualizeTrajectory(a)
# vis.plot_trajectory()

# a = Orbit_Transitions(planet = Earth, spacecraft = Condor, h_start = 518000, angle_transition = 360, n_vit = 683)
# h = a.calculate_new_height()
# deltaV = 2 * a.deltaV(518000, h)

# fuel = a.get_fuel(deltaV)
# time = a.get_time(deltaV/2)

