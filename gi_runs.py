import sys

from src import solver
from src import initial_conditions as ic
from src.composition import get_mean_molecular_mass

from math import sqrt
import matplotlib.pyplot as plt


"""
 - We will define mass loss as the fraction of the atmosphere whose particle velocity exceeds the local the 
    escape velocity.
 - The density of the atmosphere will be the mass of the BSE atmosphere divided by the number of moles of 
    the atmosphere.
 - The specific heat of the BSE atmosphere is ???
 - 
"""

# BSE composition oxide wt%, Visccher & Fegley 2013
bse_composition = {
    "SiO2": 45.40,
    'MgO': 36.76,
    'Al2O3': 4.48,
    'TiO2': 0.21,
    'Fe2O3': 0.00000,
    'FeO': 8.10,
    'CaO': 3.65,
    'Na2O': 0.349,
    'K2O': 0.031,
    'ZnO': 6.7e-3,
}

# mean molecular mass of the BSE composition
mean_molecular_mass = get_mean_molecular_mass(bse_composition) / 1000  # g/mol -> kg/mol

num_shells = 2000
r_0 = 6.4e6  # initial planet radius
mass_planet = 5.972e24  # mass earth
gamma_a = 1.4  # polytropic exponent
gamma = 1.4  # specific heat
M_a = 0.029  # kg/mol  # molecule mass, will likely have to be heavier for BSE atmosphere
T_0 = 288  # K, initial temperature (surface temperature)
P_0 = 1.01e5  # Pa, initial pressure (sea level)

s = solver.LagrangianSolver1D(
    num_shells=num_shells,
    P_0=P_0,
    T_0=T_0,
    gamma=gamma,
    gamma_a=gamma_a,
    m_a=M_a,
    r_0=r_0,
    mass_planet=mass_planet,
)
print(s.system.rho_0, s.system.c_s_0, s.system.vesc)
s.solve(timesteps=int(5e6))
