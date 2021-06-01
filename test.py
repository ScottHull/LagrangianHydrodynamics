from src import solver
from src import initial_conditions as ic

from math import sqrt
import matplotlib.pyplot as plt

num_shells = 2000
r_0 = 6.4e6  # initial planet radius
mass_planet = 5.972e24  # mass earth
gamma_a = 1.4  # polytropic exponent
gamma = 1.4  # specific heat
M_a = 0.029  # g/mol  # molecule mass
T_0 = 288
P_0 = 1.01e5

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
s.solve(timesteps=40000 * 5 * 2)
