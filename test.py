from src import solver
from src import initial_conditions as ic

from math import sqrt
import matplotlib.pyplot as plt

num_shells = 2000
r_0 = 6.4e6  # initial planet radius
mass_planet = 5.972e24  # mass earth
gamma_a = 1.4  # polytropic exponent
gamma = 1.4  # specific heat
M_a = 0.029  # g/mol  # molecule mass, will likely have to be heavier for BSE atmosphere
T_0 = 288  # K, initial temperature (surface temperature)
P_0 = 1.01e5  # Pa, initial pressure (sea level)
u_s = None  # defaults to 0.5 * u_esc in Genda and Abe 2003
outfile_dir = "test_outputs"

s = solver.LagrangianSolver1D(
    num_shells=num_shells,
    P_0=P_0,
    T_0=T_0,
    gamma=gamma,
    gamma_a=gamma_a,
    m_a=M_a,
    r_0=r_0,
    mass_planet=mass_planet,
    u_s=u_s,
    outfile_dir=outfile_dir,
)
print(s.system.rho_0, s.system.c_s_0, s.system.vesc)
s.solve(timesteps=int(5e6))