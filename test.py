from src import solver
from src import initial_conditions as ic

from math import sqrt
import matplotlib.pyplot as plt

num_shells = 2000
r_0 = 6371e3  # initial planet radius
mass_planet = 5.972 * 10 ** 24  # mass earth
r_max = r_0 + (480 * 1000)  # height of the current Earth atmosphere, m
p_0 = 92.2  # initial planet surface pressure
rho_0 = 1e-4  # initial planetary atmosphere at bottom# lambda_0 = 756.57  # escape parameter initial condition
gamma_a = 1.31  # polytropic exponent
gamma = 1.31  # specific heat
M_a = 0.018  # g/mol  # molecule mass
M_atm = 4.10 * 10 ** 18  # kg, total mass of initial atmosphere
T_0 = 2000
c_s_0 = sqrt(gamma_a * p_0 / rho_0)  # initial sound velocity at bottom of atmosphere
t0 = r_0 / c_s_0
dt = 0.01 / t0

inner_boundary_velocity = 0
inner_boundary_mass = 0
outer_boundary_pressure = 0
outer_boundary_density = 0

lambda_0 = ic.lambda_0_initial(mass_planet=mass_planet, P_0=p_0, r_0=r_0, rho_0=rho_0)
s = solver.LagrangianSolver1D(
    num_shells=num_shells,
    P_max=outer_boundary_pressure,
    rho_max=outer_boundary_density,
    P_0=p_0,
    rho_0=rho_0,
    T_0=T_0,
    gamma=gamma,
    gamma_a=gamma_a,
    lambda_0=lambda_0,
    m_a=M_a,
    m_initial=inner_boundary_mass,
    r_0=r_0,
    v_0=inner_boundary_velocity,
    timestep=dt,
    c_s_0=c_s_0,
    mass_planet=mass_planet
)
s.solve(timesteps=40000*5*2)
