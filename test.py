from src import solver
from src import initial_conditions as ic

from math import sqrt
import matplotlib.pyplot as plt

num_shells = 2000
r_0 = 3 * 3.37e6  # initial planet radius
mass_planet = 5.972 * 10 ** 24  # mass earth
# r_max = r_0 + (480 * 1000)  # height of the current Earth atmosphere, m
p_0 = 92.2  # initial planet surface pressure
rho_0 = 1e-4  # initial planetary atmosphere at bottom# lambda_0 = 756.57  # escape parameter initial condition
gamma_a = 1.31  # polytropic exponent
gamma = 0.75  # specific heat
M_a = 0.018  # g/mol  # molecule mass
M_atm = 4.10 * 10 ** 18  # kg, total mass of initial atmosphere
T_0 = 2000

inner_boundary_velocity = 0
inner_boundary_mass = 0
outer_boundary_pressure = 0
outer_boundary_density = 0

lambda_0 = ic.lambda_0_initial(mass_planet=mass_planet, P_0=p_0, r_0=r_0, rho_0=rho_0)
c_s_0 = sqrt((gamma_a * p_0) / rho_0)  # initial sound velocity at bottom of atmosphere
t0 = r_0 / c_s_0
dt = 0.01 / t0

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
    mass_planet=mass_planet,
)
# s.solve(timesteps=40000*5*2)
s.solve(timesteps=2)

fig = plt.figure(figsize=(16, 9))
ax_pressure = fig.add_subplot(221)
ax_pressure.plot(
    [p.radius for p in s.grid],
    [p.pressure for p in s.grid],
    linewidth=2.0,
    color='black'
)
ax_pressure.set_xlabel("r / r_0")
ax_pressure.set_ylabel("P / P_0")
ax_pressure.set_title("Pressure (IC)")
ax_pressure.grid()

ax_mass = fig.add_subplot(222)
ax_mass.plot(
    [p.radius for p in s.grid],
    [p.mass for p in s.grid],
    linewidth=2.0,
    color='black'
)
ax_mass.set_xlabel("r / r_0")
ax_mass.set_ylabel("m / (r_0^3 rho_0)")
ax_mass.set_title("Mass (IC)")
ax_mass.grid()

ax_velocity = fig.add_subplot(223)
ax_velocity.plot(
    [p.radius for p in s.grid],
    [p.velocity * c_s_0 / s.v_esc for p in s.grid],
    linewidth=2.0,
    color='black'
)
ax_velocity.set_xlabel("r / r_0")
ax_velocity.set_ylabel("v / v_esc")
ax_velocity.set_title("Velocity (IC)")
ax_velocity.grid()

ax_density = fig.add_subplot(224)
ax_density.plot(
    [p.radius for p in s.grid],
    [p.density for p in s.grid],
    linewidth=2.0,
    color='black'
)
ax_density.set_xlabel("r / r_0")
ax_density.set_ylabel("rho / rho_0")
ax_density.set_title("Density (IC)")
ax_density.grid()

plt.show()
