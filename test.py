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
# s.solve(timesteps=40000*5*2)
s.solve(timesteps=1)
# print([p.q for p in s.grid])
# print([p.velocity for p in s.grid])

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
    [p.velocity * s.system.c_s_0 / s.v_esc for p in s.grid],
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
