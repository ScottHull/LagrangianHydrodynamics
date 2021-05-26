from src import nondimensional
from src import initial_conditions as ic

import numpy as np
import matplotlib.pyplot as plt

r_0 = 6371e3  # initial planet radius
r_max = r_0 + (480 * 1000)  # height of the current Earth atmosphere, m
p_0 = 1.01e5  # initial planet surface pressure
rho_0 = 1.22  # initial planetary atmosphere at bottom
c_s_0 = 340  # initial sound velocity at bottom of atmosphere
lambda_0 = 756.57  # escape parameter initial condition
gamma_a = 1.4  # polytropic exponent
gamma = 1.4  # specific heat
M_a = 29  # g/mol  # molecule mass
M_atm = 4.10 * 10 ** 18  # kg, total mass of initial atmosphere

initial_conditions = []

inc = (r_max - r_0) / 1000
for r in np.arange(r_0, r_max + inc, inc):
    p_p0 = ic.pressure_initial(polytropic_exponent=gamma_a, lambda_0=lambda_0, radius=r, radius_0=r_0)
    rho_rho0 = ic.density_initial(polytropic_exponent=gamma_a, lambda_0=lambda_0, radius=r, radius_0=r_0)
    t_t0 = ic.temperature_initial(polytropic_exponent=gamma_a, lambda_0=lambda_0, radius=r, radius_0=r_0)
    p = p_p0 * p_0
    rho = rho_rho0 * rho_0
    t = t_t0 * t_t0
    initial_conditions.append((r / r_0, p_p0, rho_rho0, t_t0, p, rho, t))

fig = plt.figure(figsize=(16, 9))
fig.suptitle('Initial Conditions')
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)

ax1.plot(
    [i[0] for i in initial_conditions],
    [i[1] for i in initial_conditions],
    linewidth=2.0,
)
ax1.set_xlabel("r / r0")
ax1.set_ylabel("P / P0")
ax1.grid()

ax2.plot(
    [i[0] for i in initial_conditions],
    [i[2] for i in initial_conditions],
    linewidth=2.0,
)
ax2.set_xlabel("r / r0")
ax2.set_ylabel("rho / rho0")
ax2.grid()

ax3.plot(
    [i[0] for i in initial_conditions],
    [i[3] for i in initial_conditions],
    linewidth=2.0,
)
ax3.set_xlabel("r / r0")
ax3.set_ylabel("T / T0")
ax3.grid()

plt.show()
