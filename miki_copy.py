from math import pi, sqrt, exp
import matplotlib.pyplot as plt
import sys

n = 2001
m = 100

# parameters
gamma = 1.4
gamma_a = 1.4
Mp = 5.972e24
G = 6.67408e-11
Ma = 0.029
Rgas = 8.3
temp0 = 288
p0 = 1.01e5
rho0 = p0 * Ma / temp0 / Rgas
r0 = 6.4e6
cq = 0.75

lambda0 = sqrt((2 * G * Mp) / r0)
c0 = sqrt(gamma_a * p0 / rho0)
vesc = sqrt(2 * G * Mp / r0)
us = 0.5 * vesc
rlast = r0 / (1 - gamma_a / (gamma_a - 1) / lambda0)
t0 = r0 / c0

rr = [0 for i in range(0, n)]
rho = [0 for i in range(0, n)]
uu = [0 for i in range(0, n)]
pp = [0 for i in range(0, n)]
qq = [0 for i in range(0, n)]
mm = [0 for i in range(0, n)]
temp = [0 for i in range(0, n)]

rr_next = rr
rho_next = rho
uu_next = uu
pp_next = pp
temp_next = temp

# initial conditions
for i in range(0, n):
    rr[i] = r0 + (((rlast - r0) / n) * i)
    rho[i] = rho0 * ((gamma_a - 1) / gamma_a * lambda0 * (r0 / rr[i] - 1) + 1) ** (1 / (gamma_a - 1))
    pp[i] = p0 * ((gamma_a - 1) / gamma_a * lambda0 * (r0 / rr[i] - 1) + 1) ** (gamma_a / (gamma_a - 1))
    if i > 0:
        mm[i] = mm[i - 1] + rho[i] * 4 * pi / 3 * (rr[i] ** 3 - rr[i - 1] ** 3)
    # ideal gas
    temp[i] = pp[i] / rho[i] * Ma / Rgas
    # uu[i] = (sqrt(gamma_a * Rgas / Ma * temp[i]))
    uu[i] = 0

# inner boundary conditions
# uu[0] = us
pp[n - 1] = 0
rho[n - 1] = 0

# normalize
for i in range(0, n):
    uu[i] = uu[i] / c0
    rr[i] = rr[i] / r0
    rho[i] = rho[i] / rho0
    pp[i] = pp[i] / p0
    mm[i] = mm[i] / (rho0 * r0 ** 3)
    temp[i] = temp[i] / temp0

fig = plt.figure(figsize=(16, 9))
ax_pressure = fig.add_subplot(221)
ax_density = fig.add_subplot(222)
ax_velocity = fig.add_subplot(223)
ax_mass = fig.add_subplot(224)
ax_pressure.plot(
    rr,
    pp,
    linewidth=2.0,
    color='black'
)
ax_density.plot(
    rr,
    rho,
    linewidth=2.0,
    color='black'
)
ax_velocity.plot(
    rr,
    [i * c0 / vesc for i in uu],
    linewidth=2.0,
    color='black'
)
ax_mass.plot(
    rr,
    mm,
    linewidth=2.0,
    color='black'
)
ax_pressure.set_xlabel("r / r0")
ax_density.set_xlabel("r / r0")
ax_velocity.set_xlabel("r / r0")
ax_mass.set_xlabel("r / r0")
ax_pressure.set_ylabel("P / P0")
ax_density.set_ylabel("rho / rho0")
ax_velocity.set_ylabel("u / u_esc")
ax_mass.set_ylabel("m / (r_0^3 rho_0)")
ax_pressure.set_title("Pressure (IC)")
ax_density.set_title("Density (IC)")
ax_velocity.set_title("Velocity (IC)")
ax_mass.set_title("Mass (IC)")
ax_pressure.grid()
ax_density.grid()
ax_velocity.grid()
ax_mass.grid()

# plt.show()

# time-dependent part
dt = 0.001 / t0
tt = 0
nt = 1000
uu[0] = - (lambda0 / gamma / (rr[0] ** 2) * dt)

for l in range(0, nt):
    print("At time: ", tt)
    for i in range(0, n - 1):
        if uu[i + 1] < uu[i]:
            qq[i] = -cq * gamma * rho[i] * (uu[i + 1] - uu[i]) * (
                    sqrt(pp[i] / rho[i]) - (((gamma + 1) / 2) * (uu[i + 1] - uu[i])))
        else:
            qq[i] = 0
    qq[n - 1] = 0
    uu_next[0] = uu[0] - (lambda0 / gamma / (rr[0] ** 2) * dt)  # boundary condition
    for i in range(1, n - 1):
        # a1 = (8 * pi / gamma) * (rr[i] ** 2) * ((pp[i] - pp[i - 1] + qq[i] - qq[i - 1]) / (mm[i] - mm[i - 1]))
        # a2 = gamma_a / (gamma * (rr[i] ** 2))
        # uu_next[i] = uu[i] - (a1 * a2 * dt)
        uu_next[i] = uu[i] - ((8 * pi) / gamma * rr[i] ** 2 * (pp[i] - pp[i - 1] + qq[i] - qq[i - 1]) / (
                mm[i + 1] - mm[i]) + lambda0 / gamma / rr[i] ** 2) * dt
    uu_next[n - 1] = uu_next[n - 2]  # free-slip

    for i in range(0, n):
        rr_next[i] = rr[i] + (uu_next[i] * dt)

    for i in range(0, n - 1):
    #     a1 = 4 * pi * rho[i]
    #     a2 = (gamma * pp[i]) + ((gamma - 1) * qq[i])
    #     a3 = (((rr[i + 1] ** 2 * uu[i + 1]) - (rr[i] ** 2 * uu[i]))) * dt
    #     pp_next[i] = pp[i] - (a1 * a2 * a3)
        pp_next[i] = pp[i] - (4 * pi) * rho[i] * (gamma * pp[i] + (gamma - 1) * qq[i]) * (
                    rr[i + 1] ** 2 * uu[i + 1] - rr[i] ** 2 * uu[i]) / (mm[i + 1] - mm[i]) * dt
        rho_next[i] = (3 / (4 * pi)) * (mm[i + 1] - mm[i]) / ((rr_next[i + 1] ** 3) - (rr_next[i] ** 3))
    for i in range(0, n - 1):
        # ideal gas
        temp_next[i] = p0 * pp_next[i] / (rho0 * rho_next[i]) * Ma / Rgas / temp0

    # outer boundary condition
    pp_next[n - 1] = 0
    rho_next[n - 1] = 0

    Mloss = 0
    for i in range(0, n):
        criterion = uu[i] * c0 / sqrt(2 * G * Mp / (r0 * rr[i]))
        if criterion > 1:
            Mloss = mm[i] / mm[-1]

    tt += dt

    pp = pp_next
    rho = rho_next
    temp = temp_next
    uu = uu_next
    rr = rr_next

    print(rr)
    print(uu)

    fig = plt.figure(figsize=(16, 9))
    ax_pressure = fig.add_subplot(221)
    ax_density = fig.add_subplot(222)
    ax_velocity = fig.add_subplot(223)
    ax_mass = fig.add_subplot(224)
    ax_pressure.plot(
        rr,
        pp,
        linewidth=2.0,
        color='black'
    )
    ax_density.plot(
        rr,
        rho,
        linewidth=2.0,
        color='black'
    )
    ax_velocity.plot(
        rr,
        [i * c0 / vesc for i in uu],
        linewidth=2.0,
        color='black'
    )
    ax_mass.plot(
        rr,
        mm,
        linewidth=2.0,
        color='black'
    )
    ax_pressure.set_xlabel("r / r0")
    ax_density.set_xlabel("r / r0")
    ax_velocity.set_xlabel("r / r0")
    ax_mass.set_xlabel("r / r0")
    ax_pressure.set_ylabel("P / P0")
    ax_density.set_ylabel("rho / rho0")
    ax_velocity.set_ylabel("u / u_esc")
    ax_mass.set_ylabel("m / (r_0^3 rho_0)")
    ax_pressure.set_title("Pressure (IC)")
    ax_density.set_title("Density (IC)")
    ax_velocity.set_title("Velocity (IC)")
    ax_mass.set_title("Mass (IC)")
    ax_pressure.grid()
    ax_density.grid()
    ax_velocity.grid()
    ax_mass.grid()

    plt.show()