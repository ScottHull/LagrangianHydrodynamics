from src import solver
from src import initial_conditions as ic
import pandas as pd

import os
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
outfile_dir = "spherical_test_outputs"
output_plots_dir = "spherical_plots"
max_time = 86400  # 1 day in seconds

# s = solver.LagrangianSolver1DSpherical(
#     num_shells=num_shells,
#     P_0=P_0,
#     T_0=T_0,
#     gamma=gamma,
#     gamma_a=gamma_a,
#     m_a=M_a,
#     r_0=r_0,
#     mass_planet=mass_planet,
#     u_s=u_s,
#     outfile_dir=outfile_dir,
#     plot_separation=10000,
#     use_cfl=True,
#     show_figs=False,
#     save_figs=True,
#     output_file_interval=10000,
#     fig_save_path=output_plots_dir,
# )
# s.solve(max_time=max_time)

# read every file in the output directory and get the first row (time) as a float
# then sort the list of times
times = []
iterations = []
for f in os.listdir(outfile_dir):
    with open(os.path.join(outfile_dir, f), "r") as file:
        times.append(float(file.readline().split()[0]))
        iterations.append(int(f.split(".")[0]))
t_i = zip(times, iterations)
# sort the list of tuples by the first element (time)
t_i = list(sorted(t_i, key=lambda x: x[0]))

r_0 = 6.4e6  # initial planet radius
T_0 = 288
P_0 = 1.01e5
rho_0 = 1.223254617378986
c_s_0 = 339.99010127578305
vesc = 11160.421945428408

# create a figure with 4 subplots
fig = plt.figure(figsize=(16, 9))
ax_density = fig.add_subplot(221)
ax_pressure = fig.add_subplot(222)
ax_velocity = fig.add_subplot(223)
ax_temperature = fig.add_subplot(224)

# for each unique time, get the rows that correspond to that time
# only plot every other time to reduce the number of lines
for i in range(0, len(t_i), 2):
    time, iteration = t_i[i]
    df = pd.read_csv(outfile_dir + "/{}.csv".format(iteration), skiprows=3, header=None, index_col=None)
    # get the radius, pressure, velocity, density, and temperature
    radius = df[0].values
    pressure = df[1].values
    velocity = df[2].values
    density = df[3].values
    temperature = df[4].values

    # plot the normalized data
    ax_density.plot(radius / r_0, density / rho_0, label=f"{time:.2f} s")
    ax_pressure.plot(radius / r_0, pressure / P_0, label=f"{time:.2f} s")
    ax_velocity.plot(radius / r_0, velocity / vesc, label=f"{time:.2f} s")
    ax_temperature.plot(radius / r_0, temperature / T_0, label=f"{time:.2f} s")

# label the axes
ax_density.set_xlabel("r / r0")
ax_density.set_ylabel("rho / rho0")
ax_pressure.set_xlabel("r / r0")
ax_pressure.set_ylabel("P / P0")
ax_velocity.set_xlabel("r / r0")
ax_velocity.set_ylabel("v / vesc")
ax_temperature.set_xlabel("r / r0")
ax_temperature.set_ylabel("T / T0")

# add a grid to each subplot
for ax in fig.axes:
    ax.grid()

# add a legend to the first subplot
ax_density.legend()

plt.savefig("spherical_test.png", dpi=200)

