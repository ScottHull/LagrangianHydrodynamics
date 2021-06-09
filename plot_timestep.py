import os
import shutil
import matplotlib.pyplot as plt

from src import plot

iteration = 10
output_path = "/scratch/shull4/outputs"
fig_path = "/scratch/shull4/lagrangian_figs"

r_0 = 6.4e6  # initial planet radius
T_0 = 288
P_0 = 1.01e5
rho_0 = 1.223254617378986
c_s_0 = 339.99010127578305

if os.path.exists(fig_path):
    shutil.rmtree(fig_path)
os.mkdir(fig_path)

fig = plt.figure(figsize=(16, 9))
ax_density = fig.add_subplot(221)
ax_pressure = fig.add_subplot(222)
ax_velocity = fig.add_subplot(223)
ax_temperature = fig.add_subplot(224)

plot.plot_time(
    output_path=output_path,
    iteration=iteration,
    fig=fig,
    ax_density=ax_density,
    ax_pressure=ax_pressure,
    ax_velocity=ax_velocity,
    ax_temperature=ax_temperature,
    r_0=r_0,
    rho_0=rho_0,
    P_0=r_0,
    c_s_0=c_s_0,
    T_0=T_0,
    fig_path=fig_path
)
