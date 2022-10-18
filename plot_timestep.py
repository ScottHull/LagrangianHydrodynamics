import os
import shutil
import matplotlib.pyplot as plt

from src import plot

plt.style.use('seaborn-colorblind')

iterations_to_plot = [10, 50, 100]
output_path = "/scratch/shull4/outputs"
fig_path = "/scratch/shull4/lagrangian_figs"

r_0 = 6.4e6  # initial planet radius
T_0 = 288
P_0 = 1.01e5
rho_0 = 1.223254617378986
c_s_0 = 339.99010127578305
vesc = 11160.421945428408

if os.path.exists(fig_path):
    shutil.rmtree(fig_path)
os.mkdir(fig_path)

# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


fig = plt.figure(figsize=(16, 9))
ax_density = fig.add_subplot(221)
ax_pressure = fig.add_subplot(222)
ax_velocity = fig.add_subplot(223)
ax_temperature = fig.add_subplot(224)

for i in iterations_to_plot:
    plot.plot_time(
        output_path=output_path,
        iteration=i,
        fig=fig,
        ax_density=ax_density,
        ax_pressure=ax_pressure,
        ax_velocity=ax_velocity,
        ax_temperature=ax_temperature,
        r_0=r_0,
        rho_0=rho_0,
        P_0=r_0,
        vesc=vesc,
        T_0=T_0,
        fig_path=fig_path,
        min_x=0,
        max_x=1.005
    )

plt.savefig("lagrangian.png", format='png')

