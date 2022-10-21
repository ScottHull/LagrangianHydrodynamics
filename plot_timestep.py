import os
import shutil
import matplotlib.pyplot as plt

from src import plot

plt.style.use('seaborn-colorblind')

times_to_plot = [0.5, 1, 1.5, 2.0, 2.5, 3.0]
times_to_plot = dict(zip(times_to_plot, [False for i in times_to_plot]))
# output_path = "/scratch/shull4/outputs"
# fig_path = "/scratch/shull4/lagrangian_figs"
output_path = "spherical_test_outputs"
fig_path = ""

r_0 = 6.4e6  # initial planet radius
T_0 = 288
P_0 = 1.01e5
rho_0 = 1.223254617378986
c_s_0 = 339.99010127578305
vesc = 11160.421945428408

try:
    if os.path.exists(fig_path):
        shutil.rmtree(fig_path)
    os.mkdir(fig_path)
except OSError:
    pass

fig = plt.figure(figsize=(16, 9))
ax_density = fig.add_subplot(221)
ax_pressure = fig.add_subplot(222)
ax_velocity = fig.add_subplot(223)
ax_temperature = fig.add_subplot(224)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for f in os.listdir(output_path):
    with open(f"{output_path}/{f}", 'r') as infile:
        time = float(next(infile))
        timestep = float(next(infile))
        mass_loss_fraction = float(next(infile))
    for index, i in enumerate(times_to_plot.keys()):
        if time >= i and not times_to_plot[i]:
            times_to_plot[i] = True
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
                color=colors[index],
                min_x=1,
                max_x=1.005
            )
    if False in times_to_plot.values():
        continue
    else:
        break

plt.savefig("lagrangian.png", format='png')

