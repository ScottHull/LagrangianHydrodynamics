import os
import matplotlib.pyplot as plt

from src import plot

plt.style.use('seaborn-colorblind')

times = [0.5, 1, 1.5, 2, 3]
output_directory = "spherical_test_outputs"

# read all files in the directory, get the first row (time) as a float, then get the closest time to each time in times
# then sort the list of times
iteration_and_time = []
for f in os.listdir(output_directory):
    with open(os.path.join(output_directory, f), "r") as file:
        iteration_and_time.append((int(f.split(".")[0]), float(file.readlines()[0].split()[0])))
iteration_and_time = list(sorted(iteration_and_time, key=lambda x: x[1]))

# get the closest time to each time in times
closest_times = []
for t in times:
    closest_times.append(min(iteration_and_time, key=lambda x: abs(x[1] - t)))

closest_iterations = [x[0] for x in closest_times]

r_0 = 6.4e6  # initial planet radius
T_0 = 288
P_0 = 1.01e5
rho_0 = 1.223254617378986
c_s_0 = 339.99010127578305
vesc = 11160.421945428408

fig = plt.figure(figsize=(16, 9))
ax_density = fig.add_subplot(221)
ax_pressure = fig.add_subplot(222)
ax_velocity = fig.add_subplot(223)
ax_temperature = fig.add_subplot(224)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for index, i in enumerate(closest_iterations):
    plot.plot_time(
        output_path=output_directory,
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
        fig_path="",
        color=colors[index],
        min_x=1,
        max_x=1.005
    )

plt.savefig("lagrangian.png", format='png')
