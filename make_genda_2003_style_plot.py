import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src import plot

plt.style.use('seaborn-colorblind')

# ====================================== INPUTS ======================================

times = [
    [0.5, 1, 1.5, 2, 3],
    [3.5, 4, 4.5, 5],
    [100, 200, 400, 600, 800, 1000, 1200, 1400, 1600],
]

output_directory = "/scratch/shull4/spherical-500b073S_outputs"  # the directory containing the output files
ic_file = ic_filename = output_directory.split("/")[-1] + ".txt"  # initial conditions file

# ====================================================================================

# ======================= Collect the initial conditions =============================

ic_df = pd.read_csv(ic_filename, header=None, index_col=0)  # read the initial conditions file

r_0 = float(ic_df[1]['r_0'])  # initial planet radius
T_0 = float(ic_df[1]['T_0'])
P_0 = float(ic_df[1]['P_0'])
rho_0 = float(ic_df[1]['rho_0'])
c_s_0 = float(ic_df[1]['c_s_0'])
vesc = float(ic_df[1]['v_esc'])

# ======================= Read files and get closest iterations and times =======================

# read all files in the directory, get the first row (time) as a float, then get the closest time to each time in times
# then sort the list of times
iteration_and_time = []
for f in os.listdir(output_directory):
    with open(os.path.join(output_directory, f), "r") as file:
        iteration_and_time.append((int(f.split(".")[0]), float(file.readlines()[0].split()[0])))
iteration_and_time = list(sorted(iteration_and_time, key=lambda x: x[1]))

# get the closest time to each time in times
closest_times = []
for timeset in times:
    closest_times.append([])
    for t in timeset:
        closest_times[-1].append(min(iteration_and_time, key=lambda x: abs(x[1] - t)))


# set up the figure, where there are 4 rows and len(times) columns
fig, axes = plt.subplots(4, len(times), figsize=(16, 9))

# label the top of each column with a letter in alphabetical order
for i in range(len(times)):
    axes[0, i].set_title(chr(97 + i))

# add y labels to the leftmost column
axes[0, 0].set_ylabel(r"v / v$_{\rm esc}$")
axes[1, 0].set_ylabel(r"$ρ / ρ_{0}$")
axes[2, 0].set_ylabel(r"$P / P_{0}$")
axes[3, 0].set_ylabel(r"$T / T_{0}$")

# label the bottom of each column with r/r_0
for i in range(len(times)):
    axes[-1, i].set_xlabel(r"$r / r_{0}$")

# plot each time in each timeset in the appropriate column
# we plot velocity in row 1, density in row 2, pressure in row 3, and temperature in row 4
for i, timeset in enumerate(closest_times):
    for j, (iteration, time) in enumerate(timeset):
        df = pd.read_csv(os.path.join(output_directory, f"{iteration}.csv"), skiprows=3, header=None, index_col=0)
        radius, mass, pressure, density, velocity, temperature = df[1], df[2], df[3], df[4], df[5], df[6]
        axes[0, i].plot(
            radius / r_0, velocity / vesc, color='black', linewidth=2.0
        )
        axes[1, i].plot(
            radius / r_0, density / rho_0, color='black', linewidth=2.0
        )
        axes[2, i].plot(
            radius / r_0, pressure / P_0, color='black', linewidth=2.0
        )
        axes[3, i].plot(
            radius / r_0, temperature / T_0, color='black', linewidth=2.0
        )

        # annotate each line with the time at the maximum y value in column 1, else use the maximum y value in the column
        for index, ax in enumerate(axes[:, i]):
            if i == 0:
                # get the maximum y value of the line
                y = max(ax.lines[-1].get_ydata())
                # get the x value at the maximum y value
                x = ax.lines[-1].get_xdata()[np.argmax(ax.lines[-1].get_ydata())]
            else:
                # get the minimum y value of the line
                y = min(ax.lines[-1].get_ydata())
                # get the x value at the minimum y value
                x = ax.lines[-1].get_xdata()[np.argmin(ax.lines[-1].get_ydata())]
            # annotate the line with the time at xy
            ax.annotate(
                f"{time:.1f} s",
                xy=(x, y),
                xytext=(x, y),
                textcoords="offset points",
                xycoords="data",
                ha="center",
                va="bottom",
                fontsize=14,
                color="black",
            )

# save the figure
plt.savefig(f"{output_directory.split('/')[-1]}_genda_2003_style_plot.png", dpi=300)
