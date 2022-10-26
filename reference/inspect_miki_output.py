import os
import csv
import pandas as pd
import matplotlib.pyplot as plt

# use colorblind friendly colors
plt.style.use('seaborn-colorblind')

file_path = r"C:\Users\Scott\OneDrive\Desktop\Mt.dat"
r0 = 6.4 * 10 ** 6
r_0 = 6.4e6  # initial planet radius
T_0 = 288
P_0 = 1.01e5
rho_0 = 1.223254617378986
c_s_0 = 339.99010127578305
vesc = 11160.421945428408

if "tmp.txt" in os.listdir(os.getcwd()):
    os.remove("tmp.txt")
fixed_output = open("tmp.txt", 'w')
with open(file_path, "r") as file:
    reader = csv.reader(file, delimiter=" ")
    for line in reader:
        line = [x for x in line if x != '']
        for index, i in enumerate(line):
            if "E-00" in str(i):
                line[index] = float(str(i).replace("E-00", "e-"))
            elif "E+00" in str(i):
                line[index] = float(str(i).replace("E+00", "e+"))
        fixed_output.write(",".join([str(i) for i in line]) + "\n")

file.close()
fixed_output.close()
# os.remove(file_path)
# os.rename("tmp.txt", file_path)

df = pd.read_csv('tmp.txt', header=None)
os.remove("tmp.txt")

unique_times = df[6].unique()  # get unique times

# create a figure with 4 subplots
fig = plt.figure(figsize=(16, 9))
ax_density = fig.add_subplot(221)
ax_pressure = fig.add_subplot(222)
ax_velocity = fig.add_subplot(223)
ax_temperature = fig.add_subplot(224)

# for each unique time, get the rows that correspond to that time
# only plot every other time to reduce the number of lines
# for i in range(0, len(unique_times), 2):
for i in range(0, len(unique_times), 1):
    time = unique_times[i]
    df_time = df[df[6] == time]  # time in seconds
    # get the radius, pressure, velocity, density, and temperature
    radius = df_time[0].values
    pressure = df_time[1].values
    velocity = df_time[2].values
    density = df_time[3].values
    temperature = df_time[4].values

    # plot the normalized data
    ax_density.plot(radius / r0, density / rho_0, label=f"{time:.2f} s")
    ax_pressure.plot(radius / r0, pressure / P_0, label=f"{time:.2f} s")
    ax_velocity.plot(radius / r0, velocity / vesc, label=f"{time:.2f} s")
    ax_temperature.plot(radius / r0, temperature / T_0, label=f"{time:.2f} s")

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

plt.show()
