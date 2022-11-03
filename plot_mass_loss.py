import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# use colorblind-friendly colors
plt.style.use('seaborn-colorblind')
# output_directory = "/scratch/shull4/spherical_test_outputs"
output_directory = "/scratch/shull4/jet_test_outputs"

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
# read every file in the output directory and get the first row (time) and the third line (mass loss fraction) as floats
# then sort the list of times and mass loss fractions
times_and_mass_loss_fraction = []
for f in os.listdir(output_directory):
    with open(os.path.join(output_directory, f), "r") as file:
        lines = file.readlines()
        times_and_mass_loss_fraction.append((float(lines[0]), float(lines[2])))
    file.close()

# sort the list of tuples by the first element (time)
times_and_mass_loss_fraction = list(sorted(times_and_mass_loss_fraction, key=lambda x: x[0]))
ax.plot([x[0] for x in times_and_mass_loss_fraction], [x[1] for x in times_and_mass_loss_fraction], linewidth=2.0, label="Calculated")

# fit the curve to times and mass loss
def objective(x, a, b, c):
    return a * np.exp(-b * x) + c
x, y = np.array([i[0] for i in times_and_mass_loss_fraction if i[0] > 5]), \
       np.array([i[1] for i in times_and_mass_loss_fraction if i[0] > 5])
popt, _ = curve_fit(objective, x, y)
a, b, c = popt
# plot the fit
# y_line = objective(x, a, b, c, d, e, f)
y_line = objective(np.arange(0, 100), a, b, c)
ax.plot(np.arange(0, 100), y_line, '--', label="fit")

# plot the mass loss fraction vs time
plt.xlabel("Time (s)")
plt.ylabel("Mass Loss Fraction")
plt.title("Mass Loss Fraction vs Time")
plt.grid()
plt.savefig("mass_loss_fraction_vs_time.png", dpi=200)
