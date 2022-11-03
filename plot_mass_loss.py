import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# use colorblind-friendly colors
plt.style.use('seaborn-colorblind')
output_directory = "spherical_test_outputs"

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

# plot the mass loss fraction vs time
plt.plot([x[0] for x in times_and_mass_loss_fraction], [x[1] for x in times_and_mass_loss_fraction], linewidth=2.0)
plt.xlabel("Time (s)")
plt.ylabel("Mass Loss Fraction")
plt.title("Mass Loss Fraction vs Time")
plt.grid()
plt.savefig("mass_loss_fraction_vs_time.png", dpi=200)

# compare both given rho0 and given matm

rho0_output_directory = "/scratch/shull4/spherical_test_outputs_given_rho0"
matm_output_directory = "/scratch/shull4/spherical_test_outputs_given_matm"

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
for i in [rho0_output_directory, matm_output_directory]:

    # read every file in the output directory and get the first row (time) and the third line (mass loss fraction) as floats
    # then sort the list of times and mass loss fractions
    times_and_mass_loss_fraction = []
    for f in os.listdir(i):
        with open(os.path.join(i, f), "r") as file:
            lines = file.readlines()
            times_and_mass_loss_fraction.append((float(lines[0]), float(lines[2])))
        file.close()

# sort the list of tuples by the first element (time)
times_and_mass_loss_fraction = list(sorted(times_and_mass_loss_fraction, key=lambda x: x[0]))

# plot the mass loss fraction vs time
ax.plot([x[0] for x in times_and_mass_loss_fraction], [x[1] for x in times_and_mass_loss_fraction], linewidth=2.0, label=i.split("/")[-1])

# fit the curve to times and mass loss
def objective(x, a, b, c, d, e, f):
    return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + f
x, y = np.array([x[0] for x in times_and_mass_loss_fraction]), np.array([x[1] for x in times_and_mass_loss_fraction])
popt, _ = curve_fit(objective, x, y)
a, b, c, d, e, f = popt
# plot the fit
y_line = objective(x, a, b, c, d, e, f)
ax.plot(x, y_line, '--', label="fit")


ax.set_xlabel("Time (s)")
ax.set_ylabel("Mass Loss Fraction")
ax.set_title("Mass Loss Fraction vs Time")
ax.grid()
ax.legend()
plt.savefig("mass_loss_fraction_vs_time.png", dpi=200)
