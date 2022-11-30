import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# use colorblind-friendly colors
base_path = '/scratch/shull4/'
outputs = [
    ("jet-500b073S_outputs", "Jet 45 deg."), ("jet-500b073S_outputs_65deg", "Jet 65 deg."),
    ("spherical-500b073S_outputs", "Spherical")
]

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
for p, name in outputs:
    path = os.path.join(base_path, p)
    # read every file in the output directory and get the first row (time) and the third line (mass loss fraction) as floats
    # then sort the list of times and mass loss fractions
    times_and_mass_loss_fraction = []
    for f in os.listdir(path):
        with open(os.path.join(path, f), "r") as file:
            lines = file.readlines()
            times_and_mass_loss_fraction.append((float(lines[0]), float(lines[2])))
        file.close()

    # sort the list of tuples by the first element (time)
    times_and_mass_loss_fraction = list(sorted(times_and_mass_loss_fraction, key=lambda x: x[0]))
    ax.plot([x[0] for x in times_and_mass_loss_fraction], [x[1] for x in times_and_mass_loss_fraction], linewidth=2.0, label=name)
    # annotate the last y value in the upper right corner
    ax.annotate(f"Final Mass Loss Frac.: {times_and_mass_loss_fraction[-1][1]:.2f}", (times_and_mass_loss_fraction[-1][0], times_and_mass_loss_fraction[-1][1]))

plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Mass Loss Fraction")
plt.title("Mass Loss Fraction vs Time")
plt.grid()
plt.savefig("mass_loss_fraction_vs_time_comparison.png", dpi=200)
