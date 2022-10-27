import os
import matplotlib.pyplot as plt

# use colorblind-friendly colors
plt.style.use('seaborn-colorblind')
output_directory = "spherical_test_outputs"

# read every file in the output directory and get the first row (time) and the third line (mass loss fraction) as floats
# then sort the list of times and mass loss fractions
times_and_mass_loss_fraction = []
for f in os.listdir(output_directory):
    with open(os.path.join(output_directory, f), "r") as file:
        times_and_mass_loss_fraction.append((float(file.readline().split()[0]), float(file.readlines()[2].split()[0])))

# sort the list of tuples by the first element (time)
times_and_mass_loss_fraction = list(sorted(times_and_mass_loss_fraction, key=lambda x: x[0]))

# plot the mass loss fraction vs time
plt.plot([x[0] for x in times_and_mass_loss_fraction], [x[1] for x in times_and_mass_loss_fraction], linewidth=2.0)
plt.xlabel("Time (s)")
plt.ylabel("Mass Loss Fraction")
plt.title("Mass Loss Fraction vs Time")
plt.savefig("mass_loss_fraction_vs_time.png", dpi=200)
