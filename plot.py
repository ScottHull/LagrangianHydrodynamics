import pandas as pd
import matplotlib.pyplot as plt

mass_loss_df = pd.read_csv("lagrangian_solver.csv")

# plot mass loss fraction
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    mass_loss_df['time'],
    mass_loss_df['mass_loss'],
    linewidth=2.0,
    color='black'
)
ax.set_xlabel("Time (s)")
ax.set_ylabel("Mass Loss Fraction")
ax.set_title("Volatile Mass Loss Fraction")
ax.grid()

plt.savefig("mass_loss.png", format='png')
