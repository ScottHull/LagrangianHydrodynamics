import pandas as pd
import matplotlib.pyplot as plt


def plot_time(output_path, iteration, fig, ax_density, ax_pressure, ax_velocity, ax_temperature, r_0, rho_0, P_0, c_s_0,
              T_0, fig_path):
    f = output_path + "/{}.csv".format(iteration)
    time = None
    mass_loss_fraction = None
    with open(f, 'r') as infile:
        time = next(infile)[0]
        mass_loss_fraction = next(next(infile))[0]

    df = pd.read_csv(f, skiprows=3, header=None)
    radius, mass, pressure, density, velocity, temperature = df[1], df[2], df[3], df[4], df[5], df[6]

    # fig = plt.figure(figsize=(16, 9))
    # ax_density = fig.add_subplot(221)
    # ax_pressure = fig.add_subplot(222)
    # ax_velocity = fig.add_subplot(223)
    # ax_temperature = fig.add_subplot(224)

    norm_radius = [i / r_0 for i in radius]

    ax_density.plot(
        norm_radius,
        [i / rho_0 for i in density],
        linewidth=2.0,
        color='black'
    )
    ax_pressure.plot(
        norm_radius,
        [i / P_0 for i in pressure],
        linewidth=2.0,
        color='black'
    )
    ax_velocity.plot(
        norm_radius,
        [i / c_s_0 for i in velocity],
        linewidth=2.0,
        color='black'
    )
    ax_temperature.plot(
        norm_radius,
        [i / T_0 for i in temperature],
        linewidth=2.0,
        color='black'
    )
    ax_density.set_ylabel("$ρ / ρ_{0}$")
    ax_pressure.set_ylabel("$P / P_{0}$")
    ax_velocity.set_ylabel("$v / v_{0}$")
    ax_temperature.set_ylabel("$T / T_{0}$")

    ax_density.grid()
    ax_pressure.grid()
    ax_velocity.grid()
    ax_temperature.grid()

    fig.supxlabel("$r / r_{0}$")
    fig.suptitle("Time {} s".format(time))

    plt.savefig(fig_path + "/{}.png".format(time), format='png')
