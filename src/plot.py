import pandas as pd
import matplotlib.pyplot as plt


def __get_coord_of_max_value(radius, vals):
    max_val = max(vals)
    index = vals.index(max_val)
    r = radius[index]
    return r, max_val


def annotate(ax, time, radius, vals):
    r, m = __get_coord_of_max_value(radius=radius, vals=vals)
    ax.text(
        r,
        m + (m * 0.05),
        "t = {} s".format(round(float(time), 4)),
        rotation=90
    )


def plot_time(output_path, iteration, fig, ax_density, ax_pressure, ax_velocity, ax_temperature, r_0, rho_0, P_0, vesc,
              T_0, fig_path, min_x=0, max_x=8, color='black'):
    f = output_path + "/{}.csv".format(iteration)
    time = None
    mass_loss_fraction = None
    with open(f, 'r') as infile:
        time = next(infile)
        timestep = next(infile)
        mass_loss_fraction = next(infile)

    df = pd.read_csv(f, skiprows=3, header=None)
    radius, mass, pressure, density, velocity, temperature = df[1], df[2], df[3], df[4], df[5], df[6]

    # fig = plt.figure(figsize=(16, 9))
    # ax_density = fig.add_subplot(221)
    # ax_pressure = fig.add_subplot(222)
    # ax_velocity = fig.add_subplot(223)
    # ax_temperature = fig.add_subplot(224)

    norm_radius = [i / r_0 for i in radius]
    norm_density = [i / rho_0 for i in density]
    norm_pressure = [i / P_0 for i in pressure]
    norm_velocity = [i / vesc for i in velocity]
    norm_temperature = [i / T_0 for i in temperature]

    ax_density.plot(
        norm_radius,
        norm_density,
        linewidth=2.0,
        color=color,
        label="t = {} s".format(round(float(time), 2))
    )
    ax_pressure.plot(
        norm_radius,
        norm_pressure,
        linewidth=2.0,
        color=color
    )
    ax_velocity.plot(
        norm_radius,
        norm_velocity,
        linewidth=2.0,
        color=color
    )
    ax_temperature.plot(
        norm_radius,
        norm_temperature,
        linewidth=2.0,
        color=color
    )
    ax_density.set_ylabel("$ρ / ρ_{0}$")
    ax_pressure.set_ylabel("$P / P_{0}$")
    ax_velocity.set_ylabel("$v / v_{esc}$")
    ax_temperature.set_ylabel("$T / T_{0}$")

    for ax in [ax_density, ax_pressure, ax_velocity, ax_temperature]:
        ax.set_xlim(min_x, max_x)
        ax.grid()
    ax_density.legend()

    # annotate(ax=ax_density, time=time, radius=norm_radius, vals=norm_density)
    # annotate(ax=ax_pressure, time=time, radius=norm_radius, vals=norm_pressure)
    # annotate(ax=ax_velocity, time=time, radius=norm_radius, vals=norm_velocity)
    # annotate(ax=ax_temperature, time=time, radius=norm_radius, vals=norm_temperature)

    fig.supxlabel("$r / r_{esc}$")
    # fig.suptitle("Time {} s".format(time))

    # plt.savefig(fig_path + "/{}.png".format(time), format='png')
