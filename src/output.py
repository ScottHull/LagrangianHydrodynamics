import csv


def write_state(path, fname, time, timestep, grid, mass_fraction_loss, system):
    with open(path + "/{}.csv".format(fname), 'w') as outfile:
        outfile.write("{}\n{}\n{}\n".format(time, timestep, mass_fraction_loss))
        for p in grid:
            line = "{},{},{},{},{},{},{}\n".format(
                p.id,
                round(p.radius * system.r_0, 4),
                round(p.mass * (system.r_0 ** 3 * system.rho_0), 4),
                round(p.pressure * system.P_0, 4),
                round(p.density * system.rho_0, 4),
                round(p.velocity * system.c_s_0, 4),
                round(p.temperature * system.T_0, 4)
            )
            outfile.write(line)
    outfile.close()
