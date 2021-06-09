import csv


def write_state(path, fname, time, timestep, grid, mass_fraction_loss, system):
    with open(path + "/{}.csv".format(fname), 'w') as outfile:
        outfile.write("{}\n{}\n{}\n".format(time, timestep, mass_fraction_loss))
        for p in grid:
            line = "{},{},{},{},{},{},{}\n".format(
                p.id,
                p.radius * system.r_0,
                p.mass * (system.r_0 ** 3 * system.rho_0),
                p.pressure * system.P_0,
                p.density * system.rho_0,
                p.velocity * system.c_s_0,
                p.temperature * system.T_0
            )
            outfile.write(line)
    outfile.close()
