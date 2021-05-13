from math import sqrt


def density_nd(density, density_0):
    return density / density_0


def pressure_nd(pressure, pressure_0):
    return pressure / pressure_0


def radius_nd(radius, radius_0):
    return radius / radius_0


def mass_nd(mass, radius_0, density_0):
    return mass / ((radius_0 ** 3) * density_0)


def soundspeed_0(gamma, pressure_0, density_0):
    return sqrt(gamma * (pressure_0 / density_0))


def velocity_nd(velocity, gamma, pressure_0, density_0):
    c_s_0 = soundspeed_0(gamma=gamma, pressure_0=pressure_0, density_0=density_0)
    return velocity / c_s_0


def time_nd(time, gamma, pressure_0, density_0, radius_0):
    c_s_0 = soundspeed_0(gamma=gamma, pressure_0=pressure_0, density_0=density_0)
    return time * (c_s_0 / radius_0)
