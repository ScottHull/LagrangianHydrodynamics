

"""
Initial atmosphere conditiosn described in Section 2.3 of Genda and Abe 2003.
Assumes a hydrostatically equilibrated polytropic atmosphere with polytropic exponent gamma_a.
"""

def pressure_initial(polytropic_exponent, lambda_0, radius, radius_0):
    a1 = ((polytropic_exponent - 1) / polytropic_exponent) * lambda_0
    a2 = (radius_0 / radius) - 1
    exponent = polytropic_exponent / (polytropic_exponent - 1)
    p_p0 = (a1 * a2 + 1) ** exponent  # P / P_0
    return p_p0

def density_initial(polytropic_exponent, lambda_0, radius, radius_0):
    a1 = ((polytropic_exponent - 1) / polytropic_exponent) * lambda_0
    a2 = (radius_0 / radius) - 1
    exponent = 1 / (polytropic_exponent - 1)
    rho_rho0 = (a1 * a2 + 1) ** exponent  # rho / rho_0
    return rho_rho0

def temperature_initial(polytropic_exponent, lambda_0, radius, radius_0)
    a1 = ((polytropic_exponent - 1) / polytropic_exponent) * lambda_0
    a2 = (radius_0 / radius) - 1
    exponent = 1
    t_t0 = (a1 * a2 + 1) ** exponent  # T / T_0
    return t_t0
