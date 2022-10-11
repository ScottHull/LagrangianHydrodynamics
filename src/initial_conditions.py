from math import pi, sqrt, tan

"""
Initial atmosphere conditions described in Section 2.3 of Genda and Abe 2003.
Assumes a hydrostatically equilibrated polytropic atmosphere with polytropic exponent gamma_a.
"""


class InitialConditionsSpherical:

    def lambda_0_initial(self, mass_planet, rho_0, r_0, P_0, G=6.674 * 10 ** -11):
        return (G * mass_planet * rho_0) / (r_0 * P_0)


    def _r_last(self, polytropic_exponent, lambda_0, r_0):
        numerator = (r_0 * polytropic_exponent * lambda_0) - (r_0 * lambda_0)
        denominator = (polytropic_exponent * lambda_0) - polytropic_exponent - lambda_0
        return numerator / denominator


    def radius_initial(self, index, total_shells, polytropic_exponent, lambda_0, r_0):
        r_last = self._r_last(polytropic_exponent=polytropic_exponent, lambda_0=lambda_0, r_0=r_0)
        return r_0 + (index * ((r_last - r_0) / total_shells))


    def mass_initial(self, mass_last_index, rho_index, r_last_index, r_index):
        # print(mass_last_index, rho_index, r_last_index, r_index)
        return mass_last_index + (rho_index * (4 / 3) * pi * ((r_index ** 3) - (r_last_index ** 3)))


    def velocity_initial(self, polytropic_exponent, m_a, T_index, R_gas=8.314462):
        # return sqrt(((polytropic_exponent * R_gas) / m_a) * T_index)
        return 0.0


    def pressure_initial(self, polytropic_exponent, lambda_0, radius, radius_0, p_0):
        a1 = ((polytropic_exponent - 1) / polytropic_exponent) * lambda_0
        a2 = (radius_0 / radius) - 1
        exponent = polytropic_exponent / (polytropic_exponent - 1)
        p_p0 = (a1 * a2 + 1) ** exponent  # P / P_0
        return p_p0 * p_0


    def density_initial(self, polytropic_exponent, lambda_0, radius, radius_0, rho_0):
        a1 = ((polytropic_exponent - 1) / polytropic_exponent) * lambda_0
        a2 = (radius_0 / radius) - 1
        exponent = 1.0 / (polytropic_exponent - 1)
        rho_rho0 = (a1 * a2 + 1) ** exponent  # rho / rho_0
        return rho_rho0 * rho_0


    def temperature_initial(self, polytropic_exponent, lambda_0, radius, radius_0, T_0):
        a1 = ((polytropic_exponent - 1) / polytropic_exponent) * lambda_0
        a2 = (radius_0 / radius) - 1
        exponent = 1
        T_T0 = (a1 * a2 + 1) ** exponent  # T / T_0
        return T_T0 * T_0


class InitialConditionsJet(InitialConditionsSpherical):

    def mass_initial(self, mass_last_index, rho_index, r_last_index, r_index):
        """
        Volume of a cone: V = (1/3) * pi * r^2 * z
        where z is the height of the cone and r is the radius of the base.
        Therefore, the
        """
        # print(mass_last_index, rho_index, r_last_index, r_index)
        return mass_last_index + (rho_index * (4 / 3) * pi * ((r_index ** 3) - (r_last_index ** 3)))
