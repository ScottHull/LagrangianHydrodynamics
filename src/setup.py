import src.initial_conditions as ic
import src.nondimensional as nd

from math import sqrt


class Point:
    def __init__(self, id, q, mass, pressure, temperature, density, radius, velocity):
        self.id = id
        self.q = q
        self.mass = mass
        self.pressure = pressure
        self.temperature = temperature
        self.density = density
        self.radius = radius
        self.velocity = velocity


class System:

    def __init__(self, num_shells, gamma_a, mass_planet, r_0, rho_0, P_0, T_0, m_a, gamma, u_s):
        self.num_shells = num_shells
        self.G = 6.67408e-11

        self.gamma_a = gamma_a
        self.gamma = gamma
        self.r_0 = r_0
        self.rho_0 = rho_0
        self.P_0 = P_0
        self.T_0 = T_0
        self.m_a = m_a
        self.mass_planet = mass_planet
        self.vesc = sqrt(2 * self.G * self.mass_planet / self.r_0)
        self.lambda_0 = self.G * self.mass_planet * self.rho_0 / (self.r_0 * self.P_0)
        self.c_s_0 = sqrt(gamma_a * self.P_0 / self.rho_0)
        self.u_s = u_s  # initial shock velocity, is 0.5 * u_esc in Genda and Abe 2003
        if self.u_s is None:
            self.u_s = 0.5 * self.vesc  # is 0.5 * u_esc in Genda and Abe 2003

        self.grid = []
        self.__setup_grid()
        self.__nondimensionalize_initial()

    def __setup_grid(self):
        current = 0

        while current < self.num_shells:
            r = ic.radius_initial(index=current, lambda_0=self.lambda_0, polytropic_exponent=self.gamma_a, r_0=self.r_0,
                                  total_shells=self.num_shells)
            P = ic.pressure_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, p_0=self.P_0, radius=r,
                                    radius_0=self.r_0)
            T = ic.temperature_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, T_0=self.T_0, radius=r,
                                       radius_0=self.r_0)
            rho = ic.density_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, rho_0=self.rho_0,
                                     radius=r, radius_0=self.r_0)
            if current > 0:
                m = ic.mass_initial(mass_last_index=self.grid[current - 1].mass, r_index=r,
                                    r_last_index=self.grid[current - 1].radius, rho_index=rho)
            else:
                m = 0
            v = ic.velocity_initial(polytropic_exponent=self.gamma_a, T_index=T, m_a=self.m_a)
            q = 0
            point = Point(
                id=current,
                q=q,
                mass=m,
                pressure=P,
                temperature=T,
                density=rho,
                radius=r,
                velocity=v
            )
            self.grid.append(point)
            current += 1
        self.grid[-1].pressure = 0.0
        self.grid[-1].density = 0.0
        # self.grid[0].velocity = 0.5 * self.vesc
        self.grid[0].velocity = self.u_s

    def __nondimensionalize_initial(self):
        for p in self.grid:
            p.mass = nd.mass_nd(mass=p.mass, density_0=self.rho_0, radius_0=self.r_0)
            p.pressure = nd.pressure_nd(pressure=p.pressure, pressure_0=self.P_0)
            p.temperature = nd.temperature_nd(temperature=p.temperature, temperature_0=self.T_0)
            p.density = nd.density_nd(density=p.density, density_0=self.rho_0)
            p.radius = nd.radius_nd(radius=p.radius, radius_0=self.r_0)
            p.velocity = nd.velocity_nd(velocity=p.velocity, density_0=self.rho_0, gamma=self.gamma,
                                        pressure_0=self.P_0, c_s_0=self.c_s_0)
