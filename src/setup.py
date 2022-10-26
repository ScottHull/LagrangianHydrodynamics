import src.initial_conditions as ic
import src.nondimensional as nd

from math import sqrt
from copy import copy


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


class SphericalSystem:

    def __init__(self, num_shells, gamma_a, mass_planet, r_0, rho_0, P_0, T_0, m_a, gamma, u_s, **kwargs):
        self.num_shells = num_shells
        self.G = 6.67408e-11
        self.gamma_a = gamma_a
        self.gamma = gamma
        self.r_0 = r_0
        self.P_0 = P_0
        self.T_0 = T_0
        self.m_a = m_a
        self.mass_planet = mass_planet
        self.mass_atmosphere = kwargs.get("mass_atmosphere", None)
        self.vesc = sqrt(2 * self.G * self.mass_planet / self.r_0)
        self.rho_0 = rho_0
        self.lambda_0 = self.G * self.mass_planet * self.rho_0 / (self.r_0 * self.P_0)
        self.c_s_0 = kwargs.get("c_s_0", sqrt(self.gamma_a * self.P_0 / self.rho_0))
        self.u_s = u_s  # initial shock velocity, is 0.5 * u_esc in Genda and Abe 2003
        if self.u_s is None:
            self.u_s = 0.5 * self.vesc  # is 0.5 * u_esc in Genda and Abe 2003
        self.grid = []
        self.__setup_grid()
        if self.mass_atmosphere is not None:
            self.rho_0 = self.get_rho_0_given_atmosphere_mass(self.mass_atmosphere)
        self.__nondimensionalize_initial()

    def __setup_grid(self):
        current = 0
        s = ic.InitialConditionsSpherical()

        while current < self.num_shells:
            r = s.radius_initial(index=current, lambda_0=self.lambda_0, polytropic_exponent=self.gamma_a, r_0=self.r_0,
                                  total_shells=self.num_shells)
            P = s.pressure_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, p_0=self.P_0, radius=r,
                                    radius_0=self.r_0)
            T = s.temperature_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, T_0=self.T_0, radius=r,
                                       radius_0=self.r_0)
            rho = s.density_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, rho_0=self.rho_0,
                                     radius=r, radius_0=self.r_0)
            if current > 0:
                m = s.mass_initial(mass_last_index=self.grid[current - 1].mass, r_index=r,
                                    r_last_index=self.grid[current - 1].radius, rho_index=rho)
            else:
                m = 0
            v = s.velocity_initial(polytropic_exponent=self.gamma_a, T_index=T, m_a=self.m_a)
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
        self.grid[-1].pressure = 0.0  # vacuum boundary condition
        self.grid[-1].density = 0.0  # vacuum boundary condition
        self.grid[0].velocity = self.u_s  # the initial velocity

    def __nondimensionalize_initial(self):
        for p in self.grid:
            p.mass = nd.mass_nd(mass=p.mass, density_0=self.rho_0, radius_0=self.r_0)
            p.pressure = nd.pressure_nd(pressure=p.pressure, pressure_0=self.P_0)
            p.temperature = nd.temperature_nd(temperature=p.temperature, temperature_0=self.T_0)
            p.density = nd.density_nd(density=p.density, density_0=self.rho_0)
            p.radius = nd.radius_nd(radius=p.radius, radius_0=self.r_0)
            p.velocity = nd.velocity_nd(velocity=p.velocity, density_0=self.rho_0, gamma=self.gamma,
                                        pressure_0=self.P_0, c_s_0=self.c_s_0)

    def get_rho_0_given_atmosphere_mass(self, mass_atmosphere, increment_rho_0=0.1):
        """
        Given the mass of the atmosphere, find the initial density at the surface of the planet.
        """
        iterations = 0
        last_grid_mass = copy(self.grid[-1].mass)
        while min([self.grid[-1].mass, last_grid_mass]) < mass_atmosphere < max([self.grid[-1].mass, last_grid_mass]) or iterations == 0:
            last_grid_mass = copy(self.grid[-1].mass)
            self.lambda_0 = self.G * self.mass_planet * self.rho_0 / (self.r_0 * self.P_0)
            self.c_s_0 = sqrt(self.gamma_a * self.P_0 / self.rho_0)
            self.__setup_grid()
            if self.grid[-1].mass < mass_atmosphere:
                self.rho_0 += increment_rho_0
            else:
                self.rho_0 -= increment_rho_0
            # TODO: interpolate atmosphere mass
            # TODO: inherit new c_s_o and lambda_0 from here in main solver.py
            print("HERE", iterations, self.grid[-1].mass, mass_atmosphere)
            self.lambda_0 = self.G * self.mass_planet * self.rho_0 / (self.r_0 * self.P_0)
            self.c_s_0 = sqrt(self.gamma_a * self.P_0 / self.rho_0)
            iterations += 1
        return self.rho_0





class JetSystem:

    def __init__(self, num_shells, gamma_a, mass_planet, r_0, rho_0, P_0, T_0, m_a, gamma, u_s, jet_angle):
        self.jet_angle = jet_angle
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
        s = ic.InitialConditionsJet(jet_angle=self.jet_angle)

        while current < self.num_shells:
            r = s.radius_initial(index=current, lambda_0=self.lambda_0, polytropic_exponent=self.gamma_a, r_0=self.r_0,
                                  total_shells=self.num_shells)
            P = s.pressure_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, p_0=self.P_0, radius=r,
                                    radius_0=self.r_0)
            T = s.temperature_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, T_0=self.T_0, radius=r,
                                       radius_0=self.r_0)
            rho = s.density_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, rho_0=self.rho_0,
                                     radius=r, radius_0=self.r_0)
            if current > 0:
                m = s.mass_initial(mass_last_index=self.grid[current - 1].mass, r_index=r,
                                    r_last_index=self.grid[current - 1].radius, rho_index=rho)
            else:
                m = 0
            v = s.velocity_initial(polytropic_exponent=self.gamma_a, T_index=T, m_a=self.m_a)
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