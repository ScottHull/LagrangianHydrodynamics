import os
import sys

import src.initial_conditions as ic
import src.nondimensional as nd

from math import sqrt
from copy import copy
from scipy.interpolate import interp1d


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

    def __init__(self, num_shells, gamma_a, mass_planet, r_0, rho_0, P_0, T_0, m_a, gamma, u_s, R=8.314, **kwargs):
        self.num_shells = num_shells
        self.G = 6.67408e-11
        self.R = R
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
            self.get_rho_0_given_atmosphere_mass(self.mass_atmosphere)
            print(
                "Initial conditions found!\n\tAtmosphere mass: {} (error: {} %)\n\tInitial density: {}\n\tInitial pressure: {}\n\t"
                "Initial temperature: {}\n\tInitial radius: {}\n\tInitial velocity: {}\n\tInitial lambda: {}\n\t"
                "Initial c_s: {}\n\tMean Atmosphere Mass: {}\n\tgamma: {}\n\tgamma_a: {}".format(
                    self.grid[-1].mass, abs((self.mass_atmosphere - self.grid[-1].mass) / self.grid[-1].mass) * 100.0,
                    self.rho_0, self.P_0, self.T_0, self.r_0, self.u_s, self.lambda_0, self.c_s_0, self.m_a,
                    self.gamma, self.gamma_a)

            )
        else:
            print(
                "Initial conditions found!\n\tAtmosphere mass: {}\n\tInitial density: {}\n\tInitial pressure: {}\n\t"
                "Initial temperature: {}\n\tInitial radius: {}\n\tInitial velocity: {}\n\tInitial lambda: {}\n\t"
                "Initial c_s: {}\n\tMean Atmosphere Mass: {}\n\tgamma: {}\n\tgamma_a: {}".format(
                    self.grid[-1].mass, self.rho_0, self.P_0, self.T_0, self.r_0, self.u_s,
                    self.lambda_0, self.c_s_0, self.m_a, self.gamma, self.gamma_a)

            )
        self.__nondimensionalize_initial()
        self.write_initial_conditions_to_file(fname=kwargs.get("ic_fname", "initial_conditions.txt"))

    def __setup_grid(self):
        self.grid = []  # reset the grid
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
            if current > 0:  # first shell is the planet and contains no atmosphere mass
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

    def get_rho_0_given_atmosphere_mass(self, mass_atmosphere, increment_rho_0=0.01):
        """
        Given the mass of the atmosphere, find the initial density at the surface of the planet.
        """
        print("Iteratively solving for initial conditions that match atmosphere mass. This may take a while...")
        print("[!] Warning: P_0, lambda_0, and c_s_0 will change from initially given values.")
        iterations = 0
        last_grid_mass = copy(self.grid[-1].mass)
        last_rho_0 = copy(self.rho_0)
        # begin iterative solution for rho_0 and corresponding initial conditions & grid
        while not min([self.grid[-1].mass, last_grid_mass]) < mass_atmosphere < max(
                [self.grid[-1].mass, last_grid_mass]) or iterations == 0:
            last_grid_mass = copy(self.grid[-1].mass)
            last_rho_0 = copy(self.rho_0)
            if self.grid[-1].mass < mass_atmosphere:
                self.rho_0 += increment_rho_0
            else:
                self.rho_0 -= increment_rho_0
            # adjust other initial conditions that are dependent on density
            # ideal gas: rho_0 = (P_0 * m_a) / (R * T_0) --> P_0 = (rho_0 * R * T_0) / m_a
            self.P_0 = (self.rho_0 * self.R * self.T_0) / self.m_a  # ideal gas condition
            self.lambda_0 = self.G * self.mass_planet * self.rho_0 / (self.r_0 * self.P_0)
            self.c_s_0 = sqrt(self.gamma_a * self.P_0 / self.rho_0)
            self.__setup_grid()
            iterations += 1
        # we have a close enough solution...now interpolate to get the exact mass, re-solve initial conditions
        interp = interp1d([last_grid_mass, self.grid[-1].mass], [last_rho_0, self.rho_0])  # the interpolation function
        self.rho_0 = interp(mass_atmosphere)  # interpolate the density at the atmosphere mass
        self.P_0 = (self.rho_0 * self.R * self.T_0) / self.m_a  # ideal gas condition
        self.lambda_0 = self.G * self.mass_planet * self.rho_0 / (self.r_0 * self.P_0)
        self.c_s_0 = sqrt(self.gamma_a * self.P_0 / self.rho_0)
        self.__setup_grid()
        return self.rho_0

    def write_initial_conditions_to_file(self, fname, run_name='hydrodynamic_escape'):
        """
        Write the initial conditions to a file.
        """
        if os.path.exists(fname):
            os.remove(fname)
        with open(fname, 'w') as f:
            # write the initial conditions
            f.write(f"Initial Conditions: {run_name}\n")
            f.write("rho_0 = {}\n".format(self.rho_0))
            f.write("P_0 = {}\n".format(self.P_0))
            f.write("T_0 = {}\n".format(self.T_0))
            f.write("lambda_0 = {}\n".format(self.lambda_0))
            f.write("c_s_0 = {}\n".format(self.c_s_0))
            f.write("m_a = {}\n".format(self.m_a))
            f.write("u_s = {}\n".format(self.u_s))
            f.write("gamma = {}\n".format(self.gamma))
            f.write("gamma_a = {}\n".format(self.gamma_a))
            f.write("mass_planet = {}\n".format(self.mass_planet))
            f.write("r_0 = {}\n".format(self.r_0))
            f.write("mass_atmosphere = {}\n".format(self.grid[-1].mass))
        f.close()

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
