import os
import sys
import shutil
from math import pi, sqrt, tan
from copy import copy
import matplotlib.pyplot as plt

from src import setup, output


class LagrangianSolver1DSpherical:
    """
    1D Lagragian differential equation solver.
    Based on Richtmyer and Morton 1967 as described by Genda and Abe 2003 Appendix A.
    https://www.sciencedirect.com/science/article/pii/S0019103503001015

    m is the mass contained within radius r and is the co-moving coordinate.
    u is the velocity.
    p is the pressure.
    """

    def __init__(self, num_shells, gamma_a, r_0, P_0, T_0, m_a, gamma, mass_planet, u_s, R=8.314,
                 outfile_dir="/scratch/shull4/outputs", output_file_interval=1000, use_cfl=True, save_figs=False,
                 show_figs=False, **kwargs):
        num_shells += 1
        self.R = R
        self.P_0 = P_0
        self.rho_0 = kwargs.get("rho_0", P_0 * m_a / (T_0 * R))  # ideal gas
        self.mass_atmosphere = kwargs.get("mass_atmosphere", None)
        if self.mass_atmosphere is None:
            self.system = setup.SphericalSystem(num_shells=num_shells, gamma_a=gamma_a, mass_planet=mass_planet,
                                                r_0=r_0,
                                                rho_0=self.rho_0,
                                                P_0=self.P_0, T_0=T_0, m_a=m_a, gamma=gamma, u_s=u_s)
        else:  # iteratively solve for rho_0 given mass_atmosphere
            self.system = setup.SphericalSystem(num_shells=num_shells, gamma_a=gamma_a, mass_planet=mass_planet,
                                                r_0=r_0,
                                                rho_0=self.rho_0,
                                                P_0=self.P_0, T_0=T_0, m_a=m_a, gamma=gamma, u_s=u_s,
                                                mass_atmosphere=self.mass_atmosphere)
            self.rho_0 = self.system.rho_0
            self.P_0 = self.system.P_0
        self.grid = self.system.grid
        self.time = 0
        self.dt = 0.001 / (r_0 / self.system.c_s_0)
        self.dt_0 = copy(self.dt)
        self.output_count = 0
        self.gamma = gamma  # escape parameter, related to the ratio of gravidational energy required for escape from the planet and the thermal energy of the atmosphere
        self.q_coeff = 0.75  # numerical viscosity coefficient
        self.num_shells = num_shells
        self.mass_planet = mass_planet
        self.lambda_0 = self.system.lambda_0
        self.output_file_interval = output_file_interval
        self.outfile_dir = outfile_dir
        self.cfl_coefficient = kwargs.get("cfl_coefficient", 0.25)
        self.__use_cfl = use_cfl
        self.plot_separation = kwargs.get("plot_separation", 100)
        if self.outfile_dir in os.listdir(os.getcwd()):
            shutil.rmtree(self.outfile_dir)
        os.mkdir(self.outfile_dir)
        self.show_figs = show_figs
        self.save_figs = save_figs
        if self.save_figs:
            self.fig_save_path = kwargs.get("fig_save_path", 'plots')
            if self.fig_save_path in os.listdir(os.getcwd()):
                shutil.rmtree(self.fig_save_path)
            os.mkdir(self.fig_save_path)
        self.iteration = 0

    def solve(self, max_time: float):
        print("Solving...")
        dimensional_time = 0.0
        while dimensional_time <= max_time:
        # while self.iteration < 100001:
            if self.iteration != 0:
                dimensional_time = self.__time_dimensional()
                grid_copy = copy(self.grid)
                if self.iteration % 500 == 0:
                    print(
                        "At time {} (max: {} sec.) ({} iterations)".format(dimensional_time, max_time, self.iteration))
                self.__solve_q(grid_copy=grid_copy)
                for index, p in enumerate(grid_copy):
                    p.velocity = self.velocity(index=index)
                grid_copy[-1].velocity = grid_copy[-2].velocity
                for index, p in enumerate(grid_copy):
                    p.radius = self.radius(index=index, velocity_tplus=p.velocity)
                for index, p in enumerate(grid_copy):
                    if index < len(self.grid) - 1:
                        p.pressure = self.pressure(index=index)
                        p.density = self.density_mid_forward(index=index,
                                                             radius_tplus=p.radius,
                                                             radius_tplus_forward=grid_copy[
                                                                 index + 1].radius)
                        p.temperature = self.temperature(pressure_tplus=p.pressure,
                                                         density_tplus=p.density)
                grid_copy[-1].pressure = 0.0
                grid_copy[-1].density = 0.0
                self.grid = grid_copy
                self.__cfl_dt()
            if self.iteration % self.output_file_interval == 0 or self.iteration == 0:
                mass_loss = self.mass_loss()  # assess atmospheric mass loss
                output.write_state(
                    path=self.outfile_dir,
                    system=self.system,
                    fname=self.output_count,
                    time=dimensional_time,
                    timestep=self.iteration,
                    grid=self.grid,
                    mass_fraction_loss=mass_loss,
                )  # write output files
                self.output_count += 1  # increment output count
            self.time += self.dt
            self.plot_timestep(timestep=self.iteration)
            self.iteration += 1
        print("Finished!")

    def __time_dimensional(self):
        return self.time * (self.system.r_0 / self.system.c_s_0)

    def __cfl_dt(self):
        """
        Adjust timestep to satisfy CFL condition.
        """
        if not self.__use_cfl:
            self.dt = self.dt_0
            return self.dt
        # dt = self.dt
        # self.dt = 0.001 / (self.system.r_0 / self.system.c_s_0)
        self.dt = copy(self.dt_0)
        dt = copy(self.dt)
        for i in range(0, self.num_shells - 1):
            criterion = (self.grid[i + 1].radius - self.grid[i].radius) / self.grid[i].velocity
            criterion = abs(criterion)
            if dt > criterion:
                dt = self.cfl_coefficient * criterion
        self.dt = dt
        return self.dt

    def __solve_q(self, grid_copy):
        """
        Numerical viscosity.
        Note: this solves for index i-1/2 and i+1/2, but we will store at integer indices.
        """
        for index in range(0, self.num_shells - 1):
            if self.grid[index + 1].velocity < self.grid[index].velocity:
                self.grid[index].q = self.numerical_viscosity_mid_forward(index=index)
            else:
                self.grid[index].q = 0.0

    def mass_loss(self):
        """
        Returns the mass fraction of the atmosphere lost from the planet.
        The criterion is that the atmosphere must be moving faster than the escape velocity.
        criterion = (u * c_s) / sqrt(2 * G * M_planet) / (r_0 * r)
        """
        for p in self.grid:
            criterion = (p.velocity * self.system.c_s_0) / sqrt(
                ((2.0 * self.system.G * self.mass_planet) / (self.system.r_0 * p.radius)))
            if criterion > 1:
                mass_loss = 1.0 - (p.mass / self.grid[-1].mass)
                return mass_loss
        return 0.0

    def velocity(self, index):
        """
        Returns the gridded velocity of the atmosphere at a given grid index.
        """
        p = self.grid[index]
        if index == 0:
            return self.grid[index].velocity - (self.lambda_0 / self.gamma / (self.grid[index].radius ** 2)) * self.dt
        elif index < len(self.grid) - 1:
            m_forward = self.grid[index + 1].mass
            m_backwards = self.grid[index - 1].mass
            a1 = (8 * pi) / self.gamma
            a2 = p.radius ** 2
            a3 = (p.pressure - self.grid[index - 1].pressure + p.q - self.grid[
                index - 1].q) / (m_forward - m_backwards)
            a4 = (self.lambda_0 / (self.gamma * (p.radius ** 2)))
            return p.velocity - (((a1 * a2 * a3) + a4) * self.dt)
        # else:
        #     return self.grid[index].velocity

    def pressure(self, index):
        """
        Returns the pressure of the atmosphere at a given grid index.
        """
        p = self.grid[index]
        a1 = 4 * pi * p.density
        a2 = (self.gamma * p.pressure) + ((self.gamma - 1) * self.grid[index].q)
        a3 = (((self.grid[index + 1].radius ** 2) * self.grid[index + 1].velocity) - (
                (p.radius ** 2) * p.velocity)) / (
                     self.grid[index + 1].mass - p.mass)
        return p.pressure - ((a1 * a2 * a3) * self.dt)

    def numerical_viscosity_mid_forward(self, index):
        """
        Forward finite difference for the numerical viscosity.
        """
        p = self.grid[index]
        if index < len(self.grid) - 1:
            if (self.grid[index + 1].velocity - p.velocity) < 0:
                a1 = - self.q_coeff * self.gamma * p.density * (self.grid[index + 1].velocity - p.velocity)
                a2 = (sqrt(p.pressure / p.density) - (
                        ((self.gamma + 1) / 2) * (self.grid[index + 1].velocity - p.velocity)))
                return a1 * a2
            else:
                return 0
        return 0

    def radius(self, index, velocity_tplus):
        """
        Calculates the radial extent of the atmosphere at a given grid index.
        """
        p = self.grid[index]
        return p.radius + (velocity_tplus * self.dt)

    def density_mid_forward(self, index, radius_tplus_forward, radius_tplus):
        """
        The forward finite difference for the density.
        Note: this solves for index i-1/2 and i+1/2, but we will store at integer indices.
        """
        p = self.grid[index]
        a1 = (3 / (4 * pi))
        a2 = (self.grid[index + 1].mass - p.mass) / ((radius_tplus_forward ** 3) - (radius_tplus ** 3))
        return a1 * a2

    def temperature(self, pressure_tplus, density_tplus):
        """
        Calculates the temperature of the atmosphere at a given grid index.
        """
        return (pressure_tplus / density_tplus) * self.system.m_a / self.R

    def plot_timestep(self, timestep):
        if timestep % self.plot_separation != 0:
            return None
        fig = plt.figure(figsize=(16, 9))
        ax_pressure = fig.add_subplot(221)
        ax_density = fig.add_subplot(222)
        ax_velocity = fig.add_subplot(223)
        ax_mass = fig.add_subplot(224)

        g = self.grid
        rr = [i.radius for i in g]
        rho = [i.density for i in g]
        mm = [i.mass for i in g]
        pp = [i.pressure for i in g]
        uu = [i.velocity for i in g]

        ax_pressure.plot(
            rr,
            pp,
            linewidth=2.0,
            color='black'
        )
        ax_density.plot(
            rr,
            rho,
            linewidth=2.0,
            color='black'
        )
        ax_velocity.plot(
            rr,
            [i * self.system.c_s_0 / self.system.vesc for i in uu],
            linewidth=2.0,
            color='black'
        )
        ax_mass.plot(
            rr,
            mm,
            linewidth=2.0,
            color='black'
        )
        ax_pressure.set_xlabel("r / r0")
        ax_density.set_xlabel("r / r0")
        ax_velocity.set_xlabel("r / r0")
        ax_mass.set_xlabel("r / r0")
        ax_pressure.set_ylabel("P / P0")
        ax_density.set_ylabel("rho / rho0")
        ax_velocity.set_ylabel("u / u_esc")
        ax_mass.set_ylabel("m / (r_0^3 rho_0)")
        ax_pressure.set_title("Pressure (IC)")
        ax_density.set_title("Density (IC)")
        ax_velocity.set_title("Velocity (IC)")
        ax_mass.set_title("Mass (IC)")
        ax_pressure.grid()
        ax_density.grid()
        ax_velocity.grid()
        ax_mass.grid()
        fig.suptitle("Time: {}".format(self.__time_dimensional()))

        if self.show_figs:
            plt.show()
        if self.save_figs:
            plt.savefig("{}/{}.png".format(self.fig_save_path, timestep))


class LagrangianSolverJet(LagrangianSolver1DSpherical):
    """
    Solves for a jet as modelled by a basic cone.
    """

    def __init__(self, **kwargs):
        super(LagrangianSolverJet, self).__init__(**kwargs)
        self.theta = self.degrees_to_radians(kwargs.get('jet_angle', 45.0))  # expansion angle of the jet in degrees
        self.system = setup.JetSystem(num_shells=kwargs.get('num_shells'), gamma_a=kwargs.get('gamma_a'),
                                      mass_planet=kwargs.get('mass_planet'), r_0=kwargs.get('r_0'),
                                      rho_0=self.rho_0, P_0=kwargs.get('P_0'), T_0=kwargs.get('T_0'),
                                      m_a=kwargs.get('m_a'), gamma=kwargs.get('gamma'), u_s=kwargs.get('u_s'),
                                      jet_angle=self.theta)

    def degrees_to_radians(self, degrees):
        """
        Converts degrees to radians.
        """
        return degrees * (pi / 180)

    def velocity(self, index):
        p = self.grid[index]
        if index == 0:
            return self.grid[index].velocity - (self.lambda_0 / self.gamma / (self.grid[index].radius ** 2)) * self.dt
        elif index < len(self.grid) - 1:
            m_forward = self.grid[index + 1].mass
            m_backwards = self.grid[index - 1].mass
            a1 = 2 * (pi / self.gamma)
            a2 = (p.radius ** 2) * (tan(self.theta) ** 2)
            a3 = (p.pressure - self.grid[index - 1].pressure + p.q - self.grid[
                index - 1].q) / (m_forward - m_backwards)
            a4 = (self.lambda_0 / (self.gamma * (p.radius ** 2)))
            v = p.velocity - (((a1 * a2 * a3) + a4) * self.dt)
            return v

    def pressure(self, index):
        p = self.grid[index]
        a1 = pi * p.density * (tan(self.theta) ** 2)
        a2 = (self.gamma * p.pressure) + ((self.gamma - 1) * p.q)
        a3 = (((self.grid[index + 1].radius ** 2) * self.grid[index + 1].velocity) - (
                (p.radius ** 2) * p.velocity)) / (
                     self.grid[index + 1].mass - p.mass)
        pres = p.pressure - ((a1 * a2 * a3) * self.dt)
        return pres

    def radius(self, index, velocity_tplus):
        """
        Calculates the radial extent of the atmosphere at a given grid index.
        This is "z" in documentations, or the cone's height.
        """
        p = self.grid[index]
        return p.radius + (velocity_tplus * self.dt)

    def density_mid_forward(self, index, radius_tplus_forward, radius_tplus):
        """
        The forward finite difference for the density.
        Note: this solves for index i-1/2 and i+1/2, but we will store at integer indices.
        """
        p = self.grid[index]
        a1 = 1 / (pi * (radius_tplus ** 2) * (tan(self.theta) ** 2))
        a2 = (self.grid[index + 1].mass - p.mass) / ((radius_tplus_forward ** 3) - (radius_tplus ** 3))
        return a1 * a2
