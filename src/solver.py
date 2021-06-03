import os
import shutil
from math import pi, sqrt
from copy import copy
import matplotlib.pyplot as plt

from src import setup, output


class LagrangianSolver1D:
    """
    1D Lagragian differential equation solver.
    Based on Richtmyer and Morton 1967 as described by Genda and Abe 2003 Appendix A.
    https://www.sciencedirect.com/science/article/pii/S0019103503001015

    m is the mass constained within radius r and is the comoving coordinate.
    u is the velocity.
    p is the pressure.
    """

    def __init__(self, num_shells, gamma_a, r_0, P_0, T_0, m_a, gamma, mass_planet, R=8.314):
        num_shells += 1
        self.R = R
        rho_0 = P_0 * m_a / (T_0 * R)  # ideal gas
        self.system = setup.System(num_shells=num_shells, gamma_a=gamma_a, mass_planet=mass_planet, r_0=r_0,
                                   rho_0=rho_0,
                                   P_0=P_0, T_0=T_0, m_a=m_a, gamma=gamma)
        self.grid = self.system.grid
        self.time = 0
        self.dt = 0.001 / (r_0 / self.system.c_s_0)
        self.output_count = 0
        self.gamma = gamma
        self.q_coeff = 0.75
        self.num_shells = num_shells
        self.mass_planet = mass_planet
        self.lambda_0 = self.system.lambda_0
        self.outfile_dir = "/scratch/shull4/outputs"
        if self.outfile_dir in os.listdir(os.getcwd()):
            shutil.rmtree(self.outfile_dir)
        os.mkdir(self.outfile_dir)

    def solve(self, timesteps):
        for i in range(0, timesteps):
            print("At time {} ({}/{} steps)".format(self.__time_dimensional(), i, timesteps))
            grid_copy = copy(self.grid)
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
            if i % 100 == 0:
                mass_loss = self.mass_loss()
                output.write_state(
                    path=self.outfile_dir,
                    system=self.system,
                    fname=self.output_count,
                    time=self.time,
                    timestep=i,
                    grid=self.grid,
                    mass_fraction_loss=mass_loss,
                )
                self.output_count += 1
            self.grid = grid_copy
            self.time += self.dt
            self.__cfl_dt()
        print("Finished!")

    def __time_dimensional(self):
        return self.time * (self.system.r_0 / self.system.c_s_0)

    def __cfl_dt(self):
        dt = self.dt
        for i in range(0, self.num_shells - 1):
            criterion = (self.grid[i + 1].radius - self.grid[i].radius) / self.grid[i].velocity
            if dt > criterion:
                dt = 0.25 * ((self.grid[i + 1].radius - self.grid[i].radius) / self.grid[i].velocity)
        self.dt = dt

    def __solve_q(self, grid_copy):
        """
        Note: this solves for index i-1/2 and i+1/2, but we will store at integer indices.
        """
        for index in range(0, self.num_shells - 1):
            if self.grid[index + 1].velocity < self.grid[index].velocity:
                self.grid[index].q = self.numerical_viscosity_mid_forward(index=index)
            else:
                self.grid[index].q = 0.0

    def mass_loss(self):
        for p in self.grid:
            criterion = (p.velocity * self.system.c_s_0) / sqrt(
                ((2.0 * self.system.G * self.mass_planet) / (self.system.r_0 * p.radius)))
            if criterion > 1:
                mass_loss = 1.0 - (p.mass / self.grid[-1].mass)
                return mass_loss
        return 0.0

    def velocity(self, index):
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
        p = self.grid[index]
        a1 = 4 * pi * p.density
        a2 = (self.gamma * p.pressure) + ((self.gamma - 1) * self.grid[index].q)
        a3 = (((self.grid[index + 1].radius ** 2) * self.grid[index + 1].velocity) - (
                (p.radius ** 2) * p.velocity)) / (
                     self.grid[index + 1].mass - p.mass)
        return p.pressure - ((a1 * a2 * a3) * self.dt)

    def numerical_viscosity_mid_forward(self, index):
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
        p = self.grid[index]
        return p.radius + (velocity_tplus * self.dt)

    def density_mid_forward(self, index, radius_tplus_forward, radius_tplus):
        """
        Note: this solves for index i-1/2 and i+1/2, but we will store at integer indices.
        """
        p = self.grid[index]
        a1 = (3 / (4 * pi))
        a2 = (self.grid[index + 1].mass - p.mass) / ((radius_tplus_forward ** 3) - (radius_tplus ** 3))
        return a1 * a2

    def temperature(self, pressure_tplus, density_tplus):
        return (pressure_tplus / density_tplus) * self.system.m_a / self.R

    def plot_timestep(self, timestep, plot_separation=100, show=True):
        if timestep % plot_separation != 0:
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

        if show:
            plt.show()
