import os

from src import initial_conditions as ic
from src import nondimensional as nd
from src import setup

from math import pi, sqrt
import numpy as np
from copy import copy
import sys


class LagrangianSolver1D:
    """
    1D Lagragian differential equation solver.
    Based on Richtmyer and Morton 1967 as described by Genda and Abe 2003 Appendix A.
    https://www.sciencedirect.com/science/article/pii/S0019103503001015

    m is the mass constained within radius r and is the comoving coordinate.
    u is the velocity.
    p is the pressure.
    """

    def __init__(self, timestep, num_shells, gamma_a, lambda_0, r_0, rho_0, P_0, T_0, P_max, rho_max, m_initial, v_0,
                 m_a, gamma, c_s_0, mass_planet):
        self.system = setup.System(num_shells=num_shells, gamma_a=gamma_a, lambda_0=lambda_0, r_0=r_0, rho_0=rho_0,
                                   P_0=P_0, T_0=T_0, P_max=P_max, rho_max=rho_max, m_initial=m_initial, v_0=v_0,
                                   m_a=m_a,
                                   gamma=gamma, c_s_0=c_s_0)
        self.grid = self.system.grid
        self.time = 0
        self.dt = timestep
        self.gamma = gamma
        self.q_coeff = 0.75
        self.num_shells = num_shells
        self.lambda_0 = lambda_0
        self.c_s_0 = c_s_0
        self.mass_planet = mass_planet
        self.G = 6.674 * 10 ** -11
        self.v_esc = sqrt((2 * self.G * self.mass_planet) / self.system.rho_0)
        self.fname = "lagrangian_solver.csv"
        if self.fname in os.listdir(os.getcwd()):
            os.remove(self.fname)
        self.outfile = open(self.fname, 'w')
        self.outfile.write("time,mass_loss\n")

    def __solve_q(self, grid_copy):
        """
        Note: this solves for index i-1/2 and i+1/2, but we will store at integer indices.
        """
        for index in np.arange(0, len(self.grid) - 1, 1):
            if self.grid[index].q is None:
                self.grid[index].q = 0.0
            if self.grid[index + 1].q is None:
                self.grid[index + 1].q = 0.0
            if self.grid[index + 1].velocity < self.grid[index].velocity:
                self.grid[index].q = self.numerical_viscosity_mid_forward(index=index)
            else:
                self.grid[index].q = 0.0

    def solve(self, timesteps):
        for i in np.arange(0, timesteps, 1):
            print("At time {} ({}/{} steps)".format(self.time, i, timesteps))
            grid_copy = copy(self.grid)
            self.__solve_q(grid_copy=grid_copy)
            for index, p in enumerate(grid_copy):
                grid_copy[index].velocity = self.velocity(index=index)
                grid_copy[index].radius = self.radius(index=index, velocity_tplus=grid_copy[index].velocity)
            for index, p in enumerate(grid_copy):
                if index < len(self.grid) - 2:
                    grid_copy[index].pressure = self.pressure(index=index)
                    grid_copy[index].density = self.density_mid_forward(index=index,
                                                                        radius_tplus=grid_copy[index].radius,
                                                                        radius_tplus_forward=grid_copy[
                                                                            index + 1].radius)
            self.grid = grid_copy
            self.mass_loss()
            self.time += self.dt
        self.outfile.close()
        print("Finished!")

    def mass_loss(self):

        for p in self.grid:
            mass_loss = None
            criterion = (p.velocity * self.c_s_0) / sqrt(
                ((2 * self.G * self.mass_planet) / (self.system.r_0 * p.radius)))
            if criterion > 1:
                mass_loss = p.mass / self.grid[-1].mass
                print("MASS LOSS ", mass_loss)
            if mass_loss is not None:
                self.outfile.write("{},{}\n".format(self.time, mass_loss))

    def velocity(self, index):
        p = self.grid[index]
        if index == 0:
            return self.grid[index].velocity - (self.lambda_0 / (self.grid[index].radius ** 2)) * self.dt
        elif index < len(self.grid) - 2:
            m_forward = self.grid[index + 1].mass
            m_backwards = self.grid[index - 1].mass
            a1 = (8 * pi) / self.gamma
            a2 = (p.radius) ** 2
            a3 = (self.grid[index].pressure - self.grid[index - 1].pressure + self.grid[index].q - self.grid[
                index - 1].q) / (m_forward - m_backwards)
            a4 = (self.lambda_0 / (self.gamma * (p.radius ** 2)))

            return p.velocity - (((a1 * a2 * a3) + a4) * self.dt)
        else:
            return self.grid[index].velocity

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
        if 1 <= index < len(self.grid) - 2:
            if (self.grid[index + 1].velocity - p.velocity) < 0:
                # print(index, self.grid[index].pressure, self.grid[index].density)
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
        print((radius_tplus_forward ** 3) - (radius_tplus ** 3))
        a2 = (self.grid[index + 1].mass - p.mass) / ((radius_tplus_forward ** 3) - (radius_tplus ** 3))
        return a1 * a2
