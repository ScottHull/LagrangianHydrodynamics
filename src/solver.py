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

    def __init__(self, num_shells, gamma_a, r_0, P_0, T_0, m_a, gamma, mass_planet, R=8.3):
        num_shells += 1
        rho_0 = P_0 * m_a / (T_0 * R)  # ideal gas
        self.system = setup.System(num_shells=num_shells, gamma_a=gamma_a, mass_planet=mass_planet, r_0=r_0,
                                   rho_0=rho_0,
                                   P_0=P_0, T_0=T_0, m_a=m_a, gamma=gamma)
        self.grid = self.system.grid
        self.time = 0
        self.dt = 0.001 / r_0 / self.system.c_s_0
        self.gamma = gamma
        self.q_coeff = 0.75
        self.num_shells = num_shells
        self.mass_planet = mass_planet
        self.lambda_0 = self.system.lambda_0
        self.v_esc = sqrt((2 * self.system.G * self.mass_planet) / self.system.rho_0)
        self.fname = "lagrangian_solver.csv"
        if self.fname in os.listdir(os.getcwd()):
            os.remove(self.fname)
        self.outfile = open(self.fname, 'w')
        self.outfile.write("time,mass_loss\n")

    def solve(self, timesteps):
        for i in range(0, timesteps):
            print("At time {} ({}/{} steps)".format(self.time, i, timesteps))
            grid_copy = copy(self.grid)
            print([p.velocity for p in grid_copy])
            self.__solve_q(grid_copy=grid_copy)
            for index, p in enumerate(grid_copy):
                p.velocity = self.velocity(index=index)
                p.radius = self.radius(index=index, velocity_tplus=p.velocity)
            grid_copy[-1].velocity = grid_copy[-2].velocity
            print(grid_copy[-1].velocity, grid_copy[-2].velocity, grid_copy[-3].velocity)
            for index, p in enumerate(grid_copy):
                if index < len(self.grid) - 2:
                    p.pressure = self.pressure(index=index)
                    p.density = self.density_mid_forward(index=index,
                                                         radius_tplus=p.radius,
                                                         radius_tplus_forward=grid_copy[
                                                             index + 1].radius)
            grid_copy[-1].pressure = 0.0
            grid_copy[-1].density = 0.0
            self.grid = grid_copy
            # self.mass_loss()
            self.time += self.dt
        self.outfile.close()
        print("Finished!")

    def __solve_q(self, grid_copy):
        """
        Note: this solves for index i-1/2 and i+1/2, but we will store at integer indices.
        """
        for index in range(0, self.num_shells):
            if self.grid[index + 1].velocity < self.grid[index].velocity:
                self.grid[index].q = self.numerical_viscosity_mid_forward(index=index)
            else:
                self.grid[index].q = 0.0

    def mass_loss(self):
        for p in self.grid:
            mass_loss = None
            criterion = (p.velocity * self.system.c_s_0) / sqrt(
                ((2 * self.system.G * self.mass_planet) / (self.system.r_0 * p.radius)))
            if criterion > 1:
                mass_loss = p.mass / self.grid[-1].mass
                print("MASS LOSS ", mass_loss)
            if mass_loss is not None:
                self.outfile.write("{},{}\n".format(self.time, mass_loss))

    def velocity(self, index):
        p = self.grid[index]
        if index == 0:
            return self.grid[index].velocity - (self.lambda_0 / self.gamma / (self.grid[index].radius ** 2)) * self.dt
        elif index < len(self.grid) - 2:
            m_forward = self.grid[index + 1].mass
            m_backwards = self.grid[index - 1].mass
            a1 = (8 * pi) / self.gamma
            a2 = p.radius ** 2
            a3 = (p.pressure - self.grid[index - 1].pressure + p.q - self.grid[
                index - 1].q) / (m_forward - m_backwards)
            a4 = (self.lambda_0 / (self.gamma * (p.radius ** 2)))
            # print(self.grid[index + 1].mass, self.grid[index - 1].mass, self.gamma, p.pressure,
            #       self.grid[index - 1].pressure, p.q, self.grid[
            #           index - 1].q, self.dt, p.velocity - (((a1 * a2 * a3) + a4) * self.dt))
            print("~", index, self.grid[index + 1].mass, self.grid[index - 1].mass, p.radius, p.pressure,
                  self.grid[index - 1].pressure, p.q, self.grid[
                      index - 1].q, p.velocity, p.velocity - (((a1 * a2 * a3) + a4) * self.dt))
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
        if index < len(self.grid) - 1:
            if (self.grid[index + 1].velocity - p.velocity) < 0:
                print("***", [p.velocity for p in self.grid])
                print(index, self.grid[index + 1].velocity, p.velocity, p.density, p.pressure)
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
