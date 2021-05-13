from math import pi, sqrt


class Point:

    def __init__(self, radius, velocity, mass, pressure, density, temperature):
        self.radius = radius
        self.velocity = velocity
        self.mass = mass
        self.pressure = pressure
        self.density = density
        self.num_visc = 0
        self.temperature = temperature


class LagrangianSolver1D:
    """
    1D Lagragian differential equation solver.
    Based on Richtmyer and Morton 1967 as described by Genda and Abe 2003 Appendix A

    m is the mass constained within radius r and is the comoving coordinate.
    u is the velocity.
    p is the pressure.
    """

    def __init__(self, timestep, num_shells=1000):
        # TODO: Make a single class to describe the Lagrangian coordinates with all attributes
        self.grid = []
        self.dt = timestep
        self.gamma = 0
        self.num_visc_coeff = 0.75
        self.num_shells = num_shells
        self.lambda_0 = 0

    def construct_grid(self, inner_boundary_velocity, inner_boundary_mass, inner_boundary_radius,
                       density, temperature, pressure, shell_thickness=1):
        self.lambda_0 = self.escape_parameter(

        )
        inner_boundary = Point(
            radius=inner_boundary_radius,
            mass=inner_boundary_mass,
            velocity=inner_boundary_velocity,
            pressure=pressure,
            density=density,
            temperature=temperature
        )
        self.grid.append(inner_boundary)
        index = 1
        while index <= self.num_shells:
            radius = self.grid[index - 1] + shell_thickness
            mass = 4 * pi * (radius ** 2) * density
            p = Point(
                radius=radius,
                mass=mass,
                velocity=0,
                pressure=pressure,
                density=density,
                temperature=temperature
            )
            self.grid.append(p)
            index += 1

    def escape_parameter(self, mass_planet, initial_planet_radius, pressure_nd, density_nd):
        G = 6.674 * 10 ** -11
        # alternative form: lambda_0 = G M m_a / k T_0 r_0
        # where m_a is the mass of the molecule, k is the Boltzmann constant, and T is the
        # initial temperature at the planetary surface
        return (G * mass_planet * density_nd) / (initial_planet_radius * pressure_nd)

    def velocity(self, index):
        p = self.grid[index]
        p_mid_forward = (self.grid[index + 1].pressure - p.pressure) / 2.0
        if index == self.num_shells - 1:
            p_mid_forward = 0
        p_mid_backward = (p.pressure - self.grid[index - 1].pressure) / 2.0
        q_mid_forward = (self.grid[index + 1].num_visc - p.num_visc) / 2.0
        q_mid_backward = (p.num_visc - self.grid[index - 1].num_visc) / 2.0
        m_forward = self.grid[index + 1].mass - p.mass
        m_backwards = p.mass - self.grid[index - 1].mass

        a1 = (8 * pi) / self.gamma
        a2 = (p.radius) ** 2
        a3 = (p_mid_forward - p_mid_backward + q_mid_forward - q_mid_backward) / (m_forward - m_backwards)
        a4 = (self.lambda_0 / (self.gamma * (p.radius ** 2)))

        return p.velocity - (((a1 * a2 * a3) + a4) * self.dt)

    def pressure(self, index):
        p = self.grid[index]
        rho_mid_forward = (self.grid[index + 1].density - p.density) / 2.0
        p_mid_forward = (self.grid[index + 1].pressure - p.pressure) / 2.0
        if index == self.num_shells - 1:
            p_mid_forward = 0
        q_mid_forward = (self.grid[index + 1].num_visc - p.num_visc) / 2.0
        a1 = 4 * pi * rho_mid_forward
        a2 = (self.gamma * p_mid_forward) + ((self.gamma - 1) * q_mid_forward)
        a3 = (((self.grid[index + 1].radius ** 2) * self.grid[index + 1].velocity) - (
                (p.radius ** 2) * p.velocity)) / (
                     self.grid[index + 1].mass - p.mass)
        return p_mid_forward - ((a1 * a2 * a3) * self.dt)

    def numerical_viscosity_mid_forward(self, index):
        p = self.grid[index]
        if ((self.grid[index + 1].velocity - p.velocity)) < 0:
            rho_mid_forward = (self.grid[index + 1].density - p.density) / 2.0
            p_mid_forward = (self.grid[index + 1].pressure - p.pressure) / 2.0
            if index == self.num_shells - 1:
                p_mid_forward = 0
            a1 = - self.num_visc_coeff * self.gamma * (self.grid[index + 1].velocity - p.velocity)
            a2 = (sqrt(p_mid_forward / rho_mid_forward) - (
                    ((self.gamma + 1) / 2) * (self.grid[index + 1].velocity - p.velocity)))
            return a1 * a2
        else:
            return 0

    def radius(self, index, velocity_tplus):
        p = self.grid[index]
        return p.radius + (velocity_tplus * self.dt)

    def density_mid_forward(self, index, radius_tplus_forward, radius_tplus):
        p = self.grid[index]
        a1 = (3 / (4 * pi))
        a2 = (self.grid[index + 1].mass - p.mass) / ((radius_tplus_forward ** 3) - (radius_tplus ** 3))
        return a1 * a2
