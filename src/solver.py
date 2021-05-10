from math import pi, sqrt


class Point:

    def __init__(self):
        self.radius = 0
        self.velocity = 0
        self.mass = 0
        self.pressure = 0
        self.density = 0
        self.num_visc = 0


class LagrangianSolver1D:
    """
    1D Lagragian differential equation solver.
    Based on Richtmyer and Morton 1967 as described by Genda and Abe 2003 Appendix A

    m is the mass constained within radius r and is the comoving coordinate.
    u is the velocity.
    p is the pressure.
    """

    def __init__(self, inner_boundary, outer_boundary, timestep, num_shells=1000):
        # TODO: Make a single class to describe the Lagrangian coordinates with all attributes
        self.grid = []
        self.inner_boundary = inner_boundary
        self.outer_boundary = outer_boundary
        self.dt = timestep
        self.lambda_0 = 0
        self.gamma = 0
        self.num_visc_coeff = 0.75

    def velocity(self, index):
        p = self.grid[index]
        p_mid_forward = (self.grid[index + 1].pressure - p.pressure) / 2.0
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
