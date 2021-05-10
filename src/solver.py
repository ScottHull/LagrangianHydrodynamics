from math import pi

class Point:

    def __init__(self):
        self.radius = 0
        self.velocity = 0
        self.mass = 0
        self.pressure = 0
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
        # TODO: Make a single class to describe the Lagrangian coordinates with all
        self.grid = []
        self.inner_boundary = inner_boundary
        self.outer_boundary = outer_boundary
        self.dt = timestep
        self.lambda_0 = 0
        self.gamma = 0


    def velocity(self, index):
        p = self.grid[index]
        p_mid_forward = (self.grid[index + 1].pressure - self.grid[index].pressure) / 2.0
        p_mid_backward = (self.grid[index].pressure - self.grid[index - 1].pressure) / 2.0
        q_mid_forward = (self.grid[index + 1].num_visc - self.grid[index].num_visc) / 2.0
        q_mid_backward = (self.grid[index].num_visc - self.grid[index - 1].num_visc) / 2.0
        m_forward = self.grid[index + 1].mass - self.grid[index].mass
        m_backwards = self.grid[index].mass - self.grid[index - 1].mass

        a1 = (8 * pi) / self.gamma
        a2 = (p.radius) ** 2
        a3 = (p_mid_forward - p_mid_backward + q_mid_forward - q_mid_backward) / (m_forward - m_backwards)
        a4 = (self.lambda_0 / (self.gamma * (p.radius ** 2)))

        return p.velocity - ((a1 * a2 * a3) + a4) * self.dt



