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

    def __init__(self, num_shells, gamma_a, lambda_0, r_0, rho_0, P_0, T_0, P_max, rho_max, m_initial, v_0, m_a, gamma):
        self.num_shells = num_shells
        self.G = 6.674 * 10 ** -11

        self.gamma_a = gamma_a
        self.gamma = gamma
        self.lambda_0 = lambda_0
        self.r_0 = r_0
        self.rho_0 = rho_0
        self.P_0 = P_0
        self.T_0 = T_0
        self.v_0 = v_0
        self.m_a = m_a

        self.P_max = P_max
        self.rho_max = rho_max
        self.m_initial = m_initial

        self.grid = []
        self.__setup_grid()
        self.__nondimensionalize_initial()

    def __setup_grid(self):
        current = 1
        point_init = Point(
            id=current - 1,
            q=None,
            mass=self.m_initial,
            pressure=self.P_0,
            temperature=self.T_0,
            density=self.rho_0,
            radius=self.r_0,
            velocity=self.v_0
        )
        self.grid.append(point_init)

        while current <= self.num_shells:
            r = ic.radius_initial(index=current, lambda_0=self.lambda_0, polytropic_exponent=self.gamma_a, r_0=self.r_0,
                                  total_shells=self.num_shells)
            P = ic.pressure_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, p_0=self.P_0, radius=r,
                                    radius_0=self.r_0)
            T = ic.temperature_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, T_0=self.T_0, radius=r,
                                       radius_0=self.r_0)
            rho = ic.density_initial(polytropic_exponent=self.gamma_a, lambda_0=self.lambda_0, rho_0=self.rho_0,
                                     radius=r, radius_0=self.r_0)
            m = ic.mass_initial(mass_last_index=self.grid[current - 1].mass, r_index=r,
                                r_last_index=self.grid[current - 1].radius, rho_index=rho)
            v = ic.velocity_initial(polytropic_exponent=self.gamma_a, T_index=T, m_a=self.m_a)
            if current == self.num_shells:
                rho = 0.0
                P = 0.0
            q = None
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

    def __nondimensionalize_initial(self):
        for p in self.grid:
            p.mass = nd.mass_nd(mass=p.mass, density_0=self.rho_0, radius_0=self.rho_0)
            p.pressure = nd.pressure_nd(pressure=p.pressure, pressure_0=self.P_0)
            p.temperature = nd.temperature_nd(temperature=p.temperature, temperature_0=self.T_0)
            p.rho = nd.density_nd(density=p.density, density_0=self.rho_0)
            p.radius = nd.radius_nd(radius=p.radius, radius_0=self.r_0)
            p.velocity = nd.velocity_nd(velocity=p.velocity, density_0=self.rho_0, gamma=self.gamma,
                                        pressure_0=self.P_0)
