from src import nondimensional, solver

"""
Test of a plane-paralle, isothermal, isochoric, isobaric atmosphere.
"""

pressure = 1.01 * 10 ** 7  # Pa
temperature = 298  # K
density = 1.29  # kg/m3
velocity = 0

inner_boundary_radius = 0
inner_boundary_velocity = 0
inner_boundary_mass = 0
outer_boundary_pressure = 0  # p_(1000 + 1/2) where 1000 is the grid number


s = solver.LagrangianSolver1D(
    timestep=1
)
s.construct_grid(
    density=nondimensional.density_nd(density=density, density_0=density),
    pressure=nondimensional.pressure_nd(pressure=pressure, pressure_0=pressure),
    temperature=temperature,
    inner_boundary_mass=inner_boundary_mass,
    inner_boundary_radius=inner_boundary_radius,
    inner_boundary_velocity=inner_boundary_velocity,
\)
