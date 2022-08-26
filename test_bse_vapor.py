import sys

from src import solver
from src import initial_conditions as ic

from math import sqrt
import matplotlib.pyplot as plt


def dulong_petit_limit(N, T):
    """
    Returns the high-temperature Dulong-Petit Limit for heat capacity (c_v).
    :param N: number of moles of the atmosphere per unit mass
    :param T: temperature in K
    """
    k = 1.380649 * 10 ** -23  # Boltzmann constant, m2 kg s-2 K-1
    return 3 * N * k * T


def get_specific_heat_ratio(c_v):
    """
    For an ideal gas, c_p = c_v + R
    :param c_v: specific heat at constant volume
    :return: specific heat ratio, gamma = c_p / c_v
    """
    R = 8.314  # universal gas constant, J mol-1 K-1
    c_p = c_v + R  # specific heat at constant pressure
    return c_v / c_p


"""
 The Dulong-Petit Limit requires the average number of moles per unit mass.
 Dunite is Mg2SiO4.
 Atomic masses:
    Mg = 24.3050 g/mol
    Si = 28.0855 g/mol
    O = 15.9994 g/mol
 1 mole of dunite is (2 * 24.3050 + 1 * 28.0855 + 4 * 15.9994) g = 140.6931 g
 Therefore, the average number of moles per unit mass is 1 mol / 140.6931 g = 0.0071076691 mol/g
 Avogadro's number: 6.0221408e+23 mol/mol
 The avg. number of atoms is therefore 0.0071076691 mol/g * 6.0221408e+23 atoms/mol = 4.2803384e+21 atoms/kg
"""

num_shells = 2000
r_0 = 6.4e6  # initial planet radius
mass_planet = 5.972e24  # mass earth
# gamma_a = 1.4  # polytropic exponent
# gamma = 1.4  # specific heat ratio of the gas, gamma = c_p / c_v
M_a = 0.029  # g/mol  # molecule mass, will likely have to be heavier for BSE atmosphere
T_0 = 288  # K, initial temperature (surface temperature)
P_0 = 1.01e5  # Pa, initial pressure (sea level)
gamma = get_specific_heat_ratio(dulong_petit_limit(4.2803384e+21, T_0))  # specific heat ratio of the gas, gamma = c_p / c_v
gamma_a = gamma
print(gamma)
sys.exit()
u_s = None  # defaults to 0.5 * u_esc in Genda and Abe 2003
outfile_dir = "test_outputs"

s = solver.LagrangianSolver1D(
    num_shells=num_shells,
    P_0=P_0,
    T_0=T_0,
    gamma=gamma,
    gamma_a=gamma_a,
    m_a=M_a,
    r_0=r_0,
    mass_planet=mass_planet,
    u_s=u_s,
    outfile_dir=outfile_dir
)
print(s.system.rho_0, s.system.c_s_0, s.system.vesc)
s.solve(timesteps=int(5e6))
