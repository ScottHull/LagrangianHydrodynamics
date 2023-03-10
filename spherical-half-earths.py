from src import solver
from src import initial_conditions as ic
from src import composition

from math import sqrt
import matplotlib.pyplot as plt

bse_composition = {
    "SiO2": 45.40,
    'MgO': 36.76,
    'Al2O3': 4.48,
    'TiO2': 0.21,
    'Fe2O3': 0.00000,
    'FeO': 8.10,
    'CaO': 3.65,
    'Na2O': 0.349,
    'K2O': 0.031,
    'ZnO': 6.7e-3,
}

num_shells = 1000
r_0 = 6.4e6  # initial planet radius
mass_planet = 1.1 * 5.972e24  # mass earth
gamma = composition.get_heat_capacity_ideal_gas(molecule_type='diatomic',
                                                heat_capacity_form='cp') / composition.get_heat_capacity_ideal_gas(
    molecule_type='diatomic', heat_capacity_form='cv')  # specific heat
gamma_a = gamma  # polytropic exponent
# M_a = composition.get_mean_molecular_mass(disk_composition) / 1000  # g/mol to kg/mol
M_a = 0.036
u_s = 5199.265282773104  # defaults to 0.5 * u_esc in Genda and Abe 2003
T_0 = 3517.8334184120226  # K
# P_0, rho_0, T_0 = composition.get_P0_rho0_T0_given_S0(S_0)  # initial pressure, density, and temperature
P_0, rho_0 = composition.get_P0_and_rho0_given_T0(T_0)
# outfile_dir = "spherical_test_outputs"
# output_plots_dir = "spherical_plots"
outfile_dir = "/scratch/shull4/spherical-500b073S_outputs"
output_plots_dir = "/scratch/shull4/spherical-500b073S_plots"
max_time = 86400  # 1 day in seconds

s = solver.LagrangianSolver1DSpherical(
    num_shells=num_shells,
    P_0=P_0,
    T_0=T_0,
    rho_0=rho_0,
    gamma=gamma,
    gamma_a=gamma_a,
    m_a=M_a,
    r_0=r_0,
    mass_planet=mass_planet,
    u_s=u_s,
    outfile_dir=outfile_dir,
    plot_separation=100000,
    use_cfl=True,
    show_figs=False,
    save_figs=False,
    output_file_interval=1000,
    fig_save_path=output_plots_dir,
    run_name=outfile_dir.split("/")[-1],
)
s.solve(max_time=max_time)
