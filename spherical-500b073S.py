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
mass_planet = 5.972e24  # mass earth
gamma = composition.get_heat_capacity_ideal_gas(molecule_type='diatomic', heat_capacity_form='cp') / composition.get_heat_capacity_ideal_gas(molecule_type='diatomic', heat_capacity_form='cv')  # specific heat, gamma = cp/cv
gamma_a = gamma  # polytropic exponent
M_a = composition.get_mean_molecular_mass(bse_composition) / 1000  # g/mol  # molecule mass, will likely have to be heavier for BSE atmosphere
T_0 = 3063.18893  # K, initial temperature (surface temperature)
P_0 = 34703699057  # Pa, initial pressure (sea level)
u_s = 7086.216941  # defaults to 0.5 * u_esc in Genda and Abe 2003
rho_0 = 3656.224576
jet_angle = 45.0
# outfile_dir = "jet_test_outputs"
# output_plots_dir = "jet_plots"
outfile_dir = "/scratch/shull4/spherical-500b073S_outputs"
output_plots_dir = "/scratch/shull4/spherical-500b073S_plots"
max_time = 86400  # 1 day in seconds

s = solver.LagrangianSolver1DSpherical(
    num_shells=num_shells,
    P_0=P_0,
    T_0=T_0,
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
    output_file_interval=10000,
    fig_save_path=output_plots_dir,
    run_name=outfile_dir.split("/")[-1],
)
s.solve(max_time=max_time)
