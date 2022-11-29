from src import solver
from src import initial_conditions as ic
from src import composition

from math import sqrt
import matplotlib.pyplot as plt

# bse_composition = {
#     "SiO2": 45.40,
#     'MgO': 36.76,
#     'Al2O3': 4.48,
#     'TiO2': 0.21,
#     'Fe2O3': 0.00000,
#     'FeO': 8.10,
#     'CaO': 3.65,
#     'Na2O': 0.349,
#     'K2O': 0.031,
#     'ZnO': 6.7e-3,
# }

# disk_composition = {'O': 7.552391599681756, 'SiO2_g': 2.8254967389689796, 'Si_g': 0.0, 'Mg_g': 3.5162214342184765,
#                     'FeO_g': 4.701091252726253, 'CaO_g': 0.0, 'AlO_g': 0.01164730458854052, 'Al2O_g': 0.0,
#                     'AlO2_g': 0.0, 'Al2O2_g': 0.0, 'TiO_g': 0.0, 'TiO2_g': 0.0, 'NaO_g': 0.0, 'Na2_g': 0.0,
#                     'Na2O_g': 0.0, 'Na+': 0.0, 'KO_g': 0.0, 'K2O_g': 0.0, 'K2_g': 0.0, 'K+': 0.0, 'ZnO_g': 0.0,
#                     'e-': 0.0, 'SiO': 47.34150833750081, 'O2': 21.556803352884636, 'MgO': 1.2669515094056627,
#                     'Fe': 11.227888470024903, 'Ca': 0.0, 'Al': 0.0, 'Ti': 0.0, 'Na': 0.0, 'K': 0.0, 'ZnO': 0.0,
#                     'Zn': 0.0}

num_shells = 1000
r_0 = 6.4e6  # initial planet radius
mass_planet = 5.972e24  # mass earth
gamma = composition.get_heat_capacity_ideal_gas(molecule_type='diatomic',
                                                heat_capacity_form='cp') / composition.get_heat_capacity_ideal_gas(
    molecule_type='diatomic', heat_capacity_form='cv')  # specific heat
gamma_a = gamma  # polytropic exponent
# M_a = composition.get_mean_molecular_mass(disk_composition) / 1000  # g/mol to kg/mol
M_a = 0.036
T_0 = 3063.18893  # K, initial temperature (surface temperature)
u_s = 7086.216941  # defaults to 0.5 * u_esc in Genda and Abe 2003
P_0, rho_0 = composition.get_P0_and_rho0_given_T0(T_0)
jet_angle = 45.0
# outfile_dir = "jet_test_outputs"
# output_plots_dir = "jet_plots"
outfile_dir = "/scratch/shull4/jet-500b073S_outputs"
output_plots_dir = "/scratch/shull4/jet-500b073S_plots"
max_time = 86400  # 1 day in seconds

s = solver.LagrangianSolverJet(
    num_shells=num_shells,
    P_0=P_0,
    rho_0=rho_0,
    T_0=T_0,
    gamma=gamma,
    gamma_a=gamma_a,
    m_a=M_a,
    r_0=r_0,
    mass_planet=mass_planet,
    u_s=u_s,
    jet_angle=jet_angle,
    outfile_dir=outfile_dir,
    plot_separation=10000,
    use_cfl=True,
    show_figs=False,
    save_figs=False,
    output_file_interval=100,
    fig_save_path=output_plots_dir,
    run_name=outfile_dir.split("/")[-1],
)
s.solve(max_time=max_time)
