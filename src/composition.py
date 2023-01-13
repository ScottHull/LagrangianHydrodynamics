import re
from scipy.interpolate import interp1d
import pandas as pd


def get_molecular_mass(molecule: str):
    """
    Returns the number of moles of the given composition in weight percent.
    """
    # read in the period table
    pt = pd.read_csv("src/data/periodic_table.csv", index_col='element')
    # get the stoichiometry of the molecule
    stoich = re.findall(r'([A-Z][a-z]*)(\d*)', molecule)
    # get the sum of the masses of each atom in the molecule
    moles = 0
    for atom, count in stoich:
        if len(count) == 0:
            count = 1
        moles += int(count) * pt.loc[atom, 'atomic_mass']
    return moles


def get_mean_molecular_mass(composition: dict):
    """
    Returns the mean molecular mass of the given composition in weight percent.
    """
    # remove the gas label
    cleaned_composition = {key.split("_")[0].replace("+", "").replace("-", ""): value
                           for key, value in composition.items()}
    # remove species with 0 weight percent
    cleaned_composition = {key: value for key, value in cleaned_composition.items() if value != 0}
    # get the sum of the masses in the composition dict
    total_mass = sum(cleaned_composition.values())
    # get the number of moles in the composition dict
    total_moles = sum([(1 / get_molecular_mass(molecule)) * cleaned_composition[molecule]
                       for molecule in cleaned_composition.keys()])
    return total_mass / total_moles

def get_heat_capacity_ideal_gas(molecule_type='diatomic', heat_capacity_form='cp', R=8.31446261815324):
    """
    The heat capacity of an ideal gas, either monatomic, diatomic, or polyatomic.
    Based on https://pressbooks.online.ucf.edu/phy2048tjb/chapter/3-5-heat-capacities-of-an-ideal-gas/.
    """
    cp = 0
    cv = 0
    if molecule_type == 'monatomic':
        cp = 5 / 2
        cv = 3 / 2
    elif molecule_type == 'diatomic':
        cp = 7 / 2
        cv = 5 / 2
    elif molecule_type == 'polyatomic':
        cp = 4
        cv = 3
    cp = cp * R
    cv = cv * R
    if heat_capacity_form == 'cp':
        return cp
    elif heat_capacity_form == 'cv':
        return cv

def get_P0_and_rho0_given_T0(T_0, phase_curve_path="src/data/dunite_phase_boundary.txt"):
    """
    Using the phase curve allows for a single-state variable to be used to describe the initial conditions.
    Because temperature is the more reliable state variable from the SPH giant impact, we use it to calculate.
    Using the Hugoniot would require 2 state variables, which is not ideal.
    Note that SPH temperatures are typically underestimates.
    Using the dunite phase curve, uses the given temperature to interpolate the corresponding pressure and density.
    """
    phase_curve = pd.read_csv(phase_curve_path, skiprows=1, sep='\s+')
    # get the temperature above and below the given temperature
    T_below = phase_curve[phase_curve['T'] < T_0]['T'].iloc[0]
    T_above = phase_curve[phase_curve['T'] > T_0]['T'].iloc[-1]
    # get the corresponding pressures and densities
    P_below = phase_curve[phase_curve['T'] < T_0]['PVAP'].iloc[0]
    P_above = phase_curve[phase_curve['T'] > T_0]['PVAP'].iloc[-1]
    rho_below = phase_curve[phase_curve['T'] < T_0]['RHOVAP'].iloc[0]
    rho_above = phase_curve[phase_curve['T'] > T_0]['RHOVAP'].iloc[-1]
    # interpolate the pressure and density
    P_0 = interp1d([T_below, T_above], [P_below, P_above])(T_0)
    rho_0 = interp1d([T_below, T_above], [rho_below, rho_above])(T_0)
    return P_0 * 10 ** 9, rho_0

def get_P0_rho0_T0_given_S0(S_0, phase_curve_path="src/data/dunite_phase_boundary.txt"):
    """
    Returns P0, rho0, and T0 given S0.
    Entropy is better conserved in SPH than temperature and so it makes sense to interpolate with S.
    """
    phase_curve = pd.read_csv(phase_curve_path, skiprows=1, sep='\s+')
    entropies = phase_curve['SLIQ'].tolist() + phase_curve['SVAP'].tolist()
    pressures = phase_curve['PLIQ'].tolist() + phase_curve['PVAP'].tolist()
    densities = phase_curve['RHOLIQ'].tolist() + phase_curve['RHOVAP'].tolist()
    temperatures = phase_curve['T'].tolist()
    # get the entropy above and below the given entropy
    S_below = entropies[entropies.index(min(entropies, key=lambda x: abs(x - S_0))) + 1]
    S_above = entropies[entropies.index(min(entropies, key=lambda x: abs(x - S_0))) - 1]
    # get the corresponding pressures, densities,and temperatures
    P_below = pressures[entropies.index(S_below)]
    P_above = pressures[entropies.index(S_above)]
    rho_below = densities[entropies.index(S_below)]
    rho_above = densities[entropies.index(S_above)]
    T_below = temperatures[entropies.index(S_below)]
    T_above = temperatures[entropies.index(S_above)]
    # interpolate the pressure and density
    P_0 = interp1d([S_below, S_above], [P_below, P_above])(S_0)
    rho_0 = interp1d([S_below, S_above], [rho_below, rho_above])(S_0)
    T_0 = interp1d([S_below, S_above], [T_below, T_above])(S_0)
    return P_0 * 10 ** 9, rho_0, T_0