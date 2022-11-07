import re
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
    # get the sum of the masses in the composition dict
    total_mass = sum(composition.values())
    # get the number of moles in the composition dict
    total_moles = sum([(1 / get_molecular_mass(molecule)) * composition[molecule] for molecule in composition.keys()])
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
