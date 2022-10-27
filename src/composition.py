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
