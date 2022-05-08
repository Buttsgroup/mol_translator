# Copyright 2022 Will Gerrard, Calvin Yiu
# This file is part of autoenrich.

# autoenrich is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# autoenrich is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.

# Functions here are from utilises rdMolStandardize which originates from MolVS https://github.com/mcs07/MolVS

from typing import Union, Type, List

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def cleanup(rdmol: Type) -> Type:
    """
    Function for a general purpose cleanup of molecules based on predetermined parameters
    which can be checked through calling rdMolStandardize.CleanupParameters()

    :param rdmol: RDKit molecule object
    :return: RDKit molecule object
    """
    cleaned_rdmol = rdMolStandardize.Cleanup(rdmol)
    return cleaned_rdmol


def standardise_smiles(smiles_str: str) -> Type:
    """
    Validates SMILES string passed as argument and converts to RDKit molecule object

    :param smiles_str: str SMILES string
    :return: RDKit molecule object
    """
    return rdMolStandardize.StandardizeSmiles(smiles_str)


def normalise(rdmol: Union[str, Type], from_smiles: bool = False) -> Type:
    """
    Function for correcting functional groups and recombining charges, with the input
    being smiles string or a RDKit molecule object. Change from_smiles T/F flag depending
    on input manually, defaults to False

    :param mol: str SMILES string, or RDKit molecule object
    :param from_smiles: bool, set to true if param mol is a smiles string
    :return: RDKit molecule object
    """
    normaliser = rdMolStandardize.Normalizer()

    if from_smiles:
        return normaliser.normalize(Chem.MolFromSmiles(rdmol))

    return normaliser.normalize(rdmol)


def disconnect_metal(mol: Union[str, Type], from_smiles: bool = False) -> Type:
    """
    Function for fragmenting molecules contaning metals into charged species. Change from_smiles T/F flag depending
    on input manually, defaults to False

    Consider running rdmol_fragment_chooser after to isolate the desired structure, then rdmol_uncharger to neutralise charges

    :param mol: str SMILES string, or RDKit molecule object
    :param from_smiles: bool, set to true if param mol is a smiles string
    :return: RDKit molecule object
    """

    mol_disconnector = rdMolStandardize.MetalDisconnector()

    if from_smiles:
        return mol_disconnector.Disconnect(Chem.MolFromSmiles(mol))

    else:
        return mol_disconnector.Disconnect(mol)


def select_largest_fragment(mol: Union[str, Type], from_smiles: bool = False) -> Type:
    """
    Function for selecting the largest fragment, generally used in combination with uncharger to remove counter ions
    and neutralise charge

    Consider running rdmol_uncharger afterwards to neutralise charges

    :param mol: str SMILES string, or RDKit molecule object
    :param from_smiles: bool, set to true if param mol is a smiles string
    :return: RDKit molecule object
    """
    fragment_selector = rdMolStandardize.LargestFragmentChooser()

    if from_smiles:
        return fragment_selector.choose(Chem.MolFromSmiles(mol))

    else:
        return fragment_selector.choose(mol)


def uncharger(mol: Union[str, Type], from_smiles: bool = False) -> Type:
    """
    Function for neutralising point charges through addition of hydrogen

    :param mol: str SMILES string, or RDKit molecule object
    :param from_smiles: bool, set to true if param mol is a smiles string
    :return: RDKit molecule object
    """
    uncharger = rdMolStandardize.Uncharger()

    if from_smiles:
        return uncharger.uncharge(Chem.MolFromSmiles(mol))

    else:
        return uncharger.uncharge(mol)


def charge_parent(mol: Union[str, Type], from_smiles: bool = False) -> Type:
    """
    Function for neutralising the largest fragment of a charged species

    :param mol: str SMILES string, or RDKit molecule object
    :param from_smiles: bool, set to true if param mol is a smiles string
    :return: RDKit molecule object
    """
    if from_smiles:
        rdMolStandardize.ChargeParent(Chem.MolFromSmiles(mol))
        return mol

    else:
        rdMolStandardize.ChargeParent(mol)
        return mol


def find_resonance(mol: Union[str, Type], from_smiles: bool = False) -> List:
    """
    Function for finding viable resonance structures and returning them all in a list

    :param mol: str SMILES string, or RDKit molecule object
    :param from_smiles: bool, set to true if param mol is a smiles string
    :return: List, RDKit molecule object list
    """
    resonance = rdMolStandardize.ResonanceEnumerator()

    if from_smiles:
        return resonance.enumerate(Chem.MolFromSmiles(mol))

    else:
        return resonance.enumerate(mol)
