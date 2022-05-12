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
from typing import Type
from mol_translator.properties.nmr.nmr_ops import get_coupling_types
import os.path
import numpy as np

# Gets the necessary info for aemol to turn it into an IMP dataframe


def prep_mol(aemol: Type):
    """
    A function to prepare for the conversion of an Aemol molecule to a dataframe format. 
    Obtains all the bond/coupling information and finds the path lengths (absolute distance and bond distance)

    :param aemol: Type[Aemol], the aemol object to convert to dataframe
    :return: Aemol
    """

    aemol.get_bonds()
    aemol.get_path_lengths()
    get_coupling_types(aemol)

    return aemol


def prep_mol_nmr(aemol: Type, nmr_file: str = "", nmr_type: str = ""):
    """
    A function to prepare for the conversion of an Aemol molecule to a dataframe format with the inclusion of NMR data. 
    Obtains all the bond/coupling information and finds the path lengths (absolute distance and bond distance)
    Processes NMR file given to extract tensors/chemical shifts and coupling constants, if None is passed then fills with 0 values

    :param aemol: Type[Aemol], the aemol object to convert to dataframe
    :param nmr_file: str, file path of the correspoinding NMR datafile
    :param nmr_type: str, file type of the NMR data file, e.g. 'log', 'nmredata'
    :return: Aemol
    """

    aemol.get_bonds()
    aemol.get_path_lengths()
    get_coupling_types(aemol)
    if os.path.isfile(nmr_file):
        aemol.prop_from_file(nmr_file, 'nmr', nmr_type)
    else:
        aemol.atom_properties['shift'] = np.zeros(
            len(aemol.structure['types']), dtype=np.float64)
        aemol.pair_properties['coupling'] = np.zeros(
            (len(aemol.structure['types']), len(aemol.structure['types'])), dtype=np.float64)
    return aemol


def prep_mol_ic50(aemol: Type, ic50_file: str = "", ic50_type: str = ""):
    """
    A function to prepare for the conversion of an Aemol molecule to a dataframe format with the inclusion of ic50 data. 
    Obtains all the bond/coupling information and finds the path lengths (absolute distance and bond distance)
    Processes ic50 file given to extract ic50 values, if None is passed then fills with 0 values

    :param aemol: Type[Aemol], the aemol object to convert to dataframe
    :param ic50_file: str, file path of the correspoinding ic50 datafile
    :param ic50_type: str, file type of the ic50 data file, e.g. 'log', 'nmredata'
    :return: Aemol
    """
    aemol.get_bonds()
    aemol.get_path_lengths()
    get_coupling_types(aemol)
    if os.path.isfile(ic50_file):
        aemol.prop_from_file(ic50_file, 'ic50', ic50_type)
        aemol.atom_properties['ic50'] = np.full(
            len(aemol.structure['types']), aemol.mol_properties['ic50'], dtype=np.float64)
    else:
        print('Setting fake ic50 values', ic50_file)
        aemol.atom_properties['ic50'] = np.zeros(
            len(aemol.structure['types']), dtype=np.float64)

    return aemol


def prep_mol_ecfp4(aemol: Type):
    """
    A function to prepare for the conversion of an Aemol molecule to a dataframe format with fingerprints. 
    Obtains all the bond/coupling information and finds the path lengths (absolute distance and bond distance)
    Processes the fingerprint stored

    :param aemol: Type[Aemol], the aemol object to convert to dataframe
    :return: Aemol
    """
    aemol.get_bonds()
    aemol.get_path_lengths()
    get_coupling_types(aemol)
    aemol.atom_properties['ecfp4'] = aemol.mol_properties['ecfp4']

    return aemol


def prep_mol_mulliken(aemol: Type, mc_file: str = "", mc_format: str = ""):
    """
    A function to prepare for the conversion of an Aemol molecule to a dataframe format with the inclusion of mulliken charge data. 
    Obtains all the bond/coupling information and finds the path lengths (absolute distance and bond distance)
    Processes data file given to extract mulliken charge values, if None is passed then fills with 0 values

    :param aemol: Type[Aemol], the aemol object to convert to dataframe
    :param mc_file: str, file path of the correspoinding mulliken charge datafile
    :param mc_type: str, file type of the mulliken charge data file, e.g. 'log', 'nmredata'
    :return: Aemol
    """
    aemol.get_bonds()
    aemol.get_path_lengths()
    get_coupling_types(aemol)
    if os.path.isfile(mc_file):
        aemol.prop_from_file(mc_file, mc_format, 'mc')
    else:
        print('Setting fake mulliken charge values', mc_file)
        aemol.atom_properties['mull_chg'] = np.zeros(
            len(aemol.structure['types']), dtype=np.float64)

    return aemol
