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

from atom_descriptors import rdkit_atom_descriptors, openbabel_atom_descriptors
from pair_descriptors import rdkit_pair_descriptors, openbabel_pair_descriptors


def get_all_descriptors(aemol: Type) -> None:
    """
    Gets available descriptors from both rdkit and openbabel and stores descriptors within the aemol object

    :param aemol: Type, aemol object to obtain descriptors for
    :return None:
    """
    get_openbabel_descriptors(aemol)
    get_rdkit_descriptors(aemol)


def get_openbabel_descriptors(aemol: Type, get_atom_desc: bool = True, get_pair_desc: bool = True):
    """
    Gets available descriptors from openbabel and stores descriptors within the aemol object

    :param aemol: Type, aemol object to obtain descriptors for
    :param get_atom_desc: bool, set to True to obtain openbabel atom descriptors
    :param get_pair_desc: bool, set to True to obtain openbabel pair descriptors
    :return None:
    """
    if aemol.obmol is None:
        obmol = aemol.to_pybel()
    else:
        obmol = aemol.obmol

    if get_atom_desc:
        atom_props = openbabel_atom_descriptors.get_ob_atom_descriptors(obmol)
        aemol.atom_properties.update(atom_props)

    if get_pair_desc:
        pair_props = openbabel_pair_descriptors.get_ob_pair_descriptors(obmol)
        aemol.pair_properties.update(pair_props)


def get_rdkit_descriptors(aemol: Type, get_atom_desc: bool = True, get_pair_desc: bool = False):
    """
    Gets available descriptors from rdkit and stores descriptors within the aemol object

    :param aemol: Type, aemol object to obtain descriptors for
    :param get_atom_desc: bool, set to True to obtain rdkit atom descriptors
    :param get_pair_desc: bool, set to True to obtain rdkit pair descriptors
    :return None:
    """
    if aemol.rdmol is None:
        rdmol = aemol.to_rdkit()
    else:
        rdmol = aemol.rdmol

    if get_atom_desc:
        atom_props = rdkit_atom_descriptors.get_rd_atom_descriptors(rdmol)
        aemol.atom_properties.update(atom_props)

    if get_pair_desc:
        pair_props = rdkit_pair_descriptors.get_rd_pair_descriptors(rdmol)
        aemol.pair_properties.update(pair_props)
