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

import numpy as np
from typing import Type


def get_ob_atom_descriptors(obmol: Type) -> dict:
    """
    Gathers openbable atom features and stores the data inside a dictionary.
    Current atom descriptors included are:
        - heavy_valence, number of non-hydrogen atoms connected
        - heterovalence, number of heteroatoms connected
        - hyb, hybridisation of the atom, 1 for sp, 2 for sp2, 3 for sp3, etc...
        - partial charge, the partial charge of the atom

    :param obmol: Type, openbabel class object
    :return atom_properties: dict, dictionary containing all the atom descriptors
    """
    atom_properties = {}

    #atom_num = np.zeros(len(obmol.atoms), dtype=np.int32)
    #atom_mass = np.zeros(len(obmol.atoms), dtype=np.float64)
    heavy_valence = np.zeros(len(obmol.atoms), dtype=np.int32)
    heterovalence = np.zeros(len(obmol.atoms), dtype=np.int32)
    hyb = np.zeros(len(obmol.atoms), dtype=np.int32)
    partial_charge = np.zeros(len(obmol.atoms), dtype=np.float64)
    #formal_charge = np.zeros(len(obmol.atoms), dtype=np.int32)

    for idx, atom in enumerate(obmol.atoms):
        #atom_num[idx] = atom.atomicnum
        #atom_mass[idx] = atom.atomicmass
        heavy_valence[idx] = atom.heavy_valence
        heterovalence[idx] = atom.heterovalence
        hyb[idx] = atom.hyb
        partial_charge[idx] = atom.partialcharge
        #format_charge[idx] = atom.formalcharge

    #atom_properties['atom_num'] = atom_num
    #atom_properties['atom_mass'] = atom_mass
    atom_properties['heavy_valence'] = heavy_valence
    atom_properties['heterovalence'] = heterovalence
    atom_properties['hyb'] = hyb
    atom_properties['partial_charge'] = partial_charge
    #atom_properties['formal_charge'] = formal_charge

    return atom_properties
