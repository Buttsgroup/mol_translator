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


def get_ob_pair_descriptos(obmol: Type) -> dict:
    """
    Gathers openbable atom-pair/bond features and stores the data inside a dictionary.
    Current openbabel pair descriptors included are:
        - bond order, the bond order associated with an atom-pair connection

    :param obmol: Type, openbabel class object
    :return atom_properties: dict, dictionary containing all the atom descriptors
    """

    pair_properties = {}
    pair_properties['bond_order'] = get_bond_order_matrix_ob(obmol)
    #pair_properties['distance_matrix'] = get_distance_matrix_ob(obmol)
    #pair_properties['angle_matrix'] = get_angle_matrix_ob(obmol)

    return pair_properties


def get_bond_order_matrix_ob(obmol: Type) -> np.darray:
    """
    Generates an NxN numpy array (N = number of atoms) containing the bond order of connected bonds,
    0 is no bonds, 1 is a single bond, 2 is a double bond, 3 is a triple bonds. 

    :param obmol: Type, openbabel class object
    :return bond_matrix: np.array, NxN matrix of bond order values
    """
    bond_matrix = np.zeros(
        (len(obmol.atoms), len(obmol.atoms)), dtype=np.int32)
    for a in range(len(obmol.atoms)):
        for b in range(len(obmol.atoms)):
            if a == b:
                continue
            bond = obmol.atoms[a].OBAtom.GetBond(obmol.atoms[b].OBAtom)
            if bond is not None:
                bo = float(bond.GetBondOrder())
                bond_matrix[a][b] = bo
                bond_matrix[b][a] = bo

    return bond_matrix


def get_distance_matrix_ob(obmol: Type) -> np.darray:
    """
    Generates an NxN numpy array (N = number of atoms) containing the distance between atoms in angstroms 

    :param obmol: Type, openbabel class object
    :return distance_matrix: np.array, NxN matrix of distance values
    """
    distance_matrix = np.zeros(
        (len(obmol.atoms), len(obmol.atoms)), dtype=np.int32)
    for a in range(len(obmol.atoms)):
        for b in range(len(obmol.atoms)):
            if a == b:
                continue
            dist = obmol.atoms[a].OBAtom.GetDistance(obmol.atoms[b].OBAtom)
            distance_matrix[a][b] = dist
            distance_matrix[b][a] = dist

    return distance_matrix


def get_angle_matrix_ob(obmol: Type) -> np.darray:
    """
    Generates an NxNxN numpy array (N = number of atoms) containing the bond angle between atom triplets in degrees 

    :param obmol: Type, openbabel class object
    :return angle_matrix: np.array, NxNxN matrix of distance values
    """
    angle_matrix = np.zeros((len(obmol.atoms), len(
        obmol.atoms), len(obmol.atoms)), dtype=np.int32)
    for a in range(len(obmol.atoms)):
        for b in range(len(obmol.atoms)):
            if a == b:
                continue
            for c in range(len(obmol.atoms)):
                if a == b:
                    continue
                if b == c:
                    continue

                ang = obmol.atoms[a].OBAtom.GetAngle(
                    obmol.atoms[b].OBAtom, obmol.atoms[c].OBAtom)
                angle_matrix[a][b][c] = float(ang)
                angle_matrix[c][b][a] = float(ang)

    return angle_matrix
