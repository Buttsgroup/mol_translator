# Copyright 2020 Will Gerrard
#This file is part of autoenrich.

#autoenrich is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoenrich is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.


import numpy as np
import os


def pybmol_to_aemol(pybmol):

    type_array = np.zeros(len(mol.atoms), dtype=np.int32)
    xyz_array = np.zeros((len(mol.atoms),3), dtype=np.float64)
    conn_array = np.zeros((len(mol.atoms), len(mol.atoms)), dtype=np.int32)

    for i in range(len(mol.atoms)):
        type_array[i] = int(mol.atoms[i].atomicnum)
        xyz_array[i][0] = float(mol.atoms[i].coords[0])
		xyz_array[i][1] = float(mol.atoms[i].coords[1])
		xyz_array[i][2] = float(mol.atoms[i].coords[2])


        for j in range(len(mol.atoms)):
            if i == j:
                continue

            bond = mol.atoms[i].OBAtom.GetBond(mol.atoms[j].OBAtom)
            if bond is not None:
                conn_array[i][j] = int(bond.GetBondOrder())
                conn_array[j][i] = int(bond.GetBondOrder())

    return type_array, xyz_array, conn_array

def aemol_to_pybmol(aemol):

    # Do this a cheat way for now, can probably do this properly
