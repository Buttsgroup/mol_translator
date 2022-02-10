# Copyright 2022 Will Gerrard, Calvin Yiu
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

def get_rd_pair_descriptors(rdmol):
    pair_properties = {}
    pair_properties['bond_order_matrix'] = get_bond_order_matrix_rdkit(rdmol)

    return pair_properties

def get_bond_order_matrix_rdkit(rdmol):
    bond_matrix = np.zeros((len(rdmol.GetAtoms()), len(rdmol.GetAtoms())), dtype=np.int32)
    for a in range(len(rdmol.GetAtoms())):
        for b in range(len(rdmol.GetAtoms())):
            if a == b: continue
            bo = float(rdmol.GetBondBetweenAtoms(a,b).GetBondTypeAsDouble())
            bond_matrix[a][b] = bo
            bond_matrix[b][a] = bo

    return bond_matrix