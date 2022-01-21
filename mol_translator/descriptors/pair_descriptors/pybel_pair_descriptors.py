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

def get_pyb_pair_descriptors(pybmol):
    pair_properties = {}
    pair_properties['bond_order_matrix'] = get_bond_order_matrix_pyb(pybmol)
    #pair_properties['distance_matrix'] = get_distance_matrix_pyb(pybmol)
    #pair_properties['angle_matrix'] = get_angle_matrix_pyb(pybmol)

    return pair_properties

def get_bond_order_matrix_pyb(pybmol):
    bond_matrix = np.zeros((len(pybmol.atoms), len(pybmol.atoms)), dtype=np.int32)
    for a in range(len(pybmol.atoms)):
        for b in range(len(pybmol.atoms)):
            if a == b: continue
            bond = pybmol.atoms[a].OBAtom.GetBond(pybmol.atoms[b].OBAtom)
            if bond is not None:
                bo = float(bond.GetBondOrder())
                bond_matrix[a][b] = bo
                bond_matrix[b][a] = bo

    return bond_matrix

def get_distance_matrix_pyb(pybmol):
    distance_matrix = np.zeros((len(pybmol.atoms), len(pybmol.atoms)), dtype=np.int32)
    for a in range(len(pybmol.atoms)):
        for b in range(len(pybmol.atoms)):
            if a == b: continue
            dist = pybmol.atoms[a].OBAtom.GetDistance(pybmol.atoms[b].OBAtom)
            distance_matrix[a][b] = int(dist)
            distance_matrix[b][a] = int(dist)

    return distance_matrix

def get_angle_matrix_pyb(pybmol):
    angle_matrix = np.zeros((len(pybmol.atoms),len(pybmol.atoms), len(pybmol.atoms)), dtype=np.int32)
    for a in range(len(pybmol.atoms)):
        for b in range(len(pybmol.atoms)):
            if a == b: continue
            for c in range(len(pybmol.atoms)):
                if a == c: continue
                if b == c: continue

                ang = pybmol.atoms[a].OBAtom.GetAngle(pybmol.atoms[b].OBAtom, pybmol.atoms[c].OBAtom)
                angle_matrix[a][b][c] = float(ang)
                angle_matrix[c][b][a] = float(ang)

    return angle_matrix