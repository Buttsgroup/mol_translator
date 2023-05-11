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
import os
import openbabel.pybel as pyb

from mol_translator.structure import structure_write as strucwrt


def obmol_to_aemol(obmol):
    type_array = np.zeros(len(obmol.atoms), dtype=np.int32)
    xyz_array = np.zeros((len(obmol.atoms), 3), dtype=np.float64)
    conn_array = np.zeros((len(obmol.atoms), len(obmol.atoms)), dtype=np.int32)

    for i, atom_x in enumerate(obmol.atoms):
        type_array[i] = int(atom_x.atomicnum)
        xyz_array[i][0] = float(atom_x.coords[0])
        xyz_array[i][1] = float(atom_x.coords[1])
        xyz_array[i][2] = float(atom_x.coords[2])

        for j, atom_y in enumerate(obmol.atoms):
            if i == j:
                continue

            bond = atom_x.OBAtom.GetBond(atom_y.OBAtom)
            if bond is not None:
                conn_array[i][j] = int(bond.GetBondOrder())
                conn_array[j][i] = int(bond.GetBondOrder())

    return type_array, xyz_array, conn_array


def old_aemol_to_obmol(structure, id):
    # Do this a cheat way for now, can probably do this properly

    strucwrt.write_mol_toxyz(structure, f"{id}_tmp.xyz")
    obmol = next(pyb.readfile("xyz", f"{id}_tmp.xyz"))
    os.remove(f"{id}_tmp.xyz")

    return obmol


def aemol_to_obmol(structure):
    # Testing mol generation from scratch using openbabel
    # Generating blank pybel mol object
    obmol = pyb.ob.OBMol()

    # adding atoms based off aemol.structure['types'], adding coordinates based of aemol.structure['xyz']
    for i, atom in enumerate(structure["types"]):
        obatom = obmol.NewAtom()
        obatom.SetAtomicNum(int(atom))
        obatom.SetVector(*structure["xyz"][i])

    # adding bonds based off aemol.structure['conn']
    visited = []
    for i, bond_order_array in enumerate(structure["conn"]):
        for j, bond_order in enumerate(bond_order_array):
            if j in visited:
                continue
            elif bond_order != 0:
                obmol.AddBond(i + 1, j + 1, int(bond_order))
            else:
                continue
        visited.append(i)

    return pyb.Molecule(obmol)
