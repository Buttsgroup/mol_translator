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
from openbabel import pybel as pyb
from openbabel import openbabel as ob

from mol_translator.structure import structure_write as strucwrt

def pybmol_to_aemol(pybmol):

	type_array = np.zeros(len(pybmol.atoms), dtype=np.int32)
	xyz_array = np.zeros((len(pybmol.atoms),3), dtype=np.float64)
	conn_array = np.zeros((len(pybmol.atoms), len(pybmol.atoms)), dtype=np.int32)

	for i in range(len(pybmol.atoms)):
		type_array[i] = int(pybmol.atoms[i].atomicnum)
		xyz_array[i][0] = float(pybmol.atoms[i].coords[0])
		xyz_array[i][1] = float(pybmol.atoms[i].coords[1])
		xyz_array[i][2] = float(pybmol.atoms[i].coords[2])


		for j in range(len(pybmol.atoms)):
			if i == j:
				continue

			bond = pybmol.atoms[i].OBAtom.GetBond(pybmol.atoms[j].OBAtom)
			if bond is not None:
				conn_array[i][j] = int(bond.GetBondOrder())
				conn_array[j][i] = int(bond.GetBondOrder())

	return type_array, xyz_array, conn_array

def aemol_to_pybmol(structure, id):
	# Do this a cheat way for now, can probably do this properly

	strucwrt.write_mol_toxyz(structure, f'tmp{id}.xyz')
	pybmol = next(pyb.readfile('xyz', f'tmp{id}.xyz'))
	os.remove(f'tmp{id}.xyz')

	return pybmol
