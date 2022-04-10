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
import os
from rdkit import Chem
from rdkit.Chem import AllChem

from mol_translator.structure import structure_write as strucwrt


def rdmol_to_aemol(rdmol):

	type_array = np.zeros(rdmol.GetNumAtoms(), dtype=np.int32)
	xyz_array = np.zeros((rdmol.GetNumAtoms(),3), dtype=np.float64)
	conn_array = np.zeros((rdmol.GetNumAtoms(), rdmol.GetNumAtoms()), dtype=np.int32)

	for i, atoms in enumerate(rdmol.GetAtoms()):
		type_array[i] = atoms.GetAtomicNum()
		xyz_array[i][0] = rdmol.GetConformer().GetAtomPosition(i).x
		xyz_array[i][1] = rdmol.GetConformer().GetAtomPosition(i).y
		xyz_array[i][2] = rdmol.GetConformer().GetAtomPosition(i).z

		for j, atoms in enumerate(rdmol.GetAtoms()):
			if i == j:
				continue

			bond = rdmol.GetBondBetweenAtoms(i,j)
			if bond is not None:
				conn_array[i][j] = int(bond.GetBondTypeAsDouble())
				conn_array[j][i] = int(bond.GetBondTypeAsDouble())

	return type_array, xyz_array, conn_array

def aemol_to_rdmol(aemol, molid, sanitize, removeHs):

	strucwrt.write_mol_tosdf(aemol, f'tmp{molid}.sdf')
	molblock = open(f'tmp{molid}.sdf', 'r').read()
	rdmol = Chem.MolFromMolBlock(molblock, sanitize=sanitize, removeHs=removeHs)
	os.remove(f'tmp{molid}.sdf')

	return rdmol
