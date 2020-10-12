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
import pybel as pyb
from rdkit import Chem
from rdkit.Chem import AllChem

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

def aemol_to_pybmol(structure):
    # Do this a cheat way for now, can probably do this properly

    strucwrt.write_mol_toxyz(structure, 'tmp.xyz')
    pybmol = next(pyb.readfile('xyz', 'tmp.xyz'))
    os.remove('tmp.xyz')

    return pybmol

def rdmol_to_aemol(rdmol):

	type_array = np.zeros(rdmol.GetNumAtoms(), dtype=np.int32)
	xyz_array = np.zeros((rdmol.GetNumAtoms(),3), dtype=np.float64)
	conn_array = np.zeros((rdmol.GetNumAtoms(), rdmol.GetNumAtoms()), dtype=np.int32)

	for i, atoms in enumerate(rdmol.GetAtoms()):
	    type_array[i] = atoms.GetAtomicNum()
		xyz_array[i][0] = rdmol.GetConformer().GetAtomPosition(i).x
		xyz_array[i][1] = rdmol.GetConformer().GetAtomPosition(i).y
		xyz_array[i][2] = rdmol.GetConformer().GetAtomPosition(i).z


		for j, batoms in enumerate(rdmol.GetAtoms()):
			if i == j:
				continue

			bond = atoms.GetBonds(batoms)
			if bond is not None:
                conn_array[i][j] = int(bond.GetBondTypeAsDouble())
                conn_array[j][i] = int(bond.GetBondTypeAsDouble())

	return type_array, xyz_array, conn_array

def aemol_to_rdmol(structure):

	strucwrt.write_mol_tosdf(structure, 'tmp.sdf')
	molblock = open('tmp.sdf', 'r').read()
	rdmol = Chem.MolFromMolBlock(molblock)
	os.remove('tmp.sdf')

	return rdmol
