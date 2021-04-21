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
from openbabel import pybel as pyb
from openbabel import openbabel as ob
from mol_translator.util.periodic_table import Get_periodic_table

def pybmol_find_all_paths(pybmol, maxlen=5):
	all_paths = []
	for atom1 in range(len(pybmol.atoms)):
		for atom2 in range(len(pybmol.atoms)):
			coords1 = np.asarray(pybmol.atoms[atom1].coords)
			coords2 = np.asarray(pybmol.atoms[atom2].coords)
			if np.linalg.norm(coords1-coords2) > 2.5*maxlen:
				continue
			for length in range(1, maxlen):
				paths = pybmol_find_paths(pybmol, atom1, atom2, length)
				if len(paths) > 0:
					all_paths.extend(paths)

	return all_paths


def pybmol_find_paths(pybmol, start, end, coupling_length, path=[]):
		# append atom to start
		path = path + [start]
		# check if we have reached target atom
		if start == end:
			# if we have, return succesful path
			if len(path) == 1:
				return []
			else:
				return [path]
		# define new path
		paths = []
		# loop over neighbouring atoms
		for nbr_atom in ob.OBAtomAtomIter(pybmol.atoms[start].OBAtom):
			# get ID of neighbour
			node = nbr_atom.GetId()
			# check the neighbour is not already in the path, and that the path is not over the required length
			if node not in path and len(path) <= coupling_length:
				# get new paths for the neighbour
				newpaths = pybmol_find_paths(pybmol, node, end, coupling_length, path)
				#for each new path, check for paths of correct length
				for newpath in newpaths:
					if len(newpath) == coupling_length+1 and newpath != []:
						paths.append(newpath)

		return paths

def pybmol_get_bond_table(pybmol):

	atoms = len(pybmol.atoms)

	bond_table = np.zeros((atoms, atoms), dtype=np.int32)

	for atom1 in range(atoms):
		for atom2 in range(atom1, atoms):

			for nbr_atom in ob.OBAtomAtomIter(pybmol.atoms[atom1].OBAtom):
				check = nbr_atom.GetId()
				if atom2 != check:
					continue

				bond = pybmol.atoms[atom1].OBAtom.GetBond(nbr_atom)
				order = bond.GetBondOrder()

				bond_table[atom1][atom2] = int(order)
				bond_table[atom2][atom1] = int(order)

	return bond_table

def pybmol_get_path_lengths(pybmol, maxlen=5):
	atoms = len(pybmol.atoms)
	coupling_len = np.zeros((atoms, atoms), dtype=np.int32)

	allpaths = pybmol_find_all_paths(pybmol, maxlen=maxlen)

	for atom1 in range(atoms):
		for atom2 in range(atoms):
			if atom1 == atom2:
				shortest = 0
			else:
				shortest = maxlen+1
			for path in allpaths:
				if path[0] == atom1 and path[-1] == atom2:
					if len(path)-1 < shortest:
						shortest = len(path) - 1

			coupling_len[atom1][atom2] = shortest
			coupling_len[atom2][atom1] = shortest

	return coupling_len
