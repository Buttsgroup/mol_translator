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

from mol_translator.util.targetflag import flag_to_target
from mol_translator.util.periodic_table import Get_periodic_table

import numpy as np
import pandas as pd
from tqdm import tqdm

import sys
import collections

import pickle

def make_atom_df(aemols, progress=False, write=False):

	p_table = Get_periodic_table()

	# construct dataframes
	# atoms has: molecule_name, atom, labeled atom,
	molecule_name = [] 	# molecule name
	atom_index = []		# atom index
	typestr = []		# atom type (string)
	typeint = []		# atom type (integer)
	x = []				# x coordinate
	y = []				# y coordinate
	z = []				# z coordinate
	conns = []
	atom_props = []
	for propname in aemols[0].atom_properties.keys():
		atom_props.append([])

	if progress:
		pbar = tqdm(aemols, desc='Constructing atom dictionary', leave=False)
	else:
		pbar = aemols

	m = -1
	for aemol in pbar:
		m += 1
		# Add atom values to lists
		for t, type in enumerate(aemol.structure['types']):
			molecule_name.append(aemol.info['molid'])
			atom_index.append(t)
			typestr.append(p_table[type])
			typeint.append(type)
			x.append(aemol.structure['xyz'][t][0])
			y.append(aemol.structure['xyz'][t][1])
			z.append(aemol.structure['xyz'][t][2])
			conns.append(aemol.structure['conn'][t])
			for p, prop in enumerate(aemol.atom_properties.keys()):
				atom_props[p].append(aemol.atom_properties[prop][t])

	# Construct dataframe
	atoms = {	'molecule_name': molecule_name,
				'atom_index': atom_index,
				'typestr': typestr,
				'typeint': typeint,
				'x': x,
				'y': y,
				'z': z,
				'conn': conns
			}
	for p, propname in enumerate(aemol.atom_properties.keys()):
		atoms[propname] = atom_props[p]

	atoms = pd.DataFrame(atoms)

	if progress:
		pbar.close()

	if write:
		atoms.to_pickle('atoms.pkl')
	else:
		return atoms


def make_pair_df(aemols, progress=False, max_bond_distance=4, write=False):

	p_table = Get_periodic_table()

	# construct dataframe for pairs in molecule
	id = []				# number
	molecule_name = [] 	# molecule name
	atom_index_0 = []	# atom index for atom 1
	atom_index_1 = []	# atom index for atom 2
	dist = []			# distance between atoms
	path_len = []		# number of pairs between atoms (shortest path)
	pair_props = []
	for propname in aemols[0].pair_properties.keys():
		pair_props.append([])

	if progress:
		pbar = tqdm(aemols, desc='Constructing pairs dictionary', leave=False)
	else:
		pbar = aemols

	m = -1
	for aemol in pbar:
		m += 1

		for t, type in enumerate(aemol.structure['types']):
			for t2, type2 in enumerate(aemol.structure['types']):
				# Add pair values to lists
				molecule_name.append(aemol.info['molid'])
				atom_index_0.append(t)
				atom_index_1.append(t2)
				dist.append(np.linalg.norm(aemol.structure['xyz'][t]-aemol.structure['xyz'][t2]))
				path_len.append(aemol.structure['path_len'][t][t2])
				for p, prop in enumerate(aemol.pair_properties.keys()):
					pair_props[p].append(aemol.pair_properties[prop][t][t2])

	# Construct dataframe
	pairs = {	'molecule_name': molecule_name,
				'atom_index_0': atom_index_0,
				'atom_index_1': atom_index_1,
				'dist': dist,
				'path_len': path_len
			}
	for p, propname in enumerate(aemol.pair_properties.keys()):
		pairs[propname] = pair_props[p]

	pairs = pd.DataFrame(pairs)

	if progress:
		pbar.close()

	if write:
		pairs.to_pickle('pairs.pkl')
	else:
		return pairs
