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

from util.targetflag import flag_to_target
from util.periodic_table import Get_periodic_table

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
    atom_props = []
    for propname in aemol.atom_properties.keys():
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
            conns.append(aemol.conn[t])
            for p in enumerate(aemol.atom_properties.keys()):
                atom_props[p].append(aemol[prop][t])

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
        pickle.dump(open('atoms.pkl', 'wb'))
    else:
	    return atoms


def make_bonds_df(aemols, progress=False, big_data=False, max_bond_distance=4):

		p_table = Get_periodic_table()

		# construct dataframe for bonds in molecule
		id = []				# number
		molecule_name = [] 	# molecule name
		atom_index_0 = []	# atom index for atom 1
		atom_index_1 = []	# atom index for atom 2
		dist = []			# distance between atoms
		bnd_dist = []		# number of bonds between atoms (shortest path)

		if progress:
			pbar = tqdm(aemols, desc='Constructing bonds dictionary', leave=False)
		else:
			pbar = aemols

		m = -1
		for molrf in pbar:
			m += 1

			# If only refs supplied, import molecule from file
			if big_data:
				mol = nmrmol(molid=molrf[1])

				if molrf[2] == '':
					ftype = get_type(molrf[0])
				else:
					ftype = molrf[2]
				aemol.read_nmr(molrf[0], ftype)
			else:
				mol = molrf

			for t, type in enumerate(aemol.types):
				for t2, type2 in enumerate(aemol.types):
					if t == t2:
						continue

					# Only include pairs up to 4 bonds away
					if aemol.coupling_len[t][t2] > max_bond_distance or aemol.coupling_len[t2][t] > max_bond_distance or aemol.coupling_len[t][t2] == 0:
						continue

					# Enforce label ordering, always 1JCH, never 1JHC for example
					if type > type2:
						targetflag = str(aemol.coupling_len[t][t2]) + 'J' + p_table[type] + p_table[type2]
					else:
						targetflag = str(aemol.coupling_len[t][t2]) + 'J' + p_table[type2] + p_table[type]

					# Add pair values to lists
					molecule_name.append(aemol.molid)
					atom_index_0.append(t)
					atom_index_1.append(t2)
					dist.append(np.linalg.norm(aemol.xyz[t]-aemol.xyz[t2]))
					bnd_dist.append(aemol.coupling_len[t][t2])
					cpltype.append(targetflag)
					coupling.append(aemol.coupling[t][t2])

		# Construct dataframe
		bonds = {	'molecule_name': molecule_name,
					'atom_index_0': atom_index_0,
					'atom_index_1': atom_index_1,
					'type': cpltype,
					'dist': dist,
					'bond_dist': bnd_dist,
					'scalar_coupling_constant': coupling
				}

		bonds = pd.DataFrame(bonds)
		bonds['predicted_scalar_coupling_constant'] = 0.0
		bonds['predict'] = 0

		if progress:
			pbar.close()

		return bonds
