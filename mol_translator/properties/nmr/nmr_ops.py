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

from mol_translator.util.periodic_table import Get_periodic_table

def get_coupling_types(aemol, maxlen=6):
	p_table = Get_periodic_table()

	if 'path_len' not in aemol.structure.keys():
		aemol.get_path_lengths()

	cpl_types = []
	for t, type in enumerate(aemol.structure['types']):
		tmp_types = []
		for t2, type2 in enumerate(aemol.structure['types']):

			if type > type2:
				targetflag = str(aemol.structure['path_len'][t][t2]) + 'J' + p_table[type] + p_table[type2]
			else:
				targetflag = str(aemol.structure['path_len'][t][t2]) + 'J' + p_table[type2] + p_table[type]

			tmp_types.append(targetflag)
		cpl_types.append(tmp_types)


	aemol.pair_properties['nmr_types'] = cpl_types

def scale_chemical_shifts(aemol, scaling={6: [-1.0399, 187.136], 1: [-1.0719, 32.1254]}):
	# default scaling values are for functional: wb97xd, basis_set: 6-311g(d,p)

	for t, type in enumerate(aemol.structure['types']):
		if type not in scaling.keys():
			continue
			
		aemol.atom_properties['shift'][t] = (aemol.atom_properties['shift'][t] - scaling[type][1]) / scaling[type][0] 
		
		









#
