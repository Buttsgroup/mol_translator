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

from mol_translator.structure.structure_write import write_mol_tosdf

# Write an nmrmol object to an nmredata file
def write_nmredata(outfile, aemol, write_zeros=False,
					info=None, count_from=0, print_predicted=False):

	atoms = len(aemol.structure['types'])
	props = {}
	for label in ["predicted_shift", "shift", "shift_var"]:
		if label in aemol.atom_properties.keys():
			props[label] = aemol.atom_properties[label]
		else:
			props[label] = np.zeros(atoms, dtype=np.float64)
	for label in ["predicted_coupling", "coupling", "coupling_var"]:
		if label in aemol.pair_properties.keys():
			props[label] = aemol.pair_properties[label]
		else:
			props[label] = np.zeros((atoms,atoms), dtype=np.float64)

	if print_predicted:
		props['shift'] = props['predicted_shift']
		props['coupling'] = props['predicted_coupling']

	sdfstrings = write_mol_tosdf(aemol, "", stringsonly=True)

	lines = []
	for line in sdfstrings:
		lines.append(line)
	# assignment section
	lines.append('')
	lines.append('> <NMREDATA_ASSIGNMENT>')
	# Print chemical shifts with variance
	for i, shift, type, var in zip(range(len(aemol.structure['types'])), props['shift'], aemol.structure['types'], props['shift_var']):
		string = " {atom:<5d}, {shift:<15.8f}, {type:<5d}, {variance:<15.8f}\\".format(atom=i+count_from, shift=shift, type=type, variance=var)
		lines.append(string)

	lines.append('')
	lines.append('> <NMREDATA_J>')
	# Print couplings with variance and label
	for i in range(len(aemol.structure['types'])):
		for j in range(len(aemol.structure['types'])):
			if i >= j:
				continue
			if aemol.structure['path_len'][i][j] == 0:
				continue
			if props['coupling'][i][j] == 0 and not write_zeros:
				continue
			string = " {a1:<10d}, {a2:<10d}, {coupling:<15.8f}, {label:<10s}, {var:<15.8f}".format(a1=i+count_from,
																									a2=j+count_from,
																									coupling=props['coupling'][i][j],
																									label=aemol.pair_properties['nmr_types'][i][j],
																									var=props['coupling_var'][i][j])

			lines.append(string)

	#End section
	lines.append('')
	lines.append('$$$$')

	# Print assembled lines to output file
	with open(outfile, 'w') as f:
		for line in lines:
			print(line, file=f)
