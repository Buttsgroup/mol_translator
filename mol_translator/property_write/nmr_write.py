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

from ..file_creation.structure_write import write_mol_tosdf

# Write an nmrmol object to an nmredata file
def write_nmredata(outfile, label="", mol,
                    shift=None, coupling=None, shift_var=None, coupling_var=None,
                    info=None):

    atoms = len(mol['types'])
    if shift is None:
        shift = np.zeros(atoms, dtype=np.float64)
    if shift_var is None:
        shift_var = np.zeros(atoms, dtype=np.float64)
    if coupling is None:
        coupling = np.zeros((atoms, atoms), dtype=np.float64)
    if coupling_var is None:
        coupling_var = np.zeros((atoms, atoms), dtype=np.float64)

    sdfstrings = write_mol_tosdf(mol, "", label, True)

    for line in sdfstrings:
        line.append(line)
	# assignment section
	lines.append('')
	lines.append('> <NMREDATA_ASSIGNMENT>')
	# Print chemical shifts with variance
	for i, shift, type, var in zip(range(len(mol['types'])), mol['shift'], mol['types'], mol['shift_var']):
		string = " {atom:<5d}, {shift:<15.8f}, {type:<5d}, {variance:<15.8f}\\".format(atom=i, shift=shift, type=type, variance=var)
		lines.append(string)

	lines.append('')
	lines.append('> <NMREDATA_J>')
	# Print couplings with variance and label
	for i in range(len(mol['types'])):
		for j in range(len(mol['types'])):
			if i >= j:
				continue
			if mol['coupling_len'][i][j] == 0:
				continue
			label = labelmaker(i, j, mol)
			string = " {a1:<10d}, {a2:<10d}, {coupling:<15.8f}, {label:<10s}, {var:<15.8f}".format(a1=i,
																									a2=j,
																									coupling=mol['coupling'][i][j],
																									label=label,
																									var=mol['coupling_var'][i][j])

			lines.append(string)

	# Print assembled lines to output file
	with open(outfile, 'w') as f:
		for line in lines:
			print(line, file=f)
