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


def write_mol_toxyz(mol, outname, label="tmp"):
	periodic_table = Get_periodic_table()
	with open(outname, 'w') as f:
		print(len(mol['types']), file=f)
		print(label, file=f)

		for i in range(len(mol['types'])):
			string = "{i:<10s}\t{x:<10.6f}\t{y:<10.6f}\t{z:<10.6f}".format(i=periodic_table[mol['types'][i]],
																			x=mol['xyz'][i][0],
																			y=mol['xyz'][i][1],
																			z=mol['xyz'][i][2])
			print(string, file=f)

def write_mol_tosdf(mol, outfile, label="tmp", stringsonly=False):
	# Get periodic table
	periodic_table = Get_periodic_table()

	# Assume molecule is not chiral
	chiral = 0
	# Determine number of bonds
	bonds = 0
	for at1 in range(len(mol['types'])):
		for at2 in range(at1, len(mol['types'])):

			if mol['conn'][at1][at2] >= 1:
				bonds += 1

	atoms = len(mol['types'])

	# Start putting together file lines
	lines = []
	# print molecule name
	if '/' in outfile:
		lines.append(outfile.split('/')[-1].split('.')[0])
	else:
		lines.append(outfile.split('.')[0])
	# print file and author
	lines.append('auto-ENRICH - 2020 - ButtsGroup')
	lines.append('')

	# Structure section
	string = '{atoms:>3d}{bonds:>3d}  0  0{chiral:>3d}  0  0  0  0  0  1 V2000'.format(atoms=atoms,
																							bonds=bonds,
																							chiral=chiral)
	lines.append(string)
	# Print xyz coordinates and types
	for i, xyz in enumerate(mol['xyz']):
		string = '{x:>10.4f}{y:>10.4f}{z:>10.4f} {typechar:>3s} 0  0  0  0  0  0  0  0  0  0  0  0'.format(x=xyz[0],
																												y=xyz[1],
																												z=xyz[2],
																												typechar = periodic_table[mol['types'][i]])
		lines.append(string)
	# Print bonds and bond types
	for at1 in range(len(mol['types'])):
		for at2 in range(at1, len(mol['types'])):

			if mol['conn'][at1][at2] >= 1:
				string = '{at1:>3d}{at2:>3d}{bond:>3d}  0  0  0  0'.format(at1=at1+1,
																				at2=at2+1,
																				bond=mol['conn'][at1][at2])
				lines.append(string)
	# Terminate structure section
	lines.append('M\tEND'.format())

	if stringsonly:
		return lines
	else:
		with open(outfile, 'w') as f:
			for line in lines:
				print(line, file=f)
