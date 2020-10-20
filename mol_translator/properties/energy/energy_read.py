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


def energy_read(file, ftype, prop):
	if prop == 'scf':
		if format == 'g09':
			energy = g09_scfread(file)
		elif format == 'g16':
			energy = g16_scfread(file)
		elif format == 'orca':
			energy = orca_scfread(file)


def g09_scfread(file):
	energy = -404

	with open(file, 'r') as f:
		for line in f:
			if 'Sum of electronic and thermal Free Energies' in line:
				items = line.split()
				energy = float(items[7])

	return energy

def g16_scfread(file):
	energy = -404

	with open(file, 'r') as f:
		for line in f:
			if 'Sum of electronic and thermal Free Energies' in line:
				items = line.split()
				energy = float(items[7])

	return energy

def orca_scfread(file):
	energy = 0.0
	with open(file ,'r') as f:
		for line in f:
			if 'FINAL SINGLE POINT ENERGY' in line:
				items=line.split()
				energy = float(items[-1])

	return energy
