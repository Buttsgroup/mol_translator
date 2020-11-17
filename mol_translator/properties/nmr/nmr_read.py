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


# functions for reading nmr properties from files


# nmredata g09, g16, orca

def nmr_read(file, ftype, prop):
	if prop == 'nmr':
		if ftype == 'g09':
			shift, coupling = g09_nmrread(file)
		elif ftype == 'g16':
			shift, coupling = g16_nmrread(file)
		elif ftype == 'orca':
			shift, coupling = orca_nmrread(file)
		elif ftype == 'nmredata':
			shift, _, coupling, _ = nmredata_nmrread(file)
		else:
			print("ftype not recognised, no function to read:", ftype)
			raise ValueError('Cannout read ftype: ', ftype, ' for property: ', prop)
		return shift, coupling
	if prop == 'nmr_var':
		_, shift_var, _, coupling_var = nmredata_nmrread(file)
		return shift_var, coupling_var
	else:
		print('property not recognised or function not written yet !')
		raise ValueError('Cannot output property: ', prop)


# Read NMR inftypeion from ORCA NMR log files
def orca_nmrread(file):

	shiftswitch = False
	shifts = []
	cplswitch = False
	cpls = []
	strucswitch = False
	atoms = 0

	with open(file ,'r') as f:
		for line in f:

			if 'CARTESIAN COORDINATES (A.U.)' in line:
				strucswitch = True

			if 'CHEMICAL SHIELDING SUMMARY (ppm)' in line:
				shiftswitch = True
				cplswitch = False
				strucswitch = False

			if 'NMR SPIN-SPIN COUPLING CONSTANTS' in line:
				shiftswitch = False
				strucswitch = False
				cplswitch = True

			items = line.split()
			if len(items) == 0:
				continue

			if strucswitch and len(items) == 8 and items[0] != 'NO':
				try:
					atoms = int(items[0]) + 1
				except:
					continue

			if shiftswitch and len(items) == 4:
				try:
					int(items[0])
					float(items[2])
					float(items[3])
				except:
					continue

				shifts.append([int(items[0]), float(items[2])])

			if cplswitch and len(items) in [6, 10]:
				if 'NUCLEUS A' in line and len(items) == 10:
					a = int(items[4])
					b = int(items[9])

				if items[0] == 'Total' and len(items) == 6:
					c = float(items[5])

					cpls.append([a, b, c])

	shift = np.zeros((atoms), dtype=np.float64)
	for sh in shifts:
		shift[sh[0]] = sh[1]

	coupling = np.zeros((atoms, atoms), dtype=np.float64)
	for cp in cpls:
		coupling[cp[0]][cp[1]] = cp[2]
		coupling[cp[1]][cp[0]] = cp[2]

	return shift, coupling

def nmredata_nmrread(file):
	# Input:
	#	file: filename
	#	atoms: number of atoms in molecule (read sdf part first)

	# Returns:
	#	shift_array, array of chemical shifts (1D numpy array)
	#	shift_var, array of chemical shift variances (used in machine learning) (1D numpy array)
	#	coupling_array, array of coupling constants (2D numpy array)
	#	coupling_len, array of bond distance between atoms (2D numpy array)
	#	coupling_var, array of coupling constant variances (used in machine learning) (2D numpy array)

	atoms = 0
	with open(file, 'r') as f:
		for line in f:
			if len(line.split()) == 16:

				atoms += 1

	with open(file, 'r') as f:
		for line in f:
			if 'V2000' in line or len(line.split()) == 12:
				chkatoms = int(line.split()[0])

	# check for stupid size labelling issue
	if atoms != chkatoms:
		for i in range(1, len(str(chkatoms))):
			if atoms == int(str(chkatoms)[:-i]):
				chkatoms = atoms
				break

	assert atoms == chkatoms
	# Define empty arrays
	shift_array = np.zeros(atoms, dtype=np.float64)
	# Variance is used for machine learning
	shift_var = np.zeros(atoms, dtype=np.float64)
	coupling_array = np.zeros((atoms, atoms), dtype=np.float64)
	coupling_len = np.zeros((atoms, atoms), dtype=np.int64)
	# Variance is used for machine learning
	coupling_var = np.zeros((atoms, atoms), dtype=np.float64)

	# Go through file looking for assignment sections
	with open(file, 'r') as f:
		shift_switch = False
		cpl_switch = False
		for line in f:
			if '> <NMREDATA_ASSIGNMENT>' in line:
				shift_switch = True
			if '> <NMREDATA_J>' in line:
				shift_switch = False
				cpl_switch = True
			# If shift assignment label found, process shift rows
			if shift_switch:
				# Shift assignment row looks like this
				#  0    , -33.56610000   , 8    , 0.00000000     \
				items = line.split()
				try:
					int(items[0])
				except:
					continue
				print(items)
				shift_array[int(items[0])] = float(items[2])
				shift_var[int(items[0])] = float(items[6])
			# If coupling assignment label found, process coupling rows
			if cpl_switch:
				# Coupling row looks like this
				#  0         , 4         , -0.08615310    , 3JON      , 0.00000000
				# ['0', ',', '1', ',', '-0.26456900', ',', '5JON', ',', '0.00000000']
				items = line.split()
				try:
					int(items[0])
				except:
					continue
				length = int(items[6].strip()[0])
				coupling_array[int(items[0])][int(items[2])] = float(items[4])
				coupling_array[int(items[2])][int(items[0])] = float(items[4])
				coupling_var[int(items[0])][int(items[2])] = float(items[8])
				coupling_var[int(items[2])][int(items[0])] = float(items[8])
				coupling_len[int(items[0])][int(items[2])] = length
				coupling_len[int(items[2])][int(items[0])] = length

	return shift_array, shift_var, coupling_array, coupling_var#, coupling_len



def g09_nmrread(file):
	with open(file, 'r') as f:
		for line in f:
			if 'NAtoms=' in line:
				items = line.split()
				atomnumber = int(items[1])

	shift_array = np.zeros(atomnumber, dtype=np.float64)
	switch = False
	# Go through file to find magnetic shielding tensors
	with open(file, 'r') as f_handle:
		for line in f_handle:
			# If tensor label is found, activate switch
			if "SCF GIAO Magnetic shielding tensor (ppm)" in line:
				switch = True
			# This label comes at the end of the tensor section, so deactivate switch
			if "Fermi Contact" in line:
				switch = False

			if switch:
				# Find isotropic tensor lines
				if "Isotropic" in line:
					items = line.split()
					try:
						num = int(items[0])
					except:
						continue
					# Isotropic shielding tensor is the 5th item (0, 1, 2, 3, '4')
					shift_array[num-1] = float(items[4])

	# Define empty array for couplings
	couplings = np.zeros((atomnumber, atomnumber), dtype=float)
	# Go through file to find coupling constants
	with open(file, 'r') as f:
		switch = False
		for line in f:
			# If coupling label is found, activate switch
			if "Total nuclear spin-spin coupling J (Hz):" in line:
				switch = True
				continue
			# This label comes at the end of the coupling section, so deactivate switch
			elif "End of Minotr" in line:
				switch = False
				continue

			if switch:
				# All coupling lines contain "D", all index lines do not
				if "D" not in line:
					# Get indices for this section
					tokens = line.split()
					i_indices = np.asarray(tokens, dtype=int)
				else:
					# Assign couplings (array is diagonalised in log file, so this is fiddly)
					tokens = line.split()
					index_j = int(tokens[0]) - 1
					for i in range(len(tokens)-1):
						index_i = i_indices[i] - 1
						coupling = float(tokens[i+1].replace("D","E"))
						couplings[index_i][index_j] = coupling
						couplings[index_j][index_i] = coupling

	return shift_array, couplings

def g16_nmrread(file):
	with open(file, 'r') as f:
		for line in f:
			if 'NAtoms=' in line:
				items = line.split()
				atomnumber = int(items[1])

	shift_array = np.zeros(atomnumber, dtype=np.float64)
	switch = False
	# Go through file to find magnetic shielding tensors
	with open(file, 'r') as f_handle:
		for line in f_handle:
			# If tensor label is found, activate switch
			if "SCF GIAO Magnetic shielding tensor (ppm)" in line:
				switch = True
			# This label comes at the end of the tensor section, so deactivate switch
			if "Fermi Contact" in line:
				switch = False

			if switch:
				# Find isotropic tensor lines
				if "Isotropic" in line:
					items = line.split()
					try:
						num = int(items[0])
					except:
						continue
					# Isotropic shielding tensor is the 5th item (0, 1, 2, 3, '4')
					shift_array[num-1] = float(items[4])

	# Define empty array for couplings
	couplings = np.zeros((atomnumber, atomnumber), dtype=float)
	# Go through file to find coupling constants
	with open(file, 'r') as f:
		switch = False
		for line in f:
			# If coupling label is found, activate switch
			if "Total nuclear spin-spin coupling J (Hz):" in line:
				switch = True
				continue
			# This label comes at the end of the coupling section, so deactivate switch
			elif "End of Minotr" in line:
				switch = False
				continue

			if switch:
				# All coupling lines contain "D", all index lines do not
				if "D" not in line:
					# Get indices for this section
					tokens = line.split()
					i_indices = np.asarray(tokens, dtype=int)
				else:
					# Assign couplings (array is diagonalised in log file, so this is fiddly)
					tokens = line.split()
					index_j = int(tokens[0]) - 1
					for i in range(len(tokens)-1):
						index_i = i_indices[i] - 1
						coupling = float(tokens[i+1].replace("D","E"))
						couplings[index_i][index_j] = coupling
						couplings[index_j][index_i] = coupling

	return shift_array, couplings
