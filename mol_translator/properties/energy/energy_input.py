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

def make_orca_optin(prefs, molname, aemol, outfile):
	# Input:
	#	prefs: preferences dictionary
	#	molname: name of molecule
	#	xyz: xyz coordinates of conformer
	#	types: type list of conformer (numeric)
	#	path: path to molecule folder

	# Returns: filename for input file

	# Get preferences from prefs
	charge = prefs['mol']['charge']
	multiplicity = prefs['mol']['multiplicity']
	functional = prefs['optimisation']['functional']
	basis_set = prefs['optimisation']['basisset']
	solvent = prefs['optimisation']['solvent']
	if solvent != None:
		solventmodel = prefs['optimisation']['solventmodel']
	direct_cmd_line_opt = prefs['optimisation']['custom_cmd_line']
	processors = prefs['optimisation']['processors']
	# Get periodic table
	Periodic_table = Get_periodic_table()
	# Define instruction line for ORCA
	instr = '! ' + str(functional) + ' ' + str(basis_set) + ' TightSCF OPT miniprint'
	# Add parallel option if multiple processors requested
	if processors != 1:
		instr += ' PAL{0:<d}'.format(processors)
	# Add solvent model/solvent if requested
	if solvent != None:
		instr += ' CPCM(' + solvent + ')'
	# If direct line input specified then overwrite all of this
	if direct_cmd_line_opt:
		instr = direct_cmd_line_opt

	# Construct file strings
	strings = []
	strings.append(instr)
	strings.append('')
	strings.append("* xyz {0:<1d} {1:<1d}".format(charge, multiplicity))
	for i in range(len(aemol.structure['xyz'])):
		str_type = Periodic_table[aemol.structure['types'][i]]
		string = " {0:<2s}        {1:>10.5f}        {2:>10.5f}        {3:>10.5f}".format(str_type, aemol.structure['xyz'][i][0], aemol.structure['xyz'][i][1], aemol.structure['xyz'][i][2])
		strings.append(string)
	strings.append('*')
	strings.append('')
	strings.append('%geom')
	strings.append('     AddExtraBonds true         # switch on/off assigning bonds to atom pairs that are')
	strings.append('                                #  connected by more than <Max_Length> bonds and are less')
	strings.append('                                #  than <MaxDist> Ang. apart (default true)')
	strings.append('     AddExtraBonds_MaxLength 10 # cutoff for number of bonds connecting the two')
	strings.append('                                #  atoms (default 10)')
	strings.append('     AddExtraBonds_MaxDist 5    # cutoff for distance between two atoms (default 5 Ang.)')
	strings.append('end')
	# Write file
	with open(outfile, 'w') as f_handle:
		for string in strings:
			print(string, file=f_handle)


def make_g09_optin(prefs, molname, aemol, outfile):
## Inputs:
	#  molname      = molecule name, string. 'Progesterone'
	#  xyz          = xyz coordinates, Nx3 numpy array.
	#  type         = atom types, list(or 1D array) length N.
	#  charge       = molecule charge, integer.
	#  multiplicity = multiplcity of molecule, integer.
	#  memory       = memory value to use for gaussian, integer.
	#  processors   = number of processoers to use for gaussian, integer.
	#  instr        = instruction line for gaussian, string. 'opt=tight freq mpw1pw91/6-311g(d,p) scrf=(iefpcm,solvent=chloroform) geom=distance MaxDisk=50GB'
	Periodic_table = Get_periodic_table()

	charge = prefs['mol']['charge']
	multiplicity = prefs['mol']['multiplicity']
	memory = prefs['optimisation']['memory']
	processors = prefs['optimisation']['processors']
	opt = prefs['optimisation']['opt']
	freq = prefs['optimisation']['freq']
	functional = prefs['optimisation']['functional']
	basis_set = prefs['optimisation']['basisset']
	grid = prefs['optimisation']['grid']
	solvent = prefs['optimisation']['solvent']
	if solvent != None:
		solventmodel = prefs['optimisation']['solventmodel']
	freq = True # Without this the energies we get arent useful for boltzmann weighting


	try:
		int(memory)
		int(processors)
	except:
		return

	if freq:
		instr = 'opt=' + str(opt) + ' freq ' + str(functional) + '/' + str(basis_set) + ' integral=' + str(grid) + ' MaxDisk=50GB'
	else:
		instr = 'opt=' + str(opt) + ' ' + str(functional) + '/' + str(basis_set) + ' integral=' + str(grid) + ' MaxDisk=50GB'

	instr = 'opt=' + str(opt) + ' freq ' + str(functional) + '/' + str(basis_set) + ' integral=' + str(grid) + ' MaxDisk=50GB'

	if solvent != None:
		instr += ' scrf=(' + str(solventmodel) + ',solvent=' + str(solvent) + ')'

	comfile = outfile


	with open(comfile, 'w') as f_handle:
		strings = []
		strings.append("%Chk=OPT/{0:<1s}_OPT1".format(molname))
		strings.append("%NoSave")
		strings.append("%mem={0:<1d}GB".format(memory))
		strings.append("%NProcShared={0:<1d}".format(processors))
		strings.append("# {0:<1s}".format(instr))
		strings.append("")
		strings.append("{0:<1s} OPT".format(molname))
		strings.append("")
		strings.append("{0:<1d} {1:<1d}".format(charge, multiplicity))
		for string in strings:
			print(string, file=f_handle)
		for i in range(len(aemol.structure['xyz'])):
			str_type = Periodic_table[aemol.structure['types'][i]]
			string = " {0:<2s}        {1:>10.5f}        {2:>10.5f}        {3:>10.5f}".format(str_type, aemol.structure['xyz'][i][0], aemol.structure['xyz'][i][1], aemol.structure['xyz'][i][2])
			print(string, file=f_handle)
		for i in range(4):
			print("", file=f_handle)
		strings = []

	return comfile
