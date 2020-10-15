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

def make_orca_nmrin(prefs, molname, aemol, outfile):
	# Input:
	#	prefs: preferences dictionary
	#	molname: name of molecule
	#	xyz: xyz coordinates of conformer
	#	aemol.structure['types']: type list of conformer (numeric)
	#	path: path to molecule folder

	# Returns: input file path/name

	# Get values from preferences
	charge = prefs['mol']['charge']
	multiplicity = prefs['mol']['multiplicity']
	functional = prefs['NMR']['functional']
	basis_set = prefs['NMR']['basisset']
	aux_basis_set = prefs['NMR']['aux_basis_set']
	solvent = prefs['NMR']['solvent']
	direct_cmd_line_nmr = prefs['NMR']['custom_cmd_line']
	processors = prefs['NMR']['processors']
	# Get periodic table
	Periodic_table = Get_periodic_table()
	# Construct instruction line for ORCE
	instr = '! ' + str(functional) + ' ' + str(basis_set) + ' ' + str(aux_basis_set) +  '  TightSCF miniprint' + ' NMR '
	# Add parallel option if multiple processors requested
	if processors != 1:
		instr += ' PAL{0:<d}'.format(processors)
	# Add solvent model/solvent if requested
	if solvent != 'none':
		instr += ' CPCM(' + solvent + ')'
	# If direct line input specified then overwrite all of this
	if direct_cmd_line_nmr:
		instr = direct_cmd_line_nmr
	# Define input file path/name
	# Construct file strings
	strings = []
	strings.append(instr)
	strings.append("")
	strings.append("* xyz {0:<1d} {1:<1d}".format(charge, multiplicity))
	for i in range(len(aemol.structure['xyz'])):
		str_type = Periodic_table[aemol.structure['types'][i]]
		string = " {0:<2s}        {1:>10.6f}        {2:>10.6f}        {3:>10.6f}".format(str_type, aemol.structure['xyz'][i][0], aemol.structure['xyz'][i][1], aemol.structure['xyz'][i][2])
		strings.append(string)
	strings.append('*')
	strings.append('%eprnmr')
	# Needed for the functional we commonly use, ORCA shouted at me
	strings.append("       GIAO_2el = GIAO_2el_RIJCOSX")
	for type in prefs['NMR']['shift_nuclei']:
		strings.append("       Nuclei = all {type:<2s}".format(type=type) + '  { shift }')
	for type in prefs['NMR']['spin_nuclei']:
		strings.append("       Nuclei = all {type:<2s}".format(type=type) + '  { ssall }')
	strings.append('SpinSpinRThresh {0:<f}'.format(prefs['NMR']['spin_thresh']))
	strings.append('end')
	# Write file
	with open(outfile, 'w') as f_handle:
		for string in strings:
			print(string, file=f_handle)


def make_g09_nmrin(prefs, molname, aemol, outname):

	# Get values from preferences
	charge = prefs['mol']['charge']
	multiplicity = prefs['mol']['multiplicity']
	functional = prefs['NMR']['functional']
	basis_set = prefs['NMR']['basisset']
	solvent = prefs['NMR']['solvent']
	solventmodel = prefs['NMR']['solventmodel']
	mixed = prefs['NMR']['mixed']
	direct_cmd_line_nmr = prefs['NMR']['custom_cmd_line']
	processors = prefs['NMR']['processors']
	memory = prefs['NMR']['memory']

	Periodic_table = Get_periodic_table()
	try:
		int(memory)
		int(processors)
	except:
		return

	if mixed == "True":
		instr='nmr(giao,spinspin,mixed)' + str(functional) + '/' + str(basis_set) + ' maxdisk=50GB'
	else:
		instr='nmr(giao,spinspin)' + str(functional) + '/' + str(basis_set) + ' maxdisk=50GB'

	if solvent != 'none':
		instr += ' scrf=(' + str(solventmodel) + ',solvent=' + str(solvent) + ')'

	comfile = outname
	chkfile = outname.strip('.com') + '.chk'

	strings = []
	strings.append("%Chk={0:<1s}".format(chkfile))
	strings.append("%NoSave")
	strings.append("%mem={0:<1d}GB".format(memory))
	strings.append("NProcShared={0:<1d}".format(processors))
	strings.append("#T {0:<1s}".format(instr))
	strings.append("")
	strings.append("{0:<1s} NMR".format(molname))
	strings.append("")
	strings.append("{0:<1d} {1:<1d}".format(charge, multiplicity))
	for i in range(len(aemol.structure['xyz'])):
		str_type = Periodic_table[aemol.structure['types'][i]]
		strings.append(" {0:<2s}        {1:>10.6f}        {2:>10.6f}        {3:>10.6f}".format(str_type, aemol.structure['xyz'][i][0], aemol.structure['xyz'][i][1], aemol.structure['xyz'][i][2]))
	for i in range(4):
		strings.append("")

	with open(comfile, 'w') as f:
		for string in strings:
			print(string, file=f)

	return 0
