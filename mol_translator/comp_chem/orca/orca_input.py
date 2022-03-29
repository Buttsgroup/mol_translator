# Copyright 2022 Will Gerrard, Calvin Yiu
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

def make_orca_rootline(prefs: dict[str, str]) -> str:
    """
    Function for the automation of the root line required for setting up the gaussian calculation

    :param prefs: dictionary containing key parameters for the gaussian calculation
    :return: str root_line
    """
    if prefs['custom_cmd_line']:
        return prefs['custom_cmd_line']

    if prefs['calc_type'] == ('optimisation' or 'opt'):

        root_line = f"! {str(prefs['functional'])} {str(prefs['basis_set'])} {prefs['opt']} OPT"

        if pref['dispersion_correction']:
            root_line += f" {prefs['dispersion_correction']}"

        if prefs['freq']:
            root_line += f" FREQ"

        if prefs['solventmodel'] == 'CPCM':
    		root_line += f' CPCM( {prefs['solvent']})'

        if prefs['processors'] != 1:
            root_line += f' PAL{prefs['processors']}'

        if prefs['print_format']:
            root_line += f' {prefs['print_format']}'

        return root_line

    elif prefs['calc_type'] == ('nmr' or 'NMR'):
        # change elif assertion to switch case if RDKit updates to python 3.10

        root_line = f"! {str(prefs['functional'])} {str(prefs['basis_set'])} NMR"

    	if solvent != None:
    		root_line += f' scrf=( {str(prefs['solventmodel'])} , solvent={str(prefs['solvent'])} )'

        return root_line

    else:
        raise ValueError(f"Calculation type is set to {prefs['calc_type']}, current only 'optimisation' & 'nmr' are supported")

def write_orca_inp(prefs: dict, molname: str, aemol: Aemol, root_line: str, outfile: str) -> file:

    assert type(prefs['memory']) is int, f"{prefs['memory']} is {type(prefs['memory'])}, should be an integer value"
    assert type(prefs['processors']) is int, f"{prefs['processors']} is {type(prefs['processors'])}, should be an integer value"

	strings = []
    strings.append(f"{root_line}")
	strings.append("")
	strings.append(f"* xyz {prefs['charge']} {prefs['multiplicity']}")

    if prefs['calc_type'] == ('nmr' or 'NMR'):
        atoms_list = []

	for i in range(len(aemol.structure['xyz'])):
		atom_type = periodic_table[aemol.structure['types'][i]]
        if prefs['calc_type'] == ('nmr' or 'NMR'):
            atoms_list.append(atom_type)
		string = f" {atom_type:<2.2s}{aemol.structure['xyz'][i][0]:>18.6f}{aemol.structure['xyz'][i][1]:>18.6f}{aemol.structure['xyz'][i][2]:>18.6f}"
        strings.append(string)

    strings.append('*')

    if prefs['calc_type'] == ('optimisation' or 'opt'):
    	strings.append('')
    	strings.append('%geom')
    	strings.append('     AddExtraBonds true         # switch on/off assigning bonds to atom pairs that are')
    	strings.append('                                #  connected by more than <Max_Length> bonds and are less')
    	strings.append('                                #  than <MaxDist> Ang. apart (default true)')
    	strings.append('     AddExtraBonds_MaxLength 10 # cutoff for number of bonds connecting the two')
    	strings.append('                                #  atoms (default 10)')
    	strings.append('     AddExtraBonds_MaxDist 5    # cutoff for distance between two atoms (default 5 Ang.)')

    if prefs['calc_type'] == ('nmr' or 'NMR'):
        strings.append('%EPRNMR')

        if prefs['functional'] == 'wB97XD':
            strings.append('     GIAO_2el = GIAO_2el_RIJCOSX')

        for atom in list(set(atoms_list)):
    		strings.append(f'       Nuclei = ALL {atom:<2s} {{SHIFT, SSALL}}')

        if prefs['spin_thresh']:
            strings.append(f'SpinSpinRThresh {prefs['spin_thresh']}')

    if prefs['solventmodel'] == 'SMD':
        strings.append('%CPCM SMD TRUE')
        strings.append(f'     SMDSOLVENT "{prefs['solvent']}"')

	strings.append('end')

    with open(outfile, 'w') as f:
        for string in strings:
            print(string, file=f)





