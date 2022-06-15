# Copyright 2022 Will Gerrard, Calvin Yiu
# This file is part of autoenrich.

# autoenrich is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# autoenrich is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.

from mol_translator.util.periodic_table import Get_periodic_table
from typing import Type, Dict


def make_gaussian_rootline(prefs: Dict[str, str]) -> str:
    """
    Function for the automation of the root line required for setting up the gaussian calculation

    :param prefs: dictionary containing key parameters for the gaussian calculation
    :return: str root_line
    """

    if prefs['calc_type'] == ('optimisation' or 'opt'):
        if prefs['freq']:
            root_line = f"opt={str(prefs['opt'])} freq {str(prefs['functional'])}/{str(prefs['basis_set'])} integral={str(prefs['grid'])} MaxDisk=50GB"

        else:
            root_line = f"opt={str(prefs['opt'])} {str(prefs['functional'])}/{str(prefs['basis_set'])} integral={str(prefs['grid'])} MaxDisk=50GB"

        if prefs['solvent'] != None:
            root_line += f" scrf=( {str(prefs['solvent_model'])} , solvent={str(prefs['solvent'])} )"

        return root_line

    elif prefs['calc_type'] == ('nmr' or 'NMR'):
        # change elif assertion to switch case if RDKit updates to python 3.10
        if prefs['mixed'] == True:
            root_line = f"nmr(giao,spinspin,mixed) {str(prefs['functional'])}/{str(prefs['basis_set'])} maxdisk=50GB"

        else:
            root_line = f"nmr(giao,spinspin) {str(prefs['functional'])}/{str(prefs['basis_set'])} maxdisk=50GB"

        if prefs['solvent'] != None:
            root_line += f" scrf=( {str(prefs['solvent_model'])} , solvent={str(prefs['solvent'])} )"

        return root_line

    else:
        raise ValueError(
            f"Calculation type is set to {prefs['calc_type']}, current only 'optimisation' & 'nmr' are supported")


def write_gaussian_com(prefs: dict, molname: str, aemol: Type, root_line: str, outfile: str) -> None:
    """
    Generates gaussian input .com files based off prefs dictionary values, requires aemol object, molecule name, root line, and outfile name as input.
    Example preferences:

        >>> prefs = {
                 'charge' : 0,
                 'multiplicity : 1,
                 'calc_type' : 'opt',
                 'memory' : 12,
                 'processor' : 8,
                 'nodes' : 1,
                 'walltime' : '120:00:00'
                 'opt' : 'tight',
                 'freq' : True,
                 'functional' : 'mPW1PW',
                 'basis_set' : '6-311g(d,p)',
                 'solvent' : None,
                 'solvent_model' : None,
                 'grid' : 'ultrafine',
                 'custom_cmd_line' = False,
        }

    :param prefs: dictionary containing key parameters for the gaussian calculation
    :param molname: name of current molecule
    :param aemol: Aemol object of the molecule to process
    :param root_line: the output of make_gaussian_root function, passed as a string
    :param outfile: Name of the output file to write to

    :return: None, writes a physical file used to run gaussian calculation
    """

    assert type(
        prefs['memory']) is int, f"{prefs['memory']} is {type(prefs['memory'])}, should be an integer value"
    assert type(
        prefs['processors']) is int, f"{prefs['processors']} is {type(prefs['processors'])}, should be an integer value"

    strings = []
    if prefs['calc_type'] == ('optimisation' or 'opt'):
        strings.append(f"%Chk={molname}_OPT.chk")
    if prefs['calc_type'] == ('nmr' or 'NMR'):
        strings.append(f"%Chk={molname}_NMR.chk")
        strings.append("%NoSave")
        strings.append(f"%mem={prefs['memory']}GB")
        strings.append(f"%NProcShared={prefs['processors']}")
    if prefs['calc_type'] == ('optimisation' or 'opt'):
        strings.append(f"# {root_line}")
    if prefs['calc_type'] == ('nmr' or 'NMR'):
        strings.append(f"#T {root_line}")
        strings.append("")
    if prefs['calc_type'] == ('optimisation' or 'opt'):
        strings.append(f"{molname} OPT")
    if prefs['calc_type'] == ('NMR' or 'nmr'):
        strings.append(f"{molname} NMR")
        strings.append("")
        strings.append(f"{prefs['charge']} {prefs['multiplicity']}")

        periodic_table = Get_periodic_table()

        for i in range(len(aemol.structure['xyz'])):
            atom_type = periodic_table[aemol.structure['types'][i]]
            string = f" {atom_type:<2.2s}{aemol.structure['xyz'][i][0]:>18.6f}{aemol.structure['xyz'][i][1]:>18.6f}{aemol.structure['xyz'][i][2]:>18.6f}"
        for i in range(4):
            strings.append("")

    with open(outfile, 'w') as f:
        for string in strings:
            print(string, file=f)
