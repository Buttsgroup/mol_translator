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

import sys

from mol_translator.comp_chem.gaussian.gaussian_read import gauss_scf_read
from mol_translator.comp_chem.orca.orca_read import orca_scf_read


def energy_read(file: str, prop: str = 'scf', format: str = 'gauss') -> float:
    """
    Function for reading and processing energy data stored in the file.

    :param file: str, filepath of the file containing the data
    :param prop: str, the property to read in, defaults to 'scf'
    :param format: str, the format of the file to process, defaults to 'gauss' to read log files
    :return mc_array: float, scf value of the molecule
    """
    if prop == 'scf':
        if format == 'gauss':
            energy = gauss_scf_read(file)
        elif format == 'orca':
            energy = orca_scf_read(file)

        else:
            print('format not recognised:', format)
            sys.exit(0)
    else:
        print('Prop not recognised:', prop)
        sys.exit(0)

    return energy
