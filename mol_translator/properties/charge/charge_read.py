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

import sys
from mol_translator.comp_chem.gaussian.gaussian_read import gauss_mc_read

def charge_read(file, prop, format):
    if prop ==  'mc':
        if format == 'gauss':
            return gauss_mc_read(file)
        else:
            print('format not recognised:', format)
            sys.exit(0)
    else:
        print('Prop not recognised:', prop)
        sys.exit(0)
