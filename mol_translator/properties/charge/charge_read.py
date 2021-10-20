# Copyright 2020 Will Gerrard, Calvin Yiu
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
import numpy as np

def charge_read(file, prop, format):
    if prop ==  'mc':
        if format == 'g09':
            return g09_mcread(file)
        else:
            print('format not recognised:', format)
            sys.exit(0)
    else:
        print('Prop not recognised:', prop)
        sys.exit(0)

def g09_mcread(file):
    with open(file, 'r') as f:
        for line in f:
            if 'NAtoms=' in line:
                items = line.split()
                atomnumber = int(items[1])

    mc_array = np.zeros(atomnumber, dtype=np.float64)
    switch = False
    # Go through file to find mulliken charge
    with open(file, 'r') as f_handle:
        for line in f_handle:
            # If mulliken charge label is found, activate switch
            if "Mulliken charges:" in line:
                switch = True
            # This label comes at the end of the Mulliken section, so deactivate switch
            if "Sum of Mulliken charges" in line:
                switch = False
            # Only run this code on the first mulliken charge set found
            if switch == True:
                # Find
                items = line.split()
                if len(items) != 3:
                    continue
                # Mulliken charge is the 3rd item (0, 1, 2)
                try:
                    mc_array[int(items[0])-1] = float(items[2])
                except:
                    print(file)

    return mc_array
