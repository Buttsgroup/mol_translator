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

import numpy as np


def gauss_scf_read(file: str):
    energy = -404

    with open(file, 'r') as f:
        for line in f:
            if 'Sum of electronic and thermal Free Energies' in line:
                items = line.split()
                energy = float(items[7])

    return energy


def gauss_mc_read(file):
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


def gauss_nmr_read(file: str):
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
                        coupling = float(tokens[i+1].replace("D", "E"))
                        couplings[index_i][index_j] = coupling
                        couplings[index_j][index_i] = coupling

    return shift_array, couplings
