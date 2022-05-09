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


def orca_scf_read(file):
    """
    Reads orca .log files and extracts the energy value

    :param file: filepath of the .log file which contains the optimization logs
    :return energy: int value representing the energy
    """
    energy = 0.0
    with open(file, 'r') as f:
        for line in f:
            if 'FINAL SINGLE POINT ENERGY' in line:
                items = line.split()
                energy = float(items[-1])

    return energy


def orca_nmr_read(file: str):
    """
    Read orca .log files and extracts the nmr parameters if present

    :param file: filepath of the .log file which contains the nmr logs
    :return shift_array and coupling array: numpy array of the chemical tensors for respective atom position and atom-pair interaction
    """

    shiftswitch = False
    shifts = []
    cplswitch = False
    cpls = []
    strucswitch = False
    atoms = 0

    with open(file, 'r') as f:
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
