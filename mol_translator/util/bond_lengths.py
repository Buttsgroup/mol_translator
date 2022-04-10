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

# H, C  N , O, Cl, F, Si, P, S, Br
# Potential max average bond = C-S 2.6 Angstrom

def get_bond_lengths():
    bond_lengths = {
        'CH': [0.90, 1.20],
        'CC': [[1.30, 1.70], [1.20, 1.40], [1.10, 1.30], [1.10, 1.30]],
        'CN': [[1.30, 1.60], [1.15, 1.35], [1.00, 1.30], [1.20, 1.40]],
        'CO': [[1.30, 1.60], [1.10, 1.30]],
        'CCl': [1.60, 1.90],
        'CSi': [2.30, 2.60],
        'CP': [1.70, 2.00],
        'CS': [1.70, 2.00],
        'CBr': [1.80, 2.00],
        'NO': [1.10, 1.60],
        'SS': [1.80, 2.20],
        'SiSi': [2.30, 2.60],
        'NN': [[1.15, 1.35], [1.25, 1.45]],
    }
