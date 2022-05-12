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


def get_rd_atom_descriptors(rdmol):

    atom_properties = {}

    degree = np.zeros(len(rdmol.GetAtoms()), dtype=np.int64)
    #total_valency = np.zeros(len(rdmol.GetAtoms()), dtype=np.int64)
    #bitvector = np.zeros(len(rdmol.GetAtoms()), dtype=np.int64)
    InRing = np.full(len(rdmol.GetAtoms()), False)
    IsAromatic = np.full(len(rdmol.GetAtoms()), False)

    for idx, atom in enumerate(rdmol.GetAtoms()):
        degree[idx] = atom.GetTotalDegree()
        #total_valency[idx] = atom.GetTotalValence()
        #bitvector[idx] = atom.GetExplicitBitVectProp()
        InRing[idx] = atom.IsInRing()
        IsAromatic[idx] = atom.GetIsAromatic()

    atom_properties['degree'] = degree
    #atom_properties['total_valency'] = total_valency
    #atom_properties['bitvector'] = bitvector
    atom_properties['is_InRing'] = InRing
    atom_properties['isAromatic'] = IsAromatic

    return atom_properties
