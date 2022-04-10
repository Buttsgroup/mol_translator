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

from typing import Type
from mol_translator import Aemol
import numpy as np


def calc_pops(aemols: Type[Aemol], temp: int = 298):

    e_array = np.zeros(len(aemols), dtype=np.float64)
    for c, aemol in enumerate(aemols):
        e_array[c] = aemol.mol_properties['energy']

    kj_array = e_array * 2625.5
    min_val = np.amin(kj_array)
    rel_array = (kj_array - min_val)
    exp_array = -(rel_array*1000) / float(8.31*temp)
    exp_array = np.exp(exp_array)
    sum_val = np.sum(exp_array)

    pop_array = exp_array / sum_val

    for a, aemol in enumerate(aemols):
        aemol.mol_properties['pop'] = pop_array[a]


def boltzmann_average(aemols: Type[Aemol], pair_props: list = ['coupling'], atom_props: list = ['shift']):
    atoms = len(aemols[0].structure['types'])

    new_atom_dict = {}
    new_pair_dict = {}
    for prop in atom_props:
        new_atom_dict[prop] = np.zeros(atoms, dtype=np.float64)
    for prop in pair_props:
        new_pair_dict[prop] = np.zeros((atoms, atoms), dtype=np.float64)

    for aemol in aemols:
        for prop in atom_props:
            new_atom_dict[prop] += aemol.atom_properties[prop] * \
                aemol.mol_properties['pop']
        for prop in pair_props:
            new_pair_dict[prop] += aemol.pair_properties[prop] * \
                aemol.mol_properties['pop']

    return new_atom_dict, new_pair_dict
