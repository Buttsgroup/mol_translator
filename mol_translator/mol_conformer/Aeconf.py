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


def get_dist_array(aemols):

    for t, aemol in enumerate(aemols):
        num_atoms = len(aemol.structure['types'])
        dist_array = np.zeros((num_atoms, num_atoms), dtype=np.float64)
        for i in range(num_atoms):
            for j in range(num_atoms):
                dist_array[i][j] = np.absolute(np.linalg.norm(
                    aemol.structure['xyz'][i] - aemol.structure['xyz'][j]))

        aemol.mol_properties['dist_array'] = dist_array


def redundant_elimination(aemols, geom_threshold=0.1, e_threshold=0.1, redundant_atoms=None, achiral=False):

    elim_list = []
    e_array = np.zeros(len(aemols), dtype=np.float64)

    if redundant_atoms != None:
        redundant_atoms = redundant_atoms.split(",")
        redundant_atoms = list(map(int, redundant_atoms))

        for atoms in redundant_atoms:
            for i in range(aemols):
                dist_arrays[i][int(atoms)-1] = 0
                for k in range(atoms):
                    dist_arrays[i][k][int(atoms)-1] = 0

    with open("redundant_elimination_log.txt", 'a') as l:

        for a, aemol_a in enumerate(aemols):
            e_array[a] = aemol_a.mol_properties['energy']
            dist_array_a = aemol_a.mol_properties['dist_array']

            for b, aemol_b in enumerate(aemols):
                e_array[b] = aemol_b.mol_properties['energy']
                dist_array_b = aemol_b.mol_properties['dist_array']
                if a > b and not b in elim_list:
                    diff = np.absolute(np.sum(dist_array_a - dist_array_b))
                    energy_diff = np.absolute(e_array[a] - e_array[b]) * 2625.5

                    if diff < geom_threshold:
                        if energy_diff < e_threshold:
                            if achiral is False:
                                elim_list.append(a)
                                print(
                                    f"added mol {aemol_a.info['molid']} to eliminated_mol.txt due to geomtric similarity to {aemol_b.info['molid']}", file=l)

                            else:
                                if a - b == 1:
                                    print(
                                        f"energy difference between {aemol_a.info['molid']} & {aemol_b.info['molid']} detected but could be mirror images, please check manually", file=l)
                                    print(
                                        f"{aemol_a.info['molid']} & {aemol_b.info['molid']} energy diff & geom diff = {energy_diff} kj/mol & {diff} Angstroms", file=l)
                                else:
                                    elim_list.append(a)
                                    print(
                                        f"added mol {aemol_a.info['molid']} to eliminated_mol.txt due to geomtric similarity to {aemol_b.info['molid']} as mirror image has been found", file=l)
                        else:
                            print(
                                f"geometry threshold is passed but not energy threshold, consider changing parameters after checking {aemol_a.info['molid']} & {aemol_b.info['molid']}", file=l)
                            print(
                                f"{aemol_a.info['molid']} & {aemol_b.info['molid']} energy diff & geom diff = {energy_diff} kj/mol & {diff} Angstroms", file=l)

    elim_list = list(set(elim_list))
    removed_mol = []
    for i in range(len(elim_list)):
        if elim_list[i] == 0:
            removed_mol.append(aemols[i].info['molid'])

    with open("eliminated_molecules.txt", 'w')as f:
        for id in removed_mol:
            f.write(id + "\n")
