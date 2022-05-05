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

from typing import Type, Optional, List
from mol_translator.aemol import Aemol
import numpy as np


class Aeconf(Aemol):

    def __init__(self, aemols: List[Type[Aemol]], mol_name='conf'):
        self.mol_name = mol_name
        self.num_confs = len(aemols)
        self.aemols = aemols
        self.energy_array = np.zeros(self.num_confs, dtype=np.float64)
        self.eliminated_mols = []
        self.averaged_shift = []
        self.averaged_coupling = []

    def calc_pops(self, temp: int = 298):

        for c, aemol in enumerate(self.aemols):
            self.energy_array[c] = aemol.mol_properties['energy']

        kj_array = self.energy_array * 2625.5
        min_val = np.amin(kj_array)
        rel_array = (kj_array - min_val)
        exp_array = -(rel_array*1000) / float(8.31*temp)
        exp_array = np.exp(exp_array)
        sum_val = np.sum(exp_array)

        pop_array = exp_array / sum_val

        for a, aemol in enumerate(self.aemols):
            aemol.mol_properties['pop'] = pop_array[a]

    def boltzmann_average(self, pair_props: list = ['coupling'], atom_props: list = ['shift']):
        atoms = len(self.aemols[0].structure['types'])

        new_atom_dict = {}
        new_pair_dict = {}
        for prop in atom_props:
            new_atom_dict[prop] = np.zeros(atoms, dtype=np.float64)
        for prop in pair_props:
            new_pair_dict[prop] = np.zeros((atoms, atoms), dtype=np.float64)

        for aemol in self.aemols:
            for prop in atom_props:
                new_atom_dict[prop] += aemol.atom_properties[prop] * \
                    aemol.mol_properties['pop']
            for prop in pair_props:
                new_pair_dict[prop] += aemol.pair_properties[prop] * \
                    aemol.mol_properties['pop']

        self.averaged_shift = new_atom_dict
        self.averaged_coupling = new_pair_dict

    def get_dist_array(self):

        for t, aemol in self.aemols:
            num_atoms = len(aemol.structure['types'])
            dist_array = np.zeros((num_atoms, num_atoms), dtype=np.float64)
            for i in range(num_atoms):
                for j in range(num_atoms):
                    dist_array[i][j] = np.absolute(np.linalg.norm(
                        aemol.structure['xyz'][i] - aemol.structure['xyz'][j]))

            aemol.mol_properties['dist_array'] = dist_array

    def redundant_elimination(self, geom_threshold: float = 0.1, e_threshold: float = 0.1, redundant_atoms: Optional[list] = None, achiral: bool = False):

        if redundant_atoms != None:
            redundant_atoms = redundant_atoms.split(",")
            redundant_atoms = list(map(int, redundant_atoms))

            dist_arrays = np.zeros(len(aemols[0].structure['types']), len(
                aemols[0].structure['types']))
            for atoms in redundant_atoms:
                for i in range(aemols):
                    dist_arrays[i][int(atoms)-1] = 0
                    for k in range(atoms):
                        dist_arrays[i][k][int(atoms)-1] = 0

        with open(f"{self.mol_name}_redundant_elimination_log.txt", 'a') as l:

            for a, aemol_a in enumerate(self.aemols):
                self.energy_array[a] = aemol_a.mol_properties['energy']
                dist_array_a = aemol_a.mol_properties['dist_array']

                for b, aemol_b in enumerate(self.aemols):
                    self.energy_array[b] = aemol_b.mol_properties['energy']
                    dist_array_b = aemol_b.mol_properties['dist_array']
                    if a > b and not b in self.eliminated_mols:
                        diff = np.absolute(np.sum(dist_array_a - dist_array_b))
                        energy_diff = np.absolute(
                            self.energy_array[a] - self.energy_array[b]) * 2625.5

                        if diff < geom_threshold:
                            if energy_diff < e_threshold:
                                if achiral is False:
                                    self.eliminated_mols.append(a)
                                    print(
                                        f"added mol {aemol_a.info['molid']} to eliminated_mol.txt due to geomtric similarity to {aemol_b.info['molid']}", file=l)

                                else:
                                    if a - b == 1:
                                        print(
                                            f"energy difference between {aemol_a.info['molid']} & {aemol_b.info['molid']} detected but could be mirror images, please check manually", file=l)
                                        print(
                                            f"{aemol_a.info['molid']} & {aemol_b.info['molid']} energy diff & geom diff = {energy_diff} kj/mol & {diff} Angstroms", file=l)
                                    else:
                                        self.eliminated_mols.append(a)
                                        print(
                                            f"added mol {aemol_a.info['molid']} to eliminated_mol.txt due to geomtric similarity to {aemol_b.info['molid']} as mirror image has been found", file=l)
                            else:
                                print(
                                    f"geometry threshold is passed but not energy threshold, consider changing parameters after checking {aemol_a.info['molid']} & {aemol_b.info['molid']}", file=l)
                                print(
                                    f"{aemol_a.info['molid']} & {aemol_b.info['molid']} energy diff & geom diff = {energy_diff} kj/mol & {diff} Angstroms", file=l)

        elim_list = list(set(self.eliminated_mols))
        removed_mol = []
        for mol in elim_list:
            removed_mol.append(mol)

        with open(f"{self.mol_name}_eliminated_molecules.txt", 'w')as f:
            for id in removed_mol:
                f.write(id + "\n")
