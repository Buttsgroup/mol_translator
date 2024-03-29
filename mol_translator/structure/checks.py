# Copyright 2020 Will Gerrard, Calvin Yiu
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
import mol_translator.aemol as aemol
from mol_translator.structure.find_num_bonds import rdmol_find_num_bonds
from rdkit.Chem import rdMolTransforms as Chem
import numpy as np

def run_all_checks(aemol, post_check=False):

    if not check_valence(aemol):
        return False
    if not check_missing_H(aemol, post_check=post_check):
        return False
    if not check_bonds_are_plausible_and_atom_overlap(aemol):
        return False
    return True

def check_valence(aemol):
    '''
    Might be redundant since loading molecule into rdkit already checks this
    '''
    aemol.to_rdkit()
    num_bond_dict = rdmol_find_num_bonds(aemol.rdmol)
    for num_bonds in num_bond_dict.values():
        if len(num_bonds)>1:
            print(f"Unexpected number of valence electrons in mol, {aemol.info['molid']}")
            return False
    return True

def check_missing_H(aemol, post_check=False):
    '''
    everything reasonable should have hydrogens in it.
    Post check finds missing hydrogens
    '''
    atom_num_array = aemol.structure['types']
    atom_num_list = list(atom_num_array)
    if 1 not in atom_num_list:
        print(f"No Hs in mol, {aemol.info['molid']}")
        return False

    if post_check:
        aemol.to_rdkit()
        for atom in aemol.rdmol.GetAtoms():
            if atom.GetNumImplicitHs() != 0 :
                print(f"Hs Missing on mol {aemol.info['molid']}")
                return False

    return True

def check_bonds_are_plausible_and_atom_overlap(aemol):

    aemol.to_rdkit()
    rdmol = aemol.rdmol
    for i_idx, i_atom in enumerate(rdmol.GetAtoms()):
        for j_idx, j_atom in enumerate(rdmol.GetAtoms()):
            if i_idx == j_idx:
                continue
            bond = rdmol.GetBondBetweenAtoms(i_idx, j_idx)
            if bond is not None:
                conf = rdmol.GetConformer()
                bond_length = Chem.GetBondLength(conf, i_idx, j_idx)
                if bond_length >= 0.90 and bond_length <= 2.60:
                    continue
                else:
                    return False
            else:
                i_atom_pos = np.array([rdmol.GetConformer().GetAtomPosition(i_idx).x,
                                       rdmol.GetConformer().GetAtomPosition(i_idx).y,
                                       rdmol.GetConformer().GetAtomPosition(i_idx).z])

                j_atom_pos = np.array([rdmol.GetConformer().GetAtomPosition(j_idx).x,
                                       rdmol.GetConformer().GetAtomPosition(j_idx).y,
                                       rdmol.GetConformer().GetAtomPosition(j_idx).z])
                average_dist = np.linalg.norm(i_atom_pos - j_atom_pos)
                if average_dist < 1.4:
                    return False

    return True

    #angstrom check for possible bond lengths, likely with dictionary values, account for bond order
    # is every bond within a plausible distance range
