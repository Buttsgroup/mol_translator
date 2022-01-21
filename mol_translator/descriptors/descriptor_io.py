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

from mol_translator.descriptors.atom_descriptors import rdkit_atom_descriptors, pybel_atom_descriptors
from mol_translator.descriptors.pair_descriptors import rdkit_pair_descriptors, pybel_pair_descriptors

def get_all_descriptors(aemol):
    get_pyb_descriptors(aemol)
    get_rdkit_descriptors(aemol)

def get_pyb_descriptors(aemol):
    if aemol.pybmol is None:
        pybmol = aemol.to_pybel()
    else:
        pybmol = aemol.pybmol

    atom_props = pybel_atom_descriptors.get_pyb_atom_descriptors(pybmol)
    aemol.atom_properties.update(atom_props)

    pair_props = pybel_pair_descriptors.get_pyb_pair_descriptors(pybmol)
    aemol.pair_properties.update(pair_props)

def get_rdkit_descriptors(aemol):
    if aemol.rdmol is None:
        rdmol = aemol.to_rdkit()
    else:
        rdmol = aemol.rdmol

    atom_props = rdkit_atom_descriptors.get_rd_atom_descriptors(rdmol)
    aemol.atom_properties.update(atom_props)

    #pair_props = rdkit_pair_descriptors.get_rd_pair_descriptors(rdmol)
    #aemol.pair_properties.update(pair_props)
