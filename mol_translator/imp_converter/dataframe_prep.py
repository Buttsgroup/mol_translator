# Copyright 2020 Will Gerrard
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

from mol_translator.properties.nmr.nmr_ops import get_coupling_types
import os.path
import numpy as np

# Gets the necessary info for aemol to turn it into an IMP dataframe
def prep_mol(aemol):
    
    aemol.get_bonds()
    aemol.get_path_lengths()
    get_coupling_types(aemol)

    return aemol

def prep_mol_nmr(aemol, nmr_file="", nmr_type=""):
    
    aemol.get_bonds()
    aemol.get_path_lengths()
    get_coupling_types(aemol)
    if os.path.isfile(nmr_file):
        aemol.prop_fromfile(nmr_file, nmr_type, 'nmr')
    else:
        aemol.atom_properties['shift'] = np.zeros(len(aemol.structure['types']), dtype=np.float64)
        aemol.pair_properties['coupling'] = np.zeros((len(aemol.structure['types']), len(aemol.structure['types'])), dtype=np.float64)
    return aemol

def prep_mol_ic50(aemol, ic50_file="", ic50_type=""):
    aemol.get_bonds()
    aemol.get_path_lengths()
    get_coupling_types(aemol)
    if os.path.isfile(ic50_file):
        aemol.prop_fromfile(ic50_file, ic50_type, 'ic50')
        aemol.atom_properties['ic50'] = np.full(len(aemol.structure['types']), aemol.mol_properties['ic50'], dtype=np.float64)
    else:
        aemol.atom_properties['ic50'] = np.zeros(len(aemol.structure['types']), dtype=np.float64)
    
    return aemol
