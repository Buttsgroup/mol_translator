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

from mol_translator.aemol import aemol

import numpy as np

def read_df(atom_df, pair_df):
    
    
    molnames = atom_df.molecule_name.unique()
    
    mols = []
    for molname in molnames:
        amol = aemol(molname)
        
        ## Atom stuff
        mol_atom_df = atom_df.loc[(atom_df.molecule_name==molname)]
        atoms = len(mol_atom_df)
        xyz = np.zeros((atoms, 3), dtype=np.float64)
        types = np.zeros((atoms), dtype=np.float64)
        conn = np.zeros((atoms, atoms), dtype=np.float64)
        for idx in mol_atom_df['atom_index']:
            idx_atom_df = mol_atom_df.loc[(mol_atom_df.atom_index == idx)]
            xyz[idx][0] = idx_atom_df['x'].to_numpy().squeeze()   
            xyz[idx][1] = idx_atom_df['y'].to_numpy().squeeze()   
            xyz[idx][2] = idx_atom_df['z'].to_numpy().squeeze()   
            types[idx] = idx_atom_df['typeint'].to_numpy().squeeze()
            conn[idx] = idx_atom_df['conn'].to_numpy().squeeze()
        amol.structure['xyz'] = xyz
        amol.structure['types'] = types
        amol.structure['conn'] = conn
        
        propnames = []
        for prop in mol_atom_df.keys():
            if not prop in ['molecule_name', 'atom_index', 'typestr', 
                    'typeint', 'x', 'y', 'z', 'conn']:
                propnames.append(prop)
        
        for prop in propnames:
            amol.atom_properties[prop] = mol_atom_df[prop].to_numpy()
            
        ## Pair stuff
        mol_pair_df = pair_df.loc[(pair_df.molecule_name==molname)]
        
        propnames = []
        for prop in mol_pair_df.keys():
            if not prop in ['molecule_name', 'atom_index_0', 
                    'atom_index_1', 'dist']:
                propnames.append(prop)
        
        for prop in propnames:
            amol.pair_properties[prop] = mol_pair_df[prop].to_numpy()
        
            if prop == 'bond_dist':
                amol.pair_properties['path_len'] = amol.pair_properties['bond_dist']
                del amol.pair_properties['bond_dist']

        mols.append(amol)
        
    return mols
        







        