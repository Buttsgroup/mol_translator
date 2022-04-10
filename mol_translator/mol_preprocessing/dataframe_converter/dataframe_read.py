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

from tqdm import tqdm


def read_df(atom_df, pair_df):

    if 'bond_dist' in pair_df.keys():
        pair_df['path_len'] = pair_df['bond_dist']

    molnames = atom_df.molecule_name.unique()

    mols = []
    for molname in tqdm(molnames):
        amol = aemol(molname)

        mol_atom_df = atom_df.loc[(atom_df.molecule_name == molname)]
        mol_pair_df = pair_df.loc[(pair_df.molecule_name == molname)]

        atoms = len(mol_atom_df)
        xyz = np.zeros((atoms, 3), dtype=np.float64)
        types = np.zeros((atoms), dtype=np.int32)
        conn = []
        path_len = []
        for idx in mol_atom_df['atom_index']:
            idx_atom_df = mol_atom_df.loc[(mol_atom_df.atom_index == idx)]
            xyz[idx][0] = idx_atom_df['x'].to_numpy().squeeze()
            xyz[idx][1] = idx_atom_df['y'].to_numpy().squeeze()
            xyz[idx][2] = idx_atom_df['z'].to_numpy().squeeze()
            types[idx] = idx_atom_df['typeint'].to_numpy().squeeze()
            conn.append(*idx_atom_df['conn'].to_list())
            path = []
            for jdx in mol_atom_df['atom_index']:
                plen = mol_pair_df.loc[
                    (mol_pair_df.atom_index_0 == idx)
                    & (mol_pair_df.atom_index_1 == jdx)]['path_len'].to_list()
                if len(plen) == 1:
                    path.append(plen[0])
                else:
                    path.append(0)
            path_len.append([*path])

        amol.structure['xyz'] = xyz
        amol.structure['types'] = types
        amol.structure['conn'] = conn
        amol.structure['path_len'] = np.asarray(path_len)

        propnames = []
        for prop in mol_atom_df.keys():
            if not prop in ['molecule_name', 'atom_index', 'typestr',
                            'typeint', 'x', 'y', 'z', 'conn']:
                propnames.append(prop)
        for prop in propnames:
            amol.atom_properties[prop] = mol_atom_df[prop].to_numpy()

        propnames = []
        for prop in mol_pair_df.keys():
            if not prop in ['molecule_name', 'atom_index_0',
                            'atom_index_1', 'dist', 'path_len', 'bond_dist']:
                propnames.append(prop)
                amol.pair_properties[prop] = []

        for prop in propnames:
            for idx in mol_atom_df['atom_index']:
                prop_list = []
                for jdx in mol_atom_df['atom_index']:
                    df_prop = mol_pair_df.loc[
                        (mol_pair_df.atom_index_0 == idx)
                        & (mol_pair_df.atom_index_1 == jdx)][prop].to_list()[0]
                    prop_list.append(df_prop)
                amol.pair_properties[prop].append(prop_list)

        mols.append(amol)

    return mols
