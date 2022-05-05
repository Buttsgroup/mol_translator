import pandas as pd

import random


def df_batch_merge(folderpath: str, outfolder: str, dataset_name: str, num_batches: int = 9):
    atom_dfs = []
    pair_dfs = []

    for i in range(1, num_batches):
        atom_fp = f'{folderpath}/{dataset_name}/{dataset_name}_batch{i}_atoms_df.pkl'
        pair_fp = f'{folderpath}/{dataset_name}/{dataset_name}_batch{i}_pairs_df.pkl'

        atom_df = pd.read_pickle(atom_fp)
        pair_df = pd.read_pickle(pair_fp)

        atom_dfs.append(atom_df)
        pair_dfs.append(pair_df)

    merged_atom_df = pd.concat(atom_dfs)
    merged_pair_df = pd.concat(pair_dfs)

    merged_atom_df.to_pickle(f'{outfolder}/{dataset_name}_all_atoms_df.pkl')
    merged_pair_df.to_pickle(f'{outfolder}/{dataset_name}_all_pairs_df.pkl')


def dataset_merge(df_atom_1: str, df_pair_1: str, df_atom_2: str, df_pair_2: str, outfolder: str, dataset_name: str, shuffle: bool = True):
    df_atom_1 = pd.read_pickle(df_atom_1)
    df_pair_1 = pd.read_pickle(df_pair_1)
    df_atom_2 = pd.read_pickle(df_atom_2)
    df_pair_2 = pd.read_pickle(df_pair_2)

    atoms_df = [df_atom_1, df_atom_2]
    pairs_df = [df_pair_1, df_pair_2]

    merged_atoms_df = pd.concat(atoms_df)
    merged_pairs_df = pd.concat(pairs_df)

    if shuffle:
        mol_id = merged_atoms_df['molecule_name'].unique()
        random.shuffle(mol_id)
        merged_atoms_df = merged_atoms_df.set_index(
            'molecule_name').loc[mol_id].reset_index()
        merged_pairs_df = merged_pairs_df.set_index(
            'molecule_name').loc[mol_id].reset_index()

    atoms_outfile = f'{outfolder}/{dataset_name}_atoms_df.pkl'
    pairs_outfile = f'{outfolder}/{dataset_name}_pairs_df.pkl'

    merged_atoms_df.to_pickle(atoms_outfile)
    merged_pairs_df.to_pickle(pairs_outfile)
