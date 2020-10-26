import sys
import os
import pandas as pd
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

import mol_translator.imp_converter
from mol_translator.imp_converter import dataframe_read, dataframe_write


'''
    Write tests for every function in the imp_converter module
'''
atom_df = pd.read_pickle("./test_dataset/test_df_atom.pkl")
pair_df = pd.read_pickle("./test_dataset/test_df_pair.pkl")

def test_dataframe_read():
    test_mols = dataframe_read.read_df(atom_df, pair_df)
    assert test_mols is not None

def test_dataframe_write():
    test_mols = dataframe_read.read_df(atom_df, pair_df)
    out_atom = dataframe_write.make_atom_df(test_mols)
    out_pair = dataframe_write.make_pair_df(test_mols)

    assert atom_df.loc[0]['x'] == out_atom.loc[0]['x']
    assert pair_df.loc[0]['coupling'] == out_pair.loc[0]['coupling']
