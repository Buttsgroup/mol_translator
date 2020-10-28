import sys
import os
import pandas as pd
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

import glob

import numpy as np

import mol_translator.imp_converter
from mol_translator.imp_converter import dataframe_read, dataframe_write, dataframe_prep

import mol_translator.aemol as aemol

'''
    Write tests for every function in the imp_converter module
'''


def test_dataframe_read():
    atom_df = pd.read_pickle("tests/test_dataset/test_df_atom.pkl")
    pair_df = pd.read_pickle("tests/test_dataset/test_df_pair.pkl")
    test_mols = dataframe_read.read_df(atom_df, pair_df)
    assert test_mols is not None

def test_dataframe_readwrite():
    atom_df = pd.read_pickle("tests/test_dataset/test_df_atom.pkl")
    pair_df = pd.read_pickle("tests/test_dataset/test_df_pair.pkl")
    test_mols = dataframe_read.read_df(atom_df, pair_df)
    out_atom = dataframe_write.make_atom_df(test_mols)
    out_pair = dataframe_write.make_pair_df(test_mols)

    assert atom_df.loc[0]['x'] == out_atom.loc[0]['x']
    assert pair_df.loc[0]['coupling'] == out_pair.loc[0]['coupling']

def test_prep_mol():
    test_mol = aemol(0)
    test_mol.from_file('tests/test_mols/qm9_8.nmredata.sdf', ftype='sdf')
    test_mol = dataframe_prep.prep_mol(test_mol, nmr_file='tests/test_mols/qm9_8.nmredata.sdf', nmr_type='nmredata')
    
    assert np.array_equal(test_mol.atom_properties['shift'], np.asarray([ 46.06039042, 327.955,   3.76051871,   3.76014554,
         3.59063345,  -0.6359735 ]))
    assert test_mol.pair_properties['nmr_types'] == [['0JCC', '1JOC', '1JCH', '1JCH', '1JCH', '2JCH'], 
                                                    ['1JOC', '0JOO', '2JOH', '2JOH', '2JOH', '1JOH'], 
                                                    ['1JCH', '2JOH', '0JHH', '2JHH', '2JHH', '3JHH'], 
                                                    ['1JCH', '2JOH', '2JHH', '0JHH', '2JHH', '3JHH'], 
                                                    ['1JCH', '2JOH', '2JHH', '2JHH', '0JHH', '3JHH'], 
                                                    ['2JCH', '1JOH', '3JHH', '3JHH', '3JHH', '0JHH']]
    assert np.array_equal(test_mol.pair_properties['coupling'], np.asarray([[  0.     ,   0.     ,  87.6326 ,  87.6253 ,  90.0888 ,  -1.36995],
                                                                           [  0.     ,   0.     ,   0.     ,   0.     ,   0.     ,   0.     ],
                                                                           [ 87.6326 ,   0.     ,   0.     ,  -6.48291, -10.6098 ,   1.43077],
                                                                           [ 87.6253 ,   0.     ,  -6.48291,   0.     , -10.6153 ,   1.45233],
                                                                           [ 90.0888 ,   0.     , -10.6098 , -10.6153 ,   0.     ,  13.7865 ],
                                                                           [ -1.36995,   0.     ,   1.43077,   1.45233,  13.7865 ,   0.     ]]))
    assert np.array_equal(test_mol.structure['xyz'], np.asarray([[ 0.6607, -0.0193,  0.    ],
                                                               [-0.748 ,  0.1217,  0.    ],
                                                               [ 1.035 , -0.5428,  0.8925],
                                                               [ 1.0348, -0.5448, -0.8914],
                                                               [ 1.0819,  0.9897, -0.0012],
                                                               [-1.1318, -0.7597,  0.    ]]))
    assert np.array_equal(test_mol.structure['types'], np.asarray([6, 8, 1, 1, 1, 1]))
    assert np.array_equal(test_mol.structure['conn'], np.asarray([[0, 1, 1, 1, 1, 0],
                                                                   [1, 0, 0, 0, 0, 1],
                                                                   [1, 0, 0, 0, 0, 0],
                                                                   [1, 0, 0, 0, 0, 0],
                                                                   [1, 0, 0, 0, 0, 0],
                                                                   [0, 1, 0, 0, 0, 0]]))
    assert np.array_equal(test_mol.structure['path_len'], np.asarray([[0, 1, 1, 1, 1, 2],
                                                                       [1, 0, 2, 2, 2, 1],
                                                                       [1, 2, 0, 2, 2, 3],
                                                                       [1, 2, 2, 0, 2, 3],
                                                                       [1, 2, 2, 2, 0, 3],
                                                                       [2, 1, 3, 3, 3, 0]]))
    
def test_dataframe_write():
    
    files = glob.glob("tests/test_mols/qm9_*.nmredata.sdf")
    mols = []
    for file in files:
        f = int(file.split('_')[-1].split('.')[0])
        test_mol = aemol(f)
        test_mol.from_file(file, ftype='sdf')
        test_mol = dataframe_prep.prep_mol(test_mol, nmr_file=file, nmr_type='nmredata')
        mols.append(test_mol)
        
    atom_df = dataframe_write.make_atom_df(mols)
    pair_df = dataframe_write.make_pair_df(mols)
    
    assert np.all(pair_df.loc[(pair_df.path_len == 0)]['coupling'].to_numpy() == 0)
    assert np.all([x[2] == x[3] for x in pair_df.loc[(pair_df.path_len == 0)]['nmr_types'].to_list()])
    assert pair_df['path_len'].to_list() == [int(x[0]) for x in pair_df['nmr_types'].to_list()]
    
    mol_atoms = atom_df.loc[atom_df.molecule_name == 8]
    mol_pairs = pair_df.loc[pair_df.molecule_name == 8]
    
    assert np.array_equal(mol_atoms['shift'].to_numpy(), np.asarray([ 46.06039042, 327.955,   3.76051871,   3.76014554,
        3.59063345,  -0.6359735 ])) 
         
    assert mol_pairs['nmr_types'].to_list() == ['0JCC', '1JOC', '1JCH', '1JCH', '1JCH', '2JCH', 
                                                    '1JOC', '0JOO', '2JOH', '2JOH', '2JOH', '1JOH', 
                                                    '1JCH', '2JOH', '0JHH', '2JHH', '2JHH', '3JHH', 
                                                    '1JCH', '2JOH', '2JHH', '0JHH', '2JHH', '3JHH', 
                                                    '1JCH', '2JOH', '2JHH', '2JHH', '0JHH', '3JHH', 
                                                    '2JCH', '1JOH', '3JHH', '3JHH', '3JHH', '0JHH']
                                                    
    assert np.array_equal(mol_pairs['coupling'], np.asarray([[  0.     ,   0.     ,  87.6326 ,  87.6253 ,  90.0888 ,  -1.36995],
                                                                           [  0.     ,   0.     ,   0.     ,   0.     ,   0.     ,   0.     ],
                                                                           [ 87.6326 ,   0.     ,   0.     ,  -6.48291, -10.6098 ,   1.43077],
                                                                           [ 87.6253 ,   0.     ,  -6.48291,   0.     , -10.6153 ,   1.45233],
                                                                           [ 90.0888 ,   0.     , -10.6098 , -10.6153 ,   0.     ,  13.7865 ],
                                                                           [ -1.36995,   0.     ,   1.43077,   1.45233,  13.7865 ,   0.     ]]).flatten())
    assert np.array_equal(np.asarray([mol_atoms['x'], mol_atoms['y'], mol_atoms['z']]).transpose(), np.asarray([[ 0.6607, -0.0193,  0.    ],
                                                               [-0.748 ,  0.1217,  0.    ],
                                                               [ 1.035 , -0.5428,  0.8925],
                                                               [ 1.0348, -0.5448, -0.8914],
                                                               [ 1.0819,  0.9897, -0.0012],
                                                               [-1.1318, -0.7597,  0.    ]]))
    assert np.array_equal(mol_atoms['typeint'].to_numpy(), np.asarray([6, 8, 1, 1, 1, 1]))
    assert np.array_equal([list(x) for x in mol_atoms['conn'].to_numpy()], np.asarray([[0, 1, 1, 1, 1, 0],
                                                                   [1, 0, 0, 0, 0, 1],
                                                                   [1, 0, 0, 0, 0, 0],
                                                                   [1, 0, 0, 0, 0, 0],
                                                                   [1, 0, 0, 0, 0, 0],
                                                                   [0, 1, 0, 0, 0, 0]]))
    assert np.array_equal(mol_pairs['path_len'].to_numpy(), np.asarray([[0, 1, 1, 1, 1, 2],
                                                                       [1, 0, 2, 2, 2, 1],
                                                                       [1, 2, 0, 2, 2, 3],
                                                                       [1, 2, 2, 0, 2, 3],
                                                                       [1, 2, 2, 2, 0, 3],
                                                                       [2, 1, 3, 3, 3, 0]]).flatten())
    