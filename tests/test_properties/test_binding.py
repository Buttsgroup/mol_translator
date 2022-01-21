import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

import mol_translator.aemol as aemol
from mol_translator.imp_converter import dataframe_read, dataframe_write, dataframe_prep
import glob

import numpy as np

def test_getpchembl():
    files = glob.glob('tests/test_mols/test_pchembl/*.sdf')
    mols = []
    for file in files:
        molid = file.split('/')[-1].split('.')[0]
        test_mol = aemol(molid)
        test_mol.from_file_pyb(file, ftype='sdf')
        test_mol.from_pybel(test_mol.pybmol)
        test_mol.prop_from_file('tests/test_mols/test_pchembl/test_ic50.tsv', prop='ic50', format='tsv')

        mols.append(test_mol)

    assert mols[0].info['molid'] == "CHEMBL460962"
    assert mols[0].mol_properties['ic50'] == 8.92
    assert mols[1].info['molid'] == "CHEMBL2179452"
    assert mols[1].mol_properties['ic50'] == 7.55

'''
def test_pchemblimp():

    files = glob.glob('tests/test_mols/test_pchembl/*.sdf')
    mols = []
    for file in files:
        molid = file.split('/')[-1].split('.')[0]
        test_mol = aemol(molid)
        test_mol.from_file_pyb(file, ftype='sdf')
        test_mol.from_pybel(test_mol.pybmol)


        pmol = dataframe_prep.prep_mol_ic50(test_mol, 'tests/test_mols/test_pchembl/test_ic50.tsv', 'ic50')

        mols.append(pmol)

    atom_df = dataframe_write.make_atom_df(mols)
    pair_df = dataframe_write.make_pair_df(mols)

    assert np.array_equal(atom_df.loc[(atom_df['molecule_name'] == "CHEMBL460962")]['ic50'].to_numpy(),
                                np.full(len(atom_df.loc[(atom_df['molecule_name'] == "CHEMBL460962")]['typeint'].to_numpy()), 8.92, dtype=np.float64))
    assert np.array_equal(atom_df.loc[(atom_df['molecule_name'] == "CHEMBL2179015")]['ic50'].to_numpy(),
                                np.full(len(atom_df.loc[(atom_df['molecule_name'] == "CHEMBL2179015")]['typeint'].to_numpy()), 7.54, dtype=np.float64))
'''
