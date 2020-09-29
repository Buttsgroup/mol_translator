import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

import mol_translator.aemol as aemol


def test_basic():

    smiles = "CCCC"

    mol = aemol('testmol')
    mol.from_string(smiles)

    print(mol.structure['xyz'])
    assert False
