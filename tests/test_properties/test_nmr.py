from mol_translator.aemol import Aemol
import mol_translator.properties.nmr as nmr
import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')


'''
    Write tests for every function in the nmr module
'''


def test_getcouplingtypes():
    test_mol = Aemol(0)
    test_mol.from_file_ob('tests/test_mols/qm9_9.nmredata.sdf', ftype='sdf')
    test_mol.from_ob(test_mol.obmol)
    nmr.nmr_ops.get_coupling_types(test_mol)

    assert test_mol.pair_properties['nmr_types'] == [['0JCC', '1JCC', '2JCC', '1JCH', '1JCH', '1JCH', '3JCH'],
                                                     ['1JCC', '0JCC', '1JCC', '2JCH',
                                                      '2JCH', '2JCH', '2JCH'],
                                                     ['2JCC', '1JCC', '0JCC', '3JCH',
                                                      '3JCH', '3JCH', '1JCH'],
                                                     ['1JCH', '2JCH', '3JCH', '0JHH',
                                                      '2JHH', '2JHH', '4JHH'],
                                                     ['1JCH', '2JCH', '3JCH', '2JHH',
                                                      '0JHH', '2JHH', '4JHH'],
                                                     ['1JCH', '2JCH', '3JCH', '2JHH',
                                                      '2JHH', '0JHH', '4JHH'],
                                                     ['3JCH', '2JCH', '1JCH', '4JHH', '4JHH', '4JHH', '0JHH']]
