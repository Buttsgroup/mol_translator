import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

import mol_translator.aemol as aemol



def test_infer_bondangle():
    test_mol = aemol(0)
    test_mol.from_file('tests/test_mols/qm9_9.nmredata.sdf', ftype='sdf')
    
    test_mol.prop_infer('bond_angle')
    
    assert test_mol.trip_properties['bond_angle'].shape == (test_mol.structure['size'], test_mol.structure['size'], test_mol.structure['size'])

    
def test_infer_dihedralangle():
    test_mol = aemol(0)
    test_mol.from_file('tests/test_mols/qm9_9.nmredata.sdf', ftype='sdf')
    
    test_mol.prop_infer('dihedral_angle')
    
    assert test_mol.quad_properties['dihedral_angle'].shape == (test_mol.structure['size'], test_mol.structure['size'], test_mol.structure['size'], test_mol.structure['size'])

    
    
    
    