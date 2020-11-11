# Copyright 2020 Will Gerrard, Calvin Yiu
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

import pybel as pyb

from mol_translator.structure.pybel_converter import pybmol_to_aemol, aemol_to_pybmol
from mol_translator.structure.rdkit_converter import rdmol_to_aemol, aemol_to_rdmol
from mol_translator.structure import structure_write as strucwrt

import mol_translator.properties.property_io as prop_io

from mol_translator.structure import find_paths as pathfind

from mol_translator.properties.nmr.nmr_write import write_nmredata

from rdkit.Chem import AllChem as Chem

class aemol(object):
    """
        molecule object

        molecular information is stored in a set of dictionaries containing strings and numpy arrays

    """

    def __init__(self, molid, filepath=""):

        self.info = {   'molid': molid,
                        'filepath': filepath}

        self.structure = {  'xyz': [],
                            'types': [],
                            'conn': []}

        self.atom_properties = {}

        self.pair_properties = {}

        self.mol_properties = {'energy': -404.404}


    def from_pybel(self, pybmol):
        types, xyz, conn = pybmol_to_aemol(pybmol)
        self.structure['types'] = types
        self.structure['xyz'] = xyz
        self.structure['conn'] = conn

    def to_pybel(self):
        pybmol = aemol_to_pybmol(self.structure)

        return pybmol

    def from_rdkit(self, rdmol):
        #assumes rdmol is already 3D with Hs included
        types, xyz = rdmol_to_aemol(rdmol)
        self.structure['types'] = types
        self.structure['xyz'] = xyz

    def to_rdkit(self, removeHs=False):
        rdmol = aemol_to_rdmol(self, removeHs)

        return rdmol

    def from_file(self, file, ftype='xyz'):
        pybmol = next(pyb.readfile(ftype, file))
        self.from_pybel(pybmol)

    def from_string(self, string, stype='smi'):
        pybmol = pyb.readstring(stype, string)
        self.from_pybel(pybmol)

    def to_file_ae(self, format, filename):
        if format == 'xyz':
            strucwrt.write_mol_toxyz(self.structure, filename)
        else:
            raise_formaterror(format)

    def to_file_pyb(self, format, filename):
        pybmol = self.to_pybel()
        pybmol.write(format, filename)

    def prop_tofile(self, filename, prop='nmr', format='nmredata'):
        prop_io.prop_write(self, filename, prop, format)

    def prop_fromfile(self, filename, ftype, prop):
        prop_io.prop_read(self, filename, prop, ftype)

    def get_all_paths(self, maxlen=5):
        pybmol = self.to_pybel()
        self.structure['paths'] = pathfind.pybmol_find_all_paths(pybmol, maxlen)

    def get_bonds(self):
        pybmol = self.to_pybel()
        self.structure['conn'] = pathfind.pybmol_get_bond_table(pybmol)

    def get_path_lengths(self, maxlen=5):
        pybmol = self.to_pybel()
        self.structure['path_len'] = pathfind.pybmol_get_path_lengths(pybmol, maxlen)

    def get_pyb_fingerprint(self, fingerprint):
        pybmol = self.to_pybel()
        self.mol_properties[fingerprint] = pybmol.calcfp(fingerprint)
        # available fingerprints: ['ecfp0', 'ecfp10', 'ecfp2', 'ecfp4', 'ecfp6', 'ecfp8', 'fp2', 'fp3', 'fp4', 'maccs']

    def get_rdkit_fingerprint(self, radius=2, nBits=2048):
        rdmol = self.to_rdkit()
        fp = Chem.GetMorganFingerprintAsBitVect(rdmol,radius=radius, nBits=nBits)
        self.mol_properties['ecfp4'] = fp
