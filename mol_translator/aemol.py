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

import pybel as pyb

from .pybel_converter.pybel_converter import pybel_to_aemol, aemol_to_pybel
from .file_creation.structure_write import write_mol_toxyz

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
                            'conn': [],
                            'cpl_len': []}

        self.atom_properties = {'shift': []}

        self.bond_properties = {'coupling': []}

        self.mol_properties = {'energy': -404.404}


    def from_pybel(self, pybmol):
        types, xyz, conn = pybel_to_aemol(pybmol)

        self.structure['types'] = types
        self.structure['xyz'] = xyz
        self.structure['conn'] = conn

    def to_pybel(self):
        write_mol_toxyz(self.structure, 'tmp.xyz')
        pybmol = next(pyb.readfile('xyz', 'tmp.xyz'))
        os.remove('tmp.xyz')

        return pybmol

    def from_rdkit(self, rdmol):
        '''
            Need to write this
        '''

    def to_rdkit(self):
        '''
            Need to write this
        '''
        return rdmol

    def from_file()

    def to_file_ae(self, format, filename):
        pybmol = self.to_pybel()
        pybmol.write(format, filename)

    def to_file_pyb(self, format, filename):
        pybmol = self.to_pybel()
        pybmol.write(format, filename)
