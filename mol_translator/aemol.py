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

from .pybel_converter.pybel_converter import pybmol_to_aemol, aemol_to_pybmol
from .file_creation import structure_write as strucwrt

from .property_write import nmr_write as nmrwrit
from .property_read import nmr_read as nmrread

from .mol_operations import find_paths as pathfind

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

        self.pair_properties = {'coupling': []}

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
        '''
            Need to write this
        '''

    def to_rdkit(self):
        '''
            Need to write this
        '''
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

        if prop == 'nmr':
            if format == 'nmredata':
                nmrwrit.write_nmredata(outfile, label,
                                            self.structure,
                                            shift=self.atom_properties['shift'],
                                            coupling=self.pair_properties['coupling'],
                                            shift_var=self.atom_properties['shift_var'],
                                            coupling_var=self.atom_properties['coupling_var'])
            else:
                raise_formaterror(format)
        else:
            print('property not recognised or function not written yet !')
            raise ValueError('Cannot output property: ', prop)

    def prop_fromfile(self, filename, ftype, prop):

        if prop == 'nmr':
            if format == 'g09':
                shift, coupling = nmrread.g09_nmrread(filename)
            elif format == 'g16':
                shift, coupling = nmrread.g16_nmrread(filename)
            elif format == 'orca':
                shift, coupling = nmrread.orca_nmrread(filename)
            elif format == 'nmredata':
                shift, _, coupling, _ = nmrread.nmredata_nmrread(filename)
            else:
                raise_formaterror(format)

            self.atom_properties['shift'] = shift
            self.pair_properties['coupling'] = coupling

        elif prop == 'nmr_var':
            if format == 'nmredata':
                _, shift_var, _, coupling_var = nmrread.nmredata_nmrread(filename)
            else:
                raise_formaterror(format)

            self.atom_properties['shift_var'] = shift_var
            self.pair_properties['coupling_var'] = coupling_var

    def get_all_paths(self, maxlen=5):
        pybmol = self.to_pybel()
        self.structure['paths'] = pathfind.pybmol_find_all_paths(pybmol, maxlen)

    def get_bonds(self):
        pybmol = self.to_pybel()
        self.structure['conn'] = pathfind.pybmol_get_bond_table(pybmol)

    def get_cpl_lengths(self):
        pybmol = self.to_pybel()
        self.structure['cpl_len'] = pathfine.pybmol_get_coupling_lengths()

    def raise_formaterror(format):
        print('Format not recognised or function not written yet !')
        print('you can use to_file_pyb to create/read structure files recognised by pybel')
        raise ValueError('Cannot read/write format:', format)
