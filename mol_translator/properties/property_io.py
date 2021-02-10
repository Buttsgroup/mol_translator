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


# Functions to manage write/read of properties
from . import nmr as nmr
from . import energy as energy
from . import binding as binding

def prop_write(aemol, outfile, prop, format):
    if prop == 'nmr':
        if 'nmr_types' not in aemol.pair_properties.keys():
            aemol.pair_properties['nmr_types'] = nmr.nmr_ops.get_coupling_types(aemol)


        if format == 'nmredata':
            nmr.nmr_write.write_nmredata(outfile, aemol, write_zeros=False, count_from=0)
        else:
            raise_formaterror(format)
    else:
        print('property not recognised or function not written yet !')
        raise ValueError('Cannot output property: ', prop)


def prop_read(aemol, filename, prop, format):
    if prop == 'nmr':
        shift, coupling = nmr.nmr_read.nmr_read(filename, format, prop)
        aemol.atom_properties['shift'] = shift
        aemol.pair_properties['coupling'] = coupling
    elif prop == 'nmr_var':
        shift_var, coupling_var = nmr.nmr_read.nmr_read(filename, format, prop)
        aemol.atom_properties['shift_var'] = shift_var
        aemol.pair_properties['coupling_var'] = coupling_var
    elif prop == 'scf':
        scf = energy.energy_read.energy_read(filename, format, prop)
        aemol.mol_properties['energy'] = scf
    elif prop == 'ic50':
        ic50 = binding.binding_read.pchembl_read(filename, format, aemol.info['molid'])
        aemol.mol_properties['ic50'] = ic50
    else:
        print('property not recognised or function not written yet !')
        raise ValueError('Cannot read property: ', prop)
