# Copyright 2022 Will Gerrard, Calvin Yiu
# This file is part of autoenrich.

# autoenrich is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# autoenrich is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.


# Functions to manage write/read of properties
from .nmr import nmr_ops, nmr_read, nmr_write
from .energy import energy_read
from .binding import binding_read
from .charge import charge_read
from ..util.custom_errors import raise_formaterror
from typing import Type


def prop_write(aemol: Type, outfile: str, prop: str = 'nmr', format: str = 'nmredata') -> None:
    if prop == 'nmr':
        if 'nmr_type' not in aemol.pair_properties.keys():
            aemol.pair_properties['nmr_type'] = nmr_ops.get_coupling_types(
                aemol)

        if format == 'nmredata':
            nmr_write.write_nmredata(
                outfile, aemol, write_zeros=False, count_from=0)
        else:
            raise_formaterror(format)
    else:
        print('property not recognised or function not written yet !')
        raise ValueError(f'Cannot output property: {prop}')


def prop_read(aemol: Type, filename: str, prop: str, format: str) -> None:
    if prop == 'nmr':
        shift, coupling = nmr_read.nmr_read(filename, prop, format)
        aemol.atom_properties['shift'] = shift
        aemol.pair_properties['coupling'] = coupling
    elif prop == 'nmr_var':
        shift_var, coupling_var = nmr_read.nmr_read(filename, prop, format)
        aemol.atom_properties['shift_var'] = shift_var
        aemol.pair_properties['coupling_var'] = coupling_var
    elif prop == 'scf':
        scf = energy_read.energy_read(filename, prop, format)
        aemol.mol_properties['energy'] = scf
    elif prop == 'ic50':
        ic50 = binding_read.pchembl_read(
            filename, format, aemol.info['molid'])
        aemol.mol_properties['ic50'] = ic50
    elif prop == 'mc':
        mc = charge_read(filename, prop, format)
        aemol.atom_properties['mull_chg'] = mc
    else:
        print('property not recognised or function not written yet !')
        raise ValueError(f'Cannot read property: {prop}')
