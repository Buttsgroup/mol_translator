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

from .periodic_table import Get_periodic_table

def flag_to_target(flag):

	p_table = Get_periodic_table()
	if len(flag) == 4:
		if str(flag[0]) == 'J':
			length = int(flag[1])
		else:
			length = int(flag[0])
		atype1 = int(p_table.index(str(flag[2])))
		atype2 = int(p_table.index(str(flag[3])))

		if atype1 >= atype2:
			return [length, atype1, atype2]
		else:
			return [length, atype2, atype1]

	elif len(flag) == 3:
		atype = int(p_table.index(str(flag[0])))
		return [atype]

	else:
		raise TypeError(f'flag, {flag} not recognised, coupling flag format is <nJxy> . . .')
		raise TypeError(f'flag, {flag} not recognised, chemical shift flag format is <XCS> . . .')

def target_to_flag(target):
	p_table = Get_periodic_table()

	if len(target) == 3:
		if target[1] >= target[2]:
			flag = str(target[0]) + 'J' + str(p_table[target[1]]) + str(p_table[target[2]])
		else:
			flag = str(target[0]) + 'J' + str(p_table[target[2]]) + str(p_table[target[1]])
	elif len(target) == 1:
		flag = str(target[0]) + 'CS'
	else:
		raise TypeError(f'Error, target {target} not recognised')

	return flag
