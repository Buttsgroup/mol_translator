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

import numpy as np
from tqdm import tqdm

def get_dist_array(aemol):
    num_atoms = len(aemol.structure['types'])
    dist_array = np.zeros(num_atoms, num_atoms, dtype=np.float64)
    for i in range(num_atoms):
        for j in range(num_atoms):
            d_array[i][j] = np.absolute(np.linalg.norm(xyz_array[i] - xyz_array[j]))

    return d_array
