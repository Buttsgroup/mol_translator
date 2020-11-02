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

import pandas as pd

# Functions to read binding properties from files

# Reads pchembl value for given molecule from tsv file downloadable from ChEMBL
def pchembl_read(file, format='tsv', molid=0):
    ic50_df = pd.read_csv(file, sep='\t')
    ic50 = ic50_df.loc[(ic50_df["Molecule ChEMBL ID"] == molid)]["pChEMBL Value"].to_numpy()
    ic50 = ic50[0]
    
    return ic50
    