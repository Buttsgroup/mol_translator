# Copyright 2022 Will Gerrard, Calvin Yiu
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
from openbabel import pybel as pyb
from rdkit.Chem import AllChem
from rdkit import Chem

def rdkit_neutralise(rdmol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    match = rdmol.GetSubstructMatches(pattern)
    match_list = [y[0] for y in match]
    if len(match_list) > 0:
        for idx in match_list:
            atom = rdmol.GetAtomWithIdx(idx)
            charge = atom.GetFormalCharge()
            num_h= atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(num_h - charge)
            atom.UpdatePropertyCache()
    return rdmol

def openbabel_neutralise(obmol):
    pattern = pyb.Smarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")

    match_list = pattern.findall(obmol)
    for match in match_list:
        atom = obmol.GetAtom(match[0])
        charge = atom.GetFormalCharge()
        num_h = atom.GetImplicitHCount()
        atom.SetFormalCharge(0)
        atom.SetImplicitHCount(num_h - charge)
    return obmol
