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

def run_all_checks(aemol):

    check1 = check_valence(aemol)
    check2 = check_missing_H(aemol)
    check3 = check_bonds_are_plausible(aemol)
        
    check = check1 and check2 and check3
        
    return check


def check_valence(aemol):
    check = True
    # check valences are roughly correct
    
    return check
    
def check_missing_H(aemol):
    check = True
    
    # everything should have hydrogens in it
    
    return check
    
def check_bonds_are_plausible(aemol):
    check = True
    
    plausible_bonds = {"1": [MIN_VALUE, MAX_VALUE]}
    .
    .
    .
    
    # is every bond within a plausible distance range
    
    return True
    
def check_atoms_dont_overlap(aemol):
    
    
    