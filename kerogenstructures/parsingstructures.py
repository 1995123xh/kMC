# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 14:43:26 2019

@author: Hao
"""
"""
Output contains three variables: atoms is a character array of the atomic symbols; kIIAbm is the bond-weight matrix;
ringinfo is the lowest memeber of ring that each atom is in.
"""
import re
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

kIIA=Chem.rdmolfiles.MolFromMol2File('EFK_50A_80.mol2')
kIIAm = Chem.rdmolops.GetAdjacencyMatrix(kIIA)
Draw.MolToImage(kIIA,size=(1500, 1500),noCarbonSymbols=False,includeAtomNumbers=True)
kIIAbm=np.array(kIIAm, dtype=float)

for b in kIIA.GetBonds():       #Loop over all bonds to assign bond order values to bond order matrix
    if b.GetBondType()==Chem.rdchem.BondType(2):
        kIIAbm[b.GetBeginAtomIdx(), b.GetEndAtomIdx()]=2
        kIIAbm[b.GetEndAtomIdx(), b.GetBeginAtomIdx()]=2
    elif b.GetBondType()==Chem.rdchem.BondType(3):
        kIIAbm[b.GetBeginAtomIdx(), b.GetEndAtomIdx()]=3
        kIIAbm[b.GetEndAtomIdx(), b.GetBeginAtomIdx()]=3
    elif b.GetBondType()==Chem.rdchem.BondType(12):
        kIIAbm[b.GetBeginAtomIdx(), b.GetEndAtomIdx()]=1.5
        kIIAbm[b.GetEndAtomIdx(), b.GetBeginAtomIdx()]=1.5
#    else:
#        kIIAbm[b.GetBeginAtomIdx(), b.GetEndAtomIdx()]=100
#        kIIAbm[b.GetEndAtomIdx(), b.GetBeginAtomIdx()]=100
atoms=[]
ringinfo=np.zeros(len(kIIAm))
j=0
for a in kIIA.GetAtoms():
    atoms.append(a.GetSymbol())
    for n in range(7):
        if a.IsInRingSize(n+1):
            ringinfo[j]=n+1
            break
    j=j+1


                    
    
#import kerogen structure
