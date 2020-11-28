'''
Description: Calculate MCE-18 index based on paper \
    "Are We Opening the Door to a New Era of Medicinal Chemistry or Being Collapsed to a Chemical Singularity"
Author: Kotori Y
Date: 2020-11-28 16:40:13
LastEditors: Kotori Y
LastEditTime: 2020-11-28 19:49:49
FilePath: \MCE-18\mce18.py
AuthorMail: kotori@cbdd.me
'''

from rdkit.Chem import AllChem as Chem
from collections import Counter
from rdkit.Chem.rdmolops import GetAdjacencyMatrix


def CalculateQ1Index(mol):
    """calculate normalized quadratic index (Q1 index)

    Parameters
    ----------
    mol : rdkit.rdchem.Mol
        molecule to be calculated

    Returns
    -------
    Q1Index : float
        normalized quadratic index

    Reference
    ---------
    Balaban, Theor Chem Acc
    doi: 10.1007/BF00555695
    
    """
    matrix = GetAdjacencyMatrix(mol)
    c = Counter(matrix.sum(axis=1))
    temp = sum([(i**2 - 2*i) * v for i,v in c.items()])
    Q1Index = (temp + 2)/2
    return Q1Index


def CalculateAR(mol):
    """check the presence of an aromatic or heteroaromatic ring

    Parameters
    ----------
    mol : rdkit.rdchem.Mol
        molecule to be checked

    Returns
    -------
    AR : int, 0 or 1
        the presence of an aromatic or heteroaromatic ring (0 or 1)
    """
    smart = "a"
    aroma = Chem.MolFromSmarts(smart)
    AR = int(mol.HasSubstructMatch(aroma))
    return AR


def CalculateNAR(mol):
    """check the presence of an aliphatic or a heteroaliphatic ring

    Parameters
    ----------
    mol : rdkit.rdchem.Mol
        molecule to be checked

    Returns
    -------
    NAR : int, 0 or 1
        the presence of an aliphatic or a heteroaliphatic ring
    """
    smart = "[A;R]"
    aroma = Chem.MolFromSmarts(smart)
    NAR = int(mol.HasSubstructMatch(aroma))
    return NAR
    

if "__main__" == __name__:
    
    smiles = "C1NCCN(C2=CC3N(CC)C=C(C(=O)O)C(=O)C=3C=C2F)C1"
    mol = Chem.MolFromSmiles(smiles)

    AR = CalculateAR(mol)
    NAR = CalculateNAR(mol)
    Q1Index = CalculateQ1Index(mol)
    print("DONE")