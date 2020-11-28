'''
Description: Calculate MCE-18 index based on paper \
    "Are We Opening the Door to a New Era of Medicinal Chemistry or Being Collapsed to a Chemical Singularity"
Author: Kotori Y
Date: 2020-11-28 16:40:13
LastEditors: Kotori Y
LastEditTime: 2020-11-28 20:50:01
FilePath: \MCE-18\mce18.py
AuthorMail: kotori@cbdd.me
'''

from rdkit.Chem import AllChem as Chem
from collections import Counter
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from rdkit.Chem.rdMolDescriptors import CalcNumSpiroAtoms


class MCE18:
    """calculate the descriptor MCE-18, which can effectively 
    score molecules by novelty in terms of their cumulative sp3 complexity.
    """

    def __init__(self, mol):
        """Init

        Parameters
        ----------
        mol : rdkit.rdchem.Mol
            molecule to be calculated
        """
        self.mol = mol
        self.nC = len(
            [atom for atom in mol.GetAtoms()
                if atom.GetAtomicNum() == 6]
        )

    def _MolMatchSmarts(self, mol, smarts):
        """*internal only*
        """
        patt = Chem.MolFromSmarts(smarts)
        res = mol.GetSubstructMatches(patt)
        return res

    def CalculateQ1Index(self):
        """calculate normalized quadratic index (Q1 index)

        Returns
        -------
        Q1Index : float
            normalized quadratic index

        Reference
        ---------
        Balaban, Theor Chem Acc
        doi: 10.1007/BF00555695

        """
        matrix = GetAdjacencyMatrix(self.mol)
        c = Counter(matrix.sum(axis=1))
        temp = sum([(i**2 - 2*i) * v for i, v in c.items()])
        Q1Index = (temp + 2)/2
        return Q1Index

    def CalculateAR(self):
        """check the presence of an aromatic or heteroaromatic ring

        Returns
        -------
        AR : int, 0 or 1
            the presence of an aromatic or heteroaromatic ring (0 or 1)
        """
        smarts = "a"
        AR = bool(self._MolMatchSmarts(self.mol, smarts))
        AR = int(AR)
        return AR

    def CalculateNAR(self):
        """check the presence of an aliphatic or a heteroaliphatic ring

        Returns
        -------
        NAR : int, 0 or 1
            the presence of an aliphatic or a heteroaliphatic ring
        """
        smarts = "[A;R]"
        NAR = bool(self._MolMatchSmarts(self.mol, smarts))
        NAR = int(NAR)
        return NAR

    def CalculateCHIRAL(self):
        """check the presence of a chiral center (0 or 1)

        Returns
        -------
        CHIRAL : int, 0 or 1
            the presence of a chiral center (0 or 1)
        """
        CHIRAL = (Chem.CalcNumAtomStereoCenters(self.mol)) > 0
        CHIRAL = int(CHIRAL)
        return CHIRAL

    def CalculateSPIRO(self):
        """check the presence of a spiro center (0 or 1)

        Returns
        -------
        SPIRO : int, 0 or 1
            the presence of a spiro center (0 or 1)
        """
        SPIRO = CalcNumSpiroAtoms(self.mol) > 0
        SPIRO = int(SPIRO)
        return SPIRO

    def CalculateSP3(self):
        """calculate the portion of sp3-hybridized carbon atoms (from 0 to 1)

        Returns
        -------
        sp3 : float, from 0 to 1
            the portion of sp3-hybridized carbon atoms (from 0 to 1)
        """
        smarts = "[CX4]"
        sp3 = len(self._MolMatchSmarts(self.mol, smarts))
        sp3 = sp3/self.nC
        return sp3


if "__main__" == __name__:

    smiles = "C1NCCN(C2=CC3N(CC)C=C(C(=O)O)C(=O)C=3C=C2F)C1"
    mol = Chem.MolFromSmiles(smiles)
    
    demo = MCE18(mol)
    
    print("DONE")
