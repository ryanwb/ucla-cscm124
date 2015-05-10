"""
genotype.py

Module for genotype data

Author: Ryan Baker
"""

# Genotype class
class Genotype:

    # Encode pairs of alleles as 0s, 1s, and 2s
    HOMO_REF = 0
    HOMO_ALT = 1
    HETERO   = 2

    def __init__(self, data):
        self.data = data    # the genotype data; list of ints with values 0, 1, 2
        self.m = len(data)  # m is the number of SNPs
        for x in self.data:
            if not (x == self.HOMO_REF or x == self.HOMO_ALT or x == self.HETERO):
                raise ValueError("bad genotype data")
        return

    def __getitem__(self, x):
        return self.data[x]

    def __len__(self):
        return self.m

    def __repr__(self):
        return repr(self.data)

    def __str__(self):
        s = ""
        for x in self.data:
            s += str(x)
        return s
