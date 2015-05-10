"""
haplotype.py

Module for haplotype data

Author: Ryan Baker
"""

# Haplotype class
class Haplotype:

    # encode alleles as 0s and 1s
    REF = 0
    ALT = 1

    def __init__(self, data):
        self.data = data    # the genotype data; list of ints with values 0, 1, 2
        self.m = len(data)  # m is number of SNPs
        for x in self.data:
            if not (x == self.REF or x == self.ALT):
                raise ValueError("bad haplotype data")
        return

    def __len__(self):
        return self.m

    def __repr__(self):
        return repr(self.data)

    def __str__(self):
        s = ""
        for x in self.data:
            s += str(x)
        return s
