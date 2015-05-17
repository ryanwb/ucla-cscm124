"""
genotype.py

Module for genotype data

Author: Ryan Baker
"""

from random import randrange

# Genotype class
class Genotype:

    # Encode pairs of alleles as 0s, 1s, and 2s
    HOMO_REF = 0
    HOMO_ALT = 1
    HETERO   = 2

    # e.g. Genotype([0 2 0 1]) explicitly, or
    # can do Genotype(random=True, length=5) to generate a random genotype of length 5
    def __init__(self, data=None, random=False, length=None):
        if data is not None:
            self.data = data    # the genotype data; list of ints with values 0, 1, 2
            self.m = len(data)  # m is the number of SNPs
            for x in self.data:
                if not (x == Genotype.HOMO_REF or x == Genotype.HOMO_ALT or x == Genotype.HETERO):
                    raise ValueError("bad genotype data")
        else:
            if length is None:
                raise ValueError("must provide length for random genotype")
            self.m = length
            self.data = [randrange(3) for i in xrange(length)]
        return

    def __hash__(self):
        return hash(tuple(self.data))

    def __eq__(self, other):
        if isinstance(other, Genotype):
            return (tuple(self.data) == tuple(other.data))
        else:
            return False

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
