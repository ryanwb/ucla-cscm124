"""
haplotype.py

Module for haplotype data

Author: Ryan Baker
"""

from genotype import Genotype

# Haplotype class
class Haplotype:

    # Encode alleles as 0s and 1s
    REF = 0
    ALT = 1

    def __init__(self, data):
        self.data = data    # the genotype data; list of ints with values 0, 1, 2
        self.m = len(data)  # m is number of SNPs
        for x in self.data:
            if not (x == Haplotype.REF or x == Haplotype.ALT):
                raise ValueError("bad haplotype data")
        return

    # Generate a complementary haplotype for a given genotype
    # This is O(m), where m is the genotype length
    def complement(self, genotype):
        if self.m != len(genotype):
            raise ValueError("genotype/haplotype length mismatch")
        complement = [None] * self.m
        for x in xrange(self.m):
            if genotype[x] == Genotype.HOMO_REF or genotype[x] == Genotype.HOMO_ALT:
                if genotype[x] == Genotype.HOMO_REF and self.data[x] == Haplotype.ALT:
                    raise ValueError("non-complementable haplotype for given genotype")
                if genotype[x] == Genotype.HOMO_ALT and self.data[x] == Haplotype.REF:
                    raise ValueError("non-complementable haplotype for given genotype")
                complement[x] = self.data[x]
            else: # HETERO
                if self.data[x] == Haplotype.REF:
                    complement[x] = Haplotype.ALT
                else: # ALT
                    complement[x] = Haplotype.REF
        return Haplotype(complement)

    # Determines whether a haplotype explains a given genotype
    def explains(self, genotype):
        if self.m != len(genotype):
            raise ValueError("genotype/haplotype length mismatch")
        for x in xrange(self.m):
            if self.data[x] == Haplotype.REF and genotype[x] == Genotype.HOMO_ALT:
                return False
            if self.data[x] == Haplotype.ALT and genotype[x] == Genotype.HOMO_REF:
                return False
        return True

    # Override the hash and eq functions since we will be
    # throwing Haplotype objects into a set to count parsimony

    def __hash__(self):
        return hash(tuple(self.data))

    def __eq__(self, other):
        if isinstance(other, Haplotype):
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
