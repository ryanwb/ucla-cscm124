"""
phase.py

Module for phased haplotype data

Author: Ryan Baker
"""

from genotype import Genotype
from haplotype import Haplotype

# Phase class
class Phase:

    def __init__(self, hap_one, hap_two):
        self.haps = [hap_one, hap_two]  # haps is a list of two haplotypes
        self.m = len(hap_one)           # m is the number of SNPs
        if hap_one.m != hap_two.m:
            raise ValueError("phase data has mismatched haplotypes")
        return

    # Returns true/false based on whether this phase is congruent with a given genotype
    # O(m)
    def matches(self, genotype):
        for x in xrange(genotype.m):
            if genotype.data[x] == Genotype.HOMO_REF:
                if not (self.haps[0][x] == Haplotype.REF and self.haps[1][x] == Haplotype.REF):
                    return False
            elif genotype.data[x] == Genotype.HOMO_ALT:
                if not (self.haps[0][x] == Haplotype.ALT and self.haps[1][x] == Haplotype.ALT):
                    return False
            else: # Genotype.HETERO
                if self.haps[0][x] == self.haps[1][x]:
                    return False
        return True

    # Returns a genotype representation of this pair of haplotypes
    def to_genotype(self):
        genotype_data = [0] * self.m
        for x in xrange(self.m):
            if self.haps[0][x] == self.haps[1][x]:
                if self.haps[0][x] == Haplotype.REF:
                    genotype_data[x] = Genotype.HOMO_REF
                else:   # Haplotype.ALT
                    genotype_data[x] = Genotype.HOMO_ALT
            else:   # heterozygous case
                genotype_data[x] = Genotype.HETERO
        return Genotype(genotype_data)

    # Override equality so that the order of the haplotypes does not matter
    def __eq__(self, other):
        if isinstance(other, Phase):
            return ((self.haps[0] == other.haps[0] and self.haps[1] == other.haps[1])
                or (self.haps[0] == other.haps[1] and self.haps[1] == other.haps[0]))
        else:
            return False

    # Override hash function
    def __hash__(self):
        return hash(tuple(self.haps[0])) + hash(tuple(self.haps[1]))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getitem__(self, i):
        return self.haps[i]

    def __len__(self):
        return self.m

    def __repr__(self):
        return repr(self.haps)

    def __str__(self):
        s = ""
        for c in xrange(len(self.haps)):
            s += str(self.haps[c])
            if (c < len(self.haps) - 1):
                s += '\n'
        return s
