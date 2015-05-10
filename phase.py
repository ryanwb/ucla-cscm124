"""
phase.py

Module for phased haplotype data

Author: Ryan Baker
"""

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
            if genotype.data[x] == genotype.HOMO_REF:
                if not (self.haps[0][x] == self.haps[0].REF and self.haps[1][x] == self.haps[1].REF):
                    return False
            elif genotype.data[x] == genotype.HOMO_ALT:
                if not (self.haps[0][x] == self.haps[0].ALT and self.haps[1][x] == self.haps[1].ALT):
                    return False
            else: # genotype.HETERO
                if self.haps[0][x] == self.haps[1][x]:
                    return False
        return True

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
