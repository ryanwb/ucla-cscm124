"""
phaser.py

Here's the module that kind of handles things relating to phasing!

Author: Ryan Baker
"""

from itertools import product, combinations_with_replacement
from haplotype import Haplotype
from sets import Set
from phase import Phase

# Phaser class
class Phaser:

    def __init__(self):
        return

    # Generate all possible haplotypes of length m (SNPs)
    # Returns a list of haplotypes (length 2^m)
    def generate_haplotypes(self, m):
        haplotypes = []
        for h in product([0, 1], repeat=m):
            haplotypes.append(Haplotype(h))
        return haplotypes

    # Generate all possible combinations of two haplotypes of length m (SNPs)
    # Returns a list of length 1 + 2 + ... + 2^m, which is O( 2^(2m) )
    # Each element in the list is a Phase object (with two haplotypes)
    def generate_haplotype_combinations(self, m):
        combinations = []
        for c in combinations_with_replacement(self.generate_haplotypes(m), 2):
            combinations.append(Phase(c[0], c[1]))
        return combinations

    # Generate all possible combinations of phases
    # Each combination will have n Phase objects, each of which is m SNPs long
    # Returns a list of length O( (2^(2m))^n ) = O( 2^(2mn) )... woah.
    def generate_phasings(self, n, m):
        phases = []
        for p in product(self.generate_haplotype_combinations(m), repeat=n):
            phases.append(p)
        return phases

    # Takes a collection of phasings of n genotypes and removes those which are invalid
    def prune_phasings(self, phasings, genotypes):
        phasings_to_remove = []
        for phasing in phasings:
            if len(phasing) != len(genotypes):
                raise ValueError("number of phases must match number of genotypes")
            m = len(genotypes)
            for p in xrange(m):
                if not phasing[p].matches(genotypes[p]):
                    # bad phasing encountered. remember this phasing and break
                    phasings_to_remove.append(phasing)
                    break
        # now remove all of the bad phasings
        for bad_phasing in phasings_to_remove:
            phasings.remove(bad_phasing)

    # Computes the parsimony of a phasing
    # Input is a list of n Phase objects, representing a phasing of n genotypes
    def parsimony(self, phasing):
        # make a set of haplotypes in this phasing
        haplotype_pool = Set()
        for phase in phasing:
            haplotype_pool.add(phase[0])
            haplotype_pool.add(phase[1])
        # return the cardinality of the set of haplotypes we encountered
        return len(haplotype_pool)

    # Run the trivial brute force algorithm for phasing
    # Basic steps:
    #   1. Generate all possible phasings of n genotypes, regardless of the genotype information
    #   2. Prune this list, removing phasings that are incongruous with the genotype information
    #   3. At this point, we have all "valid" phasings... check which one has minimum parsimony
    # Note that this algorithm is O( 2^(2mn) ). No bueno... but it is guaranteed to be optimal
    def bruteforce_min_parsimony_phase(self, genotypes):
        n = len(genotypes)
        m = genotypes[0].m
        phasings = self.generate_phasings(n, m)
        self.prune_phasings(phasings, genotypes)
        parsimonies = [ self.parsimony(phasing) for phasing in phasings ]
        return phasings[parsimonies.index(min(parsimonies))]
