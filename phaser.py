"""
phaser.py

Here's the module that kind of handles things relating to phasing!

Author: Ryan Baker
"""

# TODO: Compute/check all the complexities of these algorithms!
# TODO: List vs. tuple? Switch some lists to tuples?

from itertools import product, combinations_with_replacement
from genotype import Genotype
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

    # Generates all possible phases which are valid for one specific genotype
    # This is O(2^(k-1)) where k is the number of ambiguous sites (genotype value 2)
    # Here, we only generate permutations for the ambiguous sites,
    # then "insert" each of these permutations
    # Pseudocode:
    #   generate all possible length k-1 binary strings
    #   for each possible length k-1 binary string:
    #       assign 0 to the first ambiguous site in hap_one
    #       apply the possible binary string it to hap_one
    #       hap_two = complement(hap_one)
    #       add Phase(hap_one, hap_two) to our list
    def generate_valid_haplotype_phases(self, genotype):
        ambiguous_sites = [i for i, x in enumerate(genotype) if x == genotype.HETERO]
        k = len(ambiguous_sites)
        # Make a reference list where only the ambiguous sites are unfilled
        ref_hap = [None] * len(genotype)
        for i in xrange(len(genotype)):
            if genotype[i] == Genotype.HOMO_REF:
                ref_hap[i] = Haplotype.REF
            elif genotype[i] == Genotype.HOMO_ALT:
                ref_hap[i] = Haplotype.ALT
        phases = []
        # If no ambiguous sites, return the only possibility
        if k == 0:
            hap_one = Haplotype(list(ref_hap))
            hap_two = hap_one.complement(genotype)
            phases.append(Phase(hap_one, hap_two))
            return phases
        # Arbitrarily assign to the first unambiguous site
        ref_hap[ambiguous_sites[0]] = Haplotype.REF
        # If only one ambiguous site, return the only possibility
        if k == 1:
            hap_one = Haplotype(list(ref_hap))
            hap_two = hap_one.complement(genotype)
            phases.append(Phase(hap_one, hap_two))
            return phases
        # Case where there are several ambiguous sites...
        for assignment in product([0, 1], repeat=k-1):
            for i in xrange(1,k):
                site = ambiguous_sites[i]
                ref_hap[site] = assignment[i-1]
            hap_one = Haplotype(list(ref_hap))
            hap_two = hap_one.complement(genotype)
            phases.append(Phase(hap_one, hap_two))
        return phases

    # Generate all possible valid combinations of phases
    def generate_valid_phasings(self, genotypes):
        phases = [ self.generate_valid_haplotype_phases(g) for g in genotypes ]
        # phases is a list, where each element is a list of Phase objects
        # e.g. phases[0] is a list of the possible ways to phase genotype[0]
        return list(product(*phases))

    # Generate all possible combinations of phases
    # Each combination will have n Phase objects, each of which is m SNPs long
    # Returns a list of length O( (2^(2m))^n ) = O( 2^(2mn) )... woah.
    def generate_phasings(self, n, m):
        phases = []
        for p in product(self.generate_haplotype_combinations(m), repeat=n):
            phases.append(p)
        return phases

    # Takes a collection of phasings of n genotypes and removes those which are invalid
    # This is O(n * 2^(2mn))
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
        haplotype_pool = set()
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
    # Note that this algorithm is O( n * 2^(2mn) ). No bueno... but it is guaranteed to be optimal
    def bruteforce_min_parsimony_phase(self, genotypes):
        n = len(genotypes)
        m = genotypes[0].m
        phasings = self.generate_phasings(n, m)
        self.prune_phasings(phasings, genotypes)
        parsimonies = [ self.parsimony(phasing) for phasing in phasings ]
        min_parsimony = min(parsimonies)
        return phasings[parsimonies.index(min_parsimony)], min_parsimony

    # Same as above, but we only generate the 2^(k-1) phasings for k ambiguous (heterozygous) sites
    # This is considerably faster!
    def bruteforce_min_parsimony_phase_v2(self, genotypes):
        phasings = self.generate_valid_phasings(genotypes)
        parsimonies = [ self.parsimony(phasing) for phasing in phasings ]
        min_parsimony = min(parsimonies)
        return phasings[parsimonies.index(min_parsimony)], min_parsimony