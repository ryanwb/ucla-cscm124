"""
haplophase.py

Program to run haplotype phasing!

Author: Ryan Baker
"""

from genotype import Genotype
from haplotype import Haplotype
from phase import Phase
from phaser import Phaser

def main():

    phaser = Phaser()

    print 'Genotypes:'
    genotypes = [Genotype([2, 2, 2, 2, 2, 2]), Genotype([2, 2, 2, 2, 2, 2]), Genotype([1, 1, 1, 1, 1, 1]), Genotype([2, 0, 0, 0, 0, 0])]
    for genotype in genotypes:
        print genotype

    print ''

    best_phasing, parsimony = phaser.bruteforce_min_parsimony_phase_v2(genotypes)

    print 'Minimum parsimony: %d' % parsimony
    for phase in best_phasing:
        print phase
        print ""

if __name__ == '__main__':
    main()
