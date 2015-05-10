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
    print 'Genotypes:'
    genotypes = [Genotype([1, 2, 1]), Genotype([0, 2, 2]), Genotype([2, 0, 0])]
    for genotype in genotypes:
        print genotype

    print ''

    phaser = Phaser()
    best_phasing = phaser.bruteforce_min_parsimony_phase(genotypes)

    print 'Minimum parsimony phasing:'
    for phase in best_phasing:
        print phase
        print ""

if __name__ == '__main__':
    main()
