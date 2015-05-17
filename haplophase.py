"""
haplophase.py

Program to run haplotype phasing!

Author: Ryan Baker
"""

from genotype import Genotype
from haplotype import Haplotype
from phase import Phase
from phaser import Phaser
import time

def main():

    phaser = Phaser()

    print 'Genotypes:'
    genotypes = [Genotype(random=True, length=5) for i in xrange(30)]
    for genotype in genotypes:
        print genotype

    print ''

    start_time = time.time()
    best_phasing, min_parsimony = phaser.phase_trivial_improved(genotypes)
    end_time = time.time()

    print 'Minimum parsimony: %d\nTime: %s sec' % (min_parsimony, end_time - start_time)
    for phase in best_phasing:
        print phase
        print ""

    start_time = time.time()
    greedy_phasing, greedy_parsimony = phaser.phase_greedy(genotypes)
    end_time = time.time()

    print 'Greedy parsimony: %d\nTime: %s sec' % (greedy_parsimony, end_time - start_time)
    for phase in greedy_phasing:
        print phase
        print ""

if __name__ == '__main__':
    main()
