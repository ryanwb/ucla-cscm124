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
import argparse

def main():

    parser = argparse.ArgumentParser(description='Phase genotype data into haplotypes (min parsimony)')
    parser.add_argument('-n', type=int,
                    help='number of individuals (each individual has two haplotypes)')
    parser.add_argument('-m', type=int,
                    help='number of SNPs')
    parser.add_argument("-e", "--exhaustive", action="store_true",
                    help="use exhaustive algorithm")
    parser.add_argument("-g", "--greedy", action="store_true", 
                    help="use greedy algorithm")

    args = parser.parse_args()

    phaser = Phaser()

    print 'Genotypes:'
    genotypes = [Genotype(random=True, length=args.m) for i in xrange(args.n)]
    for genotype in genotypes:
        print genotype

    print ''

    if args.greedy:
        phasing, parsimony = phaser.phase_greedy(genotypes)
    else: # if args.exhaustive
        phasing, parsimony = phaser.phase_trivial_improved(genotypes)

    print 'Parsimony: %d\n' % parsimony
    for phase in phasing:
        print phase
        print ""

if __name__ == '__main__':
    main()
