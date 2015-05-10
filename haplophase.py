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
    p = Phaser()
    for h in p.generate_haplotypes(3):
        print h
    return

if __name__ == '__main__':
    main()
