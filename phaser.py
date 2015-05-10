"""
phaser.py

Here's the module that kind of handles things relating to phasing!

Author: Ryan Baker
"""

from itertools import product
from haplotype import Haplotype

# Phaser class
class Phaser:

    def __init__(self):
        return

    # generate all possible haplotypes of length m (SNPs)
    # returns a list of haplotypes
    def generate_haplotypes(self, m):
        haplotypes = [];
        for h in product([0, 1], repeat=m):
            haplotypes.append(Haplotype(h))
        return haplotypes
