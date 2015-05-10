"""
phase.py

Module for phased haplotype data

Author: Ryan Baker
"""

# Phase class
class Phase:

    def __init__(self, hap_one, hap_two):
        self.data = [hap_one, hap_two]
        self.m = len(hap_one)
        if hap_one.m != hap_two.m:
            raise ValueError("phase data has mismatched haplotypes")
        return

    def __repr__(self):
        return repr(self.data)

    def __str__(self):
        s = ""
        for c in xrange(len(self.data)):
            s += str(self.data[c])
            if (c < len(self.data) - 1):
                s += '\n'
        return s
