"""
phaserunner.py

Contains a class that we can use which wraps the phaser class so that we can easily run the phasing algorithms

Author: Ryan Baker
"""

from genotype import Genotype
from haplotype import Haplotype
from phase import Phase
from phaser import Phaser
import time
import random

class PhaseRunner:

    def __init__(self):
        self.phaser = Phaser()
        return

    def random_phase_data(self, n, m, p):
        real_phase_data = []
        # generate a pool of haplotypes
        hap_pool = self.phaser.generate_haplotypes(m)
        # randomly use p of these haplotypes
        pool = random.sample(xrange(len(hap_pool)), p)
        for i in xrange(n):
            hap_one = hap_pool[pool[random.randrange(len(pool))]]
            hap_two = hap_pool[pool[random.randrange(len(pool))]]
            real_phase_data.append(Phase(hap_one, hap_two))
        real_parsimony = self.phaser.parsimony(real_phase_data)
        return real_phase_data, real_parsimony

    def file_phase_data(self, filename, n, m):
        real_phase_data = []
        # grab haplotype data from input file
        hapfile = open(filename)
        data = []
        # throw away first line and first two columns
        hapfile.readline()
        for i in xrange(m):
            data.append(hapfile.readline().split()[2:])
        # now parse out the haplotypes!
        for i_n in xrange(n):       # for each individual
            hap_one_data = []
            hap_two_data = []
            for i_m in xrange(m):   # for each SNP
                # arbitrarily assign 0 to the SNP that individual number 0 has at this SNP position
                ref_snp = data[i_m][0]
                hap_one_data.append(0 if data[i_m][i_n * 2] == ref_snp else 1)
                hap_two_data.append(0 if data[i_m][i_n * 2 + 1] == ref_snp else 1)
            real_phase_data.append(Phase(Haplotype(list(hap_one_data)), Haplotype(list(hap_two_data))))
        real_parsimony = self.phaser.parsimony(real_phase_data)
        return real_phase_data, real_parsimony

    def to_genotypes(self, real_phase_data):
        # now turn the haplotype data into genotypes, which we will try to phase
        return [real_phase_data[i_n].to_genotype() for i_n in xrange(len(real_phase_data))]

    def get_accuracy(self, real_phase_data, phasing):
        total_phases = len(real_phase_data)
        correct_phases = 0
        for i_n in xrange(total_phases):
            if real_phase_data[i_n] == phasing[i_n]:
                correct_phases += 1
        accuracy = float(correct_phases)/float(total_phases)
        return correct_phases, total_phases, accuracy

    def run_greedy(self, n, m, hapmapfile, random=False, p=None):
        if not random:
            real_phase_data, real_parsimony = self.file_phase_data(hapmapfile, n, m)
        else:
            real_phase_data, real_parsimony = self.random_phase_data(n, m, p)
        genotypes = self.to_genotypes(real_phase_data)
        start_time = time.time()
        phasing, parsimony = self.phaser.phase_greedy(genotypes)
        end_time = time.time()
        elapsed_time = end_time - start_time
        return elapsed_time, self.get_accuracy(real_phase_data, phasing), (parsimony, real_parsimony)

    def run_exhaustive(self, n, m, hapmapfile, random=False, p=None):
        if not random:
            real_phase_data, real_parsimony = self.file_phase_data(hapmapfile, n, m)
        else:
            real_phase_data, real_parsimony = self.random_phase_data(n, m, p)
        genotypes = self.to_genotypes(real_phase_data)
        start_time = time.time()
        phasing, parsimony = self.phaser.phase_trivial_improved(genotypes)
        end_time = time.time()
        elapsed_time = end_time - start_time
        return elapsed_time, self.get_accuracy(real_phase_data, phasing), (parsimony, real_parsimony)

    def run_hash(self, n, m, hapmapfile, random=False, p=None):
        if not random:
            real_phase_data, real_parsimony = self.file_phase_data(hapmapfile, n, m)
        else:
            real_phase_data, real_parsimony = self.random_phase_data(n, m, p)
        start_time = time.time()
        # Just use the same input file for the frequency hash table
        ref_file = open(hapmapfile)
        ref_data = []
        ref_phases = []
        # throw away first line and first two columns
        ref_file.readline()
        for line in hapfile.readlines():
            ref_data.append(line.split()[2:])
        # now parse out the haplotypes!
        for i_n in xrange(0, len(ref_data[0]), 2):       # for each haplotype pair
            hap_one_data = []
            hap_two_data = []
            for i_m in xrange(m):   # for each SNP of interest
                # arbitrarily assign 0 to the SNP that individual number 0 has at this SNP position
                ref_snp = ref_data[i_m][0]
                hap_one_data.append(0 if ref_data[i_m][i_n] == ref_snp else 1)
                hap_two_data.append(0 if ref_data[i_m][i_n + 1] == ref_snp else 1)
            ref_phases.append(Phase(Haplotype(list(hap_one_data)), Haplotype(list(hap_two_data))))
        phasing, parsimony = self.phaser.phase_hash(genotypes, ref_phases)
        end_time = time.time()
        elapsed_time = end_time - start_time
        return elapsed_time, self.get_accuracy(real_phase_data, phasing), (parsimony, real_parsimony)


if __name__ == '__main__':
    main()
