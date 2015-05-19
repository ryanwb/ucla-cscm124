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
import random
import argparse

def main():

    parser = argparse.ArgumentParser(description="Phase genotype data into haplotypes (min parsimony)")
    parser.add_argument("-n", type=int,
                    help="number of individuals (each individual has two haplotypes)")
    parser.add_argument("-m", type=int,
                    help="number of SNPs")
    parser.add_argument("-e", "--exhaustive", action="store_true",
                    help="use exhaustive algorithm")
    parser.add_argument("-g", "--greedy", action="store_true", 
                    help="use greedy algorithm")
    parser.add_argument("-x", "--hash", action="store_true",
                    help="use greedy hash lookup algorithm")
    parser.add_argument("-r", "--reffile",
                    help="reference hapmap file to fill hash table in greedy hash lookup algorithm")
    parser.add_argument("-f", "--file",
                    help="hapmap phased haplotype input data file")
    parser.add_argument("-p", type=int,
                    help="number of haplotypes in the haplotype pool (if no input data file)")
    parser.add_argument("-v", "--verbose", action="store_true",
                    help="print extra output/data")

    args = parser.parse_args()

    phaser = Phaser()

    real_phase_data = []
    if not args.file:
        # generate a pool of haplotypes
        hap_pool = phaser.generate_haplotypes(args.m)
        # randomly use p of these haplotypes
        pool = random.sample(xrange(len(hap_pool)), args.p)
        for i in xrange(args.n):
            hap_one = hap_pool[pool[random.randrange(len(pool))]]
            hap_two = hap_pool[pool[random.randrange(len(pool))]]
            real_phase_data.append(Phase(hap_one, hap_two))
    else:
        # grab haplotype data from input file
        hapfile = open(args.file)
        data = []
        # throw away first line and first two columns
        hapfile.readline()
        for i in xrange(args.m):
            data.append(hapfile.readline().split()[2:])
        # now parse out the haplotypes!
        for i_n in xrange(args.n):       # for each individual
            hap_one_data = []
            hap_two_data = []
            for i_m in xrange(args.m):   # for each SNP
                # arbitrarily assign 0 to the SNP that individual number 0 has at this SNP position
                ref_snp = data[i_m][0]
                hap_one_data.append(0 if data[i_m][i_n * 2] == ref_snp else 1)
                hap_two_data.append(0 if data[i_m][i_n * 2 + 1] == ref_snp else 1)
            real_phase_data.append(Phase(Haplotype(list(hap_one_data)), Haplotype(list(hap_two_data))))

    real_parsimony = phaser.parsimony(real_phase_data)

    # now turn the haplotype data into genotypes, which we will try to phase
    genotypes = [real_phase_data[i_n].to_genotype() for i_n in xrange(len(real_phase_data))]

    if args.verbose:
        for phase in real_phase_data:
            print phase
            print ""

    # start the timer
    start_time = time.time()

    # run the phasing algorithm!
    if args.greedy:
        phasing, parsimony = phaser.phase_greedy(genotypes)

    elif args.hash:
        ref_file = open(args.reffile)
        ref_data = []
        ref_haps = []
        # throw away first line and first two columns
        ref_file.readline()
        for line in hapfile.readlines():
            ref_data.append(line.split()[2:])
        # now parse out the haplotypes!
        for i_n in xrange(len(ref_data[0])):       # for each haplotype
            ref_hap = []
            for i_m in xrange(args.m):   # for each SNP of interest
                # arbitrarily assign 0 to the SNP that individual number 0 has at this SNP position
                ref_snp = ref_data[i_m][0]
                ref_hap.append(0 if ref_data[i_m][i_n] == ref_snp else 1)
            ref_haps.append(Haplotype(list(ref_hap)))
        phasing, parsimony = phaser.phase_hash(genotypes, ref_haps)

    else: # if args.exhaustive
        phasing, parsimony = phaser.phase_trivial_improved(genotypes)

    # stop the timer
    end_time = time.time()
    elapsed_time = end_time - start_time

    # check our accuracy
    total_phases = len(real_phase_data)
    correct_phases = 0
    for i_n in xrange(total_phases):
        if real_phase_data[i_n] == phasing[i_n]:
            correct_phases += 1
    accuracy = float(correct_phases)/float(total_phases)

    if args.verbose:
        print "Parsimony: %d\n" % parsimony
        for phase in phasing:
            print phase
            print ""

    # output will be in format:
    # elapsed_time correct_phases total_phases accuracy found_parsimony actual_parsimony
    print "%g %d %d %.5f %d %d" % (elapsed_time, correct_phases, total_phases, accuracy, parsimony, real_parsimony)

if __name__ == '__main__':
    main()
