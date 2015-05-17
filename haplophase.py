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
    parser.add_argument("-f", "--file",
                    help="hapmap phased haplotype input data file")
    parser.add_argument("-v", "--verbose", action="store_true",
                    help="print extra output/data")

    args = parser.parse_args()

    phaser = Phaser()

    # TODO: if no input file is provided, we want to create a pool of size p of haplotypes of length m,
    # and randomly choose from this pool to create our haplotype data

    '''
    print 'Genotypes:'
    genotypes = [Genotype(random=True, length=args.m) for i in xrange(args.n)]
    for genotype in genotypes:
        print genotype
    '''

    hapfile = open(args.file)
    data = []
    # throw away first line and first two columns
    hapfile.readline()
    for i in xrange(args.m):
        data.append(hapfile.readline().split()[2:])
    # now parse out the haplotypes!
    real_phase_data = []
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
    else: # if args.exhaustive
        phasing, parsimony = phaser.phase_trivial_improved(genotypes)

    # TODO: my modification to the algorithm:
    # hash the entire input data set to get haplotype frequencies
    # then for each possible haplotype, just pick the most frequent one

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
