"""
haplophase.py

Command line program to run haplotype phasing

Author: Ryan Baker
"""

from phaserunner import PhaseRunner
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
    parser.add_argument("-f", "--file",
                    help="hapmap phased haplotype input data file")
    parser.add_argument("-p", type=int,
                    help="number of haplotypes in the haplotype pool (if random haplotype pool should be used as input)")
    parser.add_argument("-v", "--verbose", action="store_true",
                    help="print extra output/data")

    args = parser.parse_args()

    pr = PhaseRunner()


    # run the phasing algorithm!
    # they return: elapsed_time, self.get_accuracy(real_phase_data, phasing), (parsimony, real_parsimony), real_phase_data, phasing
    if args.greedy:
        t, acc, pars, real_phase_data, phasing = pr.run_greedy(args.n, args.m, args.file, args.p != None, args.p)
    elif args.hash:
        t, acc, pars, real_phase_data, phasing = pr.run_hash(args.n, args.m, args.file, args.p != None, args.p)
    else: # if args.exhaustive
        t, acc, pars, real_phase_data, phasing = pr.run_exhaustive(args.n, args.m, args.file, args.p != None, args.p)

    if args.verbose:
        for phase in real_phase_data:
            print phase
            print ""
        print "Parsimony: %d\n" % pars[0]
        for phase in phasing:
            print phase
            print ""

    # output will be in format:
    # elapsed_time correct_phases/total_phases accuracy found_parsimony/actual_parsimony
    print "%g %d/%d %.5f %d/%d" % (t, acc[0], acc[1], acc[2], pars[0], pars[1])

if __name__ == '__main__':
    main()
