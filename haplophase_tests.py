"""
haplophase.py

Program to run haplotype phasing!

Author: Ryan Baker
"""

from phaserunner import PhaseRunner

def main():

    pr = PhaseRunner()
    datafile = "haps.txt"

    n = 10
    m = 4

    print pr.run_exhaustive(n, m, datafile)
    print pr.run_exhaustive(n, m, None, True, 5)

    print pr.run_greedy(n, m, datafile)
    print pr.run_greedy(n, m, None, True, 5)

    print pr.run_hash(n, m, datafile)
    print pr.run_hash(n, m, datafile, True, 5)

    # output will be in format:
    # elapsed_time correct_phases total_phases accuracy found_parsimony actual_parsimony
    # print "%g %d %d %.5f %d %d" % (elapsed_time, correct_phases, total_phases, accuracy, parsimony, real_parsimony)

if __name__ == '__main__':
    main()
