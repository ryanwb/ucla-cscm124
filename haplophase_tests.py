"""
haplophase_tests.py

Program to run haplotype phasing tests and record results

Author: Ryan Baker
"""

from phaserunner import PhaseRunner

# output will be in format:
# inputfile n m  algorithm elapsed_time correct_phases total_phases accuracy found_parsimony actual_parsimony
def out(f, n, m, algo, t, acc, pars, real_phase_data, phasing):
    if not f:
        f = "random"
    print "%s %d %d %s %g %d %d %.5f %d %d" % (f, n, m, algo, t, acc[0], acc[1], acc[2], pars[0], pars[1])

def main():
    pr = PhaseRunner()
    f = "data2.txt"
    n = 25
    p = 4
    for m in xrange(2, 4):
        t, acc, pars, real_phase_data, phasing = pr.run_exhaustive(n, m, f, True, p)
        out(f, n, m, "exhaustive", t, acc, pars, real_phase_data, phasing)
    for m in xrange(2, 16, 2):
        t, acc, pars, real_phase_data, phasing = pr.run_greedy(n, m, f, True, p)
        out(f, n, m, "greedy", t, acc, pars, real_phase_data, phasing)
    for m in xrange(2, 16, 2):
        t, acc, pars, real_phase_data, phasing = pr.run_hash(n, m, f, True, p)
        out(f, n, m, "hash", t, acc, pars, real_phase_data, phasing)

if __name__ == '__main__':
    main()
