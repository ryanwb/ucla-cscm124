"""
haplophase.py

Program to run haplotype phasing tests and record results

Author: Ryan Baker
"""

from phaserunner import PhaseRunner

# output will be in format:
# elapsed_time correct_phases/total_phases accuracy found_parsimony/actual_parsimony
def out(f, n, m, algo, t, acc, pars, real_phase_data, phasing):
    if not f:
        f = "random"
    print "%s %d %d %s %g %d/%d %.5f %d/%d" % (f, n, m, algo, t, acc[0], acc[1], acc[2], pars[0], pars[1])

def main():

    pr = PhaseRunner()
    datafiles = ["data1.txt", "data2.txt"]

    for f in datafiles:
        n = 75
        for m in xrange(1, 7):
            t, acc, pars, real_phase_data, phasing = pr.run_exhaustive(n, m, f)
            out(f, n, m, "exhaustive", t, acc, pars, real_phase_data, phasing)
        for m in xrange(1, 16):
            t, acc, pars, real_phase_data, phasing = pr.run_greedy(n, m, f)
            out(f, n, m, "greedy", t, acc, pars, real_phase_data, phasing)
        for m in xrange(1, 31, 2):
            t, acc, pars, real_phase_data, phasing = pr.run_hash(n, m, f)
            out(f, n, m, "hash", t, acc, pars, real_phase_data, phasing)
        m = 10
        for n in xrange(1, 19):
            t, acc, pars, real_phase_data, phasing = pr.run_exhaustive(n, m, f)
            out(f, n, m, "exhaustive", t, acc, pars, real_phase_data, phasing)
        for n in xrange(1, 76, 3):
            t, acc, pars, real_phase_data, phasing = pr.run_greedy(n, m, f)
            out(f, n, m, "greedy", t, acc, pars, real_phase_data, phasing)
        for n in xrange(1, 76, 3):
            t, acc, pars, real_phase_data, phasing = pr.run_hash(n, m, f)
            out(f, n, m, "hash", t, acc, pars, real_phase_data, phasing)

    f = None
    n = 50
    p = 4
    for m in xrange(1, 6):
        t, acc, pars, real_phase_data, phasing = pr.run_exhaustive(n, m, f, True, p)
        out(f, n, m, "exhaustive", t, acc, pars, real_phase_data, phasing)
    for m in xrange(1, 16, 2):
        t, acc, pars, real_phase_data, phasing = pr.run_greedy(n, m, f, True, p)
        out(f, n, m, "greedy", t, acc, pars, real_phase_data, phasing)
    for m in xrange(1, 16, 2):
        t, acc, pars, real_phase_data, phasing = pr.run_hash(n, m, f, True, p)
        out(f, n, m, "hash", t, acc, pars, real_phase_data, phasing)

if __name__ == '__main__':
    main()
