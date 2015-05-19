"""
process_haplophase_tests.py

Program to process the output from the test script

Author: Ryan Baker
"""

import matplotlib.pyplot as plt

# output will be in format:
# inputfile n m  algorithm elapsed_time correct_phases total_phases accuracy found_parsimony actual_parsimony
# def out(f, n, m, algo, t, acc, pars, real_phase_data, phasing):
#    if not f:
#        f = "random"
#    print "%s %d %d %s %g %d %d %.5f %d %d" % (f, n, m, algo, t, acc[0], acc[1], acc[2], pars[0], pars[1])

def main():
    f = open("results.txt")

    x = [list() for _ in xrange(3)]
    y = [list() for _ in xrange(3)]

    for l in f.readlines():
        data = l.split()
        f = data[0]
        n = int(data[1])
        m = int(data[2])
        algo = data[3]
        t = float(data[4])
        correct = int(data[5])
        total = int(data[6])
        acc = float(data[7])
        pars = int(data[8])
        real_pars = int(data[9])

        x_axis_value = m
        y_axis_value = pars

        if f == "random":
            if algo == "exhaustive":
                x[0].append(x_axis_value)
                y[0].append(y_axis_value)
            elif algo == "greedy":
                x[1].append(x_axis_value)
                y[1].append(y_axis_value)
            elif algo == "hash":
                x[2].append(x_axis_value)
                y[2].append(y_axis_value)

    plt.plot(x[0], y[0])
    plt.plot(x[1], y[1])
    plt.plot(x[2], y[2])

    plt.ylim(0, 16)

    plt.legend(['Exhaustive', 'Greedy', 'Hash'], loc='upper right')
    plt.title('Parsimony vs. Number of SNPs')
    plt.xlabel('Number of SNPs (m)')
    plt.ylabel('Parsimony')

    plt.show()


if __name__ == '__main__':
    main()
