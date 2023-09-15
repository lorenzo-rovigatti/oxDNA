import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from sys import stderr

def cli_parser(prog="plot_energy.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Plot oxDNA energy files")
    parser.add_argument('energy', nargs='+', help='Energy files to plot')
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The name to save the graph file to')
    parser.add_argument('-f', '--format', metavar='<histogram/trajectory/both>', nargs=1, help='Output format for the graphs.  Defaults to histogram.  Options are \"histogram\", \"trajectory\", and \"both\"')
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))    
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check
    check(["python", "numpy", "matplotlib"])

    #get file name
    energy_files = args.energy

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else:
        outfile = 'energy.png'

    #-f defines which type of graph to produce
    hist = False
    line = False
    if args.format:
        if "histogram" in args.format:
            hist = True
        if "trajectory" in args.format:
            line = True
        if "both" in args.format:
            hist = line = True
        if hist == line == False:
            raise RuntimeError("Unrecognized graph format\nAccepted formats are \"histogram\", \"trajectory\", and \"both\"")
    else:
        print("INFO: No graph format specified, defaulting to histogram", file=stderr)
        hist = True

    all_times = []
    all_energies = []
    for efile in energy_files:
        times = []
        energies = []
        with open(efile, 'r') as f:
            l = f.readline()
            while l:
                times.append(float(l.split()[0]))
                energies.append(float(l.split()[1]))
                l = f.readline()
        all_times.append(times)
        all_energies.append(energies)

    names = ["1", "2", "3"]

    if outfile and hist == True:
        if line == True:
            out = outfile[:outfile.find(".")]+"_hist"+outfile[outfile.find("."):]
        else:
            out = outfile
    
        bins = np.linspace(min([min(e) for e in all_energies]), max([max(e) for e in all_energies]), 40)
    
        artists = []
        for i,elist in enumerate(all_energies):
            a = plt.hist(elist, bins, weights=np.ones(len(elist)) / len(elist),  alpha=0.3, label=names[i], histtype=u'stepfilled', edgecolor='k')
            artists.append(a)
        plt.legend(labels=names)
        plt.xlabel("Energy per particle (SU)")
        plt.ylabel("Normalized frequency")
        if outfile:
            print("INFO: Saving histogram to {}".format(out), file=stderr)
            plt.tight_layout()
            plt.savefig(out)
        else:
            plt.show()

    #make a trajectory plot
    if outfile and line == True:
        if hist == True:
            plt.clf()
            out = outfile[:outfile.find(".")]+"_traj"+outfile[outfile.find("."):]
        else:
            out = outfile
        
        artists = []
        for tlist,elist in zip(all_times, all_energies):
            a = plt.plot(tlist, elist, alpha = 0.5)
            artists.append(a)
        plt.legend(labels=names)
        plt.xlabel("Time (SU)")
        plt.ylabel("Energy (SU)")
        if outfile:
            print("INFO: Saving line plot to {}".format(out), file=stderr)
            plt.tight_layout()
            plt.savefig(out)
        else:
            plt.show()
    
if __name__ == '__main__':
    main()    

