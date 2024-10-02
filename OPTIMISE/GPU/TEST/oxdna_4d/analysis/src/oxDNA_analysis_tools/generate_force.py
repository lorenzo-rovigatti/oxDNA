import os
import argparse
from sys import stderr
import oxpy

def cli_parser(prog="generate_force.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Create an external forces file enforcing the current base-pairing arrangement")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('configuration', type=str, nargs=1, help="The configuration to generate the forces from")
    parser.add_argument('-o', '--output', type=str, nargs=1, help='Name of the file to write the forces to. Defaults to forces.txt')
    parser.add_argument('-f', '--pairs', type=str, nargs=1, help='Name of the file to write the designed pairs list to')
    parser.add_argument('-s', '--stiff', type=float, nargs=1, help='Stiffness of the mutual trap')
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))

    args = parser.parse_args()

    #Process command line arguments:
    inputfile = args.inputfile[0]
    conf_file = args.configuration[0]

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else: 
        outfile = "forces.txt"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    if args.pairs:
        pairsfile = args.pairs[0]
    else:
        pairsfile = False

    #-s sets the stiffness of the mutual trap
    if args.stiff:
        stiff = args.stiff[0]
    else:
        stiff = 0.9

    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(inputfile)
        inp["list_type"] = "cells"
        inp["trajectory_file"] = conf_file
        inp["confs_to_analyse"] = str(1)
        inp["analysis_data_output_1"] = '{ \n name = stdout \n print_every = 1e10 \n col_1 = { \n id = my_obs \n type = hb_list \n } \n }'

        if (not inp["use_average_seq"] or inp.get_bool("use_average_seq")) and "RNA" in inp["interaction_type"]:
            print("WARNING: Sequence dependence not set for RNA model, wobble base pairs will be ignored", file=stderr)

        backend = oxpy.analysis.AnalysisBackend(inp)

        # read one conf
        backend.read_next_configuration()

        pairs = backend.config_info().get_observable_by_id("my_obs").get_output_string(backend.config_info().current_step).strip().split('\n')
        
    bonded = {}
    for p in pairs[1:]:
        p = p.split()
        bonded[int(p[0])] = int(p[1])

    lines = []
    pairlines = []
    mutual_trap_template = '{{ \ntype = mutual_trap\nparticle = {}\nstiff = {}\nr0 = 1.2\nref_particle = {}\nPBC=1\n}}\n'
    for key in sorted(bonded):
        from_particle_id = key
        to_particle_id = bonded[key]
        if from_particle_id < to_particle_id:
            if pairsfile:
                pairlines.append("{} {}\n".format(from_particle_id, to_particle_id))
            lines.append(mutual_trap_template.format(from_particle_id,stiff,to_particle_id))
            lines.append(mutual_trap_template.format(to_particle_id,stiff,from_particle_id))

    if pairsfile:
        with open(pairsfile, "w") as file:
            file.writelines(pairlines)
            print("INFO: Wrote pairs to {}".format(pairsfile), file=stderr)

    with open(outfile, "w") as file:
            file.writelines(lines)
            print("INFO: Job finished. Wrote forces to {}".format(outfile), file=stderr)

if __name__ == "__main__":
    main()