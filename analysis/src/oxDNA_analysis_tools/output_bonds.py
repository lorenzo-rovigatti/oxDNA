from typing import List, Optional
import sys
import time
import numpy as np
import argparse
from os import path
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from collections import namedtuple
from oxDNA_analysis_tools.UTILS.data_structures import TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser, get_chunk_size
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_input_parameter
import oxpy

start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "input_file",
                                              "visualize",
                                              "conversion_factor",
                                              "n_potentials",
                                              "fields"])

def parse_header(e_txt:str) -> List[str]:
    e_txt = e_txt.strip('#') # strip leading #
    e_txt = e_txt.split(',')[0] # remove time section
    e_list = e_txt.split(' ')[2:] # The first two are id1 and id2 (unless somebody writes a non-pairwise potential, then it's your problem)
    return e_list

def get_potentials(ctx) -> List[str]:
    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(ctx.input_file)
        inp["list_type"] = "cells"
        inp["trajectory_file"] = ctx.traj_info.path
        inp["log_file"] = "/dev/null"
        inp["analysis_bytes_to_skip"] = str(0)
        inp["confs_to_analyse"] = str(1)
        inp["analysis_data_output_1"] = '{ \n name = stdout \n print_every = 1e10 \n col_1 = { \n id = my_obs \n type = pair_energy \n } \n }'

        backend = oxpy.analysis.AnalysisBackend(inp)

        backend.read_next_configuration()
        e_txt = backend.config_info().get_observable_by_id("my_obs").get_output_string(backend.config_info().current_step).strip().split('\n')
        pot_names = parse_header(e_txt[0])

    return pot_names

def _resolve_fields(pot_names: List[str], field_names: Optional[List[str]]):
    """Resolve field name strings to indices into pot_names.

    Returns:
        Tuple of (indices_or_None, active_names)
        - indices_or_None: list of int indices for numpy indexing, or None for all fields
        - active_names: list of active potential names in the order requested
    """
    if field_names is None or 'all' in [f.lower() for f in field_names]:
        return None, pot_names

    # Build lookup: lowercase name -> index, with 'dh' alias for Debye-Huckel
    field_map = {}
    for i, name in enumerate(pot_names):
        field_map[name.lower()] = i
        if 'debye' in name.lower():
            field_map['dh'] = i

    indices = []
    active_names = []
    for f in field_names:
        f_lower = f.lower()
        if f_lower in field_map:
            idx = field_map[f_lower]
            if idx not in indices:
                indices.append(idx)
                active_names.append(pot_names[idx])
        else:
            log(f"Unrecognized field: '{f}'. Available fields: {', '.join(pot_names)}", level="warning")

    if not indices:
        log("No valid fields specified, using all fields.", level="warning")
        return None, pot_names

    return indices, active_names

# pragma: no cover - coverage.py cannot track execution inside oxpy.Context() if it's in a subprocess
def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    """Parallel compute function: used only when -v is the sole output mode.

    Accumulates per-nucleotide energy sums across the chunk and returns
    an ndarray of shape (nbases, n_potentials) for the callback to add into
    the pre-allocated output array.
    """
    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(ctx.input_file)
        inp["list_type"] = "cells"
        inp["trajectory_file"] = ctx.traj_info.path
        inp["log_file"] = "/dev/null"
        inp["analysis_bytes_to_skip"] = str(ctx.traj_info.idxs[chunk_id*chunk_size].offset)
        inp["confs_to_analyse"] = str(chunk_size)
        inp["analysis_data_output_1"] = '{ \n name = stdout \n print_every = 1e10 \n col_1 = { \n id = my_obs \n type = pair_energy \n } \n }'

        if (not inp["use_average_seq"] or inp.get_bool("use_average_seq")) and "RNA" in inp["interaction_type"]:
            log("Sequence dependence not set for RNA model, wobble base pairs will be ignored", level="warning")

        backend = oxpy.analysis.AnalysisBackend(inp)

        # The 9 energies in oxDNA2 are:
        # 0 fene
        # 1 bexc
        # 2 stck
        # 3 nexc
        # 4 hb
        # 5 cr_stack
        # 6 cx_stack
        # 7 Debye-Huckel <- this one is missing in oxDNA1
        # 8 total
        energies = np.zeros((ctx.top_info.nbases, ctx.n_potentials))
        while backend.read_next_configuration():
            e_txt = backend.config_info().get_observable_by_id("my_obs").get_output_string(backend.config_info().current_step).strip().split('\n')
            for e in e_txt[1:]:
                if e and e[0] != '#':
                    parts = e.split()
                    p = int(parts[0])
                    q = int(parts[1])
                    l = np.array([float(x) for x in parts[2:]]) * ctx.conversion_factor
                    if ctx.fields is not None:
                        l = l[ctx.fields]
                    energies[p] += l
                    energies[q] += l

        return energies

def output_bonds(traj_info:TrajInfo, top_info:TopInfo, inputfile:str,
                 visualize:bool=False, conversion_factor:float=1, ncpus:int=1,
                 fields:Optional[List[str]]=None, traj_view:Optional[str]=None,
                 data_file:Optional[str]=None, produce_print:bool=False,
                 units:str="oxDNA su"):
    """
        Computes the potentials in a trajectory.

        Parameters:
            traj_info (TrajInfo): Information about the trajectory.
            top_info (TopInfo): Information about the topology.
            inputfile (str): Path to the input file.
            visualize (bool): (optional) If True, accumulates mean per-nucleotide energies
                              (returned as ndarray for the caller to write as oxView JSON).
            conversion_factor (float): (optional) Conversion factor for the energies.
                              1 for oxDNA SU, 41.42 for pN nm.
            ncpus (int): (optional) Number of CPUs to use. Ignored when any data-output
                         mode is active (serial processing is required).
            fields (List[str]): (optional) Names of energy fields to include.
                              None means all fields. Use 'dh' as alias for Debye-Huckel.
            traj_view (str): (optional) Path prefix for per-frame oxView JSON files.
                              Written incrementally; one file per active field.
            data_file (str): (optional) Path to write raw pair-interaction data.
                              Written incrementally, one frame at a time.
            produce_print (bool): (optional) If True, print raw pair-interaction data
                              to stdout, one frame at a time.
            units (str): (optional) Unit label used in output file headers/keys.

        Returns:
            Tuple[np.ndarray or None, List[str]]:
                energies: mean per-nucleotide energy array of shape (nbases, n_potentials)
                          when visualize=True (already averaged over all frames), otherwise None.
                active_names: list of energy field names included in the output.
    """

    ctx = ComputeContext(traj_info, top_info, inputfile, visualize, conversion_factor, 0, None)

    # Discover available potential names from a single configuration
    pot_names = get_potentials(ctx)

    # Resolve requested fields to column indices
    field_indices, active_names = _resolve_fields(pot_names, fields)
    n_potentials = len(active_names)

    ctx = ComputeContext(traj_info, top_info, inputfile, visualize, conversion_factor, n_potentials, field_indices)

    # Determine execution path
    produce_data = produce_print or (data_file is not None)
    needs_serial = produce_data or (traj_view is not None)

    if needs_serial and ncpus > 1:
        raise RuntimeError(
            "No flags and the -t, -d, and --force_print flags require serial processing and are "
            "incompatible with -p. Please re-run without the -p flag."
        )

    if not needs_serial:
        # ── Parallel path ────────────────────────────────────────────────────
        # Used when -v is the only requested output.  compute() runs in worker
        # processes and returns per-nucleotide sums for each chunk; the callback
        # accumulates them into a single pre-allocated array.

        energies = np.zeros((top_info.nbases, n_potentials))
        def callback(i, r):
            nonlocal energies
            energies += r

        oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)
        energies /= traj_info.nconfs
        return energies, active_names

    else:
        # ── Serial path ──────────────────────────────────────────────────────
        # Used whenever raw pair data must be written/printed, or when -t is
        # requested.  The oxpy backend runs in the main thread; each frame is
        # processed and output immediately to avoid buffering the full trajectory.

        energies = np.zeros((top_info.nbases, n_potentials)) if visualize else None

        # Open data output stream (file or stdout)
        data_fh = None
        opened_data_file = False
        if data_file:
            data_fh = open(data_file, 'w')
            opened_data_file = True
        elif produce_print:
            data_fh = sys.stdout

        if data_fh:
            data_fh.write('# p1 p2 ' + ' '.join(active_names) + f' ({units})\n')

        # Open per-field file handles for -t output
        # Each entry: [fname, file_handle, is_first_frame]
        traj_fhs = {}
        if traj_view:
            for i, potential in enumerate(active_names):
                if '.json' in traj_view:
                    fname = '.'.join(traj_view.split('.')[:-1]) + "_" + potential + '.json'
                else:
                    fname = traj_view + "_" + potential + '.json'
                fh = open(fname, 'w')
                fh.write('{{\n"{} ({})": [\n'.format(potential, units))
                traj_fhs[i] = [fname, fh, True]

        # Process the entire trajectory frame by frame in the main thread
        with oxpy.Context():
            inp = oxpy.InputFile()
            inp.init_from_filename(inputfile)
            inp["list_type"] = "cells"
            inp["trajectory_file"] = traj_info.path
            inp["log_file"] = "/dev/null"
            inp["analysis_bytes_to_skip"] = str(0)
            inp["confs_to_analyse"] = str(traj_info.nconfs)
            inp["analysis_data_output_1"] = '{ \n name = stdout \n print_every = 1e10 \n col_1 = { \n id = my_obs \n type = pair_energy \n } \n }'

            if (not inp["use_average_seq"] or inp.get_bool("use_average_seq")) and "RNA" in inp["interaction_type"]:
                log("Sequence dependence not set for RNA model, wobble base pairs will be ignored", level="warning")

            backend = oxpy.analysis.AnalysisBackend(inp)

            chunk_size = get_chunk_size()
            log(f"Starting up 1 processes for {traj_info.nconfs // chunk_size} chunks")
            frame_count = 0
            while backend.read_next_configuration():
                e_txt = backend.config_info().get_observable_by_id("my_obs").get_output_string(
                    backend.config_info().current_step).strip().split('\n')

                # Parse this frame's pair interactions
                pairs = []
                comments = []
                for e in e_txt[1:]:
                    if e and e[0] != '#':
                        parts = e.split()
                        p = int(parts[0])
                        q = int(parts[1])
                        l = np.array([float(x) for x in parts[2:]]) * conversion_factor
                        if field_indices is not None:
                            l = l[field_indices]
                        pairs.append((p, q, l))
                    elif e:
                        comments.append(e)

                # Raw pair data output (file or stdout)
                if data_fh:
                    data_fh.write(e_txt[0] + '\n')
                    for comment in comments:
                        data_fh.write(comment + '\n')
                    for p, q, l in pairs:
                        if np.any(l != 0):
                            data_fh.write(f"{p} {q} " + ' '.join([str(v) for v in l]) + '\n')

                # -v: accumulate per-nucleotide sums
                if visualize:
                    for p, q, l in pairs:
                        energies[p] += l
                        energies[q] += l

                # -t: write this frame's per-nucleotide energies to the JSON file
                for i, entry in traj_fhs.items():
                    fname, fh, is_first = entry
                    frame_e = np.zeros(top_info.nbases)
                    for p, q, l in pairs:
                        frame_e[p] += l[i]
                        frame_e[q] += l[i]
                    if not is_first:
                        fh.write(',\n')
                    fh.write('[' + ', '.join([str(x) for x in frame_e]) + ']')
                    entry[2] = False  # mark as no longer the first frame

                frame_count += 1
                if frame_count % chunk_size == 0:
                    print(f"finished {frame_count}/{traj_info.nconfs}", end="\r")
            print()

        # Finalize outputs
        if opened_data_file:
            data_fh.close()
            log(f"Wrote bond data to: {data_file}")

        for i, (fname, fh, _) in traj_fhs.items():
            fh.write('\n]\n}')
            fh.close()
            log(f"Wrote oxView trajectory overlay to: {fname}")

        if visualize:
            energies /= traj_info.nconfs
        return energies, active_names

def cli_parser(prog="output_bonds.py"):
    parser = argparse.ArgumentParser(prog = prog, description="List all the interactions between nucleotides")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-v', '--view', type=str, nargs=1, dest='outfile', help='if you want instead average per-particle energy as an oxView JSON')
    parser.add_argument('-t', '--traj_view', type=str, nargs=1, dest='traj_view', help='Write an oxView trajectory json overlay.')
    parser.add_argument('-f', '--fields', type=str, nargs='+', dest='fields', help='(optional) The fields to print out // save to files.  If not specified, all fields are printed.  Recognized options are: all, fene, bexc, stck, nexc, hb, cr_stack, cx_stack, dh, total')
    parser.add_argument('-d', '--data_file', type=str, nargs=1, dest='data_file', help='(optional) If specified, the output will be written to this file instead of printed to the screen. Default if n_confs > 10')
    parser.add_argument('--force_print', dest='force_print', action='store_const', const=True, default=False, help='(optional) If specified, the output will be printed to the screen even if n_confs > 10')
    parser.add_argument('-p', '--parallel', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-u', '--units', type=str, nargs=1, dest='units', help="(optional) The units of the energy (pNnm or oxDNA)")
    parser.add_argument('-q', '--quiet', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, help="Don't print 'INFO' messages to stderr")
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    logger_settings.set_quiet(args.quiet)
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy", "oxpy"])

    traj_file = args.trajectory[0]
    inputfile = args.inputfile[0]

    top_info, traj_info  = describe(None, traj_file)

    try:
        outfile:str = args.outfile[0]
        visualize = True
    except:
        outfile = ''
        visualize = False

    traj_view = args.traj_view[0] if args.traj_view else None
    data_file = args.data_file[0] if args.data_file else None
    fields = args.fields  # list of field name strings, or None

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    if args.units:
        if args.units[0] == "pNnm":
            units = "pN nm"
            conversion_factor = 41.42
        elif args.units[0] == "oxDNA":
            units = "oxDNA su"
            conversion_factor = 1
        else:
            raise RuntimeError(f"Unrecognized units: {args.units[0]}\n Recognized options are 'pNnm' and 'oxDNA'.")
    else:
        units = "oxDNA su"
        conversion_factor = 1
        log("no units specified, assuming oxDNA su")

    # Auto-default to a data file for large trajectories when no visualization is requested
    if not visualize and not traj_view and not data_file and not args.force_print and traj_info.nconfs > 10:
        data_file = 'output_bonds.txt'
        log("Trajectory has more than 10 configurations, writing bond data to 'output_bonds.txt' instead of printing to screen. Set this filename with the -d flag.", level="warning")

    # Determine whether raw pair data should be written/printed.
    # produce_print is True when data output is needed but no file was specified.
    produce_data = data_file is not None or args.force_print or (not visualize and not traj_view)
    produce_print = produce_data and data_file is None

    needs_serial = produce_data or (traj_view is not None)
    if needs_serial and ncpus > 1:
        raise RuntimeError(
            "No flags and the -t, -d, and --force_print flags require serial processing and are "
            "incompatible with -p. Please re-run without the -p flag."
        )

    energies, pot_names = output_bonds(
        traj_info, top_info, inputfile,
        visualize=visualize,
        conversion_factor=conversion_factor,
        ncpus=ncpus,
        fields=fields,
        traj_view=traj_view,
        data_file=data_file,
        produce_print=produce_print,
        units=units
    )

    # -v output: write mean per-nucleotide energies as oxView JSON files
    if visualize:
        for i, potential in enumerate(pot_names):
            if '.json' in outfile:
                fname = '.'.join(outfile.split('.')[:-1])+"_"+potential+'.json'
            else:
                fname = outfile+"_"+potential+'.json'
            with open(fname, 'w+') as f:
                f.write("{{\n\"{} ({})\" : [".format(potential, units))
                f.write(', '.join([str(x) for x in energies[:,i]]))
                f.write("]\n}")
            log(f"Wrote oxView overlay to: {fname}")

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == "__main__":
    main()
