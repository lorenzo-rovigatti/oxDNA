import argparse
from os import path
from typing import Dict
from oxDNA_analysis_tools.UTILS.RyeReader import strand_describe
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings


# Based on vrna_db_from_ptable from ViennaRNA
# Converted by ChatGPT with a lot of help from Erik
def pairs2db(n_bases:int, pairs:Dict[int, int]) -> str:
    """
    Convert a dictionary of paired nucleotides to dot-bracket notation

    Parameters:
        n_bases (int): The total number of bases in the structure
        pairs (dict[int, int]): The IDs of paired nucleotides

    Returns:
        str: The pairs as a dot-bracket string, including any nested bases
    """

    dotbracket = ['.'] * n_bases  # Initialize dotbracket list

    bracket_open_avail = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    bracket_close_avail = ")]}>abcdefghijklmnopqrstuvwxyz"

    stack = []
    bracket_count = 0
    recheck = True

    while recheck:
        recheck = False
        for i in range(n_bases):
            # i is unpaired
            # Note that this requires pairs to be bidirectional!
            if i not in pairs.keys():
                continue

            # i is open bracket
            if pairs[i] > i:
                # Check for crossing pairs
                if len(stack) > 0 and pairs[i] > stack[-1]:
                    recheck = True
                    continue
                else:
                    stack.append(pairs[i])
                    dotbracket[i] = bracket_open_avail[bracket_count]
                    dotbracket[pairs[i]] = bracket_close_avail[bracket_count]
            if len(stack) == 0:
                continue
            if i == stack[-1]:
                # Remove pair from stack and pairdict
                pairs.pop(pairs.pop(i))
                stack.pop()

        # Repeat again with next bracket
        bracket_count += 1

        if bracket_count >= 30:
            log("Not enough bracket types available! Skipping remaining base pairs!", level=1)
            log("If you think this shouldn't be the case, you might not have a 1:1 correspondence in your pairs.", level=1)
            break

    return ''.join(dotbracket)

def cli_parser(prog="forces2pairs"):
    parser = argparse.ArgumentParser(prog = prog, description="Convert a pair file to dot-bracket notation")
    parser.add_argument('topology', type=str, help="Topology file for the structure")
    parser.add_argument('pair_file', type=str, help="The pair file to generate the dot-bracket from")
    parser.add_argument('-o', '--output', type=str, help='If set, print the dot-bracket to a file, otherwise to the screen')
    parser.add_argument('-q', '--quiet', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, help="Don't print 'INFO' messages to stderr")
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    #run system checks
    logger_settings.set_quiet(args.quiet)
    from oxDNA_analysis_tools.config import check
    check(["python"])

    # Parse CLI input
    top_file = args.topology
    _, elems = strand_describe(top_file)
    seq = ''.join([e.btype for e in elems])
    pair_file = args.pair_file

    # Convert the designed pairs into a dict
    with open(pair_file, 'r') as file:
        pairs_txt = file.readlines()

    pairs = {int(p[0]) : int(p[1]) for p in [p.split() for p in pairs_txt]}
    # pairs2db requires a bidirectional list
    for k, v in [list(x) for x in pairs.items()]:
        if v not in pairs.keys():
            pairs[v] = k
    db = pairs2db(len(elems), pairs)

    if args.output:
        out = args.output
    else:
        log("No outfile provided, printing to screen")
        print(seq)
        print(db)
        exit(0)

    with open(out, 'w+') as f:
        f.write(seq+'\n')
        f.write(db+'\n')
        log(f"Wrote dot-bracket to file {out}.")

if __name__ == '__main__':
    main()