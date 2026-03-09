"""Translate nucleotide IDs between two oxDNA topology files.

Given a reference topology and a comparison topology, maps each nucleotide ID
from the reference to the corresponding ID in the comparison file.  Strands are
matched by sequence similarity (using :mod:`difflib.SequenceMatcher`) and then
aligned base-by-base to produce the full per-nucleotide ID mapping.

Both old-style and new-style topology formats are supported in any combination,
via :func:`~oxDNA_analysis_tools.UTILS.RyeReader.strand_describe`.
"""

import argparse
from os import path
from typing import Dict
from difflib import SequenceMatcher
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from oxDNA_analysis_tools.UTILS.RyeReader import strand_describe


def IDconvert(ref_top: str, cmp_top: str,
              min_phase1_score: float = 0.75,
              min_block_size: int = 8) -> Dict[int, int]:
    """
    Translate nucleotide IDs between two oxDNA topology files.

    Strands are matched by sequence similarity and then aligned base-by-base
    using :class:`difflib.SequenceMatcher` to produce a per-nucleotide ID
    mapping.  Both old-style and new-style topology formats are supported in
    any combination.

    Parameters:
        ref_top (str): Path to the reference ``.top`` topology file.
        cmp_top (str): Path to the comparison ``.top`` topology file.
        min_phase1_score (float): Minimum similarity score (0–1) required to
            accept a strand match in phase 1.  Strands whose best match falls
            below this threshold are left for phase 2.  Default: 0.75.
        min_block_size (int): Minimum consecutive-base block length required
            in phase 2 (nicking, ligation, and circular wrap-around) matching.
            Shorter blocks are treated as spurious.  Default: 8.

    Returns:
        Dict[int, int]: Maps each matched reference nucleotide ID to the
        corresponding nucleotide ID in the comparison topology.  Nucleotides
        present in the reference but absent from the comparison are not
        included in the returned dict.
    """
    ref_system, _ = strand_describe(ref_top)
    cmp_system, _ = strand_describe(cmp_top)

    # Build chain representations in 5'→3' order using Strand.get_5_3().
    ref_strands = {s.id: s.get_5_3() for s in ref_system}
    cmp_strands = {s.id: s.get_5_3() for s in cmp_system}

    ref_total = sum(len(c) for c in ref_strands.values())
    cmp_total = sum(len(c) for c in cmp_strands.values())
    log(f"Structure 1 (reference): {ref_total} nucleotides in {len(ref_strands)} strands")
    log(f"Structure 2 (compare):   {cmp_total} nucleotides in {len(cmp_strands)} strands")

    # Phase 1: match each reference strand to the most similar strand in the compare file.
    # All pairwise scores are computed upfront and then greedily assigned from highest to
    # lowest, ensuring globally best pairs are matched first rather than letting an
    # earlier ref strand steal a match that belongs to a later one.
    all_scores = []
    for ref_id, ref_chain in ref_strands.items():
        ref_seq = ''.join(b for _, b in ref_chain)
        for cmp_id, cmp_chain in cmp_strands.items():
            score = SequenceMatcher(None, ref_seq, ''.join(b for _, b in cmp_chain)).ratio()
            all_scores.append((score, ref_id, cmp_id))

    all_scores.sort(reverse=True)
    used_ref, used_cmp = set(), set()
    strand_map, strand_scores = {}, {}
    for score, ref_id, cmp_id in all_scores:
        if ref_id in used_ref or cmp_id in used_cmp:
            continue
        if score >= min_phase1_score:
            strand_map[ref_id] = cmp_id
            strand_scores[ref_id] = score
            used_ref.add(ref_id)
            used_cmp.add(cmp_id)
    # Unmatched ref strands: no compare strand scored above threshold
    for ref_id in ref_strands:
        if ref_id not in strand_map:
            strand_map[ref_id] = None
            strand_scores[ref_id] = 0.0

    # Align matched strand pairs base-by-base to build the ref_id → cmp_id mapping
    id_map: Dict[int, int] = {}
    for ref_sid, cmp_sid in strand_map.items():
        if cmp_sid is None:
            continue  # no confident match found; nucleotides stay absent from id_map
        ref_chain = ref_strands[ref_sid]
        cmp_chain = cmp_strands[cmp_sid]
        matcher = SequenceMatcher(None, [b for _, b in ref_chain],
                                        [b for _, b in cmp_chain], autojunk=False)
        for a, b, size in matcher.get_matching_blocks():
            for i in range(size):
                id_map[ref_chain[a + i][0]] = cmp_chain[b + i][0]

    log("Strand matching (phase 1):")
    n_skipped = sum(1 for sid, score in strand_scores.items()
                    if strand_map[sid] is None and score > 0)
    for ref_sid, cmp_sid in strand_map.items():
        ref_len = len(ref_strands[ref_sid])
        cmp_len = len(cmp_strands[cmp_sid]) if cmp_sid is not None else 0
        score   = strand_scores[ref_sid]
        log(f"  ref strand {ref_sid} ({ref_len} nt) -> cmp strand {cmp_sid} ({cmp_len} nt), similarity {score:.3f}")
    if n_skipped:
        log(f"  ({n_skipped} ref strands below score threshold {min_phase1_score:.2f},"
            f" leaving their nucleotides for phase 2)")

    # Phase 2: resolve nicking, ligation, and circular strand wrap-around.
    # After 1-to-1 matching some nucleotides may be unmatched because:
    #   - Nicking:   one ref strand was split into several cmp strands.
    #   - Ligation:  several ref strands were joined into one cmp strand.
    #   - Circular:  a circular ref strand was nicked into a single linear cmp strand at a
    #                different position, rotating the sequence; the wrap-around portion
    #                appears at the start of the ref chain and end of the cmp chain.
    #   - Score threshold: a ref strand scored below the phase-1 threshold.
    # SequenceMatcher is run on only the unmatched ref nucleotides (unmatched_sub) against
    # only the unassigned cmp nucleotides (avail), which:
    #   (a) avoids monotonicity conflicts that arise when aligning full rotated chains, and
    #   (b) eliminates the need for a separate skip-already-matched guard.
    # Two constraints prevent spurious mappings:
    #   1. Only matching blocks of ≥ min_block_size consecutive bases are accepted.
    #   2. Each cmp nucleotide ID can be assigned to at most one ref nucleotide (bijection).
    used_cmp_ids = set(id_map.values())
    unmatched = ({nid for chain in ref_strands.values() for nid, _ in chain}
                 - set(id_map.keys()))
    if unmatched:
        n_before = len(unmatched)
        log(f"Phase 2 (nicking/ligation): resolving {n_before} unmatched nucleotides...")
        for ref_id, ref_chain in ref_strands.items():
            unmatched_sub = [(nid, b) for nid, b in ref_chain if nid in unmatched]
            if not unmatched_sub:
                continue
            for cmp_id, cmp_chain in cmp_strands.items():
                avail = [(nid, b) for nid, b in cmp_chain if nid not in used_cmp_ids]
                if not avail:
                    continue
                matcher = SequenceMatcher(None, [b for _, b in unmatched_sub],
                                               [b for _, b in avail], autojunk=False)
                for a, b_off, size in matcher.get_matching_blocks():
                    if size < min_block_size:
                        continue
                    for i in range(size):
                        ref_nid = unmatched_sub[a + i][0]
                        cmp_nid = avail[b_off + i][0]
                        if ref_nid in unmatched and cmp_nid not in used_cmp_ids:
                            id_map[ref_nid] = cmp_nid
                            unmatched.discard(ref_nid)
                            used_cmp_ids.add(cmp_nid)
        log(f"  resolved {n_before - len(unmatched)} nucleotides via nicking/ligation")

    deleted = ref_total - len(id_map)
    added   = cmp_total - len(set(id_map.values()))
    log(f"Deleted bases (in ref, absent in cmp): {deleted}")
    log(f"Added bases   (in cmp, absent in ref): {added}")

    return id_map


def cli_parser(prog="IDconvert") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=prog,
        description="Translate nucleotide IDs between two oxDNA topology files.")
    parser.add_argument('reference', type=str,
        help='Reference topology file whose nucleotide IDs are provided as input.')
    parser.add_argument('compare', type=str,
        help='Comparison topology file to map the IDs into.')
    parser.add_argument('-i', '--index', metavar='index_file', dest='index_file', nargs=1,
        help='File containing a space- or comma-separated list of nucleotide IDs from the '
             'reference topology to translate.  If omitted, all nucleotides in the reference are used.')
    parser.add_argument('--min-score', type=float, default=0.75, dest='min_score',
        help='Minimum similarity score (0–1) required to accept a strand match in phase 1 '
             '(default: 0.75).')
    parser.add_argument('--min-block', type=int, default=8, dest='min_block',
        help='Minimum consecutive-base block length for phase 2 nicking/ligation matching '
             '(default: 8).')
    parser.add_argument('-q', '--quiet', dest='quiet', action='store_const',
        const=True, default=False, help='Silence non-essential output.')
    return parser


def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()
    logger_settings.set_quiet(args.quiet)

    if args.index_file:
        with open(args.index_file[0]) as f:
            indexes = f.readline().replace(',', ' ').split()
        try:
            query_ids = [int(i) for i in indexes]
        except ValueError:
            raise RuntimeError("The index file must be a space- or comma-separated list of integers.")
    else:
        ref_system, _ = strand_describe(args.reference)
        query_ids = [m.id for s in ref_system for m in s.monomers]

    id_map = IDconvert(args.reference, args.compare,
                       min_phase1_score=args.min_score,
                       min_block_size=args.min_block)
    new_ids = [str(id_map.get(old_id, 'DELETED')) for old_id in query_ids]
    print(' '.join(new_ids))


if __name__ == '__main__':
    main()
