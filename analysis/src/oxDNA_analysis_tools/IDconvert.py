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


def IDconvert(ref_top: str, cmp_top: str) -> Dict[int, int]:
    """
    Translate nucleotide IDs between two oxDNA topology files.

    Strands are matched by sequence similarity and then aligned base-by-base
    using :class:`difflib.SequenceMatcher` to produce a per-nucleotide ID
    mapping.  Both old-style and new-style topology formats are supported in
    any combination.

    Parameters:
        ref_top (str): Path to the reference ``.top`` topology file.
        cmp_top (str): Path to the comparison ``.top`` topology file.

    Returns:
        Dict[int, int]: Maps each matched reference nucleotide ID to the
        corresponding nucleotide ID in the comparison topology.  Nucleotides
        present in the reference but absent from the comparison are not
        included in the returned dict.
    """
    ref_system, _ = strand_describe(ref_top)
    cmp_system, _ = strand_describe(cmp_top)

    # Build chain representations in 5'→3' order.
    # New-format strands store monomers in 5'→3' order already.
    # Old-format strands store monomers in file order (typically 3'→5'), so we
    # must traverse connectivity from the 5' end (n5 is None) following n3 links.
    def _build_chain(strand):
        if not strand.is_old():
            return [(m.id, str(m.btype)) for m in strand.monomers]
        m_by_id = {m.id: m for m in strand.monomers}
        start = next((m for m in strand.monomers if m.n5 is None), strand.monomers[0])
        chain, cur, seen = [], start, set()
        while cur is not None and cur.id not in seen:
            chain.append((cur.id, str(cur.btype)))
            seen.add(cur.id)
            cur = m_by_id.get(cur.n3)
        return chain

    def _chains(system):
        return {s.id: _build_chain(s) for s in system}

    ref_strands = _chains(ref_system)
    cmp_strands = _chains(cmp_system)

    ref_total = sum(len(c) for c in ref_strands.values())
    cmp_total = sum(len(c) for c in cmp_strands.values())
    log(f"Structure 1 (reference): {ref_total} nucleotides in {len(ref_strands)} strands")
    log(f"Structure 2 (compare):   {cmp_total} nucleotides in {len(cmp_strands)} strands")

    # Phase 1: match each reference strand to the most similar strand in the compare file.
    # Only accept a match if the similarity score meets the minimum threshold; this prevents
    # deleted or substantially-changed ref strands from spuriously "stealing" cmp strands
    # that should belong to other ref strands.
    _MIN_PHASE1_SCORE = 0.75
    used = set()
    strand_map = {}
    strand_scores = {}
    for ref_id, ref_chain in ref_strands.items():
        ref_seq = ''.join(b for _, b in ref_chain)
        best_id, best_score = None, 0.0
        for cmp_id, cmp_chain in cmp_strands.items():
            if cmp_id in used:
                continue
            score = SequenceMatcher(None, ref_seq, ''.join(b for _, b in cmp_chain)).ratio()
            if score > best_score:
                best_score, best_id = score, cmp_id
        if best_score >= _MIN_PHASE1_SCORE:
            strand_map[ref_id] = best_id
            strand_scores[ref_id] = best_score
            used.add(best_id)
        else:
            strand_map[ref_id] = None   # below threshold — likely deleted/substantially changed
            strand_scores[ref_id] = best_score

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
        log(f"  ({n_skipped} ref strands below score threshold {_MIN_PHASE1_SCORE:.2f},"
            f" leaving their nucleotides for phase 2)")

    # Phase 2: resolve nicking and ligation.
    # After 1-to-1 matching some nucleotides may be unmatched because:
    #   - Nicking:   one ref strand was split into several cmp strands.
    #   - Ligation:  several ref strands were joined into one cmp strand.
    #   - Score threshold: a ref strand scored below the phase-1 threshold because it
    #     was heavily nicked (phase 2 maps its nucleotides to the correct fragments).
    # Two constraints are enforced to prevent spurious mappings:
    #   1. Only matching blocks of at least _MIN_BLOCK_SIZE consecutive bases are accepted.
    #      Shorter blocks are likely spurious commonality between deleted strands and
    #      unrelated compare strands, not genuine nicking/ligation matches.
    #   2. Each cmp nucleotide ID can be assigned to at most one ref nucleotide (bijection).
    #      This prevents deleted-strand nucleotides from double-booking cmp slots that are
    #      already taken by legitimate phase-1 assignments.
    _MIN_BLOCK_SIZE = 8
    used_cmp_ids = set(id_map.values())   # cmp nucleotide IDs already assigned in phase 1
    unmatched = ({nid for chain in ref_strands.values() for nid, _ in chain}
                 - set(id_map.keys()))
    if unmatched:
        n_before = len(unmatched)
        log(f"Phase 2 (nicking/ligation): resolving {n_before} unmatched nucleotides...")
        for ref_id, ref_chain in ref_strands.items():
            if not any(nid in unmatched for nid, _ in ref_chain):
                continue
            for cmp_id, cmp_chain in cmp_strands.items():
                if strand_map.get(ref_id) == cmp_id:
                    continue  # already aligned in phase 1
                matcher = SequenceMatcher(None, [b for _, b in ref_chain],
                                                [b for _, b in cmp_chain], autojunk=False)
                for a, b, size in matcher.get_matching_blocks():
                    if size < _MIN_BLOCK_SIZE:
                        continue  # skip short spurious blocks
                    for i in range(size):
                        ref_nid = ref_chain[a + i][0]
                        cmp_nid = cmp_chain[b + i][0]
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

    id_map = IDconvert(args.reference, args.compare)
    new_ids = [str(id_map.get(old_id, 'DELETED')) for old_id in query_ids]
    print(' '.join(new_ids))


if __name__ == '__main__':
    main()
