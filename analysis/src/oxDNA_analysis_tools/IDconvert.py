import sys
from difflib import SequenceMatcher

# usage: oat IDconvert reference.top compare.top ids.txt
# ids.txt should contain comma-separated nucleotide IDs (from the reference file)
# Supports both old and new .top formats, in any combination.
#
# Old format header: "<N_nucs> <N_strands>"
#   Body: one line per nucleotide: "<strand_id> <base> <n3> <n5>"
#
# New format header: "<N_nucs> <N_strands> 5->3"
#   Body: one line per strand: "<sequence> id=<n> type=<DNA/RNA> circular=<bool>"
#   Nucleotide IDs are assigned sequentially (0-indexed) across strands in file order.

def _is_new_format(header_line):
    return len(header_line.strip().split()) >= 3

def _read_strands_old(filename):
    """Read old-format .top. Returns {strand_id: [(nuc_id, base), ...]} in 5'->3' order."""
    nucs = {}
    with open(filename) as f:
        f.readline()  # skip header
        for nuc_id, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            strand, base, n3, n5 = line.split()
            nucs[nuc_id] = (int(strand), base, int(n3), int(n5))

    strands = {}
    for nuc_id, (strand, base, n3, n5) in nucs.items():
        if n5 == -1:  # 5' end, follow chain toward 3'
            chain, cur = [], nuc_id
            while cur != -1:
                chain.append((cur, nucs[cur][1]))
                cur = nucs[cur][2]
            strands[strand] = chain
    return strands

def _read_strands_new(filename):
    """Read new-format .top. Returns {strand_id: [(nuc_id, base), ...]} in 5'->3' order."""
    strands = {}
    nuc_id = 0
    with open(filename) as f:
        f.readline()  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            # format: <sequence> id=<n> type=<DNA/RNA> circular=<bool>
            parts = line.split()
            seq = parts[0]
            strand_id = None
            for part in parts[1:]:
                if part.startswith("id="):
                    strand_id = int(part[3:])
                    break
            if strand_id is None:
                continue
            chain = [(nuc_id + i, base) for i, base in enumerate(seq)]
            nuc_id += len(seq)
            strands[strand_id] = chain
    return strands

def read_strands(filename):
    """Auto-detect topology format and return {strand_id: [(nuc_id, base), ...]} in 5'->3' order."""
    with open(filename) as f:
        header = f.readline()
    if _is_new_format(header):
        return _read_strands_new(filename)
    else:
        return _read_strands_old(filename)

def main():
    reference_file, compare_file, id_file = sys.argv[1], sys.argv[2], sys.argv[3]
    query_ids = [int(x) for x in open(id_file).read().split(',')]

    ref_strands = read_strands(reference_file)
    cmp_strands = read_strands(compare_file)

    # match each reference strand to the most similar strand in the compare file
    used = set()
    strand_map = {}   # ref_strand_id -> cmp_strand_id
    strand_scores = {}  # ref_strand_id -> similarity score
    for ref_id, ref_chain in ref_strands.items():
        ref_seq = ''.join(b for _, b in ref_chain)
        best_id, best_score = None, 0
        for cmp_id, cmp_chain in cmp_strands.items():
            if cmp_id in used:
                continue
            score = SequenceMatcher(None, ref_seq, ''.join(b for _, b in cmp_chain)).ratio()
            if score > best_score:
                best_score, best_id = score, cmp_id
        strand_map[ref_id] = best_id
        strand_scores[ref_id] = best_score
        used.add(best_id)

    # align matched strand pairs base-by-base to get ref_id -> cmp_id
    id_map = {}
    for ref_sid, cmp_sid in strand_map.items():
        ref_chain = ref_strands[ref_sid]
        cmp_chain = cmp_strands[cmp_sid]
        matcher = SequenceMatcher(None, [b for _, b in ref_chain],
                                        [b for _, b in cmp_chain], autojunk=False)
        for a, b, size in matcher.get_matching_blocks():
            for i in range(size):
                id_map[ref_chain[a+i][0]] = cmp_chain[b+i][0]

    # look up and print as comma-separated new IDs
    new_ids = [str(id_map.get(old_id, 'DELETED')) for old_id in query_ids]

    # summary printed to stderr so stdout remains machine-parseable
    ref_total = sum(len(chain) for chain in ref_strands.values())
    cmp_total = sum(len(chain) for chain in cmp_strands.values())
    deleted = ref_total - len(id_map)
    added   = cmp_total - len(id_map)
    print("", file=sys.stderr)
    print(f"Structure 1 (reference): {ref_total} nucleotides in {len(ref_strands)} strands", file=sys.stderr)
    print(f"Structure 2 (compare):   {cmp_total} nucleotides in {len(cmp_strands)} strands", file=sys.stderr)
    print("", file=sys.stderr)
    print(f"Deleted bases (in ref, absent in cmp): {deleted}", file=sys.stderr)
    print(f"Added bases   (in cmp, absent in ref): {added}", file=sys.stderr)
    print(f"\nStrand matching:", file=sys.stderr)
    for ref_sid, cmp_sid in strand_map.items():
        ref_len = len(ref_strands[ref_sid])
        cmp_len = len(cmp_strands[cmp_sid]) if cmp_sid is not None else 0
        score   = strand_scores[ref_sid]
        print(f"  ref strand {ref_sid} ({ref_len} nt) -> cmp strand {cmp_sid} ({cmp_len} nt), similarity {score:.3f}", file=sys.stderr)
    print("")
    print(','.join(new_ids))
    print("")
