from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Optional, Tuple, Set

RED = "\033[91m"
RESET = "\033[0m"


@dataclass
class Segment:
    kind: str  # 'paired' or 'unpaired'
    nts: List[int]
    comp_nts: Optional[List[int]] = None  # stored in descending partner order

    @property
    def size(self) -> int:
        return len(self.nts)


@dataclass
class Bead:
    bead_id: int
    nts: List[int]
    kind: str  # 'paired' or 'unpaired'
    partner_bead_id: Optional[int] = None
    partner_nts: Optional[List[int]] = None


PairMap = Dict[int, Optional[int]]
PairTypeMap = Dict[int, str]  # 'real' or 'unbound'
StrandMap = Dict[int, int]  # nt -> 0-based strand index
BaseMap = Dict[int, str]  # nt -> base string


# ── Topology I/O ──────────────────────────────────────────────────────────


def _parse_old_topology(
    lines: List[str],
) -> Tuple[int, int, List[int], List[str], List[List[int]]]:
    header = lines[0].split()
    if len(header) != 2:
        raise ValueError("Old oxDNA topology must start with: N N_strands")
    N = int(header[0])
    n_strands = int(header[1])
    body = [ln.split() for ln in lines[1:]]
    if len(body) < N:
        raise ValueError(
            f"Topology body too short: expected {N} nt lines, found {len(body)}"
        )

    raw_strand: List[int] = []
    raw_base: List[str] = []
    strand_groups: Dict[int, List[int]] = {}

    for nt_idx, parts in enumerate(body[:N]):
        if len(parts) < 2:
            raise ValueError(
                f"Invalid old topology line for nt {nt_idx}: {' '.join(parts)!r}"
            )
        sid = int(parts[0]) - 1
        base = parts[1].upper()
        raw_strand.append(sid)
        raw_base.append(base)
        strand_groups.setdefault(sid, []).append(nt_idx)

    strands = [strand_groups[s] for s in sorted(strand_groups)]
    if len(strands) != n_strands:
        raise ValueError(f"Header says {n_strands} strands, parsed {len(strands)}")
    return N, n_strands, raw_strand, raw_base, strands


def _parse_new_topology(
    lines: List[str],
) -> Tuple[int, int, List[int], List[str], List[List[int]]]:
    header = lines[0].split()
    if len(header) < 3 or header[2] != "5->3":
        raise ValueError("New oxDNA topology must start with: N N_strands 5->3")
    N = int(header[0])
    n_strands = int(header[1])
    seq_lines = [ln.split() for ln in lines[1:]]
    if len(seq_lines) < n_strands:
        raise ValueError(
            f"Topology body too short: expected {n_strands} strand lines, found {len(seq_lines)}"
        )

    raw_strand: List[int] = []
    raw_base: List[str] = []
    strands: List[List[int]] = []
    counter = 0

    def tokenize_sequence(seq: str) -> List[str]:
        out: List[str] = []
        inside = False
        token = ""
        for c in seq:
            if c == "(":
                if inside:
                    raise ValueError(f"Unbalanced parenthesis in sequence {seq!r}")
                inside = True
                token = ""
            elif c == ")":
                if not inside:
                    raise ValueError(f"Unbalanced parenthesis in sequence {seq!r}")
                inside = False
                out.append(token)
            else:
                if inside:
                    token += c
                else:
                    out.append(c)
        if inside:
            raise ValueError(f"Unbalanced parenthesis in sequence {seq!r}")
        return out

    for sid, parts in enumerate(seq_lines[:n_strands]):
        if not parts:
            raise ValueError(f"Empty strand line at strand {sid + 1}")
        seq = tokenize_sequence(parts[0])
        strand_nts: List[int] = []
        for base in seq:
            raw_strand.append(sid)
            raw_base.append(base.upper())
            strand_nts.append(counter)
            counter += 1
        strands.append(strand_nts)

    if counter != N:
        raise ValueError(f"Header says {N} nt, parsed {counter}")
    return N, n_strands, raw_strand, raw_base, strands


def read_oxdna_topology_auto(
    path: str,
) -> Tuple[StrandMap, BaseMap, List[List[int]], Dict[int, int], str]:
    """Read old or new oxDNA topology and normalize to 5'→3'.

    Returns
    -------
    strand_map : {new_nt: strand_id}
    base_map   : {new_nt: base}
    strands    : list of strands, each a list of new nt indices in 5'→3'
    old_to_new : map from the indexing convention implied by the input topology to
                 the normalized internal indexing used by the script
    topo_kind  : 'old' or 'new'
    """
    with open(path, "r", encoding="utf-8") as fh:
        lines = [
            ln.strip() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")
        ]
    if not lines:
        raise ValueError("Empty topology file")

    header = lines[0].split()
    if len(header) == 2:
        topo_kind = "old"
        N, n_strands, raw_strand, raw_base, old_strands = _parse_old_topology(lines)
        # Old format is treated as not yet normalized to 5'→3'.
        # Follow oxDNA's old->new conversion logic: reverse EACH strand, keep strand order.
        old_to_new: Dict[int, int] = {}
        strands: List[List[int]] = []
        counter = 0
        for strand in old_strands:
            new_strand: List[int] = []
            for old_nt in reversed(strand):
                old_to_new[old_nt] = counter
                new_strand.append(counter)
                counter += 1
            strands.append(new_strand)
        if counter != N:
            raise RuntimeError("Internal error while remapping old topology")

        strand_map: StrandMap = {}
        base_map: BaseMap = {}
        for old_nt in range(N):
            new_nt = old_to_new[old_nt]
            strand_map[new_nt] = raw_strand[old_nt]
            base_map[new_nt] = raw_base[old_nt]
        return strand_map, base_map, strands, old_to_new, topo_kind

    if len(header) >= 3 and header[2] == "5->3":
        topo_kind = "new"
        N, n_strands, raw_strand, raw_base, strands = _parse_new_topology(lines)
        strand_map: StrandMap = {i: raw_strand[i] for i in range(N)}
        base_map: BaseMap = {i: raw_base[i] for i in range(N)}
        old_to_new = {i: i for i in range(N)}
        return strand_map, base_map, strands, old_to_new, topo_kind

    raise ValueError(
        "Unrecognized topology format. Expected old 'N N_strands' or new 'N N_strands 5->3'."
    )


# ── Pair I/O ──────────────────────────────────────────────────────────────


def read_pairs(
    path: str, N: int, old_to_new: Optional[Dict[int, int]] = None
) -> Tuple[PairMap, PairTypeMap]:
    """Parse a pairs/contact file.

    Supported formats:
    - full symmetric:
        i
        i j
    - sparse one-entry-per-bond:
        i j
      where unbound nts are omitted and reverse j i is omitted.
    - hybrid:
      any mixture of the above, as long as it is internally consistent.

    Rules:
    - If a reciprocal pair is explicitly present, it must be consistent.
    - If it is missing, it is auto-completed.
    - Missing nts at the end are filled as unbound.
    - Contradictions raise ValueError.
    """
    pairs: PairMap = {}
    pair_types: PairTypeMap = {}

    def set_entry(i: int, j: Optional[int], kind: str, lineno: int) -> None:
        if i in pairs:
            if pairs[i] != j:
                raise ValueError(
                    f"Line {lineno}: inconsistent definition for nt {i}: "
                    f"existing partner={pairs[i]}, new partner={j}"
                )
            if pair_types[i] != kind:
                raise ValueError(
                    f"Line {lineno}: inconsistent pair type for nt {i}: "
                    f"existing type={pair_types[i]!r}, new type={kind!r}"
                )
        else:
            pairs[i] = j
            pair_types[i] = kind

    with open(path, "r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()

            if len(parts) == 1:
                i = int(parts[0])
                set_entry(i, None, "unbound", lineno)

            elif len(parts) == 2:
                i, j = map(int, parts)
                if i == j:
                    raise ValueError(
                        f"Line {lineno}: self-pairing is not allowed: {i} {j}"
                    )
                set_entry(i, j, "real", lineno)
                set_entry(j, i, "real", lineno)

            else:
                raise ValueError(
                    f"Line {lineno}: expected 1 or 2 fields, got: {line!r}"
                )

    for i in range(N):
        if i not in pairs:
            pairs[i] = None
            pair_types[i] = "unbound"

    if old_to_new is not None:
        remapped_pairs: PairMap = {}
        remapped_types: PairTypeMap = {}
        for old_i in range(N):
            new_i = old_to_new[old_i]
            old_j = pairs.get(old_i)
            remapped_pairs[new_i] = old_to_new[old_j] if old_j is not None else None
            remapped_types[new_i] = pair_types.get(old_i, "unbound")
        pairs = remapped_pairs
        pair_types = remapped_types

    for i in range(N):
        j = pairs[i]
        if j is None:
            continue
        if j < 0 or j >= N:
            raise ValueError(f"Partner {j} of nucleotide {i} out of range 0..{N-1}")
        if pairs[j] != i:
            raise ValueError(f"Pairing is not symmetric: {i}->{j}, but {j}->{pairs[j]}")
        if pair_types[j] != pair_types[i]:
            raise ValueError(
                f"Pair type mismatch: {i}->{j} is {pair_types[i]!r}, but {j}->{i} is {pair_types[j]!r}"
            )

    return pairs, pair_types


# ── Small helpers ─────────────────────────────────────────────────────────


def same_strand_nts(nts: List[int], strand_map: Optional[StrandMap]) -> bool:
    if strand_map is None or not nts:
        return True
    sid = strand_map[nts[0]]
    return all(strand_map[nt] == sid for nt in nts)


def same_strand_pair(a: int, b: int, strand_map: Optional[StrandMap]) -> bool:
    if strand_map is None:
        return True
    return strand_map[a] == strand_map[b]


def corresponding_key(val: int, dictionary: Dict[int, List[int]]) -> Optional[int]:
    for k, v in dictionary.items():
        if val in v:
            return k
    return None


def recompute_loop_size_dic(loop_dic: Dict[int, List[int]]) -> Dict[int, List[int]]:
    out: Dict[int, List[int]] = {}
    for k, v in loop_dic.items():
        out.setdefault(len(v), []).append(k)
    return out


# ── Segmentation ──────────────────────────────────────────────────────────


def build_initial_segments(
    pairs: PairMap,
    strand_map: Optional[StrandMap] = None,
) -> Tuple[Dict[int, List[int]], Dict[int, List[int]], Dict[int, List[int]]]:
    """Build initial paired groups and unpaired loops.

    A paired segment can continue only if:
    - current nt continues on the same strand,
    - current partner continues on the same partner strand,
    - partner progression is consecutive (p == prev_p - 1).
    """
    n = len(pairs)
    group_dic: Dict[int, List[int]] = {}
    group_comp_dic: Dict[int, List[int]] = {}
    loop_dic: Dict[int, List[int]] = {}

    prev_state: Optional[str] = None
    current_group_key: Optional[int] = None
    current_loop_key: Optional[int] = None

    for i in range(n):
        p = pairs.get(i)
        state = "paired" if p is not None else "unpaired"

        if state == "paired":
            continue_group = False
            if prev_state == "paired":
                prev_p = pairs[i - 1]
                same_main = same_strand_pair(i, i - 1, strand_map)
                same_comp = prev_p is not None and same_strand_pair(
                    p, prev_p, strand_map
                )
                if prev_p is not None and same_main and same_comp and p == prev_p - 1:
                    continue_group = True

            if continue_group:
                assert current_group_key is not None
                group_dic[current_group_key].append(i)
                group_comp_dic[current_group_key].append(p)
            else:
                current_group_key = i
                group_dic[current_group_key] = [i]
                group_comp_dic[current_group_key] = [p]
            current_loop_key = None

        else:
            continue_loop = prev_state == "unpaired" and same_strand_pair(
                i, i - 1, strand_map
            )
            if continue_loop:
                assert current_loop_key is not None
                loop_dic[current_loop_key].append(i)
            else:
                current_loop_key = i
                loop_dic[current_loop_key] = [i]
            current_group_key = None

        prev_state = state

    return group_dic, group_comp_dic, loop_dic


# ── Pre-processing / local cleanup ────────────────────────────────────────


def merge_symmetric_small_loops(
    group_dic: Dict[int, List[int]],
    group_comp_dic: Dict[int, List[int]],
    loop_dic: Dict[int, List[int]],
    loop_len: int,
    strand_map: Optional[StrandMap] = None,
) -> bool:
    changed = False
    loop_size_dic = recompute_loop_size_dic(loop_dic)
    for k in list(loop_size_dic.get(loop_len, [])):
        if k not in loop_dic:
            continue
        left_anchor = k - 1
        right_anchor = k + loop_len
        kgp = corresponding_key(left_anchor, group_dic)
        kgd = corresponding_key(right_anchor, group_dic)
        if kgp is None or kgd is None:
            continue

        kcp = corresponding_key(kgp, group_comp_dic)
        kcd = corresponding_key(kgd, group_comp_dic)
        if kcp is None or kcd is None:
            continue

        klc = kcp - loop_len
        if klc not in loop_dic:
            continue

        expected_next_comp = group_comp_dic[kgp][-1] - (loop_len + 1)
        if expected_next_comp not in group_comp_dic[kgd]:
            continue

        # Keep merges strand-safe on both sides.
        left_merge_nts = group_dic[kgp] + loop_dic[k] + group_dic[kgd]
        right_merge_nts = group_dic[kcd] + loop_dic[klc] + group_dic[kcp]
        if not same_strand_nts(left_merge_nts, strand_map):
            continue
        if not same_strand_nts(right_merge_nts, strand_map):
            continue

        group_dic[kgp] = left_merge_nts
        group_comp_dic[kgp] = group_comp_dic[kgp] + loop_dic[klc] + group_comp_dic[kgd]
        group_dic.pop(kgd)
        group_comp_dic.pop(kgd)

        group_dic[kcd] = right_merge_nts
        group_comp_dic[kcd] = group_comp_dic[kcd] + loop_dic[k] + group_comp_dic[kcp]
        group_dic.pop(kcp)
        group_comp_dic.pop(kcp)

        loop_dic.pop(k)
        loop_dic.pop(klc)
        changed = True
    return changed


def attach_single_loops_with_complement_support(
    group_dic: Dict[int, List[int]],
    group_comp_dic: Dict[int, List[int]],
    loop_dic: Dict[int, List[int]],
    strand_map: Optional[StrandMap] = None,
) -> bool:
    changed = False
    loop_size_dic = recompute_loop_size_dic(loop_dic)
    for k in list(loop_size_dic.get(1, [])):
        if k not in loop_dic:
            continue

        kgp = corresponding_key(k - 1, group_dic)
        kgd = corresponding_key(k + 1, group_dic)
        if kgp is None or kgd is None:
            continue
        kcp = corresponding_key(kgp, group_comp_dic)
        kcd = corresponding_key(kgd, group_comp_dic)
        if kcp is None or kcd is None:
            continue

        kclp = kcp - 1
        kcld = group_comp_dic[kgd][0] + 1

        if (
            kclp in loop_dic
            and corresponding_key(kclp + 1, group_dic) == kcp
            and k in loop_dic
        ):
            left_new = group_dic[kgp] + [k]
            right_new = loop_dic[kclp] + group_dic[kcp]
            if not same_strand_nts(left_new, strand_map) or not same_strand_nts(
                right_new, strand_map
            ):
                continue
            group_dic[kgp] = left_new
            group_comp_dic[kgp] = group_comp_dic[kgp] + loop_dic[kclp][::-1]

            group_dic[kclp] = right_new
            group_comp_dic[kclp] = [k] + group_comp_dic[kcp]

            group_dic.pop(kcp)
            group_comp_dic.pop(kcp)
            loop_dic.pop(k)
            loop_dic.pop(kclp)
            changed = True
            continue

        if (
            kcld in loop_dic
            and corresponding_key(kcld - 1, group_dic) == kcd
            and k in loop_dic
        ):
            left_new = [k] + group_dic[kgd]
            right_new = group_dic[kcd] + loop_dic[kcld]
            if not same_strand_nts(left_new, strand_map) or not same_strand_nts(
                right_new, strand_map
            ):
                continue
            group_dic[k] = left_new
            group_comp_dic[k] = loop_dic[kcld][::-1] + group_comp_dic[kgd]

            group_dic[kcd] = right_new
            group_comp_dic[kcd] = group_comp_dic[kcd] + [k]

            group_dic.pop(kgd)
            group_comp_dic.pop(kgd)
            loop_dic.pop(k)
            loop_dic.pop(kcld)
            changed = True
    return changed


def absorb_residual_single_loops(
    group_dic: Dict[int, List[int]],
    loop_dic: Dict[int, List[int]],
    strand_map: Optional[StrandMap] = None,
) -> bool:
    changed = False
    loop_size_dic = recompute_loop_size_dic(loop_dic)
    for k in list(loop_size_dic.get(1, [])):
        if k not in loop_dic:
            continue
        kgp = corresponding_key(k - 1, group_dic)
        if kgp is not None:
            new_nts = group_dic[kgp] + loop_dic[k]
            if not same_strand_nts(new_nts, strand_map):
                continue
            group_dic[kgp] = new_nts
            loop_dic.pop(k)
            changed = True
    return changed


def _merge_nt_into_loops(
    nt: int,
    loop_dic: Dict[int, List[int]],
    strand_map: Optional[StrandMap] = None,
    prefer: str = "left",
) -> bool:
    left_key = corresponding_key(nt - 1, loop_dic)
    right_key = corresponding_key(nt + 1, loop_dic)

    if left_key is not None and not same_strand_pair(nt, nt - 1, strand_map):
        left_key = None
    if right_key is not None and not same_strand_pair(nt, nt + 1, strand_map):
        right_key = None

    if left_key is not None and right_key is not None and left_key != right_key:
        merged = loop_dic[left_key] + [nt] + loop_dic[right_key]
        if not same_strand_nts(merged, strand_map):
            return False
        loop_dic[left_key] = merged
        loop_dic.pop(right_key)
        return True

    if prefer == "left":
        if left_key is not None:
            merged = loop_dic[left_key] + [nt]
            if not same_strand_nts(merged, strand_map):
                return False
            loop_dic[left_key] = merged
            return True
        if right_key is not None:
            merged = [nt] + loop_dic[right_key]
            if not same_strand_nts(merged, strand_map):
                return False
            loop_dic[nt] = merged
            loop_dic.pop(right_key)
            return True
    else:
        if right_key is not None:
            merged = [nt] + loop_dic[right_key]
            if not same_strand_nts(merged, strand_map):
                return False
            loop_dic[nt] = merged
            loop_dic.pop(right_key)
            return True
        if left_key is not None:
            merged = loop_dic[left_key] + [nt]
            if not same_strand_nts(merged, strand_map):
                return False
            loop_dic[left_key] = merged
            return True

    loop_dic[nt] = [nt]
    return True


def convert_paired_singletons_to_loops(
    group_dic: Dict[int, List[int]],
    group_comp_dic: Dict[int, List[int]],
    loop_dic: Dict[int, List[int]],
    strand_map: Optional[StrandMap] = None,
) -> bool:
    """Conservative fallback for paired singletons.

    Convert only when both sides are already bracketed by loops on the same strand.
    """
    changed = False
    processed: Set[int] = set()

    for k in list(sorted(group_dic.keys())):
        if k not in group_dic or k in processed:
            continue
        if len(group_dic[k]) != 1:
            continue
        comp_anchor = group_comp_dic[k][0]
        partner_key = corresponding_key(comp_anchor, group_dic)
        if partner_key is None or partner_key == k or partner_key in processed:
            continue
        if len(group_dic[partner_key]) != 1:
            continue

        nt_left = group_dic[k][0]
        nt_right = group_dic[partner_key][0]

        left_has_loop_neighbor = (
            corresponding_key(nt_left - 1, loop_dic) is not None
            and same_strand_pair(nt_left, nt_left - 1, strand_map)
        ) or (
            corresponding_key(nt_left + 1, loop_dic) is not None
            and same_strand_pair(nt_left, nt_left + 1, strand_map)
        )
        right_has_loop_neighbor = (
            corresponding_key(nt_right - 1, loop_dic) is not None
            and same_strand_pair(nt_right, nt_right - 1, strand_map)
        ) or (
            corresponding_key(nt_right + 1, loop_dic) is not None
            and same_strand_pair(nt_right, nt_right + 1, strand_map)
        )
        if not (left_has_loop_neighbor and right_has_loop_neighbor):
            continue

        group_dic.pop(k)
        group_comp_dic.pop(k)
        group_dic.pop(partner_key)
        group_comp_dic.pop(partner_key)

        if not _merge_nt_into_loops(nt_left, loop_dic, strand_map, prefer="left"):
            continue
        if not _merge_nt_into_loops(nt_right, loop_dic, strand_map, prefer="right"):
            continue

        processed.add(k)
        processed.add(partner_key)
        changed = True

    return changed


def finalize_segments(
    group_dic: Dict[int, List[int]],
    group_comp_dic: Dict[int, List[int]],
    loop_dic: Dict[int, List[int]],
) -> List[Segment]:
    used_paired_keys: Set[int] = set()
    segments: List[Segment] = []
    all_keys = sorted(set(group_dic) | set(loop_dic))
    for k in all_keys:
        if k in loop_dic:
            segments.append(Segment(kind="unpaired", nts=loop_dic[k]))
            continue
        if k in group_dic and k not in used_paired_keys:
            partner_key = corresponding_key(group_comp_dic[k][-1], group_dic)
            used_paired_keys.add(k)
            if partner_key is not None:
                used_paired_keys.add(partner_key)
            segments.append(
                Segment(kind="paired", nts=group_dic[k], comp_nts=group_comp_dic[k])
            )
    return segments


# ── Partitioning ──────────────────────────────────────────────────────────


@lru_cache(None)
def partition_length(length: int) -> Optional[Tuple[int, ...]]:
    if length == 0:
        return ()
    if length < 2:
        return None

    candidates = []
    for k in (3, 2, 4):
        tail = partition_length(length - k)
        if tail is not None:
            candidates.append((k,) + tail)
    if not candidates:
        return None

    def score(part: Tuple[int, ...]) -> Tuple[int, int, int, int, int, Tuple[int, ...]]:
        non3 = part.count(2) + part.count(4)
        return (
            -part.count(3),  # maximize 3-nt beads
            non3,  # then minimize non-3 beads
            len(part),  # then fewer total beads
            part.count(2),  # then fewer 2s
            part.count(4),  # then fewer 4s
            part,
        )

    return min(candidates, key=score)


def split_by_partition(seq: List[int], part: Tuple[int, ...]) -> List[List[int]]:
    out: List[List[int]] = []
    pos = 0
    for size in part:
        out.append(seq[pos : pos + size])
        pos += size
    return out


def check_single_real_partner(
    block_a: List[int],
    block_b: List[int],
    pairs: PairMap,
    pair_types: PairTypeMap,
) -> bool:
    aset = set(block_a)
    bset = set(block_b)
    real_partners_a = {
        pairs[i] for i in block_a if pairs[i] is not None and pair_types[i] == "real"
    }
    if not real_partners_a.issubset(bset):
        return False
    real_partners_b = {
        pairs[j] for j in block_b if pairs[j] is not None and pair_types[j] == "real"
    }
    if not real_partners_b.issubset(aset):
        return False
    return True


def check_internal_order(
    block_a: List[int], block_b: List[int], pairs: PairMap
) -> bool:
    bset = set(block_b)
    mapped = [
        pairs[i] for i in sorted(block_a) if pairs[i] is not None and pairs[i] in bset
    ]
    return mapped == sorted(mapped, reverse=True)


def paired_blocks_to_beads(
    segment: Segment,
    pairs: PairMap,
    pair_types: PairTypeMap,
) -> Tuple[List[List[int]], List[List[int]]]:
    assert segment.comp_nts is not None
    left = segment.nts
    right_desc = segment.comp_nts
    L = len(left)

    part = partition_length(L)
    if part is None:
        raise ValueError(
            f"Cannot partition paired segment {left} <-> {right_desc} of length {L} using only 2/3/4 nt beads."
        )

    left_blocks = split_by_partition(left, part)
    right_blocks_desc = split_by_partition(right_desc, part)
    right_blocks = [sorted(blk) for blk in right_blocks_desc[::-1]]

    # Pair left block i with mirrored right block.
    for a, b in zip(left_blocks, right_blocks[::-1]):
        if not check_single_real_partner(a, b, pairs, pair_types):
            raise ValueError(
                f"Invalid paired bead candidate: {a} would not have a unique real partner bead in {b}."
            )
        if not check_internal_order(a, b, pairs):
            raise ValueError(
                f"Order mismatch between candidate paired beads {a} and {b}."
            )

    return [sorted(b) for b in left_blocks], right_blocks


def unpaired_blocks_to_beads(segment: Segment) -> List[List[int]]:
    L = len(segment.nts)
    part = partition_length(L)
    if part is None:
        raise ValueError(
            f"Cannot partition unpaired segment of length {L} using only 2/3/4 nt beads."
        )
    return [sorted(b) for b in split_by_partition(segment.nts, part)]


# ── Bead construction / repair ────────────────────────────────────────────


def build_beads(
    segments: List[Segment], pairs: PairMap, pair_types: PairTypeMap
) -> List[Bead]:
    beads: List[Bead] = []
    bead_id = 0

    for seg in segments:
        if seg.kind == "unpaired":
            for blk in unpaired_blocks_to_beads(seg):
                beads.append(Bead(bead_id=bead_id, nts=blk, kind="unpaired"))
                bead_id += 1
        else:
            left_blocks, right_blocks = paired_blocks_to_beads(seg, pairs, pair_types)
            left_ids = list(range(bead_id, bead_id + len(left_blocks)))
            right_ids = list(
                range(
                    bead_id + len(left_blocks),
                    bead_id + len(left_blocks) + len(right_blocks),
                )
            )

            for local_idx, blk in enumerate(left_blocks):
                partner_idx = len(right_blocks) - 1 - local_idx
                beads.append(
                    Bead(
                        bead_id=left_ids[local_idx],
                        nts=blk,
                        kind="paired",
                        partner_bead_id=right_ids[partner_idx],
                        partner_nts=right_blocks[partner_idx],
                    )
                )

            for local_idx, blk in enumerate(right_blocks):
                partner_idx = len(left_blocks) - 1 - local_idx
                beads.append(
                    Bead(
                        bead_id=right_ids[local_idx],
                        nts=blk,
                        kind="paired",
                        partner_bead_id=left_ids[partner_idx],
                        partner_nts=left_blocks[partner_idx],
                    )
                )

            bead_id += len(left_blocks) + len(right_blocks)

    return beads


def build_nt_to_bead(beads: List[Bead]) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for idx, bead in enumerate(beads):
        for nt in bead.nts:
            out[nt] = idx
    return out


def unique_real_partner_bead_idx(
    bead: Bead,
    nt_to_bead: Dict[int, int],
    pairs: PairMap,
    pair_types: PairTypeMap,
) -> Optional[int]:
    targets = {
        nt_to_bead[pairs[nt]]
        for nt in bead.nts
        if pairs.get(nt) is not None
        and pair_types.get(nt) == "real"
        and pairs[nt] in nt_to_bead
    }
    if len(targets) == 1:
        return next(iter(targets))
    return None


def steal_nt_from_bead(beads: List[Bead], src_idx: int, nt: int) -> None:
    beads[src_idx].nts = sorted(x for x in beads[src_idx].nts if x != nt)


def expand_singleton_beads(
    beads: List[Bead],
    pairs: PairMap,
    pair_types: PairTypeMap,
    strand_map: Optional[StrandMap] = None,
) -> bool:
    changed = False
    nt_to_bead = build_nt_to_bead(beads)
    for idx, bead in enumerate(beads):
        if len(bead.nts) != 1:
            continue
        nt = bead.nts[0]
        if pairs.get(nt) is None:
            continue

        candidate_nts = [nt]
        used: Set[int] = {nt}
        target_strand = strand_map[nt] if strand_map is not None else None
        left = nt - 1
        right = nt + 1
        while len(candidate_nts) < 4:
            took = False
            if (
                left in nt_to_bead
                and left not in used
                and pairs.get(left) is None
                and (strand_map is None or strand_map[left] == target_strand)
            ):
                candidate_nts.insert(0, left)
                used.add(left)
                left -= 1
                took = True
            else:
                left -= 1 if left in nt_to_bead else 10**9

            if len(candidate_nts) >= 4:
                break
            if (
                right in nt_to_bead
                and right not in used
                and pairs.get(right) is None
                and (strand_map is None or strand_map[right] == target_strand)
            ):
                candidate_nts.append(right)
                used.add(right)
                right += 1
                took = True
            else:
                right += 1 if right in nt_to_bead else 10**9

            if not took:
                break

        if len(candidate_nts) <= 1 or not same_strand_nts(candidate_nts, strand_map):
            continue

        for other_nt in list(candidate_nts):
            if other_nt == nt:
                continue
            src_idx = nt_to_bead[other_nt]
            if src_idx == idx:
                continue
            steal_nt_from_bead(beads, src_idx, other_nt)
            changed = True

        beads[idx].nts = sorted(candidate_nts)
        changed = True
        nt_to_bead = build_nt_to_bead(beads)

    beads[:] = [b for b in beads if b.nts]
    return changed


def recompute_partners(
    beads: List[Bead], pairs: PairMap, pair_types: PairTypeMap
) -> None:
    nt_to_bead = build_nt_to_bead(beads)
    for bead in beads:
        bead.partner_bead_id = None
        bead.partner_nts = None

    candidate_map: Dict[int, Optional[int]] = {}
    for idx, bead in enumerate(beads):
        candidate_map[idx] = unique_real_partner_bead_idx(
            bead, nt_to_bead, pairs, pair_types
        )

    for idx, partner_idx in candidate_map.items():
        if partner_idx is None:
            continue
        if candidate_map.get(partner_idx) != idx:
            continue
        beads[idx].partner_bead_id = partner_idx
        beads[idx].partner_nts = sorted(beads[partner_idx].nts)
        beads[idx].kind = "paired"

    for bead in beads:
        if bead.partner_bead_id is None:
            bead.kind = "unpaired"


def repair_noncontiguous_beads(
    beads: List[Bead],
    strand_map: Optional[StrandMap] = None,
) -> bool:
    changed = False
    i = 0
    while i < len(beads):
        bead = beads[i]
        is_noncontig = any(
            bead.nts[k + 1] != bead.nts[k] + 1 for k in range(len(bead.nts) - 1)
        )
        if not is_noncontig:
            i += 1
            continue

        cluster = [i]
        union = sorted(beads[i].nts)
        sizes = [len(beads[i].nts)]
        target_strand = strand_map[union[0]] if strand_map is not None else None
        j = i + 1
        while j < len(beads):
            # Never repair across strand boundaries.
            if strand_map is not None and any(
                strand_map[nt] != target_strand for nt in beads[j].nts
            ):
                break
            new_union = sorted(set(union) | set(beads[j].nts))
            cluster.append(j)
            union = new_union
            sizes.append(len(beads[j].nts))
            j += 1
            if union[-1] - union[0] + 1 == len(union):
                break

        if union[-1] - union[0] + 1 != len(union):
            i += 1
            continue
        if sum(sizes) != len(union):
            i += 1
            continue
        if not same_strand_nts(union, strand_map):
            i += 1
            continue

        pos = 0
        for bead_idx, size in zip(cluster, sizes):
            new_nts = union[pos : pos + size]
            if beads[bead_idx].nts != new_nts:
                beads[bead_idx].nts = new_nts
                changed = True
            pos += size
        i = cluster[-1] + 1
    return changed


def repair_beads(
    beads: List[Bead],
    pairs: PairMap,
    pair_types: PairTypeMap,
    strand_map: Optional[StrandMap] = None,
) -> List[Bead]:
    changed = True
    while changed:
        changed = False
        if expand_singleton_beads(beads, pairs, pair_types, strand_map):
            changed = True
        recompute_partners(beads, pairs, pair_types)
        if repair_noncontiguous_beads(beads, strand_map):
            changed = True
            recompute_partners(beads, pairs, pair_types)

    beads[:] = [b for b in beads if b.nts]
    beads.sort(key=lambda b: b.nts[0])
    recompute_partners(beads, pairs, pair_types)

    # Reassign sequential bead ids after sorting.
    for i, b in enumerate(beads):
        b.bead_id = i
    nt_to_bead = build_nt_to_bead(beads)
    recompute_partners(beads, pairs, pair_types)
    for bead in beads:
        if bead.partner_bead_id is not None:
            # candidate_map stores positional index, so renumber using any real partner nt.
            real_nt = next(
                (
                    nt
                    for nt in bead.nts
                    if pairs.get(nt) is not None and pair_types.get(nt) == "real"
                ),
                None,
            )
            if real_nt is not None:
                bead.partner_bead_id = nt_to_bead[pairs[real_nt]]
                bead.partner_nts = sorted(beads[bead.partner_bead_id].nts)
    return beads


# ── Validation / reporting ────────────────────────────────────────────────


def colorize_nt(nt: int, pair_types: PairTypeMap) -> str:
    status = pair_types.get(nt, "unbound")
    text = str(nt)
    if status == "unbound":
        return f"{RED}{text}{RESET}"
    return text


def fmt_nts(nts: List[int], pair_types: PairTypeMap) -> str:
    return "[" + ", ".join(colorize_nt(nt, pair_types) for nt in nts) + "]"


def print_beads(beads: List[Bead], pair_types: PairTypeMap) -> None:
    print(f"{'BEAD':<8} | {'NUCLEOTIDES':<35} | PARTNER")
    print("-" * 95)
    for bead in beads:
        nt_txt = fmt_nts(bead.nts, pair_types)
        partner_txt = (
            "UNBOUND"
            if bead.partner_bead_id is None
            else f"bead {bead.partner_bead_id:<3} {fmt_nts(bead.partner_nts or [], pair_types)}"
        )
        print(f"bead {bead.bead_id:<3} | {nt_txt:<35} | {partner_txt}")


def bead_is_contiguous(bead: Bead) -> bool:
    return all(bead.nts[i + 1] == bead.nts[i] + 1 for i in range(len(bead.nts) - 1))


def bead_is_single_strand(bead: Bead, strand_map: Optional[StrandMap]) -> bool:
    return same_strand_nts(bead.nts, strand_map)


def compute_warnings(
    beads: List[Bead],
    pairs: PairMap,
    pair_types: PairTypeMap,
    strand_map: Optional[StrandMap] = None,
) -> List[str]:
    warnings: List[str] = []

    singletons = [bead.bead_id for bead in beads if len(bead.nts) == 1]
    if singletons:
        preview = ", ".join(map(str, singletons[:10]))
        suffix = " ..." if len(singletons) > 10 else ""
        warnings.append(f"{len(singletons)} singleton bead(s): {preview}{suffix}")

    noncontig = [bead.bead_id for bead in beads if not bead_is_contiguous(bead)]
    if noncontig:
        preview = ", ".join(map(str, noncontig[:10]))
        suffix = " ..." if len(noncontig) > 10 else ""
        warnings.append(f"{len(noncontig)} non-contiguous bead(s): {preview}{suffix}")

    if strand_map is not None:
        cross_strand = [
            bead.bead_id
            for bead in beads
            if not bead_is_single_strand(bead, strand_map)
        ]
        if cross_strand:
            preview = ", ".join(map(str, cross_strand[:10]))
            suffix = " ..." if len(cross_strand) > 10 else ""
            warnings.append(
                f"{len(cross_strand)} cross-strand bead(s): {preview}{suffix}"
            )

    paired_in_unbound = []
    for bead in beads:
        if bead.partner_bead_id is None:
            bad_nts = [
                nt
                for nt in bead.nts
                if pairs.get(nt) is not None and pair_types.get(nt) == "real"
            ]
            if bad_nts:
                paired_in_unbound.append((bead.bead_id, bad_nts))
    if paired_in_unbound:
        preview = "; ".join(f"bead {bid}: {nts}" for bid, nts in paired_in_unbound[:5])
        suffix = " ..." if len(paired_in_unbound) > 5 else ""
        warnings.append(
            f"{len(paired_in_unbound)} bead(s) marked UNBOUND still contain real-paired nucleotide(s): {preview}{suffix}"
        )

    nt_to_bead = build_nt_to_bead(beads)
    nonunique_real = []
    for bead in beads:
        targets = {
            nt_to_bead[pairs[nt]]
            for nt in bead.nts
            if pairs.get(nt) is not None
            and pair_types.get(nt) == "real"
            and pairs[nt] in nt_to_bead
        }
        if len(targets) > 1:
            nonunique_real.append((bead.bead_id, sorted(targets)))
    if nonunique_real:
        preview = "; ".join(
            f"bead {bid} -> {targets}" for bid, targets in nonunique_real[:5]
        )
        suffix = " ..." if len(nonunique_real) > 5 else ""
        warnings.append(
            f"{len(nonunique_real)} bead(s) violate unique real-partner closure: {preview}{suffix}"
        )

    return warnings


def compute_summary(beads: List[Bead]) -> Dict[str, object]:
    size_counts = {2: 0, 3: 0, 4: 0}
    exact_counts = {(2, 2): 0, (3, 3): 0, (4, 4): 0}
    mismatch_counts = {(2, 3): 0, (3, 4): 0, (2, 4): 0}

    for bead in beads:
        size = len(bead.nts)
        if size in size_counts:
            size_counts[size] += 1

    seen_pairs: Set[Tuple[int, int]] = set()
    for bead in beads:
        if bead.partner_bead_id is None:
            continue
        a, b = bead.bead_id, bead.partner_bead_id
        key = (min(a, b), max(a, b))
        if key in seen_pairs:
            continue
        seen_pairs.add(key)
        sa = len(beads[a].nts)
        sb = len(beads[b].nts)
        pair = tuple(sorted((sa, sb)))
        if pair in exact_counts:
            exact_counts[pair] += 1
        elif pair in mismatch_counts:
            mismatch_counts[pair] += 1

    return {
        "num_beads": len(beads),
        "num_paired_pairs": len(seen_pairs),
        "size_counts": size_counts,
        "exact_counts": exact_counts,
        "mismatch_counts": mismatch_counts,
    }


def print_report(
    beads: List[Bead],
    pairs: PairMap,
    pair_types: PairTypeMap,
    strand_map: Optional[StrandMap] = None,
) -> None:
    warnings = compute_warnings(beads, pairs, pair_types, strand_map)
    summary = compute_summary(beads)

    print("-" * 95)
    print("REPORT")
    print("-" * 95)
    print(f"Total beads: {summary['num_beads']}")
    print(f"Paired bead pairs: {summary['num_paired_pairs']}")
    print(
        "Bead sizes: "
        f"2 nt = {summary['size_counts'][2]}, "
        f"3 nt = {summary['size_counts'][3]}, "
        f"4 nt = {summary['size_counts'][4]}"
    )
    print(
        "Dimensional matches: "
        f"2-2 = {summary['exact_counts'][(2, 2)]}, "
        f"3-3 = {summary['exact_counts'][(3, 3)]}, "
        f"4-4 = {summary['exact_counts'][(4, 4)]}"
    )
    print(
        "Dimensional mismatches: "
        f"2-3 = {summary['mismatch_counts'][(2, 3)]}, "
        f"3-4 = {summary['mismatch_counts'][(3, 4)]}, "
        f"2-4 = {summary['mismatch_counts'][(2, 4)]}"
    )
    if warnings:
        print("Warnings:")
        for w in warnings:
            print(f"  - {w}")


# ── Output for ANNaMo ─────────────────────────────────────────────────────


def bead_sequence(bead: Bead, base_map: BaseMap) -> str:
    return "".join(base_map[nt] for nt in bead.nts)


def build_strand_bead_lists(
    beads: List[Bead], strands: List[List[int]], strand_map: StrandMap
) -> List[List[Bead]]:
    nt_to_order: Dict[int, int] = {}
    for sid, strand_nts in enumerate(strands):
        for pos, nt in enumerate(strand_nts):
            nt_to_order[nt] = pos

    strand_beads: List[List[Bead]] = [[] for _ in range(len(strands))]
    for bead in beads:
        if not same_strand_nts(bead.nts, strand_map):
            raise ValueError(
                f"Cross-strand bead encountered during JSON export: bead {bead.bead_id} {bead.nts}"
            )
        sid = strand_map[bead.nts[0]]
        strand_beads[sid].append(bead)

    for sid in range(len(strands)):
        strand_beads[sid].sort(key=lambda b: nt_to_order[b.nts[0]])
    return strand_beads


def write_annamo_json(
    beads: List[Bead],
    strands: List[List[int]],
    base_map: BaseMap,
    strand_map: StrandMap,
    path: str,
    template: Optional[dict] = None,
) -> None:
    strand_beads = build_strand_bead_lists(beads, strands, strand_map)
    strand_seqs = [[bead_sequence(b, base_map) for b in sb] for sb in strand_beads]

    out = template.copy() if template else {}
    out["strands"] = strand_seqs
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"Written ANNaMo JSON : {path}")


def write_native_bonds(beads: List[Bead], path: str) -> None:
    seen: Set[Tuple[int, int]] = set()
    lines: List[str] = []
    for bead in beads:
        if bead.partner_bead_id is None:
            continue
        a, b = bead.bead_id, bead.partner_bead_id
        key = (min(a, b), max(a, b))
        if key not in seen:
            seen.add(key)
            lines.append(f"{key[0]} {key[1]}")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# bead_i bead_j (0-indexed, sorted pairs)\n")
        fh.write("\n".join(lines) + "\n")
    print(f"Written native bonds: {path}  ({len(lines)} pairs)")


# ── Main ──────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Automatic bead division from nucleotide pairing list, with support for "
            "old and new oxDNA topology formats and multi-strand systems."
        )
    )
    parser.add_argument(
        "pairs_file",
        help=(
            "Pairing file. Per line: 'i' (unbound) or 'i j' (paired). "
            "The file may be complete or sparse: unbound nts can be omitted, "
            "and each pair can be specified once (i j) or twice (i j and j i). "
            "Missing entries are treated as unbound and pairs are made symmetric automatically."
        ),
    )
    parser.add_argument(
        "--topology",
        "-t",
        metavar="TOP",
        default=None,
        help=(
            "oxDNA topology file (old or new format). Required for strand-aware bead division "
            "and for JSON output (including nucleotide/strand information). "
            "If not provided, only the pairing-based division is performed and "
            "native bonds can be written without strand or sequence details."
        ),
    )
    parser.add_argument(
        "--json-out",
        metavar="FILE",
        default=None,
        help="Output path for ANNaMo JSON. Requires --topology.",
    )
    parser.add_argument(
        "--bonds-out",
        metavar="FILE",
        default=None,
        help="Output path for native bead-bead bond pairs.",
    )
    parser.add_argument(
        "--json-template",
        metavar="FILE",
        default=None,
        help="Optional JSON template merged into the ANNaMo output. The 'strands' field is overwritten.",
    )
    args = parser.parse_args()

    strand_map: Optional[StrandMap] = None
    base_map: Optional[BaseMap] = None
    strands: Optional[List[List[int]]] = None
    old_to_new: Optional[Dict[int, int]] = None
    topo_kind: Optional[str] = None
    N: Optional[int] = None

    if args.topology:
        strand_map, base_map, strands, old_to_new, topo_kind = read_oxdna_topology_auto(
            args.topology
        )
        N = len(strand_map)
    else:
        # Fall back to inferring N from the pairs file only.
        max_idx = -1
        with open(args.pairs_file, "r", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                max_idx = max(max_idx, int(parts[0]))
                if len(parts) >= 2:
                    max_idx = max(max_idx, int(parts[1]))
        if max_idx < 0:
            raise ValueError("Empty pairs file")
        N = max_idx + 1

    assert N is not None
    pairs, pair_types = read_pairs(args.pairs_file, N, old_to_new=old_to_new)

    group_dic, group_comp_dic, loop_dic = build_initial_segments(pairs, strand_map)

    while True:
        changed = False
        changed |= merge_symmetric_small_loops(
            group_dic, group_comp_dic, loop_dic, 1, strand_map
        )
        changed |= merge_symmetric_small_loops(
            group_dic, group_comp_dic, loop_dic, 2, strand_map
        )
        changed |= attach_single_loops_with_complement_support(
            group_dic, group_comp_dic, loop_dic, strand_map
        )
        changed |= absorb_residual_single_loops(group_dic, loop_dic, strand_map)
        changed |= convert_paired_singletons_to_loops(
            group_dic, group_comp_dic, loop_dic, strand_map
        )
        if not changed:
            break

    segments = finalize_segments(group_dic, group_comp_dic, loop_dic)
    beads = build_beads(segments, pairs, pair_types)
    beads = repair_beads(beads, pairs, pair_types, strand_map)

    print_beads(beads, pair_types)
    print_report(beads, pairs, pair_types, strand_map)

    if args.bonds_out:
        write_native_bonds(beads, args.bonds_out)

    if args.json_out:
        if not args.topology:
            print("ERROR: --json-out requires --topology.", file=sys.stderr)
            sys.exit(1)
        assert base_map is not None and strands is not None and strand_map is not None
        template: Optional[dict] = None
        if args.json_template:
            with open(args.json_template, "r", encoding="utf-8") as fh:
                template = json.load(fh)
        write_annamo_json(beads, strands, base_map, strand_map, args.json_out, template)


if __name__ == "__main__":
    main()
