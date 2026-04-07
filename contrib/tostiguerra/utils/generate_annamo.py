import importlib
import json
import sys

import numpy as np

# Import the raw matrix calculator
from ftg_SL import interaction_matrix as compute_matrix_raw


def load_terminal_penalties(material):
    if material.upper() == "DNA":
        mod = importlib.import_module("DNA_SL")
    elif material.upper() == "RNA":
        mod = importlib.import_module("RNA_22")
    else:
        raise ValueError("Material must be DNA or RNA")
    return mod.dH_terminal_penalties, mod.dS_terminal_penalties


def write_topology(full_bead_list, filename="topology.top"):
    """
    Write the topology file in the required CUSTOM format:

    N_beads N_strands
    idx type n_neighbors
    neighbor_1 [neighbor_2]
    ...
    (empty line between strands)
    """

    # 1. Identify the indices of 'x' to calculate lengths
    # The list is structured as ['x', seq1..., 'x', seq2..., 'x']
    x_idx = [i for i, e in enumerate(full_bead_list) if e == "x"]

    n_beads_in_each_strand = [
        x_idx[i + 1] - x_idx[i] - 1 for i in range(len(x_idx) - 1)
    ]

    number_of_strands = len(n_beads_in_each_strand)
    total_beads = sum(n_beads_in_each_strand)

    print(
        f"Writing topology (Custom Format): {number_of_strands} strands, {total_beads} beads."
    )

    with open(filename, "w") as topo_file:
        # Header
        topo_file.write("{} {}\n".format(total_beads, number_of_strands))

        # Physical sequential index (0, 1, 2...)
        idx = 0

        # Index to scroll through the full list (includes 'x' to determine TYPE)
        # Start from 1 because index 0 is always 'x'
        current_list_idx = 1

        for strand_id, n in enumerate(n_beads_in_each_strand):

            for i in range(n):
                # The TYPE is the index in the full list (skipping x for counting, but including them in indices)
                # In practice: current_list_idx points to the current bead in the list with x
                bead_type = current_list_idx

                # Logic for neighbors
                if n == 1:
                    # Single-bead strand: no neighbors
                    topo_file.write("{} {} 0\n".format(idx, bead_type))

                elif i == 0:
                    # Start of strand: 1 neighbor (the next)
                    topo_file.write("{} {} 1\n".format(idx, bead_type))
                    topo_file.write("{}\n".format(idx + 1))

                elif i == n - 1:
                    # End of strand: 1 neighbor (the previous)
                    topo_file.write("{} {} 1\n".format(idx, bead_type))
                    topo_file.write("{}\n".format(idx - 1))

                else:
                    # Internal strand: 2 neighbors (previous and next)
                    topo_file.write("{} {} 2\n".format(idx, bead_type))
                    topo_file.write("{} {}\n".format(idx - 1, idx + 1))

                idx += 1
                current_list_idx += 1

            # End of strand:
            # 1. Add empty line as separator
            topo_file.write("\n")
            # 2. Skip the 'x' in the global list to align the type of the next strand
            current_list_idx += 1

    print(f"Topology written to: {filename}")


def main(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)

    material = data.get("material", "DNA")
    salt = data.get("salt_concentration", 1.0)
    strands = validate_and_normalize_strands(data, material)

    print("--- ANNaMo Generator (Preserving Type Gaps) ---")

    # --- 1. Prepare Flat List (with 'x' separators) ---
    # This list determines the matrix indices and types
    full_bead_list = ["x"]

    # Info for corrections: list_index -> (strand_id, strand_len)
    bead_info = {}

    current_idx = 1  # Index 0 is the first 'x'

    for s_id, strand in enumerate(strands):
        s_len = len(strand)
        for bead in strand:
            full_bead_list.append(bead)
            bead_info[current_idx] = (s_id, s_len)
            current_idx += 1

        # Add separator between strands
        full_bead_list.append("x")
        current_idx += 1

    # --- 2. Compute Raw Matrix ---
    print("Computing raw interaction matrix...")
    H_raw, S_raw = compute_matrix_raw(full_bead_list, material, salt)
    dH_term, dS_term = load_terminal_penalties(material)

    # --- 3. Apply Corrections and Write ---
    matrix_lines = []

    num_entries = len(full_bead_list)

    for i in range(num_entries):
        # If it's an 'x', it's not an interacting bead, skip the row
        if full_bead_list[i] == "x":
            continue

        for j in range(i + 1, num_entries):
            # If j is an 'x', skip the column
            if full_bead_list[j] == "x":
                continue

            val_H = H_raw[i][j]
            val_S = S_raw[i][j]

            dH_init, dS_init = 0.0, 0.0
            dH_AT, dS_AT = 0.0, 0.0

            s_id_i, len_i = bead_info[i]
            s_id_j, len_j = bead_info[j]

            # A. Initiation
            if s_id_i != s_id_j:
                min_len = min(len_i, len_j)
                if min_len > 0:
                    dH_init = 0.2 / min_len
                    dS_init = -5.7 / min_len

            # B. Terminal Penalties
            prev_is_x = full_bead_list[i - 1] == "x"
            next_is_x = full_bead_list[i + 1] == "x"
            j_prev_is_x = full_bead_list[j - 1] == "x"
            j_next_is_x = full_bead_list[j + 1] == "x"

            if prev_is_x and j_next_is_x:
                term_pair = full_bead_list[i][0] + full_bead_list[j][-1]
                if term_pair in dH_term:
                    dH_AT = dH_term[term_pair]
                    dS_AT = dS_term[term_pair]

            elif next_is_x and j_prev_is_x:
                term_pair = full_bead_list[j][0] + full_bead_list[i][-1]
                if term_pair in dH_term:
                    dH_AT = dH_term[term_pair]
                    dS_AT = dS_term[term_pair]

            # C. Writing
            final_H = np.round(val_H + dH_init + dH_AT, 1)
            final_S = np.round(val_S + dS_init + dS_AT, 1)

            matrix_lines.append(f"dH[{i}][{j}] = {final_H}\n")
            matrix_lines.append(f"dS[{i}][{j}] = {final_S}\n")

    with open("dHdS_matrix.dat", "w") as fout:
        fout.writelines(matrix_lines)

    print("Interaction matrix written to: dHdS_matrix.dat")

    # We pass the complete list (with the x's) to assign the correct types
    write_topology(full_bead_list, "topology.top")


def validate_and_normalize_strands(data, material):
    if "strands" not in data:
        raise ValueError("Input JSON must contain a 'strands' field")

    strands = data["strands"]
    if not isinstance(strands, list) or len(strands) == 0:
        raise ValueError("'strands' must be a non-empty list")

    material_u = material.upper()
    if material_u not in {"DNA", "RNA"}:
        raise ValueError("Material must be DNA or RNA")

    allowed = {"A", "C", "G", "T"} if material_u == "DNA" else {"A", "C", "G", "U"}
    replacement_from = "U" if material_u == "DNA" else "T"
    replacement_to = "T" if material_u == "DNA" else "U"

    replaced_count = 0
    normalized = []

    for s_idx, strand in enumerate(strands):
        # 1. AUTO-COARSE-GRAINING
        # If the user provides a raw sequence string ("AGCT..."), we split it into beads of 3 nt.
        if isinstance(strand, str):
            bead_size = 3
            # Division in chunks (the last one will take the remainder)
            auto_beads = [
                strand[i : i + bead_size] for i in range(0, len(strand), bead_size)
            ]

            print(
                f"WARNING: Strand {s_idx} provided as raw sequence. "
                f"Automatically coarse-grained into {len(auto_beads)} beads "
                f"(target size: {bead_size} nt).",
                file=sys.stderr,
            )
            # We replace the string with the list of beads
            strand = auto_beads

        elif not isinstance(strand, (list, tuple)):
            raise ValueError(f"Strand {s_idx} must be a string or list/tuple of beads")

        # 2. NORMALIZATION AND VALIDATION (on list of beads)
        strand_seq = []
        for bead in strand:
            if not isinstance(bead, str):
                raise ValueError(f"Strand {s_idx} contains a non-string bead: {bead}")

            bead_u = bead.upper()

            # Replacement T <-> U if necessary
            if replacement_from in bead_u:
                bead_u = bead_u.replace(replacement_from, replacement_to)
                replaced_count += bead.upper().count(replacement_from)

            # Character validation
            for ch in bead_u:
                if ch not in allowed:
                    raise ValueError(
                        f"Invalid base '{ch}' in strand {s_idx} bead '{bead_u}' for material {material_u}"
                    )

            strand_seq.append(bead_u)

        normalized.append(strand_seq)

    if replaced_count > 0:
        print(
            f"WARNING: Replaced {replaced_count} '{replacement_from}' base(s) with '{replacement_to}' for {material_u}.",
            file=sys.stderr,
        )

    return normalized


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_annamo.py input.json")
    else:
        main(sys.argv[1])
