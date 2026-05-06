"""
Minimal mmCIF reader/writer for oxDNA analysis tools.

Handles coordinate files: _atom_site, _entity, _entity_poly, _entity_src_gen.
No external dependencies beyond the Python standard library.
"""
from __future__ import annotations
from typing import List, Dict, Tuple, TextIO


# ── Tokenizer ─────────────────────────────────────────────────────────────────

def tokenize_cif(text: str) -> List[str]:
    """Return a flat list of CIF string tokens.

    Handles bare tokens, single-quoted strings, double-quoted strings,
    semicolon-delimited multi-line values, and #-comments.
    """
    tokens: List[str] = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c in ' \t\r\n':
            i += 1
        elif c == '#':                              # comment to EOL
            while i < n and text[i] != '\n':
                i += 1
        elif c == ';' and (i == 0 or text[i - 1] == '\n'):   # semicolon block
            i += 1
            start = i
            while i < n:
                if text[i] == '\n' and i + 1 < n and text[i + 1] == ';':
                    break
                i += 1
            tokens.append(text[start:i].strip())
            i += 2                                  # skip \n;
        elif c == "'":                              # single-quoted string
            i += 1
            start = i
            while i < n:
                if text[i] == "'" and (i + 1 >= n or text[i + 1] in ' \t\r\n'):
                    break
                i += 1
            tokens.append(text[start:i])
            i += 1
        elif c == '"':                              # double-quoted string
            i += 1
            start = i
            while i < n and text[i] != '"':
                i += 1
            tokens.append(text[start:i])
            i += 1
        else:                                       # bare token
            start = i
            while i < n and text[i] not in ' \t\r\n':
                i += 1
            tokens.append(text[start:i])
    return tokens


# ── Reader ────────────────────────────────────────────────────────────────────

def parse_atom_site(cif_text: str) -> List[Dict[str, str]]:
    """
    Extract _atom_site rows from mmCIF text.

    Returns a list of dicts keyed by column name without the '_atom_site.'
    prefix (e.g. ``'Cartn_x'``, ``'label_atom_id'``).

    Raises ValueError if no _atom_site loop_ is present.
    """
    tokens = tokenize_cif(cif_text)
    n = len(tokens)
    i = 0
    columns: List[str] = []
    data_start = -1

    while i < n:
        tok = tokens[i]
        if tok.lower() == 'loop_':
            i += 1
            loop_cols: List[str] = []
            while i < n and tokens[i].startswith('_'):
                loop_cols.append(tokens[i])
                i += 1
            if any(c.startswith('_atom_site.') for c in loop_cols):
                columns = [
                    c[len('_atom_site.'):] if c.startswith('_atom_site.') else c
                    for c in loop_cols
                ]
                data_start = i
                break
            # skip data rows of a non-_atom_site loop
            while i < n and not tokens[i].startswith('_') \
                    and tokens[i].lower() != 'loop_' \
                    and not tokens[i].lower().startswith('data_'):
                i += 1
        else:
            i += 1

    if data_start == -1:
        raise ValueError("No _atom_site loop_ block found in mmCIF data")

    ncols = len(columns)
    rows: List[Dict[str, str]] = []
    i = data_start
    while i + ncols <= n:
        if tokens[i].lower() in ('loop_', 'stop_') \
                or tokens[i].startswith('_') \
                or tokens[i].lower().startswith('data_'):
            break
        rows.append({columns[j]: tokens[i + j] for j in range(ncols)})
        i += ncols

    return rows


# ── Internal helpers ──────────────────────────────────────────────────────────

def _quote_if_needed(s: str) -> str:
    """Wrap *s* in double quotes only when necessary for valid CIF output."""
    if not s or s in ('.', '?'):
        return s
    # Must quote if starts with a reserved char or contains whitespace
    if s[0] in ('"', "'", '$', ';', '_', '[') or any(c in s for c in ' \t\n\r'):
        escaped = s.replace('"', r'\"')
        return f'"{escaped}"'
    return s


def _element_from_name(atom_name: str) -> str:
    """Derive the element symbol from an atom name (e.g. ``"C4'"`` → ``'C'``)."""
    for c in atom_name:
        if c.isalpha():
            return c.upper()
    return 'X'


def _next_chain_id(chain_id: str) -> str:
    """Advance the chain identifier through A-Z, a-z, 0-9, cycling at overflow."""
    if chain_id == 'Z':
        return 'a'
    if chain_id == 'z':
        return '1'
    if chain_id == '9':
        return 'A'
    return chr(ord(chain_id) + 1)


# ── Writer ────────────────────────────────────────────────────────────────────

def write_data_block_header(out: TextIO, block_name: str = 'oat_output') -> None:
    """Write the ``data_<block_name>`` line."""
    print(f'data_{block_name}', file=out)
    print('#', file=out)


def write_atom_site_block(
    all_chain_data: List[Tuple[str, List[List[Dict]], int]],
    out: TextIO,
) -> None:
    """
    Write the ``_atom_site`` loop_ block.

    Parameters
    ----------
    all_chain_data:
        List of ``(chain_id, strand_pdb, entity_id)`` tuples.
        *strand_pdb* is a ``List[List[Dict]]`` as produced by
        ``_build_nucleic_strand_atoms`` — outer list over residues, inner list
        over atoms, each atom a dict with keys ``'name'``, ``'residue_name'``,
        ``'pos'``, ``'bfactor'``.
    out:
        Open writable file handle.
    """
    columns = [
        'group_PDB', 'id', 'type_symbol',
        'label_atom_id', 'label_comp_id',
        'label_asym_id', 'label_entity_id', 'label_seq_id',
        'Cartn_x', 'Cartn_y', 'Cartn_z',
        'B_iso_or_equiv',
        'auth_asym_id', 'auth_seq_id',
    ]

    print('loop_', file=out)
    for col in columns:
        print(f'_atom_site.{col}', file=out)

    atom_id = 1
    for chain_id, strand_pdb, entity_id in all_chain_data:
        for seq_id, residue_atoms in enumerate(strand_pdb, 1):
            for atom in residue_atoms:
                name = atom['name']
                resn = atom['residue_name']
                x, y, z = atom['pos']
                bfac = atom['bfactor']
                elem = _element_from_name(name)
                print(
                    f"ATOM {atom_id} {elem} {_quote_if_needed(name)} "
                    f"{_quote_if_needed(resn)} {chain_id} {entity_id} {seq_id} "
                    f"{x:.3f} {y:.3f} {z:.3f} "
                    f"{bfac:.2f} "
                    f"{chain_id} {seq_id}",
                    file=out,
                )
                atom_id += 1
    print('#', file=out)


def write_entity_block(
    entity_rows: List[Tuple[int, str]],  # (entity_id, description)
    out: TextIO,
) -> None:
    """Write the ``_entity`` loop_ block."""
    print('loop_', file=out)
    for item in ['id', 'type', 'src_method', 'pdbx_description']:
        print(f'_entity.{item}', file=out)
    for eid, desc in entity_rows:
        print(f"{eid} polymer syn {_quote_if_needed('oxDNA coarse-grained model: ' + desc)}", file=out)
    print('#', file=out)


def write_entity_poly_block(
    entity_poly_rows: List[Tuple[int, str, str]],  # (entity_id, poly_type, sequence)
    out: TextIO,
) -> None:
    """Write the ``_entity_poly`` loop_ block."""
    print('loop_', file=out)
    for item in ['entity_id', 'type', 'pdbx_seq_one_letter_code', 'nstd_linkage', 'nstd_monomer']:
        print(f'_entity_poly.{item}', file=out)
    for eid, ptype, seq in entity_poly_rows:
        print(f"{eid} {_quote_if_needed(ptype)} {_quote_if_needed(seq)} no no", file=out)
    print('#', file=out)


