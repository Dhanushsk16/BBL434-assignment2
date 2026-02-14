"""
main.py — Sequence Aligner with Affine Gap Penalties.

A single-file CLI tool that:
    1. Parses FASTA sequences from input1.fasta & input2.fasta
    2. Reads scoring parameters from parameter.txt
    3. Performs local alignment  (Smith-Waterman) and
       global alignment (Needleman-Wunsch) using three DP
       matrices (M, Ix, Iy) for affine gap handling
    4. Outputs results in MARKX 0 format to terminal and output.txt
"""

import os
import math

# ======================================================================
#  Sentinel
# ======================================================================
NEG_INF = -math.inf


# ======================================================================
#  PARSING
# ======================================================================

def parse_fasta(filepath: str) -> str:
    """
    Parse a FASTA file and return the concatenated sequence as an
    uppercase string.  Header lines (starting with '>') are ignored.
    """
    sequence_parts = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            sequence_parts.append(line.upper())
    return "".join(sequence_parts)


def parse_parameters(filepath: str) -> dict:
    """
    Parse the parameter configuration file.

    Expected format (colon-delimited):
        match: 2
        mismatch: -1
        gap_open: -12
        gap_extend: -1

    Returns a dict with integer values for each key.
    """
    params = {}
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            key, value = line.split(":")
            params[key.strip()] = int(value.strip())
    return params


# ======================================================================
#  SMITH-WATERMAN ALIGNER  (Affine Gap Penalties)
# ======================================================================

def _score(a: str, b: str, match: int, mismatch: int) -> int:
    """Return the substitution score for two characters."""
    return match if a == b else mismatch


def smith_waterman(seq1: str, seq2: str, match: int, mismatch: int,
                   gap_open: int, gap_extend: int) -> dict:
    """
    Perform Smith-Waterman local alignment with affine gap penalties.

    Three DP matrices:
        M[i][j]  — best score ending with seq1[i-1] aligned to seq2[j-1]
        Ix[i][j] — best score ending with a gap in seq2 (deletion)
        Iy[i][j] — best score ending with a gap in seq1 (insertion)

    Recurrences:
        Ix[i][j] = max(M[i-1][j] + gap_open,  Ix[i-1][j] + gap_extend)
        Iy[i][j] = max(M[i][j-1] + gap_open,  Iy[i][j-1] + gap_extend)
        M[i][j]  = max(0,
                       M[i-1][j-1]  + score(seq1[i-1], seq2[j-1]),
                       Ix[i-1][j-1] + score(seq1[i-1], seq2[j-1]),
                       Iy[i-1][j-1] + score(seq1[i-1], seq2[j-1]))

    Total gap penalty = gap_open + (L - 1) * gap_extend

    Returns a dict with:
        score, seq1_aligned, seq2_aligned,
        seq1_start, seq1_end, seq2_start, seq2_end   (all 1-based)
    """
    n = len(seq1)
    m = len(seq2)

    # Initialise matrices
    M  = [[0]       * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Iy = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    # Traceback pointers for each matrix
    # trace_M:  0 = stop, 1 = from M diag, 2 = from Ix diag, 3 = from Iy diag
    # trace_Ix: 0 = opened from M, 1 = extended from Ix
    # trace_Iy: 0 = opened from M, 1 = extended from Iy
    trace_M  = [[0] * (m + 1) for _ in range(n + 1)]
    trace_Ix = [[0] * (m + 1) for _ in range(n + 1)]
    trace_Iy = [[0] * (m + 1) for _ in range(n + 1)]

    best_score = 0
    best_i, best_j = 0, 0
    best_matrix = 'M'

    # ---- Fill -----------------------------------------------------------
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Ix: gap in seq2 (consume seq1[i-1], skip seq2)
            open_ix   = M[i - 1][j] + gap_open
            extend_ix = Ix[i - 1][j] + gap_extend
            if open_ix >= extend_ix:
                Ix[i][j] = open_ix
                trace_Ix[i][j] = 0       # opened from M
            else:
                Ix[i][j] = extend_ix
                trace_Ix[i][j] = 1       # extended from Ix

            # Iy: gap in seq1 (consume seq2[j-1], skip seq1)
            open_iy   = M[i][j - 1] + gap_open
            extend_iy = Iy[i][j - 1] + gap_extend
            if open_iy >= extend_iy:
                Iy[i][j] = open_iy
                trace_Iy[i][j] = 0       # opened from M
            else:
                Iy[i][j] = extend_iy
                trace_Iy[i][j] = 1       # extended from Iy

            # M: match / mismatch — diagonal from ALL three matrices
            s = _score(seq1[i - 1], seq2[j - 1], match, mismatch)
            from_M  = M[i - 1][j - 1]  + s
            from_Ix = Ix[i - 1][j - 1] + s
            from_Iy = Iy[i - 1][j - 1] + s

            best_m = max(0, from_M, from_Ix, from_Iy)
            M[i][j] = best_m

            # Traceback pointer for M
            # Prefer Ix/Iy sources over M when tied (longer alignment paths)
            if best_m == 0:
                trace_M[i][j] = 0         # local: stop
            elif best_m == from_Ix:
                trace_M[i][j] = 2         # from Ix diagonal
            elif best_m == from_Iy:
                trace_M[i][j] = 3         # from Iy diagonal
            elif best_m == from_M:
                trace_M[i][j] = 1         # from M diagonal

            # Track global max across ALL three matrices
            # Use >= to prefer later (longer) alignments when scores tie
            for val, mat in [(M[i][j], 'M'), (Ix[i][j], 'Ix'), (Iy[i][j], 'Iy')]:
                if val >= best_score and val > 0:
                    best_score = val
                    best_i, best_j = i, j
                    best_matrix = mat

    # ---- Traceback ------------------------------------------------------
    aligned1, aligned2 = [], []
    i, j = best_i, best_j
    current_matrix = best_matrix

    while i > 0 and j > 0:
        if current_matrix == 'M':
            if M[i][j] <= 0:
                break
            direction = trace_M[i][j]
            if direction == 0:
                break
            # Diagonal move — align seq1[i-1] with seq2[j-1]
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            if direction == 1:
                current_matrix = 'M'
            elif direction == 2:
                current_matrix = 'Ix'
            else:
                current_matrix = 'Iy'
            i -= 1
            j -= 1

        elif current_matrix == 'Ix':
            # Gap in seq2 — align seq1[i-1] to '-'
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            if trace_Ix[i][j] == 1:
                current_matrix = 'Ix'     # extended
            else:
                current_matrix = 'M'      # opened from M
            i -= 1

        elif current_matrix == 'Iy':
            # Gap in seq1 — align '-' to seq2[j-1]
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            if trace_Iy[i][j] == 1:
                current_matrix = 'Iy'     # extended
            else:
                current_matrix = 'M'      # opened from M
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return {
        "score":        best_score,
        "seq1_aligned": "".join(aligned1),
        "seq2_aligned": "".join(aligned2),
        "seq1_start":   i + 1,
        "seq1_end":     best_i,
        "seq2_start":   j + 1,
        "seq2_end":     best_j,
    }


# ======================================================================
#  NEEDLEMAN-WUNSCH ALIGNER  (Global, Affine Gap Penalties)
# ======================================================================

def needleman_wunsch(seq1: str, seq2: str, match: int, mismatch: int,
                     gap_open: int, gap_extend: int) -> dict:
    """
    Perform Needleman-Wunsch global alignment with affine gap penalties.

    Three DP matrices (Gotoh formulation):
        M[i][j]  — best score ending with seq1[i-1] aligned to seq2[j-1]
        Ix[i][j] — best score ending with a gap in seq2 (deletion)
        Iy[i][j] — best score ending with a gap in seq1 (insertion)

    Recurrences:
        Ix[i][j] = max(M[i-1][j] + gap_open,  Ix[i-1][j] + gap_extend)
        Iy[i][j] = max(M[i][j-1] + gap_open,  Iy[i][j-1] + gap_extend)
        M[i][j]  = max(M[i-1][j-1]  + score(seq1[i-1], seq2[j-1]),
                       Ix[i-1][j-1] + score(seq1[i-1], seq2[j-1]),
                       Iy[i-1][j-1] + score(seq1[i-1], seq2[j-1]))

    Unlike Smith-Waterman, there is NO floor at 0 — scores may go negative.
    Traceback starts from the bottom-right corner (n, m) and ends at (0, 0).

    Returns a dict with:
        score, seq1_aligned, seq2_aligned,
        seq1_start, seq1_end, seq2_start, seq2_end   (all 1-based)
    """
    n = len(seq1)
    m = len(seq2)

    # Initialise matrices
    M  = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Iy = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0  # Base case: empty alignment

    # Initialise first row (gaps in seq1)
    for j in range(1, m + 1):
        Iy[0][j] = gap_open + (j - 1) * gap_extend
        M[0][j]  = Iy[0][j]

    # Initialise first column (gaps in seq2)
    for i in range(1, n + 1):
        Ix[i][0] = gap_open + (i - 1) * gap_extend
        M[i][0]  = Ix[i][0]

    # Traceback pointers for each matrix
    trace_M  = [[0] * (m + 1) for _ in range(n + 1)]
    trace_Ix = [[0] * (m + 1) for _ in range(n + 1)]
    trace_Iy = [[0] * (m + 1) for _ in range(n + 1)]

    # ---- Fill -----------------------------------------------------------
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Ix: gap in seq2
            open_ix   = M[i - 1][j] + gap_open
            extend_ix = Ix[i - 1][j] + gap_extend
            if open_ix >= extend_ix:
                Ix[i][j] = open_ix
                trace_Ix[i][j] = 0
            else:
                Ix[i][j] = extend_ix
                trace_Ix[i][j] = 1

            # Iy: gap in seq1
            open_iy   = M[i][j - 1] + gap_open
            extend_iy = Iy[i][j - 1] + gap_extend
            if open_iy >= extend_iy:
                Iy[i][j] = open_iy
                trace_Iy[i][j] = 0
            else:
                Iy[i][j] = extend_iy
                trace_Iy[i][j] = 1

            # M: match / mismatch — diagonal from ALL three matrices
            s = _score(seq1[i - 1], seq2[j - 1], match, mismatch)
            from_M  = M[i - 1][j - 1]  + s
            from_Ix = Ix[i - 1][j - 1] + s
            from_Iy = Iy[i - 1][j - 1] + s

            best_m = max(from_M, from_Ix, from_Iy)  # No floor at 0
            M[i][j] = best_m

            # Traceback pointer for M
            if best_m == from_Ix:
                trace_M[i][j] = 2
            elif best_m == from_Iy:
                trace_M[i][j] = 3
            else:
                trace_M[i][j] = 1

    # ---- Determine best ending matrix at (n, m) ------------------------
    best_score = max(M[n][m], Ix[n][m], Iy[n][m])
    if best_score == Ix[n][m]:
        current_matrix = 'Ix'
    elif best_score == Iy[n][m]:
        current_matrix = 'Iy'
    else:
        current_matrix = 'M'

    # ---- Traceback from (n, m) to (0, 0) --------------------------------
    aligned1, aligned2 = [], []
    i, j = n, m

    while i > 0 or j > 0:
        if current_matrix == 'M':
            if i == 0 or j == 0:
                break
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            direction = trace_M[i][j]
            if direction == 1:
                current_matrix = 'M'
            elif direction == 2:
                current_matrix = 'Ix'
            else:
                current_matrix = 'Iy'
            i -= 1
            j -= 1

        elif current_matrix == 'Ix':
            if i == 0:
                break
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            if trace_Ix[i][j] == 1:
                current_matrix = 'Ix'
            else:
                current_matrix = 'M'
            i -= 1

        elif current_matrix == 'Iy':
            if j == 0:
                break
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            if trace_Iy[i][j] == 1:
                current_matrix = 'Iy'
            else:
                current_matrix = 'M'
            j -= 1

    # Handle remaining characters
    while i > 0:
        aligned1.append(seq1[i - 1])
        aligned2.append('-')
        i -= 1
    while j > 0:
        aligned1.append('-')
        aligned2.append(seq2[j - 1])
        j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return {
        "score":        best_score,
        "seq1_aligned": "".join(aligned1),
        "seq2_aligned": "".join(aligned2),
        "seq1_start":   1,
        "seq1_end":     n,
        "seq2_start":   1,
        "seq2_end":     m,
    }


# ======================================================================
#  MARKX 0 FORMATTER
# ======================================================================

def _percent_identity(seq1_aligned: str, seq2_aligned: str) -> float:
    """Return the percentage of positions that are identical matches."""
    length = len(seq1_aligned)
    if length == 0:
        return 0.0
    matches = sum(1 for a, b in zip(seq1_aligned, seq2_aligned) if a == b)
    return (matches / length) * 100.0


def format_markx0(seq1_aligned: str, seq2_aligned: str, score: int,
                  seq1_start: int, seq2_start: int,
                  title: str = "Local Alignment (Smith-Waterman, Affine Gap Penalties)") -> str:
    """
    Format an alignment result in MARKX 0 style.

    • Header with score, identity, gaps, and alignment length.
    • 60-character alignment blocks.
    • '|' for matches, ' ' for mismatches / gaps.
    • 1-based coordinate markers on each line.
    """
    identity = _percent_identity(seq1_aligned, seq2_aligned)
    length   = len(seq1_aligned)
    matches  = sum(1 for a, b in zip(seq1_aligned, seq2_aligned) if a == b)
    gaps     = seq1_aligned.count('-') + seq2_aligned.count('-')

    lines = []
    lines.append("=" * 60)
    lines.append(title)
    lines.append("=" * 60)
    lines.append("")
    lines.append(f"  Alignment Score : {score}")
    lines.append(f"  Identity        : {matches}/{length} ({identity:.1f}%)")
    lines.append(f"  Gaps            : {gaps}/{length}")
    lines.append(f"  Alignment Length: {length}")
    lines.append("")
    lines.append("-" * 60)
    lines.append("")

    # Build middle (match) line
    middle_str = "".join('|' if a == b else ' '
                         for a, b in zip(seq1_aligned, seq2_aligned))

    block_size = 60
    pos1 = seq1_start
    pos2 = seq2_start

    for start in range(0, length, block_size):
        end = min(start + block_size, length)

        chunk1 = seq1_aligned[start:end]
        chunk_m = middle_str[start:end]
        chunk2 = seq2_aligned[start:end]

        non_gap1 = sum(1 for ch in chunk1 if ch != '-')
        non_gap2 = sum(1 for ch in chunk2 if ch != '-')

        end_pos1 = pos1 + non_gap1 - 1
        end_pos2 = pos2 + non_gap2 - 1

        label_width = 8
        lines.append(f"Seq1 {pos1:<{label_width}} {chunk1} {end_pos1}")
        lines.append(f"     {'':<{label_width}} {chunk_m}")
        lines.append(f"Seq2 {pos2:<{label_width}} {chunk2} {end_pos2}")
        lines.append("")

        pos1 = end_pos1 + 1
        pos2 = end_pos2 + 1

    lines.append("-" * 60)
    return "\n".join(lines)


def write_output(text: str, filepath: str) -> None:
    """Write the formatted report to a file."""
    with open(filepath, "w") as f:
        f.write(text)


# ======================================================================
#  MAIN
# ======================================================================

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))

    fasta1_path = os.path.join(base_dir, "input1.fasta")
    fasta2_path = os.path.join(base_dir, "input2.fasta")
    param_path  = os.path.join(base_dir, "parameter.txt")
    output_path = os.path.join(base_dir, "output.txt")

    # 1. Parse inputs
    seq1   = parse_fasta(fasta1_path)
    seq2   = parse_fasta(fasta2_path)
    params = parse_parameters(param_path)

    print(f"Sequence 1 length: {len(seq1)}")
    print(f"Sequence 2 length: {len(seq2)}")
    print(f"Parameters       : {params}")
    print()

    # 2. Local Alignment (Smith-Waterman)
    local_result = smith_waterman(
        seq1, seq2,
        match=params["match"],
        mismatch=params["mismatch"],
        gap_open=params["gap_open"],
        gap_extend=params["gap_extend"],
    )

    local_report = format_markx0(
        local_result["seq1_aligned"],
        local_result["seq2_aligned"],
        local_result["score"],
        local_result["seq1_start"],
        local_result["seq2_start"],
        title="Local Alignment (Smith-Waterman, Affine Gap Penalties)",
    )

    # 3. Global Alignment (Needleman-Wunsch)
    global_result = needleman_wunsch(
        seq1, seq2,
        match=params["match"],
        mismatch=params["mismatch"],
        gap_open=params["gap_open"],
        gap_extend=params["gap_extend"],
    )

    global_report = format_markx0(
        global_result["seq1_aligned"],
        global_result["seq2_aligned"],
        global_result["score"],
        global_result["seq1_start"],
        global_result["seq2_start"],
        title="Global Alignment (Needleman-Wunsch, Affine Gap Penalties)",
    )

    # 4. Combine and output
    full_report = local_report + "\n\n" + global_report

    print(full_report)

    write_output(full_report, output_path)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
