# Sequence Aligner — Smith-Waterman & Needleman-Wunsch with Affine Gap Penalties

A single-file Python CLI tool that performs both **local** (Smith-Waterman) and **global** (Needleman-Wunsch) pairwise sequence alignment using **affine gap penalties**.

---

## Introduction

This tool performs both **local** and **global** pairwise DNA sequence alignment:

- **Local Alignment (Smith-Waterman)** — Finds the highest-scoring sub-region shared between two sequences. Ideal for detecting conserved motifs, functional domains, or homologous regions within otherwise divergent sequences.
- **Global Alignment (Needleman-Wunsch)** — Aligns the two sequences end-to-end, accounting for every character. Best suited when the sequences are of similar length and expected to be related across their full span.

Both algorithms use an **affine gap model**, which distinguishes between *opening* a new gap (expensive) and *extending* an existing gap (cheaper). This produces biologically realistic alignments because real insertions and deletions tend to occur in contiguous stretches.

## Algorithm

### Smith-Waterman (Local Alignment)

The algorithm builds a dynamic-programming scoring matrix where every cell `(i, j)` represents the best alignment score ending at position `i` in Sequence 1 and `j` in Sequence 2. A **local** constraint forces any cell to be at least 0, allowing the alignment to start and end anywhere.

- **Traceback** begins at the cell with the **highest score** in the entire matrix and stops when a cell with value 0 is reached.
- This means only the best-matching sub-region is reported.

### Needleman-Wunsch (Global Alignment)

Like Smith-Waterman, this algorithm fills a scoring matrix — but with two key differences:

1. **No floor at 0** — scores can go negative, because every character in both sequences must be included in the alignment.
2. **Border initialization** — the first row and column are initialized with cumulative gap penalties (not zeros), representing the cost of aligning one sequence entirely to gaps.

- **Traceback** starts at the **bottom-right corner** `(n, m)` and walks back to `(0, 0)`, guaranteeing a full end-to-end alignment.

### Local vs. Global — When to Use Which

| Aspect | Local (Smith-Waterman) | Global (Needleman-Wunsch) |
|--------|----------------------|--------------------------|
| **Use case** | Finding conserved sub-regions | Aligning full-length sequences |
| **Score floor** | 0 (alignment can restart) | None (scores may go negative) |
| **Traceback start** | Highest-scoring cell anywhere | Bottom-right corner (n, m) |
| **Best for** | Divergent sequences with shared motifs | Similar-length, closely related sequences |

### Affine vs. Linear Gaps

| Model | Penalty Formula | Behaviour |
|-------|----------------|-----------|
| **Linear** | `G × L` | Every gap position costs the same — tends to scatter many small gaps. |
| **Affine** | `G_open + (L − 1) × G_extend` | Opening a gap is expensive; extending is cheap — produces fewer, longer gaps. |

Both algorithms use the affine model with **three DP matrices** (Gotoh formulation):

- **M[i][j]** — best score ending with a match / mismatch.
- **Ix[i][j]** — best score ending with a gap in Sequence X.
- **Iy[i][j]** — best score ending with a gap in Sequence Y.

**Recurrences (same for both algorithms):**

```
Ix[i][j] = max(M[i-1][j] + gap_open,  Ix[i-1][j] + gap_extend)
Iy[i][j] = max(M[i][j-1] + gap_open,  Iy[i][j-1] + gap_extend)
M[i][j]  = max(M[i-1][j-1]  + s,
               Ix[i-1][j-1] + s,
               Iy[i-1][j-1] + s)   // plus max(0, ...) for local only
```


## Project Structure

```
alignment/
├── main.py          # All-in-one: parser, both aligners, formatter, CLI
├── input1.fasta     # First DNA sequence (editable)
├── input2.fasta     # Second DNA sequence (editable)
├── parameter.txt    # Scoring parameters (editable)
├── output.txt       # Generated alignment report
├── README.md
└── LICENSE
```

## How to Use

### 1. Change Sequences

Open `input1.fasta` or `input2.fasta` and paste your DNA sequences in standard FASTA format:

```
>My_Sequence_Header
ACGTACGTACGT
GCTAGCTAGCTA
```

### 2. Tune Scoring Parameters

Edit `parameter.txt` to adjust sensitivity:

```
match: 2
mismatch: -1
gap_open: -12
gap_extend: -1
```

- **Higher `gap_open` magnitude** → fewer, longer gaps.
- **Lower `gap_open` magnitude** → more, shorter gaps.

### 3. Run the Tool

```bash
python main.py
```

Results are:
- Printed directly to the terminal.
- Saved automatically to `output.txt` (both local and global alignments).

## Installation

**Requirements:** Python 3.7+ (no external packages needed).

```bash
# Clone the repository
git clone https://github.com/Dhanushsk16/BBL434-assignment2.git
cd BBL434-assignment2

# Run
python main.py
```

## Output Format (MARKX 0)

Running `python main.py` produces **two alignment reports** (both printed to the terminal and saved to `output.txt`):

1. **Local Alignment** — Smith-Waterman result
2. **Global Alignment** — Needleman-Wunsch result

Each report includes:

- A **header** with the alignment score, percent identity, and gap count.
- **60-character alignment blocks** with:
  - `|` for matches
  - ` ` (space) for mismatches or gaps
  - Start and end coordinate markers on each line.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
