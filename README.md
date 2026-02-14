# Sequence Aligner — Smith-Waterman & Needleman-Wunsch with Affine Gap Penalties

A single-file Python CLI tool that performs both **local** (Smith-Waterman) and **global** (Needleman-Wunsch) pairwise sequence alignment using **affine gap penalties**.

---

## Introduction

This tool finds the best local alignment between two DNA sequences. Unlike global alignment (Needleman-Wunsch), local alignment identifies the most similar sub-regions, making it ideal for detecting conserved motifs and functional domains within divergent sequences.

The aligner uses an **affine gap model**, which distinguishes between *opening* a new gap (expensive) and *extending* an existing gap (cheaper). This produces biologically realistic alignments because real insertions and deletions tend to occur in contiguous stretches.

## Algorithm

### Smith-Waterman (Local Alignment)

The algorithm builds a dynamic-programming scoring matrix. Every cell `(i, j)` represents the best alignment score ending at position `i` in Sequence 1 and position `j` in Sequence 2. A **local** constraint forces any cell to be at least 0, effectively allowing the alignment to start anywhere.

### Affine vs. Linear Gaps

| Model | Penalty Formula | Behaviour |
|-------|----------------|-----------|
| **Linear** | `G × L` | Every gap position costs the same — tends to scatter many small gaps. |
| **Affine** | `G_open + (L − 1) × G_extend` | Opening a gap is expensive; extending is cheap — produces fewer, longer gaps. |

The affine model uses **three matrices**:

- **M** — best score ending with a match / mismatch.
- **Ix** — best score ending with a gap in Sequence X.
- **Iy** — best score ending with a gap in Sequence Y.

## Project Structure

```
alignment/
├── main.py          # All-in-one: parser, aligner, formatter, and CLI
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
- Saved automatically to `output.txt`.

## Installation

**Requirements:** Python 3.7+ (no external packages needed).

```bash
# Clone the repository
git clone <repo-url>
cd alignment

# Run
python main.py
```

## Output Format (MARKX 0)

The report includes:

- A **header** with the alignment score, percent identity, and gap count.
- **60-character alignment blocks** with:
  - `|` for matches
  - ` ` (space) for mismatches or gaps
  - Start and end coordinate markers on each line.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
