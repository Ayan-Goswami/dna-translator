# DNA Translator v2.0

[![CI — Lint & Biological Tests](https://github.com/Ayan-Goswami/dna-translator/actions/workflows/ci.yml/badge.svg)](https://github.com/Ayan-Goswami/dna-translator/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://dna-translator-v2.streamlit.app/)

A Python bioinformatics package, command-line tool, and interactive Streamlit dashboard for 6-frame nucleotide translation, coordinate-aware Open Reading Frame (ORF) discovery, *in silico* mutation impact simulation, and protein physicochemical profiling.

## What Changed from v1.0 to v2.0

| Feature | v1.0 | v2.0 |
| :--- | :--- | :--- |
| **Translation Engine** | Frame `+1`, Standard genetic code only | All **6 reading frames** (`+1/+2/+3`, `-1/-2/-3`) across **NCBI Genetic Code Tables** (1, 2, 11) |
| **ORF Discovery** | Raw end-to-end translation | Coordinate-aware **6-frame ORF finder** mapped to 1-based forward reference coordinates |
| **Variant Analysis** | Not available | **Wild-Type vs. Mutant simulator** classifying Silent, Missense, Nonsense, Readthrough, In-frame Indel, and Frameshift variants with HGVS protein notation |
| **Protein Profiling** | Amino acid sequence string | ExPASy-equivalent **ProtParam metrics** (MW, pI, GRAVY, Instability Index, charge at pH 7.0) + **Kyte-Doolittle hydropathy curves** |
| **Interface & Quality** | Single standalone script | Modular `src/` package, **Typer CLI** (CSV/JSON/FASTA export), **Streamlit + Plotly web app**, `pytest` suite, and GitHub Actions CI |

## Quick Start

### 1. Local Installation

```bash
git clone (https://github.com/Ayan-Goswami/dna-translator.git)
cd dna-translator
pip install -e .
```

### 2. Launch the Interactive Web Dashboard

```bash
streamlit run app.py
```

The dashboard includes two interactive modules pre-loaded with human beta-globin (*HBB*, RefSeq `NM_000518.4`):
* **6-Frame & ORF Explorer:** Upload or paste FASTA/raw DNA or RNA sequences, inspect interactive 6-frame ORF genome tracks, profile predicted proteins, and export results as `.csv` or `.fasta`.
* **Mutation Impact Simulator:** Compare a wild-type coding sequence against a mutant sequence (with one-click demos for Sickle Cell `p.Glu7Val`, Beta-Thalassemia `p.Glu7*`, and Frameshift deletions), inspect physicochemical delta cards, overlay Kyte-Doolittle hydropathy curves, and view highlighted residue alignments.

### 3. Command-Line Interface (Batch Workflows)

```bash
# Scan all 6 reading frames, extract ORFs, and export reports
dna-translator orfs data/hbb_coding.fasta --min-aa 30 --csv orfs.csv --fasta proteins.fasta

# Translate across all 6 reading frames (supports NCBI table IDs: 1, 2, 11)
dna-translator translate data/hbb_coding.fasta --table 1

# Compare wild-type vs. mutant CDS (e.g., Sickle Cell c.20A>T / p.Glu7Val)
dna-translator mutate data/hbb_coding.fasta data/hbb_sickle.fasta
```

Example `mutate` output on human *HBB* (`NM_000518.4`) vs. Sickle Cell (`c.20A>T`):

```text
Variant Classification: Missense
WT Length:  147 aa  (HBB_coding_sequence_RefSeq_NM_000518.4)
Mut Length: 147 aa  (HBB_Sickle_p.Glu7Val)  [Delta: +0 aa]

Residue Changes:
  • p.Glu7Val  |  Negatively charged -> Hydrophobic
```

## Project Architecture

```text
dna-translator/
├── .github/workflows/
│   └── ci.yml                 # Automated ruff linting + pytest CI pipeline
├── data/
│   ├── hbb_coding.fasta       # Human HBB wild-type CDS (RefSeq NM_000518.4)
│   └── hbb_sickle.fasta       # Clinical Sickle Cell mutant (c.20A>T, p.Glu7Val)
├── src/dna_translator/
│   ├── core.py                # Sequence validation, 6-frame translation & ORF finder
│   ├── properties.py          # Biopython ProtParam metrics & Kyte-Doolittle profile
│   ├── mutations.py           # Variant classification, HGVS notation & polarity shifts
│   ├── viz.py                 # Interactive Plotly ORF tracks & hydropathy curves
│   └── cli.py                 # Typer CLI for batch FASTA processing
├── tests/
│   └── test_core.py           # 11 biological unit tests (forward/reverse ORFs, RNA, Table 2, variants)
├── app.py                     # Streamlit interactive dashboard
├── pyproject.toml             # Modern Python packaging & CLI entry point
└── requirements.txt           # Cloud deployment dependencies
```

## Running Tests & Linting

```bash
pytest -v
ruff check src tests app.py
```

## Author

**Ayan Goswami** — M.Sc. Biotechnology | Bioinformatics & Healthcare Data Analytics  
[GitHub](https://github.com/Ayan-Goswami)
