DNA Translator
A beginner-friendly Python command-line bioinformatics tool for translating DNA coding sequences into protein sequences and exploring core concepts in molecular genetics.
<p align="center">
  <strong>DNA → Codons → Protein</strong>
</p>

Overview
DNA Translator is an educational project that applies fundamental molecular biology concepts through Python. It provides command-line operations for sequence translation, open reading frame (ORF) detection, codon-table inspection, single-nucleotide substitution analysis, and FASTA export.

The project uses the standard genetic code and includes a human beta-globin (HBB) coding sequence as an example dataset.
Features
- DNA translation — translate a coding DNA sequence into a one-letter amino acid sequence.
- ORF detection — search three forward reading frames or all six reading frames for complete start-to-stop ORFs.
- Codon table — generate and display the 64-codon standard genetic code.
- Mutation analysis — classify possible single-base substitutions as silent, missense, nonsense, or stop-codon changes.
- FASTA export — save translated protein sequences for downstream educational analyses.
- Input validation — reject empty sequences and characters outside A, C, G, and T.
- Automated tests — use Python's built-in unittest framework to check core behavior.
Requirements
- Python 3.9 or newer recommended
- No third-party packages required
Getting started
1. Clone the repository
git clone https://github.com/Ayan-Goswami/dna-translator.git
cd dna-translator
Alternatively, download the repository as a ZIP file from GitHub and extract it.
2. Show the codon table
python dna_translator.py codons
3. Translate a test sequence
python dna_translator.py translate --sequence ATGGTGTAA
Expected protein:
MV
The sequence contains ATG (methionine), GTG (valine), and TAA (a stop codon). The stop codon terminates translation and is not counted as an amino acid.
Usage
Run commands from the repository directory.
Translate a sequence
python dna_translator.py translate --sequence ATGGTGTAA
Include a * marker when a stop codon is reached:
python dna_translator.py translate --sequence ATGGTGTAA --include-stop
Read DNA from a FASTA file
python dna_translator.py translate --file hbb_coding.fasta
Find open reading frames
Scan all six reading frames:
python dna_translator.py orf --file hbb_coding.fasta
Scan only the three forward frames:
python dna_translator.py orf --file hbb_coding.fasta --forward-only
Set the minimum ORF protein length:
python dna_translator.py orf --file hbb_coding.fasta --min-aa 20
The ORF report includes strand, frame, start and end positions, codon count, and protein length. The codon count includes the terminating stop codon; the protein length excludes it. Coordinates are one-based and inclusive on the strand being scanned. For reverse-strand ORFs, coordinates refer to the reverse-complement sequence.
Analyze single-base substitutions
python dna_translator.py mutations --file hbb_coding.fasta
Print details for each possible substitution:
python dna_translator.py mutations --file hbb_coding.fasta --details
The input is treated as a coding sequence starting at its first base. Incomplete trailing bases are excluded from mutation analysis.
Export a protein sequence in FASTA format
python dna_translator.py translate --file hbb_coding.fasta --fasta hbb_protein.fasta --header human_beta_globin
Example dataset: human beta-globin (HBB)
The included hbb_coding.fasta file is intended for practicing coding-sequence translation and reading-frame analysis.

The reference protein begins with:
MVHLTPEEKSAVTALWGKV
For the 444-base coding sequence in the project reference, the expected translated protein has 147 amino acids, excluding the terminal stop codon (148 codons including the stop).
Run the tests
python -m unittest -v
The tests cover codon-table size, start/stop codons, DNA validation, reverse-complement generation, translation, ORF detection, substitution-count arithmetic, and FASTA output.
How it works
1. Validate input. Normalize case, remove whitespace, and allow only A, C, G, and T.
2. Build the genetic code. Generate codons from the four bases and map them to standard amino acid symbols.
3. Translate codons. Read the sequence in triplets and stop at a termination codon.
4. Search ORFs. Examine each selected reading frame for start-to-stop regions.
5. Analyze substitutions. Change one nucleotide at a time and compare the amino acid encoded by the affected codon.
6. Export results. Write translated protein sequences to FASTA when requested.
Scope and limitations
This is an educational project, not a clinically validated variant-analysis pipeline. It assumes the standard genetic code and simplified ORF rules. Alternative genetic codes, splice-aware transcript processing, expression evidence, and clinical interpretation are outside the current scope. A predicted ORF alone does not establish that a sequence is expressed or produces a functional protein.

Mutation categories are computational descriptions of codon-level consequences and should not be interpreted as a clinical assessment of a genetic variant.

Project structure

dna-translator
├── dna_translator.py
├── test_dna_translator.py
├── hbb_coding.fasta
├── README.md
├── .gitignore
└── LICENSE

Learning objectives
- Python functions, loops, dictionaries, and command-line arguments
- DNA sequence validation and string manipulation
- Codons, translation, start/stop signals, and reading frames
- Reverse-complement logic
- Basic mutation classification
- File input/output and FASTA formatting
- Unit testing and reproducible checks
  
License
See the included LICENSE file.
Built as a learning project at the intersection of Biotechnology, Molecular Biology, Python, and Bioinformatics
