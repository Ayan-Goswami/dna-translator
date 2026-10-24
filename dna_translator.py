"""Beginner-friendly DNA translator and small genetics toolkit.

Features:
- Translate a coding DNA sequence into a protein.
- Search all six reading frames for complete ATG-to-stop ORFs.
- Build and display the standard genetic code.
- Compare all possible single-base substitutions at each coding position.
- Write protein sequences in FASTA format.

This educational tool assumes the standard genetic code and DNA input (A, C, G, T).
It is not a clinical variant interpretation tool.
"""

import argparse
from pathlib import Path
import sys

BASES = "TCAG"
AMINO_ACIDS = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)


def build_codon_table():
    """Build the 64-codon standard genetic code in T-C-A-G order."""
    table = {}
    index = 0
    for first in BASES:
        for second in BASES:
            for third in BASES:
                table[first + second + third] = AMINO_ACIDS[index]
                index += 1
    return table


CODON_TABLE = build_codon_table()
STOP_CODONS = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def clean_dna(sequence):
    """Uppercase DNA and remove whitespace; reject characters outside A/C/G/T."""
    sequence = "".join(sequence.upper().split())
    if not sequence:
        raise ValueError("DNA sequence cannot be empty.")
    invalid = sorted(set(sequence) - set("ACGT"))
    if invalid:
        raise ValueError("Invalid DNA character(s): " + ", ".join(invalid))
    return sequence


def read_sequence(sequence_arg=None, file_arg=None):
    """Read DNA from --sequence or a text/FASTA file."""
    if sequence_arg and file_arg:
        raise ValueError("Use either --sequence or --file, not both.")
    if file_arg:
        content = Path(file_arg).read_text(encoding="utf-8")
        # Ignore FASTA header lines beginning with >
        lines = [line for line in content.splitlines() if not line.startswith(">")]
        raw = "".join(lines)
    elif sequence_arg:
        raw = sequence_arg
    else:
        raw = input("Paste DNA sequence: ")
    return clean_dna(raw)


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    return sequence.translate(COMPLEMENT)[::-1]


def translate_coding_sequence(sequence, include_stop=False):
    """Translate from the first base in triplets; stop terminates translation."""
    sequence = clean_dna(sequence)
    protein = []
    complete_length = len(sequence) - (len(sequence) % 3)
    for position in range(0, complete_length, 3):
        codon = sequence[position:position + 3]
        amino_acid = CODON_TABLE[codon]
        if amino_acid == "*":
            if include_stop:
                protein.append("*")
            break
        protein.append(amino_acid)
    return "".join(protein)


def find_orfs(sequence, both_strands=True, min_aa=1):
    """Find complete ATG-to-stop ORFs in 3 or 6 reading frames.

    Every ATG is considered as a possible start. Coordinates are 1-based and
    inclusive on the strand being scanned. Reverse-strand coordinates are
    reported in reverse-complement coordinates and labelled with a minus sign.
    """
    sequence = clean_dna(sequence)
    strands = [("+", sequence)]
    if both_strands:
        strands.append(("-", reverse_complement(sequence)))

    found = []
    for strand_label, strand_sequence in strands:
        for frame in range(3):
            starts = []
            for pos in range(frame, len(strand_sequence) - 2, 3):
                codon = strand_sequence[pos:pos + 3]
                if codon == "ATG":
                    starts.append(pos)
                if codon in STOP_CODONS:
                    for start in starts:
                        coding = strand_sequence[start:pos + 3]
                        aa_count = (len(coding) // 3) - 1
                        if aa_count >= min_aa:
                            found.append({
                                "strand": strand_label,
                                "frame": frame + 1,
                                "start": start + 1,
                                "end": pos + 3,
                                "codons_including_stop": len(coding) // 3,
                                "protein_length": aa_count,
                                "protein": translate_coding_sequence(coding),
                                "dna": coding,
                            })
                    starts = []

    return sorted(
        found,
        key=lambda item: (item["protein_length"], item["codons_including_stop"]),
        reverse=True,
    )


def codon_table_rows():
    """Return codons and amino acids in alphabetical codon order."""
    return [(codon, CODON_TABLE[codon]) for codon in sorted(CODON_TABLE)]


def mutation_analysis(sequence):
    """Classify every possible single-base substitution in complete codons.

    The sequence is treated as a coding sequence starting at its first base.
    The final 1-2 bases, if any, are ignored because they do not form a codon.
    Categories are silent, missense, and nonsense. Changes from a stop codon
    are reported as stop_lost/stop_changed rather than folded into these totals.
    """
    sequence = clean_dna(sequence)
    complete_length = len(sequence) - (len(sequence) % 3)
    coding = sequence[:complete_length]
    results = []
    counts = {"silent": 0, "missense": 0, "nonsense": 0, "stop_changed": 0}

    for pos, original_base in enumerate(coding):
        codon_start = (pos // 3) * 3
        original_codon = coding[codon_start:codon_start + 3]
        original_aa = CODON_TABLE[original_codon]

        for new_base in "ACGT":
            if new_base == original_base:
                continue
            mutated_codon_list = list(original_codon)
            mutated_codon_list[pos % 3] = new_base
            mutated_codon = "".join(mutated_codon_list)
            mutated_aa = CODON_TABLE[mutated_codon]

            if original_aa == "*":
                category = "stop_changed"
            elif mutated_aa == "*":
                category = "nonsense"
            elif mutated_aa == original_aa:
                category = "silent"
            else:
                category = "missense"

            counts[category] += 1
            results.append({
                "position": pos + 1,
                "original_base": original_base,
                "new_base": new_base,
                "original_codon": original_codon,
                "mutated_codon": mutated_codon,
                "original_aa": original_aa,
                "mutated_aa": mutated_aa,
                "category": category,
            })

    return counts, results, complete_length


def write_fasta(protein, output_path, header="translated_protein", line_width=60):
    """Write a protein sequence to a FASTA file."""
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(f">{header}\n")
        for start in range(0, len(protein), line_width):
            handle.write(protein[start:start + line_width] + "\n")


def command_translate(args):
    dna = read_sequence(args.sequence, args.file)
    print(f"DNA length: {len(dna)} bases")
    if args.reverse:
        dna = reverse_complement(dna)
        print("Using reverse-complement strand.")
    protein = translate_coding_sequence(dna, include_stop=args.include_stop)
    print(f"Protein length: {len(protein.replace('*', ''))} amino acids")
    print("Protein:")
    print(protein or "(empty protein; no complete amino-acid codons before a stop)")
    if args.fasta:
        fasta_protein = protein.rstrip("*")
        write_fasta(fasta_protein, args.fasta, header=args.header)
        print(f"FASTA saved to: {args.fasta}")


def command_orfs(args):
    dna = read_sequence(args.sequence, args.file)
    orfs = find_orfs(dna, both_strands=not args.forward_only, min_aa=args.min_aa)
    if not orfs:
        print("No complete ATG-to-stop ORFs found for the selected settings.")
        return
    print(f"Complete ORFs found: {len(orfs)}")
    print("Rank | Strand | Frame | Start | End | Codons incl. stop | Protein aa")
    for rank, item in enumerate(orfs[:args.limit], start=1):
        print(
            f"{rank:>4} | {item['strand']:>6} | {item['frame']:>5} | "
            f"{item['start']:>5} | {item['end']:>3} | "
            f"{item['codons_including_stop']:>16} | {item['protein_length']:>10}"
        )
    best = orfs[0]
    print(f"\nLongest ORF protein (strand {best['strand']}, frame {best['frame']}):")
    print(best["protein"])
    if args.fasta:
        write_fasta(best["protein"], args.fasta, header=args.header)
        print(f"FASTA saved to: {args.fasta}")


def command_codons(argparse):
    print("Codon  Amino acid")
    print("-----  -----------")
    for codon, aa in codon_table_rows():
        label = "STOP" if aa == "*" else aa
        print(f"{codon}    {label}")
    synonymous = {}
    for codon, aa in CODON_TABLE.items():
        if aa != "*":
            synonymous.setdefault(aa, []).append(codon)
    six_codon = sorted(aa for aa, codons in synonymous.items() if len(codons) == 6)
    print("\nAmino acids encoded by six codons:", ", ".join(six_codon))


def command_mutations(args):
    dna = read_sequence(args.sequence, args.file)
    counts, results, complete_length = mutation_analysis(dna)
    total = len(results)
    print(f"Input length: {len(dna)} bases")
    print(f"Bases included in complete codons: {complete_length}")
    print(f"Possible single-base substitutions tested: {total}")
    print("Category totals:")
    for category in ("silent", "missense", "nonsense", "stop_changed"):
        number = counts[category]
        percentage = (number / total * 100) if total else 0
        print(f"  {category:12s} {number:5d} ({percentage:5.1f}%)")
    if args.details:
        print("\nPosition  Change  Codon       Amino acid  Classification")
        for item in results:
            old_aa = "STOP" if item["original_aa"] == "*" else item["original_aa"]
            new_aa = "STOP" if item["mutated_aa"] == "*" else item["mutated_aa"]
            print(
                f"{item['position']:>8}  {item['original_base']}>{item['new_base']}     "
                f"{item['original_codon']}->{item['mutated_codon']}  "
                f"{old_aa:>3}->{new_aa:<3}       {item['category']}"
            )


def build_parser():
    parser = argparse.ArgumentParser(
        description="DNA translator: translation, six-frame ORFs, codons, and mutation analysis."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    translate_parser = subparsers.add_parser("translate", help="Translate a coding DNA sequence.")
    source = translate_parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sequence", help="DNA sequence, quoted if it contains spaces.")
    source.add_argument("--file", help="Text or FASTA file containing DNA.")
    translate_parser.add_argument("--reverse", action="store_true", help="Translate the reverse complement instead.")
    translate_parser.add_argument("--include-stop", action="store_true", help="Include '*' if a stop codon is reached.")
    translate_parser.add_argument("--fasta", help="Save the protein sequence to this FASTA file.")
    translate_parser.add_argument("--header", default="translated_protein", help="FASTA header (default: translated_protein).")
    translate_parser.set_defaults(func=command_translate)

    orf_parser = subparsers.add_parser("orf", help="Find complete ORFs in three or six frames.")
    source = orf_parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sequence", help="DNA sequence.")
    source.add_argument("--file", help="Text or FASTA file containing DNA.")
    orf_parser.add_argument("--forward-only", action="store_true", help="Scan only the three forward frames.")
    orf_parser.add_argument("--min-aa", type=int, default=1, help="Minimum protein length (default: 1).")
    orf_parser.add_argument("--limit", type=int, default=20, help="Maximum ORFs to list (default: 20).")
    orf_parser.add_argument("--fasta", help="Save the longest ORF protein to this FASTA file.")
    orf_parser.add_argument("--header", default="longest_orf_protein", help="FASTA header.")
    orf_parser.set_defaults(func=command_orfs)

    codon_parser = subparsers.add_parser("codons", help="Print the 64-codon standard genetic code.")
    codon_parser.set_defaults(func=command_codons)

    mutation_parser = subparsers.add_parser("mutations", help="Classify all single-base substitutions in a coding sequence.")
    source = mutation_parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sequence", help="Coding DNA sequence.")
    source.add_argument("--file", help="Text or FASTA file containing coding DNA.")
    mutation_parser.add_argument("--details", action="store_true", help="Print every tested substitution.")
    mutation_parser.set_defaults(func=command_mutations)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if hasattr(args, "min_aa") and args.min_aa < 1:
        parser.error("--min-aa must be at least 1.")
    if hasattr(args, "limit") and args.limit < 1:
        parser.error("--limit must be at least 1.")
    try:
        args.func(args)
    except (ValueError, OSError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

