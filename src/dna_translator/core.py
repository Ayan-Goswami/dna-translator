from dataclasses import dataclass
from Bio.Data import CodonTable
from Bio.Seq import Seq


@dataclass(frozen=True)
class ORF:
    """Represents a single Open Reading Frame with 1-based genomic coordinates."""
    frame: str
    strand: str
    start: int          # 1-based inclusive start coordinate on forward strand
    end: int            # 1-based inclusive end coordinate on forward strand
    nt_length: int
    aa_length: int
    dna_seq: str
    protein_seq: str


def clean_dna(raw_seq: str) -> str:
    """Normalize input sequence: uppercase, strip whitespace, convert RNA U to T."""
    cleaned = "".join(raw_seq.split()).upper().replace("U", "T")
    if not cleaned:
        raise ValueError("Input sequence is empty.")
    valid_bases = set("ATGCNRYKMSWBDHV")  # Standard IUPAC nucleotide codes
    invalid = set(cleaned) - valid_bases
    if invalid:
        raise ValueError(f"Invalid nucleotide characters found: {sorted(invalid)}")
    return cleaned


def _trim_to_codons(seq_str: str) -> str:
    """Trim 1-2 trailing bases so Biopython does not emit Partial Codon warnings."""
    remainder = len(seq_str) % 3
    return seq_str[:-remainder] if remainder else seq_str


def translate_six_frames(dna_seq: str, table: int = 1) -> dict[str, str]:
    """Translate a DNA sequence across all 3 forward and 3 reverse-complement frames."""
    cleaned = clean_dna(dna_seq)
    fwd = Seq(cleaned)
    rev = fwd.reverse_complement()

    frames: dict[str, str] = {}
    for offset in range(3):
        fwd_slice = _trim_to_codons(str(fwd[offset:]))
        rev_slice = _trim_to_codons(str(rev[offset:]))
        frames[f"+{offset + 1}"] = str(Seq(fwd_slice).translate(table=table))
        frames[f"-{offset + 1}"] = str(Seq(rev_slice).translate(table=table))
    return frames


def find_orfs(
    dna_seq: str,
    min_aa_length: int = 30,
    table: int = 1,
    allow_alt_starts: bool = False,
) -> list[ORF]:
    """Scan all 6 frames and return ORFs sorted by amino acid length (longest first)."""
    cleaned = clean_dna(dna_seq)
    seq_len = len(cleaned)
    codon_table = CodonTable.unambiguous_dna_by_id[table]

    stop_codons = set(codon_table.stop_codons)
    start_codons = set(codon_table.start_codons) if allow_alt_starts else {"ATG"}

    strands = [
        ("+", cleaned),
        ("-", str(Seq(cleaned).reverse_complement())),
    ]

    orfs: list[ORF] = []
    for strand, strand_seq in strands:
        for offset in range(3):
            frame_label = f"{strand}{offset + 1}"
            orf_start_idx: int | None = None

            for pos in range(offset, seq_len - 2, 3):
                codon = strand_seq[pos : pos + 3]

                # Open a new ORF at the earliest valid start codon
                if orf_start_idx is None and codon in start_codons:
                    orf_start_idx = pos
                    continue

                # Close the ORF when the first in-frame stop codon hits
                if orf_start_idx is not None and codon in stop_codons:
                    orf_end_exclusive = pos + 3
                    orf_dna = strand_seq[orf_start_idx:orf_end_exclusive]
                    protein = str(Seq(orf_dna).translate(table=table, to_stop=True))

                    if len(protein) >= min_aa_length:
                        if strand == "+":
                            genomic_start = orf_start_idx + 1
                            genomic_end = orf_end_exclusive
                        else:
                            # Map reverse-strand slice back to 1-based forward coordinates
                            genomic_start = seq_len - orf_end_exclusive + 1
                            genomic_end = seq_len - orf_start_idx

                        orfs.append(
                            ORF(
                                frame=frame_label,
                                strand=strand,
                                start=genomic_start,
                                end=genomic_end,
                                nt_length=len(orf_dna),
                                aa_length=len(protein),
                                dna_seq=orf_dna,
                                protein_seq=protein,
                            )
                        )
                    orf_start_idx = None

    return sorted(orfs, key=lambda o: o.aa_length, reverse=True)