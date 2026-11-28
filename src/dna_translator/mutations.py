from dataclasses import dataclass

from Bio.Seq import Seq
from Bio.SeqUtils import seq3

from dna_translator.core import _trim_to_codons, clean_dna

AA_CLASSES = {
    "A": "Hydrophobic", "V": "Hydrophobic", "I": "Hydrophobic",
    "L": "Hydrophobic", "M": "Hydrophobic", "F": "Aromatic/Hydrophobic",
    "W": "Aromatic/Hydrophobic", "Y": "Aromatic/Polar", "P": "Special (Proline)",
    "G": "Special (Glycine)", "C": "Special (Cysteine)",
    "S": "Polar uncharged", "T": "Polar uncharged", "N": "Polar uncharged",
    "Q": "Polar uncharged",
    "D": "Negatively charged", "E": "Negatively charged",
    "K": "Positively charged", "R": "Positively charged", "H": "Positively charged",
}


@dataclass(frozen=True)
class ResidueChange:
    position: int           # 1-based amino acid index
    wt_aa: str
    mut_aa: str
    hgvs_p: str             # e.g., p.Glu7Val
    property_shift: str     # e.g., Negatively charged -> Hydrophobic


@dataclass(frozen=True)
class MutationReport:
    variant_type: str       # Identical, Synonymous, Missense, Nonsense, Readthrough, Frameshift, In-frame Indel
    wt_protein: str
    mut_protein: str
    length_delta_aa: int
    changes: list[ResidueChange]


def _translate_cds(dna: str, table: int) -> str:
    """Translate a coding sequence up to the first stop codon."""
    trimmed = _trim_to_codons(dna)
    return str(Seq(trimmed).translate(table=table, to_stop=True))


def analyze_mutation(wt_dna: str, mut_dna: str, table: int = 1) -> MutationReport:
    """Compare a wild-type CDS against a mutant CDS and classify the functional impact."""
    wt = clean_dna(wt_dna)
    mut = clean_dna(mut_dna)

    wt_prot = _translate_cds(wt, table)
    mut_prot = _translate_cds(mut, table)
    len_delta = len(mut_prot) - len(wt_prot)

    if wt == mut:
        return MutationReport("Identical", wt_prot, mut_prot, 0, [])

    nt_diff = len(mut) - len(wt)
    if nt_diff % 3 != 0:
        # Find the first altered residue caused by the frameshift
        changes = _diff_residues(wt_prot, mut_prot, stop_at_first=True)
        return MutationReport("Frameshift", wt_prot, mut_prot, len_delta, changes)

    if wt_prot == mut_prot:
        return MutationReport("Synonymous (Silent)", wt_prot, mut_prot, 0, [])

    if nt_diff != 0:
        changes = _diff_residues(wt_prot, mut_prot, stop_at_first=True)
        return MutationReport("In-frame Indel", wt_prot, mut_prot, len_delta, changes)

    if len(mut_prot) < len(wt_prot):
        stop_pos = len(mut_prot) + 1
        wt_aa = wt_prot[len(mut_prot)]
        change = ResidueChange(
            position=stop_pos,
            wt_aa=wt_aa,
            mut_aa="*",
            hgvs_p=f"p.{seq3(wt_aa)}{stop_pos}*",
            property_shift=f"{AA_CLASSES.get(wt_aa, 'Unknown')} -> Premature Stop",
        )
        return MutationReport("Nonsense (Truncating)", wt_prot, mut_prot, len_delta, [change])

    if len(mut_prot) > len(wt_prot):
        return MutationReport("Readthrough (Stop-Loss)", wt_prot, mut_prot, len_delta, [])

    changes = _diff_residues(wt_prot, mut_prot, stop_at_first=False)
    return MutationReport("Missense", wt_prot, mut_prot, 0, changes)


def _diff_residues(wt_prot: str, mut_prot: str, stop_at_first: bool) -> list[ResidueChange]:
    """Identify mismatched amino acid positions and biochemical class shifts."""
    changes: list[ResidueChange] = []
    for idx, (w, m) in enumerate(zip(wt_prot, mut_prot), start=1):
        if w != m:
            w_cls = AA_CLASSES.get(w, "Unknown")
            m_cls = AA_CLASSES.get(m, "Unknown")
            shift = "Conservative" if w_cls == m_cls else f"{w_cls} -> {m_cls}"
            changes.append(
                ResidueChange(
                    position=idx,
                    wt_aa=w,
                    mut_aa=m,
                    hgvs_p=f"p.{seq3(w)}{idx}{seq3(m)}",
                    property_shift=shift,
                )
            )
            if stop_at_first:
                break
    return changes