from dataclasses import dataclass

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class ProteinProfile:
    length: int
    molecular_weight_kda: float
    isoelectric_point: float
    gravy: float
    instability_index: float
    stability_label: str
    aromaticity: float
    charge_at_ph7: float
    dropped_ambiguous_count: int


def profile_protein(protein_seq: str) -> ProteinProfile:
    """Compute ExPASy-equivalent physicochemical properties for a protein sequence."""
    stripped = protein_seq.upper().replace("*", "")
    clean_chars = [aa for aa in stripped if aa in STANDARD_AA]
    dropped = len(stripped) - len(clean_chars)

    if not clean_chars:
        raise ValueError("No standard amino acid residues available for profiling.")

    clean_seq = "".join(clean_chars)
    pa = ProteinAnalysis(clean_seq)
    instability = pa.instability_index()

    return ProteinProfile(
        length=len(clean_seq),
        molecular_weight_kda=round(pa.molecular_weight() / 1000.0, 2),
        isoelectric_point=round(pa.isoelectric_point(), 2),
        gravy=round(pa.gravy(), 3),
        instability_index=round(instability, 2),
        stability_label="Stable" if instability < 40.0 else "Unstable",
        aromaticity=round(pa.aromaticity(), 3),
        charge_at_ph7=round(pa.charge_at_pH(7.0), 2),
        dropped_ambiguous_count=dropped,
    )


def hydropathy_profile(protein_seq: str, window: int = 9) -> tuple[list[int], list[float]]:
    """Calculate sliding-window Kyte-Doolittle hydropathy values for plotting."""
    clean_seq = "".join(aa for aa in protein_seq.upper() if aa in STANDARD_AA)
    if len(clean_seq) < window:
        raise ValueError(
            f"Protein length ({len(clean_seq)}) must be at least the window size ({window})."
        )
    pa = ProteinAnalysis(clean_seq)
    scores = pa.protein_scale(param_dict=kd, window=window)
    # Center each score on its middle residue (1-based indexing)
    half = window // 2
    positions = list(range(1 + half, 1 + half + len(scores)))
    return positions, [round(s, 3) for s in scores]