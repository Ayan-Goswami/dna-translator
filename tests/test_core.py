from dna_translator.core import find_orfs, translate_six_frames
from dna_translator.mutations import analyze_mutation
from dna_translator.properties import profile_protein

# First 10 codons of human HBB (Beta-Globin) + TAA stop codon
HBB_WT = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTTAA"
# Sickle Cell mutation: codon 7 GAG (Glu) -> GTG (Val)
HBB_SICKLE = "ATGGTGCACCTGACTCCTGTGGAGAAGTCTTAA"
# Nonsense mutation: codon 7 GAG (Glu) -> TAG (Stop)
HBB_NONSENSE = "ATGGTGCACCTGACTCCTTACGAGAAGTCTTAA".replace("TAC", "TAG")


def test_six_frame_keys():
    frames = translate_six_frames(HBB_WT)
    assert set(frames.keys()) == {"+1", "+2", "+3", "-1", "-2", "-3"}
    assert frames["+1"] == "MVHLTPEEKS*"


def test_orf_finder_forward():
    orfs = find_orfs(HBB_WT, min_aa_length=5)
    assert len(orfs) == 1
    assert orfs[0].frame == "+1"
    assert orfs[0].start == 1
    assert orfs[0].end == 33
    assert orfs[0].protein_seq == "MVHLTPEEKS"


def test_sickle_cell_missense():
    report = analyze_mutation(HBB_WT, HBB_SICKLE)
    assert report.variant_type == "Missense"
    assert len(report.changes) == 1
    assert report.changes[0].hgvs_p == "p.Glu7Val"
    assert "Negatively charged -> Hydrophobic" in report.changes[0].property_shift


def test_nonsense_truncation():
    report = analyze_mutation(HBB_WT, HBB_NONSENSE)
    assert report.variant_type == "Nonsense (Truncating)"
    assert report.changes[0].hgvs_p == "p.Glu7*"
    assert report.length_delta_aa == -4


def test_protein_profile():
    prof = profile_protein("MVHLTPEEKS*")
    assert prof.length == 10
    assert prof.molecular_weight_kda > 0