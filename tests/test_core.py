import pytest

from dna_translator.core import clean_dna, find_orfs, translate_six_frames
from dna_translator.mutations import analyze_mutation
from dna_translator.properties import profile_protein

# First 10 codons of human HBB (Beta-Globin) + TAA stop codon
HBB_WT = "ATGGTGCACCTGACTCCTGAGGAGAAGTCTTAA"
# Sickle Cell mutation: codon 7 GAG (Glu) -> GTG (Val)
HBB_SICKLE = "ATGGTGCACCTGACTCCTGTGGAGAAGTCTTAA"
# Nonsense mutation: codon 7 GAG (Glu) -> TAG (Stop)
HBB_NONSENSE = "ATGGTGCACCTGACTCCTTACGAGAAGTCTTAA".replace("TAC", "TAG")
# Synonymous (Silent) mutation: codon 7 GAG (Glu) -> GAA (Glu)
HBB_SILENT = "ATGGTGCACCTGACTCCTGAAGAGAAGTCTTAA"
# Frameshift mutation: 1-bp deletion at position 20 (index 19)
HBB_FRAMESHIFT = HBB_WT[:19] + HBB_WT[20:]


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


def test_orf_finder_reverse_strand_coordinates():
    # Reverse complement of ATG-GAA-AAA-GCT-GTT-TAA (MEKAV*) is TTAAACAGCTTTTTCCAT
    # Flanked by 2 bases on the left (CC) and 1 on the right (G) -> forward coords 3..20
    seq = "CCTTAAACAGCTTTTTCCATG"
    orfs = find_orfs(seq, min_aa_length=5)
    assert len(orfs) == 1
    assert orfs[0].frame == "-2"
    assert orfs[0].strand == "-"
    assert orfs[0].start == 3
    assert orfs[0].end == 20
    assert orfs[0].protein_seq == "MEKAV"


def test_rna_input_normalization():
    rna_hbb = HBB_WT.replace("T", "U")
    assert translate_six_frames(rna_hbb)["+1"] == "MVHLTPEEKS*"


def test_vertebrate_mitochondrial_table():
    # In Table 1, TGA is Stop (*). In Table 2 (Vert Mito), TGA is Trp (W) and AGA is Stop (*).
    seq = "ATGTGAGAAAGA"
    assert translate_six_frames(seq, table=1)["+1"] == "M*ER"
    assert translate_six_frames(seq, table=2)["+1"] == "MWE*"


def test_invalid_characters_rejected():
    with pytest.raises(ValueError, match="Invalid nucleotide characters"):
        clean_dna("ATGXYZ")


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


def test_synonymous_mutation():
    report = analyze_mutation(HBB_WT, HBB_SILENT)
    assert report.variant_type == "Synonymous (Silent)"
    assert report.changes == []


def test_frameshift_mutation():
    report = analyze_mutation(HBB_WT, HBB_FRAMESHIFT)
    assert report.variant_type == "Frameshift"
    assert len(report.changes) == 1
    assert report.changes[0].position == 7


def test_protein_profile():
    prof = profile_protein("MVHLTPEEKS*")
    assert prof.length == 10
    assert prof.molecular_weight_kda > 0