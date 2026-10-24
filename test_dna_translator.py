"""Validation tests for dna_translator.py; run with: python -m unittest -v"""
import os
import tempfile
import unittest

from dna_translator import (
    CODON_TABLE,
    clean_dna,
    find_orfs,
    mutation_analysis,
    reverse_complement,
    translate_coding_sequence,
    write_fasta,
)


class DnaTranslatorTests(unittest.TestCase):
    def test_codon_table_has_64_entries(self):
        self.assertEqual(len(CODON_TABLE), 64)

    def test_start_and_stop_codons(self):
        self.assertEqual(CODON_TABLE["ATG"], "M")
        self.assertEqual(CODON_TABLE["TAA"], "*")
        self.assertEqual(CODON_TABLE["TAG"], "*")
        self.assertEqual(CODON_TABLE["TGA"], "*")

    def test_clean_dna_removes_whitespace_and_uppercases(self):
        self.assertEqual(clean_dna(" atg\n ccc "), "ATGCCC")

    def test_invalid_dna_is_rejected(self):
        with self.assertRaises(ValueError):
            clean_dna("ATGB")

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATGCC"), "GGCAT")

    def test_translation_stops_before_stop_codon(self):
        self.assertEqual(translate_coding_sequence("ATGGTTTAA"), "MV")

    def test_translation_can_include_stop_marker(self):
        self.assertEqual(translate_coding_sequence("ATGGTTTAA", include_stop=True), "MV*")

    def test_translation_ignores_incomplete_trailing_bases(self):
        self.assertEqual(translate_coding_sequence("ATGG"), "M")

    def test_orf_finds_expected_forward_frame(self):
        result = find_orfs("CCCATGAAATAG", both_strands=False)
        self.assertTrue(result)
        self.assertEqual(result[0]["strand"], "+")
        self.assertEqual(result[0]["frame"], 1)
        self.assertEqual(result[0]["protein"], "MK")
        self.assertEqual(result[0]["codons_including_stop"], 3)

    def test_mutation_count_is_three_per_complete_coding_base(self):
        counts, results, coding_length = mutation_analysis("ATGGCTTAA")
        self.assertEqual(coding_length, 9)
        self.assertEqual(len(results), 27)
        self.assertEqual(sum(counts.values()), 27)

    def test_fasta_export(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "protein.fasta")
            write_fasta("MVHLT", path, header="test_protein")
            with open(path, encoding="utf-8") as handle:
                content = handle.read()
            self.assertEqual(content, ">test_protein\nMVHLT\n")


if __name__ == "__main__":
    unittest.main()

