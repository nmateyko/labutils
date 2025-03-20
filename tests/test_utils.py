import unittest
from labutils.sequence_analysis.utils import (
    read_fastq,
    read_fasta,
    revcomp,
    revcomp_read,
    levenshtein,
    hamming_dist,
    is_DNA,
)
from io import StringIO

class TestUtils(unittest.TestCase):

    def test_read_fastq(self):
        # Create a mock FASTQ file with a valid sequence entry
        fastq_data = "@SEQ_ID\nGATTACA\n+\n!!''*((((***\n"
        with StringIO(fastq_data) as fastqfile:
            result = list(read_fastq(fastqfile))
            self.assertEqual(len(result), 1)
            header, seq, qual = result[0]
            self.assertEqual(header, "@SEQ_ID")
            self.assertEqual(seq, "GATTACA")
            self.assertEqual(qual, "!!''*((((***")

        # Test with invalid FASTQ entry (missing a line)
        incomplete_fastq = "@SEQ_ID\nGATTACA\n+\n"
        with StringIO(incomplete_fastq) as fastqfile:
            with self.assertRaises(EOFError):
                list(read_fastq(fastqfile))

        # Test with incorrect headers
        invalid_header_fastq = "@SEQ_ID\nGATTACA\n-\n!!''*((((***\n"
        with StringIO(invalid_header_fastq) as fastqfile:
            with self.assertRaises(ValueError):
                list(read_fastq(fastqfile))

    def test_revcomp(self):
        # Valid sequence
        self.assertEqual(revcomp("ACTG"), "CAGT")
        self.assertEqual(revcomp("NNNN"), "NNNN")
        self.assertEqual(revcomp(""), "")  # Edge case: empty sequence
        
        # Invalid sequence containing characters not in ACTGN
        with self.assertRaises(ValueError):
            revcomp("ACXG")

    def test_levenshtein(self):
        # Simple cases
        self.assertEqual(levenshtein("kitten", "sitting"), 3)
        self.assertEqual(levenshtein("flaw", "lawn"), 2)
        self.assertEqual(levenshtein("gumbo", "gambol"), 2)

        # Edge cases
        self.assertEqual(levenshtein("", ""), 0)  # Both empty
        self.assertEqual(levenshtein("abc", ""), 3)  # One empty
        self.assertEqual(levenshtein("", "abc"), 3)  # One empty reversed
        self.assertEqual(levenshtein("abc", "abc"), 0)  # Identical strings

    def test_hamming_dist(self):
        # Valid cases
        self.assertEqual(hamming_dist("ACTG", "ACTG"), 0)  # Same sequence
        self.assertEqual(hamming_dist("ACTG", "ACGG"), 1)
        self.assertEqual(hamming_dist("GATTACA", "GACTATA"), 2)

        # Edge case: Different lengths (should raise an error)
        with self.assertRaises(ValueError):
            hamming_dist("ACTG", "ACT")

    def test_is_DNA(self):
        # Valid cases
        self.assertTrue(is_DNA("ACTG"))
        self.assertTrue(is_DNA("NNNN"))
        self.assertTrue(is_DNA("ACGTN"))
        self.assertTrue(is_DNA(""))

        # Invalid cases
        self.assertFalse(is_DNA("ACXG"))
        self.assertFalse(is_DNA("ACGTX"))
        self.assertFalse(is_DNA("ACGTN "))
        self.assertFalse(is_DNA("ACGTN\n"))

    def test_read_fasta(self):
        # Valid FASTA input
        fasta_data = """>seq1
        ATCGATCG
        >seq2
        GGGTTTCCC
        """
        with StringIO(fasta_data) as fasta_file:
            result = list(read_fasta(fasta_file))
            self.assertEqual(len(result), 2)
            self.assertEqual(result[0], ("seq1", "ATCGATCG"))
            self.assertEqual(result[1], ("seq2", "GGGTTTCCC"))

        # Missing header
        invalid_fasta = """ATCGATCG"""
        with StringIO(invalid_fasta) as fasta_file:
            with self.assertRaises(ValueError):
                list(read_fasta(fasta_file))

        # Header with no sequence
        empty_seq_fasta = """>seq1
        >seq2
        GGGTTTCCC
        """
        with StringIO(empty_seq_fasta) as fasta_file:
            with self.assertRaises(ValueError):
                list(read_fasta(fasta_file))

    def test_revcomp_read(self):
        read = ("@SEQ_ID", "GATTACA", "!!''*(((")
        rev_read = revcomp_read(read)
        self.assertEqual(rev_read, ("@SEQ_ID", "TGTAATC", "(((*''!!"))

if __name__ == "__main__":
    unittest.main()