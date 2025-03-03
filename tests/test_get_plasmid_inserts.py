import unittest
from labutils.sequence_analysis.get_plasmid_inserts import get_align_pos, get_insert_pos, extract_insert_and_flanks, extract_insert_fastq
from labutils.sequence_analysis.utils import revcomp

class TestGetPlasmidInserts(unittest.TestCase):

    def test_get_align_pos(self):
        # Test for perfect match
        seq = "ATCGATCGATCG"
        kmer = "ATCG"
        dist = 0
        result = get_align_pos(seq, kmer, dist)
        self.assertEqual(result, (0, 3))  # Expect match at start
        
        # Test for mismatch within threshold
        kmer = "ATCA"
        dist = 1
        result = get_align_pos(seq, kmer, dist)
        self.assertEqual(result, (0, 2))  # Expect best match at start

        # Test for no match within distance threshold
        dist = 0
        result = get_align_pos(seq, "AAAA", dist)
        self.assertEqual(result, -1)

    def test_get_insert_pos(self):
        read = ("header", "ATTTCCCCCCGTGTGTGTGTAAAAATT", "quality_placeholder")
        up = "CCCCC"
        down = "AAAAA"
        t = 1

        # Test insert detection with forward orientation
        result = get_insert_pos(read, up, down, t)
        self.assertIsNotNone(result)
        self.assertEqual(result["up"], (4, 8))
        self.assertEqual(result["down"], (20, 24))
        self.assertFalse(result["rc"])

        # Test insert detection with reverse complement orientation
        read_rc = ("header", revcomp(read[1]), "quality_placeholder"[::-1])
        result = get_insert_pos(read_rc, up, down, t)
        self.assertIsNotNone(result)
        self.assertTrue(result["rc"])

        # Test with no match
        no_upstream_read = ("header", "GGGGGGGGGGGGGG", "quality_placeholder")
        result = get_insert_pos(no_upstream_read, up, down, t)
        self.assertIsNone(result)

    def test_extract_insert_and_flanks(self):
        read = ("header", "ATCGATCGATCGTACGTAGC", "ABCDEFGHIJKLMNOPQRST")
        pos_dict = {
            "up": (0, 4),
            "down": (14, 18),
            "rc": False
        }

        # Test with forward orientation
        upstream, insert, downstream = extract_insert_and_flanks(read, pos_dict)
        self.assertEqual(upstream[0], "ATCGA")
        self.assertEqual(insert[0], "TCGATCGTA")
        self.assertEqual(downstream[0], "CGTAG")
        self.assertEqual(upstream[1], "ABCDE")
        self.assertEqual(insert[1], "FGHIJKLMN")
        self.assertEqual(downstream[1], "OPQRS")

        # Test with reverse complement orientation
        pos_dict["rc"] = True
        upstream, insert, downstream = extract_insert_and_flanks(read, pos_dict)
        self.assertEqual(upstream[0], revcomp("GTAGC"))
        self.assertEqual(insert[0], revcomp("CGATCGTAC"))
        self.assertEqual(downstream[0], revcomp("TCGAT"))
        self.assertEqual(upstream[1], "TSRQP")
        self.assertEqual(insert[1], "ONMLKJIHG")
        self.assertEqual(downstream[1], "FEDCB")

    def test_extract_insert_fastq(self):
        read = ("header", "ATCGATCGATCGTACGTAGC", "ABCDEFGHIJKLMNOPQRST")
        pos_dict = {
            "up": (0, 4),
            "down": (14, 18),
            "rc": False
        }

        # Test with forward orientation
        header, seq, qual = extract_insert_fastq(read, pos_dict)
        self.assertEqual(header, "header insert_pos=5-14 rc=False")
        self.assertEqual(seq, "TCGATCGTA")
        self.assertEqual(qual, "FGHIJKLMN")

        # Test with reverse complement orientation
        pos_dict["rc"] = True
        header, seq, qual = extract_insert_fastq(read, pos_dict)
        self.assertEqual(header, "header insert_pos=5-14 rc=True")
        self.assertEqual(seq, revcomp("CGATCGTAC"))
        self.assertEqual(qual, "ONMLKJIHG")

if __name__ == "__main__":
    unittest.main()