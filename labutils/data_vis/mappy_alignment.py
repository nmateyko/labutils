import mappy as mp
from labutils.sequence_analysis.utils import read_fasta
from labutils.data_vis.alignment import AlignmentProcessor

class Minimap2AlignmentProcessor(AlignmentProcessor):
    def __init__(self, ref_seq_file):
        self.aligner = mp.Aligner(ref_seq_file, preset="map-ont", extra_flags=0x4000000)
        with open(ref_seq_file) as f:
            _, self.ref_seq = next(read_fasta(f))

    def align(self, seq):
        return next(self.aligner.map(seq), None)

    def encode_alignment(self, alignment):
        ref_length = len(self.ref_seq)
        alignment_row = [4] * ref_length  # Default to padding
        insertions = []

        if alignment is None:
            return alignment_row, insertions

        ref_start = alignment.r_st
        ref_pos = ref_start
        
        cigar_map = {
            7: 0,  # Match "="
            8: 1,  # Mismatch "X"
            1: 2,  # Insertion "I"
            2: 3,  # Deletion "D"
        }

        for length, op in alignment.cigar:
            if op in {7, 8, 2}:  # Match, Mismatch, Deletion
                alignment_row[ref_pos:ref_pos + length] = [cigar_map[op]] * length
                ref_pos += length
            elif op == 1:  # Insertion
                insertions.append((ref_pos, length))
        
        return alignment_row, insertions