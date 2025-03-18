import fastcluster
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio.Align import PairwiseAligner
from labutils.sequence_analysis.utils import revcomp
from matplotlib.patches import Rectangle

class AlignmentProcessor:
    def __init__(self, ref_seq):
        self.ref_seq = ref_seq

    def align(self, seq):
        raise NotImplementedError

    def encode_alignment(self, alignment):
        raise NotImplementedError
    
    def process_seqs(self, seqs):
        encoded_alignments = []
        insertions_list = []
        for seq in seqs:
            alignment = self.align(seq)
            encoded_alignment, insertions = self.encode_alignment(alignment)
            encoded_alignments.append(encoded_alignment)  # Convert to list for consistency
            insertions_list.append(insertions)
        return encoded_alignments, insertions_list

class BiopythonAlignmentProcessor(AlignmentProcessor):
    def __init__(self, ref_seq, match_score=2, mismatch_score=-1, open_gap_score=-1,
                 extend_gap_score=-1, try_revcomp=True, local=True):
        super().__init__(ref_seq)
        self.aligner = PairwiseAligner()
        if local:
            self.aligner.mode = 'local'
        else:
            self.aligner.mode = 'global'
        self.aligner.match_score = match_score
        self.aligner.mismatch_score = mismatch_score
        self.aligner.open_gap_score = open_gap_score
        self.aligner.extend_gap_score = extend_gap_score
        self.try_revcomp = try_revcomp
        self.local = local

    def align(self, seq):
        alignments = self.aligner.align(self.ref_seq, seq)
        best_alignment = next(alignments)

        if self.try_revcomp:
            seq_rc = revcomp(seq)
            alignments_rc = self.aligner.align(self.ref_seq, seq_rc)
            best_alignment_rc = next(alignments_rc)
            if best_alignment_rc.score > best_alignment.score:
                best_alignment = best_alignment_rc

        return best_alignment

    def encode_alignment(self, alignment):
        ref_aligned, query_aligned = (alignment[0], alignment[1])
        ref_length = len(self.ref_seq)
        
        alignment_row = [4] * ref_length  # Default to padding
        insertions = []
        
        # For local alignment, full reference sequence isn't returned,
        # so determine padding from the alignment positions
        if self.local:
            lead_padding = alignment.aligned[0][0][0]
            trail_padding = ref_length - alignment.aligned[0][-1][1]
            # Use the full returned aligned sequences when iterating below
            query_chunk = query_aligned
            ref_chunk = ref_aligned
        
        # For global alignment, determine padding directly from the alignment strings
        else:
            lead_padding = 0
            for x in query_aligned:
                if x == '-':
                    lead_padding += 1
                else:
                    break
            trail_padding = 0
            for x in reversed(query_aligned):
                if x == '-':
                    trail_padding += 1
                else:
                    break
            # Only iterate over the aligned portion below
            query_chunk = query_aligned[lead_padding:-trail_padding]
            ref_chunk = ref_aligned[lead_padding:-trail_padding]

        # Encode alignment
        ref_pos = lead_padding
        prev_insertion = False # Keep track of consecutive insertions
        for q_char, r_char in zip(query_chunk, ref_chunk):
            if r_char != '-':
                if q_char == '-':
                    alignment_row[ref_pos] = 3 # Deletion
                elif q_char == r_char:
                    alignment_row[ref_pos] = 0 # Match
                else:
                    alignment_row[ref_pos] = 1 # Mismatch
                ref_pos += 1
                prev_insertion = False
            else: 
                if not prev_insertion: # Start new insertion entry
                    insertions.append((ref_pos, 1))
                    prev_insertion = True
                else: # Increment length of current insertion
                    insertions[-1] = (insertions[-1][0], insertions[-1][1] + 1)

        return alignment_row, insertions

def plot_alignments(encoded_alignments, insertions_info, show_insertions=True, color_map=None):
    '''
    Plot the encoded alignments with insertions highlighted.
    
    Args:
        encoded_alignments (list): List of encoded alignments, where each alignment is a list of integers
            representing the alignment at each reference position. The integers should correspond to the
            following:
                0: Match
                1: Mismatch
                3: Deletion
                4: Padding
        insertions_info (list): List of lists, where each sublist contains tuples of (ref_pos, length) for
            each insertion in the corresponding alignment.
        show_insertions (bool): Whether to highlight insertions in the plot.
        color_map (list): List of 5 colors to use for the alignment plot. The colors should correspond to
            the following:
                0: Match
                1: Mismatch
                2: Insertion
                3: Deletion
                4: Padding
            If None, the default color map is used.
            
    '''
    if color_map is None:
        color_map = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#e0e0e0"]
    elif len(color_map) != 5:
        raise ValueError("color_map must have 5 colors.")
    color_labels = ["Match", "Mismatch", "Insertion", "Deletion", "Padding"]

    ref_length = len(encoded_alignments[0])
    alignment_df = pd.DataFrame(encoded_alignments, columns=range(ref_length))

    g = sns.clustermap(
        alignment_df, cmap=color_map, xticklabels=100, yticklabels=False,
        row_cluster=True, col_cluster=False, dendrogram_ratio=(0, 0), cbar_pos=None
    )

    if show_insertions:
        ax = g.ax_heatmap
        n_rows = alignment_df.shape[0]
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        fig_width = x_max - x_min
        fig_height = y_max - y_min
        row_height = fig_height / n_rows
        col_width = fig_width / ref_length

        reordered_row_idx = g.dendrogram_row.reordered_ind
        idx_map = {j: i for i, j in enumerate(reordered_row_idx)}

        for i in alignment_df.index:
            insertions = insertions_info[i]
            reordered_row_idx = idx_map[i] # Find the new row index after clustering

            for ref_pos, length in insertions:
                rect_width = length * col_width
                rect_start_x = ref_pos - rect_width / 2 # Center the rectangle on the reference position
                rect_start_y = reordered_row_idx + 1 # +1 because rectangles are drawn upwards
                rect = Rectangle((rect_start_x, rect_start_y), rect_width, row_height, 
                                 color=color_map[2], linewidth=0)
                ax.add_patch(rect)

        legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(color_map, color_labels)]
    else:
        legend_patches = [mpatches.Patch(color=color, label=label) for color, label in zip(color_map, color_labels) if label != "Insertion"]

    plt.legend(handles=legend_patches, loc="upper right", bbox_to_anchor=(1.2, 1.0), title="Alignment")
    g.ax_heatmap.yaxis.set_label_position('left')
    plt.xlabel("Reference Position (bp)")
    plt.ylabel("Reads")
    plt.show()