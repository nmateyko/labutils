# This script extracts inserts from raw whole plasmid Nanopore reads.

import argparse
import edlib
from labutils.sequence_analysis.utils import read_fastq, revcomp, is_DNA


def get_align_pos(seq, kmer, dist):
    '''Return the start and end index of the best alignment of kmer in seq.
        Parameters:
            seq (str): A sequence to search in.
            kmer (str): A sequence to search for.
            dist (int): Maximum distance for alignment.
        Returns:
            ((int, int)): Start and end indices of the first best alignment of
            kmer in seq, or -1 if no alignment found.
    '''
    result = edlib.align(kmer, seq, mode="HW", task="locations", k=dist)
    return result["locations"][0] if result["locations"] else -1


def get_insert_pos(read, up, down, t):
    '''Return position of an insert sequence flanked by upstream and downstream sequence.
        Parameters:
            read ((str, str, str)): A (header, sequence, quality) tuple.
            up (str): Sequence upstream of the insert.
            down (str): Sequence downstream of the insert.
            t (int): Distance threshold for matching up/downstream.
        Returns:
            A dictionary with the following keys:
                "up": (int, int) - Start and end indices of the upstream sequence.
                "down": (int, int) - Start and end indices of the downstream sequence.
                "rc": bool - True if the insert is in reverse complement orientation.
    '''
    _, seq, _ = read

    up_match = get_align_pos(seq, up, t)
    rc = False
    if up_match == -1:
        # Try reverse complement if not found in forward orientation
        seq = revcomp(seq)
        up_match = get_align_pos(seq, up, t)
        if up_match == -1:
            return None
        rc = True

    seq_remain = seq[up_match[1]:] # only search downstream of up_match!
    down_match = get_align_pos(seq_remain, down, t)
    if down_match == -1:
        return None
    
    # account for offset from up_match
    down_match = (down_match[0] + up_match[1], down_match[1] + up_match[1])
    return {"up": up_match, "down": down_match, "rc": rc}
        
def extract_insert_and_flanks(read, pos_dict):
    '''Extract insert sequence from read based on position dictionary.
        Parameters:
            read ((str, str, str)): A (header, sequence, quality) tuple.
            pos_dict (dict): A dictionary with the following keys:
                "up": (int, int) - Start and end indices of the upstream sequence.
                "down": (int, int) - Start and end indices of the downstream sequence.
                "rc": bool - True if the insert is in reverse complement orientation.
        Returns:
            A tuple of three (sequence, quality) tuples for upstream, insert, and downstream sequences.
    '''
    _, seq, qual = read
    if pos_dict['rc']:
        seq = revcomp(seq)
        qual = qual[::-1]
    up_start, up_end = pos_dict["up"]
    down_start, down_end = pos_dict["down"]
    upstream = (seq[up_start:up_end + 1], qual[up_start:up_end + 1])
    insert = (seq[up_end + 1:down_start], qual[up_end + 1:down_start])
    downstream = (seq[down_start:down_end + 1], qual[down_start:down_end + 1])
    return upstream, insert, downstream

def extract_insert_fastq(read, pos_dict):
    '''Extract insert sequence from read based on position dictionary.
        Parameters:
            read ((str, str, str)): A (header, sequence, quality) tuple.
            pos_dict (dict): A dictionary with the following keys:
                "up": (int, int) - Start and end indices of the upstream sequence.
                "down": (int, int) - Start and end indices of the downstream sequence.
                "rc": bool - True if the insert is in reverse complement orientation.
        Returns:
            A (header, sequence, quality) tuple for the insert sequence. Header
            is modified to include the position of the insert in the original read
            and the orientation of the insert.
    '''
    header, seq, qual = read
    if pos_dict['rc']:
        seq = revcomp(seq)
        qual = qual[::-1]
    up_start, up_end = pos_dict["up"]
    down_start, down_end = pos_dict["down"]
    insert = (seq[up_end + 1:down_start], qual[up_end + 1:down_start])
    insert_header = f'{header} insert_pos={up_end + 1}-{down_start} rc={pos_dict["rc"]}'
    return insert_header, insert[0], insert[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="path to fastq file (unzipped)")
    parser.add_argument("-u", "--upstream", required=True, type=str, help="sequence upstream of insert")
    parser.add_argument("-d", "--downstream", required=True, type=str, help="sequence downstream of insert")
    parser.add_argument("-t", "--threshold", default=None, type=int, help="threshold (maximum distance) for upstream and downstream matching")
    parser.add_argument("-n", "--no-empty", action='store_true', help="don't save matches with no insert in table")
    parser.add_argument("-r", "--revcomp", action='store_true', help="also print reverse complement of inserts in table")
    parser.add_argument("-q", "--quality", action='store_true', help="also print quality of inserts and upstream/downstream in table")
    parser.add_argument("-o", "--output", default="inserts", type=str, help="file name for output")
    args = parser.parse_args()

    upstream = args.upstream.upper()
    downstream = args.downstream.upper()

    if not is_DNA(upstream) or not is_DNA(downstream):
        raise ValueError("Upstream and downstream sequences must only contain ACTGN")

    if args.threshold is None:
        threshold = round(min(len(upstream), len(downstream)) * 0.05)
        print(f'Using default distance threshold of {threshold}')
    elif args.threshold < 0:
        raise ValueError("Threshold must be a non-negative integer")
    else:
        threshold = args.threshold

    # read in sequences from fastq
    with open(args.input, 'r') as f:
        fastq_reader = read_fastq(f)
        reads = [read for read in fastq_reader]
    
    # extract inserts
    extracted = []
    fastq_inserts = []
    for read in reads:
        # concat read to itself in case plasmid linearized within insert
        concat_read = (read[0], read[1] + read[1], read[2] + read[2]) 
        pos_dict = get_insert_pos(concat_read, upstream, downstream, t=threshold)
        extract = extract_insert_and_flanks(concat_read, pos_dict) if pos_dict else None
        if extract:
            extracted.append(extract)
            fastq_inserts.append(extract_insert_fastq(read, pos_dict))

    print(f'Number of reads: {len(reads)}')
    print(f'Number of matches found: {len(extracted)}')

    if args.no_empty:
        extracted = [x for x in extracted if len(x[1][0]) > 0]
        print(f'Number of non-empty inserts found: {len(extracted)}')
    
    # write sequences to table
    insert_width = max([len(insert[0]) for up, insert, down in extracted])
    up_width = max([len(up[0]) for up, insert, down in extracted])
    down_width = max([len(down[0]) for up, insert, down in extracted])
    table = [f"{'Upstream':<{up_width}}\t{'Insert':<{insert_width}}\t{'Downstream':<{down_width}}"]
    for upstream, insert, downstream in extracted:
        table.append(f"{upstream[0]:<{up_width}}\t{insert[0]:<{insert_width}}\t{downstream[0]:<{down_width}}")
        if args.quality:
            table.append(f"{upstream[1]:<{up_width}}\t{insert[1]:<{insert_width}}\t{downstream[1]:<{down_width}}\n")
        if args.revcomp:
            table.append(f"{'':<{up_width}}\t{revcomp(insert[0]):<{insert_width}}\t{'':<{down_width}}")
            if args.quality:
                table.append(f"{'':<{up_width}}\t{insert[1][::-1]:<{insert_width}}\t{'':<{down_width}}")
            table.append('')

    print(f'Saving inserts table to {args.output}.txt')
    with open(f'{args.output}.txt', 'w') as f:
        f.write('\n'.join(table))

    print(f'Saving fastq of inserts to {args.output}.fastq')
    with open(f'{args.output}.fastq', 'w') as f:
        for header, seq, qual in fastq_inserts:
            f.write(f'{header}\n{seq}\n+\n{qual}\n')

if __name__ == '__main__':
    main()
