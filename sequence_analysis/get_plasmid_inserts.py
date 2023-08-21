# This script extracts inserts from Plasmidsaurus raw whole plasmid Nanopore reads.

import argparse
import pickle
from sequence_analysis.utils import levenshtein, hamming_dist, read_fastq, revcomp


def all_kmers(seq, k):
    '''Return all kmers of length k in seq.'''
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]


def first_index_of_kmer(seq, kmer, dist, dist_func):
    '''Return the index of the first kmer in seq that matches
    the kmer argument within a distance of the dist argument.'''
    k = len(kmer)
    kmers = all_kmers(seq, k)
    for i, s in enumerate(kmers):
        if dist_func(kmer, s) <= dist:
            return i
    return -1


def best_index_of_kmer(seq, kmer, dist, dist_func):
    '''Return the index of the kmer in seq that most closely matches
    the kmer argument. If no match within dist, then return -1.'''
    k = len(kmer)
    kmers = all_kmers(seq, k)
    lowest_distance = dist + 1
    best_index = -1
    for i, s in enumerate(kmers):
        distance = dist_func(kmer, s)
        if distance <= dist and distance < lowest_distance:
            lowest_distance = distance
            best_index = i
    return best_index


def get_insert(read, up, down, length, t, dist_func):
    '''Return an insert sequence flanked by upstream and downstream sequence.
        Parameters:
            read ((str, str, str)): A (header, sequence, quality) tuple.
            up (str): Sequence upstream of the insert.
            down (str): Sequence downstream of the insert.
            length (int): Expected maximum length of the insert.
            t (int): Distance threshold for matching up/downstream.
            dist_func (returns int): Distance function for matching up/downstream.
        Returns:
            A tuple containing three (sequence, quality) tuples in the order of
            upstream match, insert, and downstream match, or None if no match found.
    '''
    _, seq, qual = read
    # concatenate each sequence to itself in case the read starts
    # in the middle of the insert site
    seq = seq + seq
    qual = qual + qual
    
    # get first instance of upstream sequence and the sequence downstream of it
    index = first_index_of_kmer(seq, up, t, dist_func)
    if index >= 0:
        insert_region = seq[index:index + length * 2 + len(up) + len(down)]
        insert_qual = qual[index:index + length * 2 + len(up) + len(down)]
    else:
        index_rev = first_index_of_kmer(revcomp(seq), up, t, dist_func)
        if index_rev >= 0:
            insert_region = revcomp(seq)[index_rev:index_rev + length * 2 + len(up) + len(down)]
            insert_qual = qual[::-1][index_rev:index_rev + length * 2 + len(up) + len(down)]
        else:
            return None
    # find best match of upstream and downstream sequence and extract insert
    up_idx = best_index_of_kmer(insert_region, up, t, dist_func)
    down_idx = best_index_of_kmer(insert_region, down, t, dist_func)
    if down_idx >= 0:
        up_match = (insert_region[up_idx:up_idx + len(up)], insert_qual[up_idx:up_idx + len(up)])
        insert_match = (insert_region[up_idx + len(up): down_idx], insert_qual[up_idx + len(up): down_idx])
        down_match = (insert_region[down_idx:down_idx + len(down)], insert_qual[down_idx:down_idx + len(down)])
        return (up_match, insert_match, down_match)
    else:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="path to fastq file (unzipped)")
    parser.add_argument("-u", "--upstream", required=True, type=str, help="sequence upstream of insert")
    parser.add_argument("-d", "--downstream", required=True, type=str, help="sequence downstream of insert")
    parser.add_argument("-l", "--length", required=True, type=int, help="expected maximum insert length")
    parser.add_argument("-t", "--threshold", default=8, type=int, help="threshold (maximum distance) for upstream and downstream matching")
    parser.add_argument("-f", "--dist-func", default='levenshtein', choices = ['levenshtein', 'hamming'], type=str, help="distance function for upstream and downstream matching")
    parser.add_argument("-n", "--no-empty", action='store_true', help="don't save matches with no insert")
    parser.add_argument("-r", "--revcomp", action='store_true', help="also print reverse complement of inserts")
    parser.add_argument("-q", "--quality", action='store_true', help="also print quality of inserts and upstream/downstream")
    parser.add_argument("-o", "--output", default="inserts", type=str, help="file name for output")
    args = parser.parse_args()

    upstream = args.upstream.upper()
    downstream = args.downstream.upper()

    if args.dist_func == 'levenshtein':
        try:
            import Levenshtein
            dist_func = Levenshtein.distance
        except ModuleNotFoundError:
            print("Levenshtein module (https://github.com/maxbachmann/Levenshtein) not found. "
                    "A much slower implementation of Levenshtein distance will be used.")
            dist_func = levenshtein
    else:
        dist_func = hamming_dist

    # read in sequences from fastq
    with open(args.input, 'r') as f:
        fastq_reader = read_fastq(f)
        reads = [read for read in fastq_reader]
    
    # extract inserts
    inserts = []
    for read in reads:
        insert = get_insert(read, upstream, downstream, length=args.length,
                           t=args.threshold, dist_func=dist_func)
        if insert:
            inserts.append(insert)

    print(f'Number of reads: {len(reads)}')
    print(f'Number of matches found: {len(inserts)}')

    if args.no_empty:
        inserts = [insert for insert in inserts if len(insert[1][0]) > 0]
        print(f'Number of non-empty inserts found: {len(inserts)}')
    
    # write sequences to table
    insert_width = max([len(insert[0]) for up, insert, down in inserts])
    up_width = len(upstream)
    down_width = len(downstream)
    table = [f"{'Upstream':<{up_width}}\t{'Insert':<{insert_width}}\t{'Downstream':<{down_width}}"]
    for upstream, insert, downstream in inserts:
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

    print(f'Saving list of inserts to {args.output}.pkl')
    with open(f'{args.output}.pkl', 'wb') as f:
        pickle.dump(inserts, f)

if __name__ == '__main__':
    main()
