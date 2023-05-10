# This script extracts inserts from Plasmidsaurus raw whole plasmid Nanopore reads.

import argparse
import pickle
from itertools import islice

def levenshtein(s1, s2):
    '''
    Calculates Levenshtein distance between two sequences. Much slower than the Levenshtein module function.
    Source: https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python.

        Parameters:
            s1 (str): First sequence for distance calculation.
            s2 (str): Second sequence for distance calculation.

        Returns:
            (int): Levenshtein distance between s1 and s2.
    '''
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]

def hamming_dist(a, b):
  return sum(i != j for i, j in zip(a, b))

# stolen from https://gist.github.com/jakebiesinger/759018/1b7d6bd6967780a8bbae743760c37885bdf86467
def readFastq(fastqfile):
    "parse a fastq-formatted file, yielding a (header, sequence, quality) tuple"
    fastqiter = (l.strip('\n') for l in fastqfile)  # strip trailing newlines 
    fastqiter = filter(lambda l: l, fastqiter)  # skip blank lines
    while True:
        fqlines = list(islice(fastqiter, 4))
        if len(fqlines) == 4:
            header1,seq,header2,qual = fqlines
        elif len(fqlines) == 0:
            return
        else:
            raise EOFError("Failed to parse four lines from fastq file!")

        if header1.startswith('@') and header2.startswith('+'):
            yield header1[1:], seq, qual
        else:
            raise ValueError("Invalid header lines: %s and %s for seq %s" % (header1, header2, seq))

def rev_comp(seq):
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(comp[base] for base in seq[::-1])

def all_kmers(seq, k):
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]

def first_index_of_kmer(seq, kmer, dist, dist_func):
  k = len(kmer)
  kmers = all_kmers(seq, k)
  for i, s in enumerate(kmers):
    if dist_func(kmer, s) <= dist:
      return i
  return -1

def best_index_of_kmer(seq, kmer, dist, dist_func):
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
    parsed_fastq = readFastq(f)
    reads = [seq for _, seq, _ in parsed_fastq]

  # concatenate each sequence to itself in case the read starts in the middle of the insert site
  reads_concat = [i + i for i in reads]

  # get first instance of upstream sequence and the sequence downstream of it
  insert_regions = []
  for seq in reads_concat:
    index = first_index_of_kmer(seq, upstream, args.threshold, dist_func)
    if index >= 0:
      insert_regions.append(seq[index:index + args.length * 2 + len(upstream) + len(downstream)])
    else:
      index_rev = first_index_of_kmer(rev_comp(seq), upstream, args.threshold, dist_func)
      if index_rev >= 0:
        insert_regions.append(rev_comp(seq)[index_rev:index_rev + args.length * 2 + len(upstream) + len(downstream)])

  # find best match of upstream and downstream sequence and extract insert
  inserts = []
  for seq in insert_regions:
    up_idx = best_index_of_kmer(seq, upstream, args.threshold, dist_func)
    down_idx = best_index_of_kmer(seq, downstream, args.threshold, dist_func)
    if down_idx >= 0:
      inserts.append((seq[up_idx:up_idx + len(upstream)], seq[up_idx + len(upstream): down_idx], seq[down_idx:down_idx + len(downstream)]))

  print(f'Number of matches found: {len(inserts)}')

  if args.no_empty:
    inserts = [insert for insert in inserts if len(insert[1]) > 0]
    print(f'Number of non-empty inserts found: {len(inserts)}')

  insert_width = max([len(insert) for up, insert, down in inserts])
  up_width = len(upstream)
  down_width = len(downstream)
  table = [f"{'Upstream':<{up_width}}\t{'Insert':<{insert_width}}\t{'Downstream':<{down_width}}"]
  for upstream, insert, downstream in inserts:
    table.append(f"{upstream:<{up_width}}\t{insert:<{insert_width}}\t{downstream:<{down_width}}")
    if args.revcomp:
      table.append(f"{'':<{up_width}}\t{rev_comp(insert):<{insert_width}}\t{'':<{down_width}}\n")

  print(f'Saving inserts table to {args.output}.txt')
  with open(f'{args.output}.txt', 'w') as f:
    f.write('\n'.join(table))

  print(f'Saving list of inserts to {args.output}.pkl')
  with open(f'{args.output}.pkl', 'wb') as f:
    pickle.dump(inserts, f)

if __name__ == '__main__':
  main()
