from itertools import islice

# stolen from https://gist.github.com/jakebiesinger/759018/1b7d6bd6967780a8bbae743760c37885bdf86467
def read_fastq(fastqfile, skip_blank=True):
    '''Parse a fastq-formatted file, yielding a (header, sequence, quality) tuple'''
    fastqiter = (l.strip('\n') for l in fastqfile)  # strip trailing newlines
    if skip_blank:
        fastqiter = filter(lambda l: l, fastqiter)  # skip blank lines
    while True:
        fqlines = list(islice(fastqiter, 4))
        if len(fqlines) == 4:
            header1, seq, header2, qual = fqlines
        elif len(fqlines) == 0:
            return
        else:
            raise EOFError("Failed to parse four lines from fastq file!")

        if header1.startswith('@') and header2.startswith('+'):
            yield header1, seq, qual
        else:
            raise ValueError("Invalid header lines: %s and %s for seq %s" % (header1, header2, seq))
        
def read_fasta(fastafile):
    """
    Reads a FASTA file and yields (header, sequence) tuples.

    Args:
        filepath (str): Path to the FASTA file.

    Yields:
        tuple: (header, sequence) where:
            - header (str) is the FASTA header without the '>'.
            - sequence (str) is the full sequence for that header.

    Raises:
        ValueError: If the file is not properly formatted.
    """
    header = None
    sequence = []

    for line in fastafile:
        line = line.strip()
        if not line:
            continue  # Ignore empty lines

        if line.startswith(">"):
            # If we already have a sequence, yield it
            if header:
                if not sequence:
                    raise ValueError(f"FASTA header '{header}' has no sequence.")
                yield (header, "".join(sequence))

            # Start a new sequence
            header = line[1:].strip()
            sequence = []
        else:
            if header is None:
                raise ValueError("FASTA file must start with a header (line beginning with '>').")
            sequence.append(line)

    # Yield the last sequence if present
    if header:
        if not sequence:
            raise ValueError(f"FASTA header '{header}' has no sequence.")
        yield (header, "".join(sequence))

COMP_TABLE = str.maketrans("ACTGN", "TGACN")

def revcomp(seq):
    '''Return the reverse complement of a DNA sequence'''
    if not set(seq).issubset({'A', 'C', 'G', 'T', 'N'}):
        raise ValueError(f"Sequence ({seq}) must only contain ACTGN")
    return seq.translate(COMP_TABLE)[::-1]

def revcomp_read(read):
    header, seq, qual = read
    return header, revcomp(seq), qual[::-1]

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
    if len(a) != len(b):
        raise ValueError("Sequences must be of equal length")
    return sum(i != j for i, j in zip(a, b))

def is_DNA(seq):
    return set(seq).issubset({'A', 'C', 'G', 'T', 'N'})
