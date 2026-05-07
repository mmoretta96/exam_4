import sys


def validate_sequence(sequence, k):
    """Return True if a sequence is long enough and does not contain numbers."""
    # Sequences shorter than k cannot produce a complete k-mer.
    if len(sequence) < k:
        return False

    # Remove lead/trail whitespace and convert sequence to uppercase.
    sequence = sequence.strip().upper()

    for nucleotide in sequence:
        if nucleotide not in 'ACGT':
            return False
    return True


def update_kmer_count(kmer_data, kmer, next_char=None):
    """Update the count for one k-mer and the character that follows."""
    # Create new dict entry for each unique k-mer.
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}

    kmer_data[kmer]['count'] += 1

    # Track how often each following char appears after k-mer.
    if next_char is not None:  # Run only when there is a next char.
        if next_char not in kmer_data[kmer]['next_chars']:
            kmer_data[kmer]['next_chars'][next_char] = 0
        kmer_data[kmer]['next_chars'][next_char] += 1

    return kmer_data


def count_kmers_with_context(sequence, k):
    """Count kmers in one sequence along with each following character."""
    kmer_data = {}

    # Check sequence one char at a time to capture overlapping kmers.
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if i + k < len(sequence):
            next_char = sequence[i+k]
        else:
            next_char = None

        kmer_data = update_kmer_count(kmer_data, kmer, next_char)

    return kmer_data


def write_results_to_file(kmer_data, output_filename):
    """Write kmer following character summary to an output file."""
    # Sort kmers.
    sorted_kmers = sorted(kmer_data.keys())

    with open(output_filename, 'w') as f:
        for kmer in sorted_kmers:
            next_chars = kmer_data[kmer]['next_chars']

            next_char_str = " ".join(
                f"{char}:{freq}"
                for char, freq in sorted(next_chars.items())
            )

            if next_char_str:
                f.write(f"{kmer} {kmer_data[kmer]['count']} {next_char_str}\n")
            else:
                f.write(f"{kmer} {kmer_data[kmer]['count']}\n")


def main():
    """Read cmd line args, analyze sequences, and write results to a file."""
    # Read cmd line args.
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]

    print(f"Reading sequences from {sequence_file}...")
    all_kmer_data = {}

    with open(sequence_file, 'r') as f:
        for sequence in f:
            sequence = sequence.strip()

            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence")
                continue

            kmer_data = count_kmers_with_context(sequence, k)

            for kmer, data in kmer_data.items():
                if kmer not in all_kmer_data:
                    all_kmer_data[kmer] = {'count': 0, 'next_chars': {}}

                all_kmer_data[kmer]['count'] += data['count']

                for next_char, freq in data['next_chars'].items():
                    if next_char not in all_kmer_data[kmer]['next_chars']:
                        all_kmer_data[kmer]['next_chars'][next_char] = 0
                    all_kmer_data[kmer]['next_chars'][next_char] += freq

    write_results_to_file(all_kmer_data, output_file)


if __name__ == '__main__':
    main()

