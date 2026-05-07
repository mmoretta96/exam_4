import sys
"""Return True if a sequence is long enough and does not contain numbers."""
def validate_sequence(sequence, k):
#sequences shorter than k cannot produce a complete k-mer.
    if len(sequence) < k:
        return False
#treat numeric chars as invalid sequence
    for nucleotide in sequence:
        if nucleotide in '1234567890':
            return False
    return True

"""Update the count for one k-mer and the character that follows."""
def update_kmer_count(kmer_data, kmer, next_char):
#create new dict entry for each unique k-mer
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 1, 'next_chars': {}}
    
    kmer_data[kmer]['count'] += 1
#track how often each following char appears after k-mer    
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0
    kmer_data[kmer]['next_chars'][next_char] += 1

    return kmer_data
"""Count kmers in one sequences along with each following character"""
def count_kmers_with_context(sequence, k):
    kmer_data = {}
#check sequence one char at a time to capture overlapping kmers    
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        next_char = sequence[i+k]
        
        kmer_data = update_kmer_count(kmer_data, kmer, next_char)
    
    return kmer_data

"""Write kmer following character summary to an output file"""
def write_results_to_file(kmer_data, output_filename):
#sort kmers
    sorted_kmers = sorted(kmer_data.keys())
    
    with open(output_filename, 'w') as f:
        for kmer in sorted_kmers:
            next_chars = kmer_data[kmer]['next_chars']
            
            next_char_str = " ".join(
                f"{char}:{freq}" 
                for char, freq in sorted(next_chars.items())
            )
            
            f.write(f"{kmer} {next_char_str}\n")

"""Read cmd line args, analyze sequences, and write results to a file"""
def main():
#read cmd line args 
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]
    
    print(f"Reading sequences from {sequence_file}...")

    with open(sequence_file, 'r') as f:
        for sequence in f:
            sequence = sequence.strip()

            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence")
                continue
            
            kmer_data = count_kmers_with_context(sequence, k) 
            
            write_results_to_file(kmer_data, output_file)

if __name__ == '__main__':
    main()
