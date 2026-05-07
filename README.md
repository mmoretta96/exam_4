This repo contains .py script for analyzing DNA sequence fragments.
Goal is to count each k-mer and record the frequency of the character that appears after each k-mer.
# UPDATED Genome K-mer Analyzer 

This project analyzes DNA sequence fragments from a text file. For a user-specified value of `k`, the script counts each overlapping k-mer and records how often each character appears immediately after that k-mer.

## Files

`kmer_analyzer.py`: Python script for analyzing sequence fragments.
`tests/test_kmer_analyzer.py`: Pytest test file for the script functions.
`tests/conftest.py`: Test configuration so pytest can import the script.
`sample_data/fragments.txt`: Example input file.
`sample_data/expected_k2.csv`: Example expected output for `k = 2`.

## How to Run the Script

Use the command line format:

```bash
python3 kmer_analyzer.py input_file k output_file

## AI Use Statement
I used AI to help better understand assignment requirements, troubleshoot Git/GitHub workflow, and help design better pytest tests when identifiny bugs in the starter script.


