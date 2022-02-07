"""Test module constants."""
from pathlib import Path

TEST_FOLDER = Path(__file__).parent
dna_seq = str(TEST_FOLDER.joinpath("dna_seq.txt"))
fasta_file = str(TEST_FOLDER.joinpath("fasta_dummy.fasta"))
fastq_file = str(TEST_FOLDER.joinpath("fastq_dummy.fastq"))
kmers_seqs = str(TEST_FOLDER.joinpath("kmer_seqs.txt"))

result1_file = str(TEST_FOLDER.joinpath('result1.txt'))
result2_file = str(TEST_FOLDER.joinpath('result2.csv'))
result3_file = str(TEST_FOLDER.joinpath('result3.tsv'))
wrong_path = str(TEST_FOLDER.joinpath('result4.pdf'))

