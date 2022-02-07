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

# if __name__ == '__main__':
#     print(TEST_FOLDER, 'jhk', type(TEST_FOLDER))
#     print(dna_seq, fastq_file, fasta_file, kmers_seqs, result3_file, result2_file, result1_file, wrong_path)
#     print(type(str(dna_seq)), 'hsdcbs', str(dna_seq))
#     print(f'hello{dna_seq}')
#     print(f'hagjhv{str(dna_seq)}')