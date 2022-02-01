from sequence_assembly import MSA, Assembly
import pytest
import os

from itertools import combinations

final_seq = 'CTTTATTTTCACCACAAACAGAGATTAAAGAAAGTGTTGGATTAAAAGATGGTGTTAAAGAGTAAAAATTGACTTATTATACTCCTGAATACGAAACCAAAGATAC\
TGATATCTTGGCAGCATTCCGAGTAACTCCTCAACCTGGAGTTCCACCCGAAGAAGCAGGGGCCGCGGTAGCTGCCGAATCTTCTACCGGTACA'


class TestMSA:
    """ Test class for MSA class in sequence assembly module."""

    def test_pairwise_alignment(self):
        test_seq1 = 'CTGGTAGCAGTAGCGTTTAAAGCGC'
        test_seq2 = 'GCAGTAGCGTTTAAAGCGCAATGCA'
        test_seq3 = 'ctttgagcgtttaaagcagtgccgg'

        score, identity, aligned_seq = MSA().pairwise_alignment(seq1=test_seq1, seq2=test_seq2)

        assert score == 19
        assert identity == 100.0
        assert aligned_seq == 'GCAGTAGCGTTTAAAGCGC'

        score2, identity2, aligned_seq2 = MSA().pairwise_alignment(seq1=test_seq2, seq2=test_seq3)

        assert score2 == 12
        assert identity2 == 100.0
        assert aligned_seq2 == 'AGCGTTTAAAGC'

    def test_perform_mas(self):
        test_seqs = ['CTGGTAGCAGTAGCGTTTAAAGCGC',
                     'GCAGTAGCGTTTAAAGCGCAATGCA',
                     'GAAATGCTATTCAATGCATTGATGC']

        pairs = [(i, j) for i, j in combinations(test_seqs, 2)]
        scores = MSA().perform_msa(sequences=test_seqs)

        assert len(scores) == 3
        assert list(scores.keys()) == pairs
        assert scores[pairs[1]] == (100.0, 3, 'GCA')


class TestAssembly:
    """ Test class for MSA class in sequence assembly module."""

    def test_de_novo_assembly(self):
        input_file = 'kmer_seqs.txt'
        out_file = 'outfile.txt'
        wrong_path = 'outfile.png'

        obj = Assembly(sequences=input_file, output_path=out_file)
        out_seq = obj.assembled_sequence

        assert out_seq == final_seq
        assert os.path.exists(out_file) is True  # checks if output file exists
        os.remove(out_file)

        # check if raises error for wrong file format
        with pytest.raises(ValueError, match=r".* not supported .*"):
            Assembly(sequences=input_file, output_path=wrong_path)
