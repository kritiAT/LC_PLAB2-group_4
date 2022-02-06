""" Tests for sequence_assembly module. """

from sequence_assembly import MSA, Assembly, FastaTools
import pytest
import os

from itertools import combinations

check_seq1 = 'CTTTATTTTCACCACAAACAGAGATTAAAGAAAGTGTTGGATTAAAAGATGGTGTTAAAGAGTAAAAATTGACTTATTATACTCCTGAATACGAAACCAAAGATAC\
TGATATCTTGGCAGCATTCCGAGTAACTCCTCAACCTGGAGTTCCACCCGAAGAAGCAGGGGCCGCGGTAGCTGCCGAATCTTCTACCGGTACA'


class TestFastaTools:
    """ Test class for FastaTools class in sequence assembly module. """

    def test_fasta_list(self):
        """ Checks if input FASTA files are read correctly. """
        check = 'GGTTAAGATGATCAATTAACAA'
        fasta = FastaTools().fasta_list('fasta_dummy.fasta')

        assert len(fasta) == 2
        assert fasta[1] == check

    def test_fastq_list(self):
        """ Checks if input FASTQ files are read correctly. """
        check = 'ATATCTAGGGCATAT'
        fastq = FastaTools().fastq_list('fastq_dummy.fastq')

        assert len(fastq) == 4
        assert fastq[2] == check


class TestMSA:
    """ Test class for MSA class in sequence assembly module. """

    def test_pairwise_alignment(self):
        """ Checks weather pairwise alignment algorithm correctly aligns identical sequences. """
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

    def test_perform_msa(self):
        """ Tests if multiple sequence alignment functions as expected. """
        test_seqs = ['CTGGTAGCAGTAGCGTTTAAAGCGC',
                     'GCAGTAGCGTTTAAAGCGCAATGCA',
                     'GAAATGCTATTCAATGCATTGATGC']

        pairs = [(i, j) for i, j in combinations(test_seqs, 2)]
        scores = MSA().perform_msa(sequences=test_seqs)

        assert len(scores) == 3
        assert list(scores.keys()) == pairs
        assert scores[pairs[1]] == (100.0, 3, 'GCA')


class TestAssembly:
    """ Test class for Assembly class in sequence assembly module. """

    def test_de_novo_assembly(self):
        """ Checks weather the reconstructed sequence is correctly assembled.
        Tests if input files are accepted extensions are accepted.
        Tests if output saved to file with correct format. """
        input_file1 = 'kmer_seqs.txt'
        input_file2 = 'fastq_dummy.fastq'
        input_file3 = 'fasta_dummy.fasta'
        check_seq2 = 'CCCATATCTAGGGCATATTTACTCCCGTATCA'
        check_seq3 = 'AAAAGAGTCTAAAGGTTAAGATGATCAATTAACAA'
        out_file = 'outfile.txt'
        wrong_path = 'file.png'

        obj = Assembly(sequences=input_file1, output_path=out_file)
        out_seq = obj.assembled_sequence

        assert out_seq == check_seq1
        assert os.path.exists(out_file) is True  # checks if output file exists
        os.remove(out_file)

        obj2 = Assembly(sequences=input_file2)
        seq2 = obj2.assembled_sequence

        assert seq2 == check_seq2

        obj3 = Assembly(sequences=input_file3)
        seq3 = obj3.assembled_sequence

        assert seq3 == check_seq3

        # check if raises error for wrong file format
        with pytest.raises(ValueError, match=r".* not supported .*"):
            Assembly(sequences=input_file1, output_path=wrong_path)

        with pytest.raises(ValueError, match=r".* not supported .*"):
            Assembly(sequences=wrong_path)
