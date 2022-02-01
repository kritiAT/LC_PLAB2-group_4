import os
import pytest

from gene_finder import Transcribe


class TestTranscribe:
    """ Test class for Transcribe class in gene finder module."""

    def test_gene_sequences(self):
        dna_sequence = 'dna_seq.txt'

        obj = Transcribe(dna_sequence, reverse=True, threshold=3, output_path='result1.txt')
        result1 = obj.genes
        pos, seq = list(result1['reverse'].items())[1]

        assert pos == (89, 107)
        assert seq == 'AUGCUUUGGUUUCUAUGA'
        assert os.path.exists('result1.txt') is True  # checks if output file exists
        os.remove('result1.txt')

        obj2 = Transcribe(dna_sequence, reverse=False, threshold=3, output_path='result2.csv')
        result2 = obj2.genes

        assert result2['reverse'] is None
        assert os.path.exists('result2.csv') is True  # checks if output file exists
        os.remove('result2.csv')

        obj3 = Transcribe(dna_sequence, reverse=True, threshold=25, output_path='result3.tsv')
        result3 = obj3.genes

        assert len(result3['forward']) == 0
        assert os.path.exists('result3.tsv') is True  # checks if output file exists
        os.remove('result3.tsv')

        # check if raises error for wrong file format
        with pytest.raises(ValueError, match=r".* not supported .*"):
            Transcribe(dna_sequence, output_path='result4.pdf')

        with pytest.raises(ValueError, match=r".* not supported .*"):
            Transcribe(dna='test_seq.png')

    def test_store_genes(self):
        dna_sequence = 'dna_seq.txt'

        with pytest.raises(ValueError, match=r".* not supported .*"):
            Transcribe(dna_sequence).store_genes(path='result5.png')
