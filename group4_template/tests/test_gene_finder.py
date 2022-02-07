""" Test module for gene_finder module. """

import os
import pytest
from .constants import *

from group4.gene_finder import Transcribe, Translate


class TestTranscribe:
    """ Test class for Transcribe class in gene finder module."""

    def test_gene_sequences(self):
        obj = Transcribe(dna_seq, reverse=True, threshold=3, output_path=result1_file)
        result1 = obj.genes
        pos, seq = list(result1['reverse'].items())[1]

        assert pos == (89, 107)
        assert seq == 'AUGCUUUGGUUUCUAUGA'
        assert os.path.exists(result1_file) is True  # checks if output file exists
        os.remove(result1_file)

        obj2 = Transcribe(dna_seq, reverse=False, threshold=3, output_path=result2_file)
        result2 = obj2.genes

        assert result2['reverse'] is None
        assert os.path.exists(result2_file) is True  # checks if output file exists
        os.remove(result2_file)

        obj3 = Transcribe(dna_seq, reverse=True, threshold=25, output_path=result3_file)
        result3 = obj3.genes

        assert len(result3['forward']) == 0
        assert os.path.exists(result3_file) is True  # checks if output file exists
        os.remove(result3_file)

        # check if raises error for wrong file format
        with pytest.raises(ValueError, match=r".* not supported .*"):
            Transcribe(dna_seq, output_path=wrong_path)

        with pytest.raises(ValueError, match=r".* not supported .*"):
            Transcribe(dna=wrong_path)

    def test_store_genes(self):
        with pytest.raises(ValueError, match=r".* not supported .*"):
            Transcribe(dna_seq).store_genes(path=wrong_path)


class TestTranslate:
    """ Test class for Translate class in gene finder module."""

    def test_translation(self):
        obj = Translate(dna=dna_seq, reverse=True, threshold=3, output_path=result2_file)
        table = obj.proteins_table
        proteins = obj.proteins
        testlist = ['SLTTVL', 'LWFL']

        assert table.shape == (4, 3)
        assert [proteins[1], proteins[3]] == testlist
        assert os.path.exists(result2_file) is True  # checks if output file exists
        os.remove(result2_file)

        with pytest.raises(AssertionError, match=r".* should be greater .*"):
            Translate(dna=dna_seq, reverse=True, threshold=-3)
