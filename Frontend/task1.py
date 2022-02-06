import numpy as np
from itertools import combinations
from copy import deepcopy


class FastaTools:
    """ Tools for dealing with fasta files. """
    def __init__(self):
        pass

    @staticmethod
    def fasta_list(filepath: str):
        """
        Extracts the sequences from a FASTA file
        Args:
            filepath: (str) path to file

        Returns:
            (list) list of sequences

        """
        with open(filepath, 'r') as data:
            file = data.read()
            sequences = file.split('>')[1:]
            list_seqs = []
            for sequence in sequences:
                seq = sequence.partition('\n')[2].replace('\n', '')
                list_seqs.append(seq)
            return list_seqs

    @staticmethod
    def fastq_list(filepath: str):
        """
        Extracts the sequences from a FASTQ file
        Args:
            filepath: (str) path to file

        Returns:
            (list) list of sequences

        """
        with open(filepath, 'r') as data:
            file = data.read()
            records = file.split('@')[1:]
            return [record.split('\n')[1] for record in records]


class MSA:
    """ Tools for performing sequence alignment. """

    def __init__(self):
        # Scoring matrix for aligning identical fragments nucleotide sequence.
        self.scoring_matrix = {'AA': 1, 'AT': -4, 'AG': -4, 'AC': -4, 'TT': 1,
                               'GT': -4, 'CT': -4, 'GG': 1, 'CG': -4, 'CC': 1}

    def pairwise_alignment(self, seq1: str, seq2: str, print_output: bool = False) -> tuple:
        """
        Performs pairwise sequence alignment using local alignment method.
        Smith-Waterman algorithm is used with some modifications to align identical sequence ends.
        Args:
            seq1: (str) sequence 1
            seq2: (str) sequence 2
            print_output: (bool) option to print formatted output. Default False

        Returns:
            tuple(int, float, str)
            alignment_score: (int) highest score of alignmnet
            seq_identity: (float) percentage of identical matches
            aligned_seq1: (str) matched fragment on sequence

        """
        seq1 = seq1.upper()
        seq2 = seq2.upper()

        n = len(seq1)
        m = len(seq2)

        # INITIALIZE MATRIX
        mat = np.zeros((n + 1, m + 1), dtype=int)
        for j in range(m + 1):  # first row with gap penalties
            mat[0, j] = j * -4
        for i in range(n + 1):  # first column with gap penalties
            mat[i, 0] = i * -4

        # FILL MATRIX
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                pair = ''.join(sorted((seq1[i - 1], seq2[j - 1])))
                diag = mat[i - 1, j - 1] + self.scoring_matrix[pair]
                up = mat[i - 1, j] + -4
                left = mat[i, j - 1] + -4
                mat[i, j] = max(diag, up, left, 0)  # replacing all negative values with zero

        alignment_score = np.amax(mat)

        # TRACEBACK ALIGNMENT
        alseq1 = ''
        alseq2 = ''  # storing the gaps in alignment

        max_index = np.where(mat == alignment_score)
        i, j = max_index[0][0], max_index[1][0]  # starting from max value
        while mat[i, j] != 0:
            key = ''.join(sorted((seq1[i - 1], seq2[j - 1])))
            diag_score = mat[i - 1, j - 1] + self.scoring_matrix[key]
            up_score = mat[i - 1, j] + -4
            left_score = mat[i, j - 1] + -4
            # finding the direction of pointer from which the score was taken
            if mat[i, j] == diag_score:
                alseq1 += seq1[i - 1]
                alseq2 += seq2[j - 1]
                i -= 1
                j -= 1
            elif mat[i, j] == up_score:
                alseq1 += seq1[i - 1]
                alseq2 += '-'
                i -= 1
            elif mat[i, j] == left_score:
                alseq1 += '-'
                alseq2 += seq2[j - 1]
                j -= 1
            else:
                raise Exception('This should never happen')

        # correct order of string since we traced backwards
        aligned_seq1 = alseq1[::-1]
        aligned_seq2 = alseq2[::-1]

        # for matching base pair at the start of a sequence
        # for a match preceding the alignment at start of any one sequence
        pos1 = seq1.find(aligned_seq1) - 1
        pos2 = seq2.find(aligned_seq2) - 1

        if seq1[pos1] == seq2[pos2]:
            aligned_seq1 = seq1[pos1] + aligned_seq1
            aligned_seq2 = seq2[pos2] + aligned_seq2
            alignment_score += 1

        # Alignment statistics
        matches = ''
        id_count = 0  # sequence identity count
        ss_count = 0  # sequence similarity count
        for aa1, aa2 in zip(aligned_seq1, aligned_seq2):
            pair = ''.join(sorted((aa1, aa2)))
            if aa1 == '-' or aa2 == '-':
                matches += ' '
            elif aa1 == aa2:
                id_count += 1
                matches += '|'
            elif self.scoring_matrix[pair] >= 0:
                ss_count += 1
                matches += ':'
            else:
                matches += ' '

        seq_identity = round(id_count / len(aligned_seq1) * 100, 2)
        seq_similarity = round(ss_count / len(aligned_seq1) * 100, 2)

        # FORMATTED OUTPUT
        if print_output is True:
            print(f'Alignment Score = {alignment_score}')
            print(f'Sequence Identity = {seq_identity}%')
            print(f'Sequence Similarity = {seq_similarity}%')

            for i in range(0, len(aligned_seq1) - 1, 80):
                print(aligned_seq1[i: i + 80])
                print(matches[i: i + 80])
                print(aligned_seq2[i: i + 80])
                print()

        return alignment_score, seq_identity, aligned_seq1

    def perform_msa(self, sequences: list) -> dict:
        """
        Performs multiple sequence alignment of given sequences using pairwise local alignment method.
        Args:
            sequences: (list) list of sequences

        Returns:
            pairwise_scores: (dict) scores (values) of all possible pairs of sequences (keys)

        """
        all_pairs = [(i, j) for i, j in combinations(sequences, 2)]
        pairwise_scores = {}
        # perform pairwise alignment of all sequences
        for seq1, seq2 in all_pairs:
            score, seq_id, alignment = self.pairwise_alignment(seq1=seq1, seq2=seq2)
            pairwise_scores[(seq1, seq2)] = (seq_id, score, alignment)

        return pairwise_scores


class Assembly(FastaTools, MSA):
    """ Performs reconstruction of nucleotide sequence from kmers/ short reads. """

    def __init__(self, sequences: str, output_path: str = None):
        FastaTools.__init__(self)
        MSA.__init__(self)
        self._check_input_path(sequences)
        self.seqs = self._read_seq_file(sequences)

        if output_path is not None:
            self._check_path(output_path)

        self.assembled_sequence = self.de_novo_assembly(sequences=self.seqs, output_path=output_path)

    @staticmethod
    def _check_input_path(path: str) -> None:
        """
        Checks if the input path format is valid.
        Args:
            path: path: (str) path to input file

        Raises:
            ValueError if output file format not in ['txt', 'fasta', 'fastq']

        """
        file_format = path.split('.')[-1]
        if file_format not in ['txt', 'fasta', 'fastq']:
            raise ValueError(f'Format "{file_format}" is not supported for input (supported formats: txt, fasta, fastq)')

    @staticmethod
    def _check_path(path: str) -> None:
        """
        Checks if the output path format is valid.
        Args:
            path: (str) path to output file

        Raises:
            ValueError if output file format not in ['txt']

        """
        file_format = path.split('.')[-1]
        if file_format not in ['txt']:
            raise ValueError(f'Format "{file_format}" is not supported (supported formats: txt)')

    def _read_seq_file(self, file_path: str) -> list:
        """
        Reads file containing list of nucleotide sequences
        Args:
            file_path: (str) path to file

        Returns:
            (list) list of sequences

        """
        file_format = file_path.split('.')[-1]
        if file_format == 'txt':
            with open(file_path, 'r') as file:
                return [line.strip().upper() for line in file.readlines()]
        elif file_format == 'fasta':
            return self.fasta_list(file_path)
        elif file_format == 'fastq':
            return self.fastq_list(file_path)

    @staticmethod
    def _find_overlapping_ends(seq1: str, seq2: str, aligned_seq: str) -> bool:
        """
        Finds if aligned sequence is present at ends of both sequence.
        Args:
            seq1: (str) sequence 1
            seq2: (str) sequence 2
            aligned_seq: (str) fragment of aligned sequence

        Returns:
            (bool) False if aligned sequence present in middle of at least one sequence

        """
        if seq1.startswith(aligned_seq) and seq2.endswith(aligned_seq):
            return True
        elif seq2.startswith(aligned_seq) and seq1.endswith(aligned_seq):
            return True
        else:
            return False

    @staticmethod
    def _merge_overlapping_ends(seq1: str, seq2: str, aligned_seq: str) -> str:
        """
        Merges overlapping ends of two sequences.
        Args:
            seq1: (str) sequence 1
            seq2: (str) sequence 2
            aligned_seq: (str) fragment of sequence overlapping

        Returns:
            (str) merged sequence

        """
        ind1, ind2 = seq1.find(aligned_seq), seq2.find(aligned_seq)
        # assert that aligned sequence is present in both sequences
        assert ind1 != -1
        assert ind2 != -1

        if ind1 > ind2:
            return seq1[:ind1]+seq2
        elif ind2 > ind1:
            return seq2[:ind2]+seq1

    def de_novo_assembly(self, sequences: list, output_path: str = None) -> str:
        """
        Performs De Novo sequence assembly using Greedy algorithm.
        Args:
            sequences: (list) list of kmers/ short reads
            output_path: (str) path to output file to store assembled sequence

        Returns:
            (str) assembled sequence

        """
        seqs = deepcopy(sequences)
        # joining sequence ends to return one long continuous sequence (contig)
        while len(seqs) > 1:
            # perform msa
            pairwise_scores = self.perform_msa(sequences=seqs)
            # sort pairwise alignments based on scores
            sorted_scores = sorted(pairwise_scores.items(), key=lambda x: x[1], reverse=True)

            for scores in sorted_scores:
                seq1, seq2 = scores[0]
                alignment = scores[1][2]

                # find overlapping ends
                if self._find_overlapping_ends(seq1=seq1, seq2=seq2, aligned_seq=alignment):
                    merged_seq = self._merge_overlapping_ends(seq1=seq1, seq2=seq2, aligned_seq=alignment)
                    seqs.remove(seq1)
                    seqs.remove(seq2)
                    seqs.append(merged_seq)
                    break
                else:  # if no overlapping ends, check next highest scoring pairwise alignment
                    pass

        if output_path is not None:
            self._check_path(output_path)
            with open(output_path, 'w') as file:
                file.write(f'{seqs[0]}\n')

        return seqs[0]
