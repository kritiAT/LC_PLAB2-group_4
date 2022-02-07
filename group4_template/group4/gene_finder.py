import pandas as pd
from group4_template.group4.constants import *


class Transcribe:
    """ Class for performing transcription process on a dna sequence """

    def __init__(self, dna: str, reverse: bool = True, threshold: int = 20, output_path: str = None):
        if '.' in dna:
            self._check_input_path(dna)
            self.dna = self._read_dna_file(dna)
        else:
            self.dna = dna
        self.reverse = reverse
        self.threshold = threshold

        self.cdna = self.complementary_dna(self.dna)  # Complementary dna strand
        self.f_mrna = self.transcribe(self.dna)  # Forward strand mrna
        self.r_mrna = self.transcribe(self.cdna)  # Reverse strand mrna

        self.genes = self.gene_sequences()
        self.genes_data = self.__genes_data()

        if output_path is not None:
            self._check_path(output_path)
            self.store_genes(path=output_path)

    @staticmethod
    def _check_input_path(path: str) -> None:
        """
        Checks if the input path format is valid.
        Args:
            path: (str) path to input file

        Raises:
            ValueError if output file format not in ['txt']

        """
        file_format = path.split('.')[-1]
        if file_format not in ['txt']:
            raise ValueError(f'Format "{file_format}" is not supported for input (supported formats: txt)')

    @staticmethod
    def _read_dna_file(file_path: str) -> str:
        """
        Reads text file containing a dna sequence.
        Args:
            file_path: (str) path to file

        Returns:
            (str) dna sequence

        """
        with open(file_path, 'r') as file:
            return file.read().strip()

    @staticmethod
    def _check_path(path: str) -> None:
        """
        Checks if the output path format is valid.
        Args:
            path: (str) path to output file

        Raises:
            ValueError if output file format not in ['txt', 'csv', 'tsv']

        """
        file_format = path.split('.')[-1]
        if file_format not in ['txt', 'csv', 'tsv']:
            raise ValueError(f'Format "{file_format}" is not supported for output (supported formats: txt, csv, tsv)')

    @staticmethod
    def complementary_dna(dna_seq: str) -> str:
        """
        The function returns the complementary strand of the input dna sequence in 5 to 3 direction.
        Args:
            dna_seq: (str) dna sequence

        Returns:
            (str) 5'-3' complementary strand

        """
        return dna_seq.upper().replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g").upper()[::-1]

    @staticmethod
    def transcribe(dna_seq: str) -> str:
        """
        The function returns the mrna strand of the input dna sequence.
        Args:
            dna_seq: (str) dna sequence

        Returns:
            (str) mRNA strand

        """
        return dna_seq.upper().replace("T", "U")

    @staticmethod
    def gene_finder(mrna: str) -> set:
        """
        Finds positions of the longest protein coding sequences in each open reading frame (ORF).
        Args:
            mrna: (str) mRNA strand

        Returns:
            all_positions: (set[tuple]) tuples of coding sequence positions on given sequence for each ORF.

        """
        all_positions = set()
        for n in range(3):
            frame = mrna[n:]
            store_end = 0  # for storing the positions of stop codons in one frame
            for i in range(0, len(frame), 3):
                codon1 = frame[i: i + 3]
                if codon1 == 'AUG':
                    start = i
                    dna_start = start + n  # position on the dna strand regardless of frame
                    for j in range(start, len(frame), 3):
                        codon2 = frame[j: j + 3]
                        if codon2 in ['UAA', 'UAG', 'UGA']:
                            end = j + 3
                            dna_end = end + n
                            if store_end != end:
                                all_positions.add((dna_start, dna_end))  # longest sequence from one stop codon
                                store_end = end
                            break
        return all_positions

    @staticmethod
    def _filter_sequence(positions: set, threshold: int) -> list:
        """
        Filters out short sequences with length less than the threshold.
        Args:
            positions: (set[tuple]) tuples of coding sequence positions.
            threshold: (int) minimum length of amino-acid sequence (excluding start and stop codon)

        Returns:
            (list[tuple])  tuples of coding sequence positions with length at least equal to the threshold.

        """
        if not threshold >= 1:
            raise AssertionError('Threshold should be greater than or equal to 1')
        return list(filter(lambda x: abs(x[1] - x[0]) >= int(threshold + 2) * 3, positions))

    def gene_sequences(self, print_output: bool = False) -> dict:
        """
        Performs exon splicing. Extracts coding sequences from mRNA strand
        The function returns list of all the ORF sequences.
        Args:
            print_output: (bool) option to print formatted output. Default False

        Returns:
            genes: (dict) spliced exon sequences with positions on corresponding strand

        """
        # find position of genes
        gene_pos_dna = self.gene_finder(self.f_mrna)
        gene_pos_cdna = self.gene_finder(self.r_mrna)

        # filter out short sequences
        gene_pos_dna = self._filter_sequence(gene_pos_dna, self.threshold)
        gene_pos_cdna = self._filter_sequence(gene_pos_cdna, self.threshold)

        # splice the genes
        genes_dna = {(start, end): self.f_mrna[start: end] for start, end in gene_pos_dna}
        genes_cdna = {(start, end): self.r_mrna[start: end] for start, end in gene_pos_cdna}

        genes = {'forward': None, 'reverse': None}

        if self.reverse is True:
            genes['forward'] = genes_dna
            genes['reverse'] = genes_cdna

        else:
            genes['forward'] = genes_dna

        if print_output is True:
            print('Gene sequences in forward strand: ')
            for pos, seq in genes['forward'].items():
                print(f'{pos} : {seq} ')
            if self.reverse is True:
                print('Gene sequences in reverse strand: ')
                for pos, seq in genes['reverse'].items():
                    print(f'{pos} : {seq} ')

        return genes

    def __genes_data(self):
        """ Returns pandas dataframe with spliced exons information """
        genes_data = pd.DataFrame(columns=['strand', 'position', 'mrna'])

        strand = ['forward' for _ in range(len(self.genes['forward']))]
        position = [f'{i}-{j}' for i, j in self.genes['forward'].keys()]
        seqs = [i for i in self.genes['forward'].values()]

        if self.reverse is True:
            strand.extend(['reverse' for _ in range(len(self.genes['reverse']))])
            position.extend([f'{i}-{j}' for i, j in self.genes['reverse'].keys()])
            seqs.extend([i for i in self.genes['reverse'].values()])

            genes_data['strand'] = strand
            genes_data['position'] = position
            genes_data['mrna'] = seqs

        else:
            genes_data['strand'] = strand
            genes_data['position'] = position
            genes_data['mrna'] = seqs

        return genes_data

    def store_genes(self, path: str) -> None:
        """
        Writes spliced exon information in output file
        Args:
            path: (str) path to output path

        Raises:
            ValueError if output file format not in ['txt', 'csv', 'tsv']

        """
        file_format = path.split('.')[-1]

        if file_format == 'csv':
            self.genes_data.to_csv(path, sep=',', header=False, index=False)
        elif file_format == 'tsv':
            self.genes_data.to_csv(path, sep='\t', header=False, index=False)
        elif file_format == 'txt':
            self.genes_data.to_csv(path, sep=' ', header=False, index=False)
        else:
            raise ValueError(f'Format "{file_format}" is not supported (supported formats: csv, tsv, txt)')

        return None


class Translate(Transcribe):
    """ Class to translate mRNA sequence to amino-acid sequence."""

    def __init__(self, dna: str = None, reverse: bool = True, threshold: int = 20, output_path: str = None):
        self.codon_table = codon_table

        super().__init__(dna=dna, reverse=reverse, threshold=threshold)
        self.proteins = None
        self.proteins_table = None
        self._translation()  # perform translation of all mrna sequences
        if output_path is not None:
            self._check_path(output_path)
            self.store_proteins(path=output_path)

    def translate(self, mrna: str) -> str:
        """
        The function returns amino-acid sequence of mRNA strand
        Args:
            mrna: (str) mRNA sequence

        Returns:
            protein: (str) amino-acid sequence

        """
        protein = ''
        for i in range(0, len(mrna), 3):
            codon = mrna[i:i + 3]
            if codon == 'AUG':  # exclude start codon from amino-acid sequence
                pass
            elif self.codon_table[codon] == '*':  # stop when stop codon is read
                break
            else:
                protein += self.codon_table[codon]
        return protein

    def _translation(self):
        """ Wrapper function to perform translation of orf genes. """
        self.genes_data['amino_acid_sequence'] = self.genes_data.apply(lambda x: self.translate(x.mrna), axis=1)
        self.proteins = list(self.genes_data['amino_acid_sequence'])
        self.proteins_table = self.genes_data.drop(columns='mrna')

    def store_proteins(self, path: str) -> None:
        """
        Writes proteins information in output file
        Args:
            path: (str) path to output path

        Raises:
            ValueError if output file format not in ['txt', 'csv', 'tsv']

        """
        file_format = path.split('.')[-1]

        if file_format == 'csv':
            self.proteins_table.to_csv(path, sep=',', header=False, index=False)
        elif file_format == 'tsv':
            self.proteins_table.to_csv(path, sep='\t', header=False, index=False)
        elif file_format == 'txt':
            self.proteins_table.to_csv(path, sep=' ', header=False, index=False)
        else:
            raise ValueError(f'Format "{file_format}" is not supported (supported formats: csv, tsv, txt)')

        return None
