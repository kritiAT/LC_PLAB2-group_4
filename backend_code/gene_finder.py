# for protein prediction
import pandas as pd


class Transcribe:

    def __init__(self, dna: str, reverse: bool = True, threshold: int = 20, output_path: str = None):
        self._check_input_path(dna)
        self.dna = self._read_dna_file(dna)
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
        file_format = path.split('.')[-1]
        if file_format not in ['txt']:
            raise ValueError(f'Format "{file_format}" is not supported for input (supported formats: txt)')

    @staticmethod
    def _read_dna_file(file_path: str):
        with open(file_path, 'r') as file:
            return file.read().strip()

    @staticmethod
    def _check_path(path: str) -> None:
        file_format = path.split('.')[-1]
        if file_format not in ['txt', 'csv', 'tsv']:
            raise ValueError(f'Format "{file_format}" is not supported for output (supported formats: txt, csv, tsv)')

    @staticmethod
    def complementary_dna(dna_seq):
        """The function returns the complementary strand of the input dna sequence in 5 to 3 direction."""
        return dna_seq.upper().replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g").upper()[::-1]

    @staticmethod
    def transcribe(dna_seq):
        """The function returns the mrna strand of the input dna sequence."""
        return dna_seq.upper().replace("T", "U")

    @staticmethod
    def gene_finder(mrna):
        """The function finds position of longest protein coding sequence in each open reading frame.
        Returns set of tuples of coding sequence positions on given sequence for each ORF."""
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
    def _filter_sequence(positions, threshold):
        """The function returns a dictionary of positions of genes with minimum length of amino acid sequence
        atleast equal to the threshold."""
        return list(filter(lambda x: abs(x[1] - x[0]) >= int(threshold + 2) * 3, positions))

    def gene_sequences(self, print_output: bool = False):
        """The function returns list of all the ORF sequences."""
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
        genes_data = pd.DataFrame(columns=['strand', 'position', 'sequence'])

        strand = ['forward' for _ in range(len(self.genes['forward']))]
        position = [f'{i}-{j}' for i, j in self.genes['forward'].keys()]
        seqs = [i for i in self.genes['forward'].values()]

        if self.reverse is True:
            strand.extend(['reverse' for _ in range(len(self.genes['reverse']))])
            position.extend([f'{i}-{j}' for i, j in self.genes['reverse'].keys()])
            seqs.extend([i for i in self.genes['reverse'].values()])

            genes_data['strand'] = strand
            genes_data['position'] = position
            genes_data['sequence'] = seqs

        else:
            genes_data['strand'] = strand
            genes_data['position'] = position
            genes_data['sequence'] = seqs

        return genes_data

    def store_genes(self, path: str):
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


# Add function to translate the mrna sequence to amino acid sequence
# using the codon table
# Test mRNA seq =
# 'AUGGUCGUCCCCUGUUGGGAUGUACUAUUAAACCUAAAUUGGGGUUAUCCGCUAAAAACUACGGUAGGGCAUGUUAUGAAUGUCUUCGUGGUGGACUUGAUUUUACCAAAGAUGAUGAAAACGUGA'
# A wrapper function could also be added to take list of mRNA sequences
# and return list of corresponding amino-acid sequence
class Translate:
    """ Class to translate mRNA sequence to amino-acid sequence."""

    # initialise codon table
    # translation function
    pass

