""" Script for a working command line interface (CLI) """

from group4.sequence_assembly import MSA, Assembly
from group4.gene_finder import Transcribe, Translate
from group4.Blast import Blast_orfs

import click


@click.group()
def protein_prediction():
    pass


@protein_prediction.command('pairwise')
@click.argument('sequence1')
@click.argument('sequence2')
@click.option('-p', '--show', is_flag=True, default=False, help='Option to print output')
def pw_alignment(sequence1, sequence2, show):
    """
    Performs pairwise local alignment of given sequences.
    Args:
        sequence1: (str) sequence 1
        sequence2: (str) sequence 2
        show: (bool) option to print output to standout

    """
    obj = MSA()
    obj.pairwise_alignment(seq1=sequence1, seq2=sequence2, print_output=show)
    return None


@protein_prediction.command('assemble')
@click.argument('sequences', type=click.Path(exists=True))
@click.option('-o', '--output', default=None, help='Option to store output')
@click.option('-p', '--show', is_flag=True, default=False, help='Option to print output')
def assemble_seq(sequences, output, show):
    """
    Performs De Novo sequence assembly
    Args:
        sequences: (str) input file with sequences
        output: (str) output to output file
        show: (bool) option to print output to standout

    """
    obj = Assembly(sequences=sequences, output_path=output)
    seq = obj.assembled_sequence
    if show is True:
        click.echo(seq)
    return None


@protein_prediction.command('transcribe')
@click.argument('dna_sequence', type=click.Path(exists=True))
@click.option('-r', '--reverse', is_flag=True, default=True, help='Use reverse strand also')
@click.option('-t', '--threshold', type=int, default=20, help='Threshold for amino-acid sequence length')
@click.option('-o', '--output', default=None, help='Option to store output')
@click.option('-p', '--show', is_flag=True, default=False, help='Option to print output')
def transcription(dna_sequence, reverse, threshold, output, show):
    """
    Performs process of transcription on a given dna sequence.
    Args:
        dna_sequence: (str) output to file containing dna sequence
        reverse: (bool) option to use complementary strand of dna also
        threshold: (int) minimum length of amino-acid sequence (excluding start and stop codon)
        output: (str) output to output file
        show: (bool) option to print output to standout

    """
    obj = Transcribe(dna=dna_sequence, reverse=reverse, threshold=threshold, output_path=output)
    genes_table = obj.genes_data

    if show is True:
        click.echo(genes_table)


@protein_prediction.command('translate')
@click.argument('dna_sequence', type=click.Path(exists=True))
@click.option('-r', '--reverse', is_flag=True, default=True, help='Use reverse strand also')
@click.option('-t', '--threshold', type=int, default=20, help='Threshold for aminoacid sequence length')
@click.option('-o', '--output', default=None, help='Option to store output')
@click.option('-p', '--show', is_flag=True, default=False, help='Option to print output')
def translation(dna_sequence, reverse, threshold, output, show):
    """
    Performs process of translation on a given dna sequence.
    Args:
        dna_sequence: (str) output to file containing dna sequence
        reverse: (bool) option to use complementary strand of dna also
        threshold: (int) minimum length of amino-acid sequence (excluding start and stop codon)
        output: (str) output to output file
        show: (bool) option to print output to standout

    """
    obj = Translate(dna=dna_sequence, reverse=reverse, threshold=threshold, output_path=output)
    protein_table = obj.proteins_table

    if show is True:
        click.echo(protein_table)


@protein_prediction.command('predict')
@click.argument('orflist')
@click.option('-k', '--keep_files', default=True, help="Decide if to store the file results or erase")
@click.option('-f', '--filename', default="temporary", help="Give a file name where the results will be saved.")
@click.option('-t', '--file_type', default="html",
              help="Decide in which format to save the results.Allowed format are html or zip")
def predict(orflist, filename, file_type, keep_files):
    """
    Given a list of orfs it blasts them against pdb database and return a dictionary with the orf as key and the list of
    alignments as value.
    Args:
        orflist: (list) list of sequences to be aligned
        filename: (str) file where the results will be stored if keep_file is in default mode
        file_type: (str) type of file to save
        keep_files: (bool) If False will erase the saved file results, else it will keep them

    """
    list_results = Blast_orfs(orflist, filename, file_type, keep_files)
    print(list_results)


@protein_prediction.command('predict_pro')
@click.argument('sequences', type=click.Path(exists=True))
@click.option('-o', '--output', default=None, help='Option to store output')
@click.option('-p', '--show', is_flag=True, default=False, help='Option to print formatted output')
@click.option('-r', '--reverse', is_flag=True, default=True, help='Use reverse strand also')
@click.option('-n', '--threshold', type=int, default=20, help='Threshold for amino-acid sequence length')
@click.option('-k', '--keep_files', default=True, help="Decide if to store the file results or erase")
@click.option('-f', '--filename', default="temporary", help="Give a file name where the results will be saved.")
@click.option('-t', '--file_type', default="html", help="Decide in which format to save the results.Allowed format are html or zip")
def predict_pro(sequences, reverse, threshold, show, output, filename, file_type, keep_files):
    """
    Performs Likelihood Protein Identification via k-mer Genetic Sequence Assembly.
    Args:
        sequences: (str) input file with sequences
        reverse: (bool) option to use complementary strand of dna also
        threshold: (int) minimum length of amino-acid sequence (excluding start and stop codon)
        show: (bool) option to print formatted output to standout
        output: (str) output to output file
        filename: (str) file where the results will be stored if keep_file is in default mode
        file_type: (str) format type of file
        keep_files: (bool) If False will erase the saved file results, else it will keep them

    """
    assembled_dna = Assembly(sequences=sequences).assembled_sequence
    obj = Translate(dna=assembled_dna, reverse=reverse, threshold=threshold)
    protein_list = obj.proteins
    protein_table = obj.proteins_table
    protein_dict = Blast_orfs(protein_list, filename, file_type, keep_files)
    protein_table['predicted_proteins'] = protein_table['amino_acid_sequence'].apply(lambda x: protein_dict[x])
    if show is True:
        click.echo(f'Assembled sequence : \n {assembled_dna}\n')
        click.echo(f'mRNA sequence (forward) : \n {obj.f_mrna}\n')
        click.echo(f'mRNA sequence (reverse) : \n {obj.r_mrna}\n')
        click.echo(protein_table)

    if output is not None:
        file_format = output.split('.')[-1]
        if file_format == 'csv':
            protein_table.to_csv(output, sep=',', header=False, index=False)
        elif file_format == 'tsv':
            protein_table.to_csv(output, sep='\t', header=False, index=False)
        elif file_format == 'txt':
            protein_table.to_csv(output, sep=' ', header=False, index=False)
        else:
            raise ValueError(f'Format "{file_format}" is not supported (supported formats: csv, tsv, txt)')


if __name__ == '__main__':
    protein_prediction()
