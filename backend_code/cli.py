from sequence_assembly import MSA, Assembly
from gene_finder import Transcribe
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
        sequences: (str) path to file with sequences
        output: (str) path to output file
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
        dna_sequence: (str) path to file containing dna sequence
        reverse: (bool) option to use complementary strand of dna also
        threshold: (int) minimum length of amino-acid sequence (excluding start and stop codon)
        output: (str) path to output file
        show: (bool) option to print output to standout

    """
    obj = Transcribe(dna=dna_sequence, reverse=reverse, threshold=threshold, output_path=output)
    genes_table = obj.genes_data

    if show is True:
        click.echo(genes_table)


if __name__ == '__main__':
    protein_prediction()
