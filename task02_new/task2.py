import re
import click
@click.group()
def main():
    """Entry method."""
    pass
@click.command()
@click.argument('dnaseq')
@click.argument('frame')

def AA_prediction(dnaseq:str,frame:int):
    '''this transcrib the genetic seq in mrna seq then predict possible aminoacid sequence from generated mrna sequence'''
    dna_new = dnaseq.replace("\n", "")
    transcript = dna_new.maketrans('ATGC', 'UACG')
    # print(transcript)
    rna = dna_new.translate(transcript)
    RNA_Codons = {
        # 'M' - START, '*' - STOP
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "UGU": "C", "UGC": "C",
        "GAU": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "UUU": "F", "UUC": "F",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "CAU": "H", "CAC": "H",
        "AUA": "I", "AUU": "I", "AUC": "I",
        "AAA": "K", "AAG": "K",
        "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUG": "M",
        "AAU": "N", "AAC": "N",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UGG": "W",
        "UAU": "Y", "UAC": "Y",
        "UAA": "*", "UAG": "*", "UGA": "*"
    }
    for i in range(frame,len(rna), 3):
        codon1 = rna[i:i + 3]
        if codon1 == 'AUG':
            position1 = i
            for j in range(position1, len(rna), 3):
                codon2 = rna[j:j + 3]
                if codon2 in ['UAA', 'UAG', 'UGA']:
                    position2 = j
                    yield (position2 - position1 + 3, rna[position1:position2 + 3])
                    break


dnase='''AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA
TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC
ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG
CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA
GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC
AGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTG
AAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT
GACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTT'''

main.add_command(AA_prediction)
if __name__ == "__main__":
    main()

#print(list(AA_prediction(dnase, 0)))
