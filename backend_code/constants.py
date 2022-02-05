""" Constants for the package. """

# Scoring matrix for aligning identical fragments nucleotide sequence.
scoring_matrix = {'AA': 1, 'AT': -4, 'AG': -4, 'AC': -4, 'TT': 1,
                  'GT': -4, 'CT': -4, 'GG': 1, 'CG': -4, 'CC': 1}

# Codon table for translation
codon_table = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S", "UCC": "S",
               "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y", "UGU": "C", "UGC": "C",
               "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L", "CCU": "P",
               "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H", "CAC": "H", "CAA": "Q",
               "CAG": "Q", "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AUU": "I",
               "AUC": "I", "AUA": "I", "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T",
               "ACG": "T", "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
               "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V", "GUA": "V",
               "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A", "GAU": "D",
               "GAC": "D", "GAA": "E", "GAG": "E", "GGU": "G", "GGC": "G", "GGA": "G",
               "GGG": "G", "UAG": "*", "UAA": "*", "UGA": "*"}
