import numpy as np

from itertools import combinations
from copy import deepcopy
from main import dna_sequence

# dna_sequence= 'CTTTATTTTCACCACAAACAGAGATTAAAGAAAGTGTTGGATTAAAAGATGGTGTTAAAGAGTAAAAATTGACTTATTATACTCCTGAATACGAAACCAAAGATACTGATATCTTGGCAGCATTCCGAGTAACTCCTCAACCTGGAGTTCCACCCGAAGAAGCAGGGGCCGCGGTAGCTGCCGAATCTTCTACCGGTACATGGACAAGTGTGTGGACCGATGGACTTACCAGCCTTGATCGTTACAAAGGGCGATGCTACCACATCGAGCCCGTTCCTGGAGAAGAAGAGCAATATATTGCTTATGTAGCTTACCCCTTAGACCTTTTTGAAGAAGGTTCTGTTACTAACATGTTTACTTCCATTGTAGGTAATGTATTTGGGTTCAAAGCCCTACGCGCTCTACGTCTGGAAGATTTGCGAATCCCTCCGGCTTATATTAAAACCTTCCAAGGCCCGCCTCATGGCATCCAGGTTGAGAGAGATAAATTGAACAAATATGGTCGTCCCCTGTTGGGATGTACTATTAAACCTAAATTGGGGTTATCCGCTAAAAACTACGGTAGGGCATGTTATGAATGTCTTCGTGGTGGACTTGATTTTACCAAAGATGATGAAAACGTGAACTCCCAACCGTTTATGCGTTGGAGAGATCGATTCTTATTTGTAGCAAAAA'


def kmers(seq,k):
        kmer_dic={}
        kmers_list = []
        for i in range(len(seq)-(k-1)):
                k_mer=seq[i:i+k]
                kmers_list.append(k_mer)
                if k_mer in kmer_dic:
                        kmer_dic[k_mer]+=1
                else:
                        kmer_dic[k_mer]=1
        return kmers_list

dna_kmers = kmers(dna_sequence, 50)
dna_kmers_random = [dna_kmers[0], dna_kmers[15], dna_kmers[25], dna_kmers[45], dna_kmers[58],
                    dna_kmers[68], dna_kmers[80], dna_kmers[105], dna_kmers[135], dna_kmers[150]]

dna_kmers_random

def read_blosum(file):
    with open(file, 'r') as f:
        matrix = f.readlines()
    del matrix[-1]
    score_dic = {}  # dictionary to store scores of all the pairs
    cols = matrix[0].split()
    for line in range(1, len(matrix)):
        entries = matrix[line].split()
        row = entries[0]  # nucleotide
        scores = entries[1:]  # scores
        for score, col in zip(scores, cols):
            key = ''.join(sorted((row, col)))
            score_dic[key] = int(score)
    return score_dic


def pairwise_alignment(seq1, seq2, w, scoring, print_output=False):
    """Smith-Waterman local alignment"""
    n = len(seq1)
    m = len(seq2)

    # INITIALIZE MATRIX
    mat = np.zeros((n + 1, m + 1), dtype=int)
    for j in range(m + 1):  # first row
        mat[0, j] = j * w
    for i in range(n + 1):  # first column
        mat[i, 0] = i * w

    # FILL MATRIX
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pair = ''.join(sorted((seq1[i - 1], seq2[j - 1])))
            diag = mat[i - 1, j - 1] + scoring[pair]
            up = mat[i - 1, j] + w
            left = mat[i, j - 1] + w
            mat[i, j] = max(diag, up, left, 0)  # replacing all negative values with zero

    alignment_score = np.amax(mat)

    # TRACEBACK ALIGNMENT
    alseq1 = ''
    alseq2 = ''  # storing the gaps in alignment

    max_index = np.where(mat == alignment_score)
    i, j = max_index[0][0], max_index[1][0]  # starting from max value
    while mat[i, j] != 0:
        key = ''.join(sorted((seq1[i - 1], seq2[j - 1])))
        diag_score = mat[i - 1, j - 1] + scoring[key]
        up_score = mat[i - 1, j] + w
        left_score = mat[i, j - 1] + w
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
        elif scoring[pair] >= 0:
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

scoring_matrix = {'AA': 1, 'AT': -4, 'AG': -4, 'AC': -4, 'TT': 1,
                               'GT': -4, 'CT': -4, 'GG': 1, 'CG': -4, 'CC': 1}
pairwise_alignment(seq1=dna_kmers_random[0], seq2=dna_kmers_random[1], w=-4, scoring=scoring_matrix, print_output=False)


def find_overlapping_ends(seq1, seq2, aligned_seq):
    """Find if aligned sequence is present at ends of both sequence.
    Returns false if present in middle or if not present."""
    if seq1.startswith(aligned_seq) and seq2.endswith(aligned_seq):
        return True
    elif seq2.startswith(aligned_seq) and seq1.endswith(aligned_seq):
        return True
    else:
        return False

def merge_overlapping_ends(seq1, seq2, aligned_seq):
    """aligned seq should be present in both and ovelap at ends"""
    ind1, ind2 = seq1.find(aligned_seq), seq2.find(aligned_seq)
    assert ind1 != -1
    assert ind2 != -1
    if ind1 > ind2:
        return seq1[:ind1]+seq2
    elif ind2 > ind1:
        return seq2[:ind2]+seq1
    
def perform_msa(sequences, score_mat):
    all_pairs = [(i,j) for i, j in combinations(sequences, 2)]
    pairwise_scores = {}
    # perform pairwise alignment of all sequences
    for seq1, seq2 in all_pairs:
        score, seq_id, alignment = pairwise_alignment(seq1=seq1, seq2=seq2, w=-4, scoring=score_mat)
        pairwise_scores[(seq1, seq2)] = (seq_id, score, alignment)

    return pairwise_scores


def seq_assembly(sequences, score_mat):
    # all sequences in multiple sequence alignment
    # all_seqs = set()
    # for s1, s2 in pairwise_scores.keys():
    #     all_seqs.update([s1, s2])

    seqs = deepcopy(sequences)
    while len(seqs) > 1:
        # perform msa
        pairwise_scores = perform_msa(sequences=seqs, score_mat=score_mat)
        # sort pairwise alignments based on scores
        sorted_scores = sorted(pairwise_scores.items(), key=lambda x: x[1], reverse=True)
        for scores in sorted_scores:
            seq1, seq2 = scores[0]
            alignment = scores[1][2]

            # find overlapping ends
            if find_overlapping_ends(seq1=seq1, seq2=seq2, aligned_seq=alignment):
                merged_seq = merge_overlapping_ends(seq1=seq1, seq2=seq2, aligned_seq=alignment)
                seqs.remove(seq1)
                seqs.remove(seq2)
                seqs.append(merged_seq)
                break
            else:  # if no overlapping ends, check next highest scoring pairwise alignment
                pass
    return seqs[0]

#seq_assembly(sequences=dna_kmers_random, score_mat=scoring_matrix)

