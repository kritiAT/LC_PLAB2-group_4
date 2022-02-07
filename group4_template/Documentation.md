# Group4 Package Documentation

The package performs likelihood protein identification via k-mer genetic sequence assembly.

The ``group4`` package a pipeline for assembling fixed-length DNA strand reads into a complete sequence and predict which protein the completed strand may represent.
The package is used to:
1. Perform de novo sequence assembly on a set of DNA strands (k-mers).
2. Transcribe and translate the assembled DNA sequence.
3. Predict protein the assembled DNA strand may represent.
4. For having an interactive experience with a web application which allows one to upload a series of DNA strand reads and obtain a list protein predictions.

## De Novo Sequence Assembly

De novo sequence assembly is a technique for assembling shorter fragments of DNA into a longer one through the use of multiple alignments. 
The ``sequence_assembly`` module performs de novo sequence assembly by assembling short fragments of nucleotide sequence to result one continuous sequence.
Greedy algorithm is used for de novo sequence assembly.
A modified version of Smith-Waterman algorithm was used for multiple sequence alignment. Pairwise local alignment for two sequences can also be done using the module.

### Usage Guidelines

#### 1. Input requirements
The input should be a file containing a list of short fragments/ kmers which spans the whole genome or a continuous part of genome. 
The fragments should have overlapping ends. There **should not** be any gaps or missing fragments in the short reads.
- Accepted file formats : txt(each read on new line), fasta, fastq

#### 2. Module usage 
The ``sequence_assembly`` module can be imported from the package to use pairwise alignment and sequence assembly methods.
The module has class ``FastaTools`` for handling fasta files.
Some examples:
- _Perform pairwise sequence alignment_
  - alignment_score, seq_identity, aligned_seq = ``MSA().pairwise_alignment(seq1=sequence1, seq2=sequence2, print_output=True)``
  - CLI : ``group4 pairwise seq1 seq2 -p``
- _Store assembled sequence in a file_
  - Object = ``Assembly(sequences=input_file, output_path=output_file)``
  - CLI : ``group4 assemble input_file -o output_file -p``

## Central Dogma
The genetic sequence undergoes the process of transcription and translation to form proteins.
The ``gene_finder`` module performs transcription on the assembled sequence to result mRNA strands which are translated to a set of possible amino-acids sequences.
Options to consider complementary (reverse) strand of dna and to set threshold for minimum length of resulting amino-acid sequences are available.

### Usage Guidelines

#### 1. Input requirements
The input should be a nucleotide sequence without any gaps.
A text file with the sequence or the sequence (of type string) can be used.

#### 2. Module usage 
The ``gene_finder`` module can be imported from the package to perform transcription and translation of a dna sequence.
Some examples:
- _Perform transcription of dna sequence_
  - Object = ``Transcribe(dna=dna_sequence, reverse=boolean, threshold=threshold_value, output_path=output_file)``
  - ``Object.genes_data`` returns a dataframe of mRNA sequences with the position on corresponding strand.
  - CLI : ``group4 transcribe dna_sequence -r -t 50 -o output_file -p``
- _Perform translation of dna sequence_
  - Object = ``Translate(dna=dna_sequence, reverse=boolean, threshold=threshold_value, output_path=output_file)``
  - ``Object.proteins_table`` results in a dataframe with all possible amino-acid sequences.
  - CLI : ``group4 translate dna_sequence -r -t 50 -o output_file -p``

