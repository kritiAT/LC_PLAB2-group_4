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


## Protein Prediction

The aminoacidic sequences are provided as a list of strings which are used to 
build a query which will be submitted to blast server for an alignment against sequences in 
a choice database (default is PDB). 

In the request is possible to specify the database, 
the blast server to use, the email address for contact purposes and other blast parameters.
The Request module is used to submit the query. 

The result response will be output by the function *query_Blast*. 
The Request ID (rid) will be extracted by the function *extract_attribute* and it will used to
check the request status by the function *check_request_status*.

Every time the function will wait 60 sec to assess the status and when it will become ready it will
breake and the function *get_results* will return the response object. The response object will be used to
extract the html (default) or json results. In the case of json results, they will be downloaded as zip file.
These files will be red by the function *html_reader* or *zip_reader* will take as input the string path of the downloaded files
and will extract the result IDs of the predicted proteins, which will be returned as a list.

The wrapper function *Blast_orfs* will loop into a list of compiled aminoacidic sequence and
will return a dictionary where the keys are the sequences and the values are the list of IDs.


## GUI

Using flask to build a web application to process fixed-length DNA strands.

The upload function only allows .fasta .fastq .txt files.

After uploading your file of DNA strands and submitting the file:

1. Depending on the length of your sequence, you need to wait until the assembly is finished.

2. You will also need to wait for BLAST tool to analyze your sequence and give you the predicted proteins list.
 
### Run Flask

`$ python main.py`

In flask, Default port is `5000`

Index page:  `http://127.0.0.1:5000/`

Upload page: `http://127.0.0.1:5000/upload`


### Run with Docker



`$ docker build . -t plab2:latest`

`$ docker run --name plab2_test -p 5000:5000 -d plab2:latest`



### References

Offical Website

- [Flask](http://flask.pocoo.org/)

Tutorial

- [Flask Overview](https://www.slideshare.net/maxcnunes1/flask-python-16299282)
- [In Flask we trust](http://igordavydenko.com/talks/ua-pycon-2012.pdf)

[Wiki Page](https://github.com/tsungtwu/flask-example/wiki)

