# Protein Identification via k-mer Genetic Sequence Assembly

### Take a list of DNA strands as inputs, and perform multiple sequence alignment to reconstruct the genetic sequence.
### The compiled genetic sequence is then transcribed to its corresponding mRNA sequence.ORFs are then identified  to predict possible amino acid sequences.
### The generated amino acid sequences are then aligned against a protein database using blast server. 
### Finally, the list of possible proteins is returned.

The aminoacidic sequences are provided as a list of strings which are used to 
build a query which will be submitted to blast server for an alignment against sequences in 
a choice database (default is PDB). 

In the request is possible to specify the database, 
the blast server to use, the email adress for contact purposes and other blast parameters.
The Request module is used to sumbit the query. 

The result response will be output by the function *query_Blast*. 
The Request ID (rid) will be extracted by the function *extract_attribute* and it will used to
check the request status by the function *check_request_status*.

Every time the function will wait 60 sec to assess the status and when it will become ready it will
breake and the function *get_results* will return the response object
