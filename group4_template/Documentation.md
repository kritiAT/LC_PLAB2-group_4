# Group4 Package Documentation



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

