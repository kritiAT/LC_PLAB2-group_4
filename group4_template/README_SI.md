## Task 4 - GUI

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

