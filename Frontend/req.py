import os

from pathlib import Path
from flask import Flask, flash, render_template, request, redirect, url_for

home_dir = Path.home()
PROJECT_DIR = home_dir.joinpath(".project3")
DATA_DIR = PROJECT_DIR.joinpath("data")
UPLOAD_FOLDER = os.path.join(DATA_DIR, 'uploads')
#faapath = os.path.join(UPLOAD_FOLDER, 'file.txt')
#textpath= os.path.join(UPLOAD_FOLDER, 'file.txt')


# url endpoint
URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?"

# Requests
PUT_Request = "CMD=put&"
GET_Request = "CMD=get&"

# PUT Parameters
#Query = "QUERY="+toy_seq+"&"   # take from Faiza
Program = "PROGRAM=blastp&"
Database_PDB = "DATABASE=pdb"

# GET Parameters
RID = "RID ="
GET_query_head = URL+ GET_Request + "RID="