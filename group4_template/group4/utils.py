import os
from pathlib import Path

home_dir = Path.home()
PROJECT_DIR = home_dir.joinpath(".project3", "Group4")
DATA_DIR = PROJECT_DIR.joinpath("data")
UPLOAD_FOLDER = os.path.join(DATA_DIR, 'uploads')
DATA_CACHE = os.path.join(PROJECT_DIR, "data")
os.makedirs(DATA_CACHE, exist_ok=True)


# url endpoint
URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?"

# Requests
PUT_Request = "CMD=put&"
GET_Request = "CMD=get&"

# PUT Parameters
# Query = "QUERY="+toy_seq+"&"
Program = "PROGRAM=blastp&"
Database_PDB = "DATABASE=pdb"

# GET Parameters
RID = "RID ="
GET_query_head = URL + GET_Request + "RID="