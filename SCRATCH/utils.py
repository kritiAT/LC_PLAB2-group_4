# import logging
# logger = logging.getLogger(__name__)

from pathlib import Path
import os
home_dir = str(Path.home())
PROJECT_DIR = os.path.join(home_dir, ".Project03", "Group4")
DATA_CACHE = os.path.join(PROJECT_DIR, "data")
os.makedirs(DATA_CACHE, exist_ok=True)

# LOG_DIR = os.path.join(PROJECT_DIR, "logs")
# os.makedirs(LOG_DIR, exist_ok=True)
# # Logging Configuration
# LOG_FILE_PATH = os.path.join(LOG_DIR, "Group4_log.log")
# logging.basicConfig(filename=LOG_FILE_PATH,
#                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')



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

