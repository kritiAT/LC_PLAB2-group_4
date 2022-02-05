import re
import requests
import time
from req import * # URL, PUT_Request, GET_Request, Program, Database_PDB, RID, GET_query_head, url_request_head
import pandas as pd

ORF1 = "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLG"
OrfList =[ORF1]
#ORF2 = "QGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGD"
#OrfList =[ORF1,ORF2]

def query_Blast (seq, program=Program, database = Database_PDB, filters="", email=""):
    
    """ Submit a PUT query to Blast and return the response object"""
    
    Query = "QUERY="+ seq + "&"
    PUT_query = URL + PUT_Request + Query + Program + Database_PDB + filters + email 
    print (PUT_query) 
    # Submit the request to the BLAST site
    p = requests.put(PUT_query) 
    return p
    

def extract_attribute(p, rid_attr):
    
    """ will go line by line through a file, searching for mentions of the rid_attribute and extract/return the value 
    associated with it.
    """
    put_response = p.text
    for line in put_response.splitlines():
        if rid_attr in line:
            attribute_value = re.sub(rid_attr, "", line)
            attribute_value="".join(attribute_value.split())
            return attribute_value

def write_file(filename, s):
    
    """Given the results of a request (s) it write them in a file with the desired filename"""
    with open(filename, "w") as file_handle:
        file_handle.write(s.text)
        print ("The results have been saved in the current folder")
    
        
def check_request_status(rid):
    
    """ Check the status ID and when ready download the results as html file"""
    
    url_submit = URL +'CMD=Get&ORMAT_OBJECT=SearchInfo&RID=' + rid  
    # Submit the request of the status:
    submit_request = requests.get(url_submit) # original: put --> Maybe is get??
    query_status = extract_attribute(submit_request, "Status=")
    return query_status


def get_results (rid):
    
    """Once the status of the request is ready, it is possible to get the results of the query"""

    GET_query = GET_query_head + rid    # /CMD=Get&RID=6NTRF1YLO8
    s = requests.get(GET_query)
    if not s.ok:
        #logger.error(f"{GET_query} returned bad status code: {s.status_code}")
        print (f"{GET_query} returned bad status code: {s.status_code}") 
        return None
    else:
        print (GET_query)
        return s
    
def is_html(output_file):
    
    """Assess that the output file is in html format."""
    
    assert output_file[-4:]== "html", print ("Save the results in html file format!")
    
def html_reader(file):
    
    """Read an html file and returns the list of accessions"""
    
    with open (file,"r") as f:
        list_accessions =[]
        for line in f.readlines():
            if "<tr id=" in line:
                i = line.index("<tr id=")
                i2 = line.index ("ind=")
                accession = (line [i+12:i2-2])
                list_accessions.append(accession)
    return list_accessions
    

def Blast_sequence (seq, filename):
    
    """seq= sequence to blast, filename= name of the output file, n = numerical identifier"""
    
    p = query_Blast (seq)                       # Submit the request to the BLAST site
    rid = extract_attribute (p, RID)       # Extract the request id that will be used for next steps
    while True:
        Status = check_request_status(rid)
        if Status == 'READY': 
            print("The staus is ready")
            break
        else:
            print("Will wait and check in 60 seconds")  # debug msg
            time.sleep(60)                              # Do not poll for any single RID more often than once a minute. 
    s = get_results (rid)
    print ("The sequence was succesfully aligned!")
    is_html(filename)
    write_file(filename, s)
    list_matches = html_reader(filename)
    return list_matches
    
def Blast_orfs (OrfList, filename):
    
    """Given a list of orfs it blast them against pdb database and return a dictionary with the orf as key and the list of
    alignments as value"""
    
    n=0
    dict_matches={}
    for orf in OrfList:
        n+=1
        file_name = "orf"+ str(n) + "_" + filename
        list_matches = Blast_sequence (orf, file_name)
        dict_matches[orf] = list_matches
        time.sleep(10)                 # Do not contact the server more often than once every 10 seconds.
    print ("The job is complete for all the sequences!")
    return dict_matches





#df_t
