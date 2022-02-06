import os
import re
import time
import requests
from group4_template.group4.utils import *  # URL, PUT_Request, GET_Request, Program, Database_PDB, RID, GET_query_head, url_request_head, DATA_CACHE


# ok
# Debug/auxiliary functions:

def _check_query(seq, program=Program, database=Database_PDB, filters="", email=""):
    """ Function for printing the query to submit. Use for debug purpose"""

    Query = "QUERY=" + seq + "&"
    PUT_query = URL + PUT_Request + Query + program + database + filters + email
    print(PUT_query)
    return PUT_query


def _clean_cache(file=None):
    """ Delete all the cached files in the cache directory. Is possible to specify which file to delete selectively"""

    for f in os.listdir(DATA_CACHE):
        if file == None:
            os.remove(os.path.join(DATA_CACHE, f))
        else:
            os.remove(os.path.join(DATA_CACHE, file))
    return


################################################################################################################

# Main functions:

def query_Blast(seq, program=Program, database=Database_PDB, filters="", email="", debug_mode=False):
    """ Submit a PUT query to Blast and return the response object"""

    Query = "QUERY=" + seq + "&"
    PUT_query = URL + PUT_Request + Query + Program + Database_PDB + filters + email
    if debug_mode == True:
        print(PUT_query)
    p = requests.put(PUT_query)  # Submit the request to the BLAST site
    return p


def extract_attribute(p, rid_attr=RID):
    """ will go line by line through a file, searching for mentions of the rid_attribute and extract/return the value
    associated with it.
    """
    put_response = p.text  # get html/text file
    for line in put_response.splitlines():
        if rid_attr in line:
            attribute_value = re.sub(rid_attr, "", line)
            attribute_value = "".join(attribute_value.split())
            return attribute_value


def check_request_status(rid):
    """ Check the status ID and when ready download the results as html file"""

    url_submit = URL + 'CMD=Get&FORMAT_OBJECT=SearchInfo&RID=' + rid
    # Submit the request of the status:
    submit_request = requests.get(url_submit)  # original: put --> Maybe is get??
    query_status = extract_attribute(submit_request, "Status=")
    return query_status


def get_results(rid, response_format="html"):
    """Once the status of the request is ready, it is possible to get the results of the query in html or json.zip format"""

    GET_query = GET_query_head + rid  # /CMD=get&RID=rid"
    if response_format == "json":
        GET_query = GET_query_head + rid + "&FORMAT_TYPE=JSON2"  # /CMD=get&RID=rid&FORMAT_TYPE=JSON2""   
    s = requests.get(GET_query)
    if not s.ok:
        print(f"{GET_query} returned bad status code: {s.status_code}")
        return
    else:
        print(GET_query)
        return s


def download_html_file(s, filename="temporary.html"):
    """Given the result of a request (s) with html content, save it in a file stored in the home directory."""

    assert filename.lower().endswith('.html'), print("save with html file extension")
    path = os.path.join(DATA_CACHE, filename)

    with open(path, "w") as file_handle:
        file_handle.write(s.text)
        print(f"The html results have been saved in the home folder as {filename}")
    return path


def download_zip_json_file(s, filename="temporary.zip"):
    """Given the results of a request (s) with zipped json content, save it in a file stored in the home directory."""

    path = os.path.join(DATA_CACHE, filename)

    with open(path, 'wb') as file_handle:  # write the zip file
        file_handle.write(s.content)
        print(f"The zip/json results have been saved in the home folder as {filename}")
    return path


def html_reader(filepath):
    """Read a html file and returns the list of accessions IDs"""

    assert filepath.lower().endswith('.html'), print("accept only html files")
    with open(filepath, "r") as f:
        list_accessions = []
        for line in f.readlines():
            if "<tr id=" in line:
                i = line.index("<tr id=")
                i2 = line.index("ind=")
                accession = (line[i + 12:i2 - 2])
                list_accessions.append(accession)
    return list_accessions


def zip_reader(filepath):
    """Extract a zip json file and returns the list of accessions IDs"""

    from zipfile import ZipFile
    import json
    assert filepath.lower().endswith('.zip'), print("accept only zip files")

    list_accessions = []
    zf = ZipFile(filepath)  # convert the zip file to a python object
    for filename in zf.namelist():
        if "_" in filename:  # Take only the file with the results
            with zf.open(filename) as f:
                data = f.read()
                d = json.loads(data)
                list_of_hits = d["BlastOutput2"]['report']['results']['search']['hits']
                for e in list_of_hits:
                    ID = e['description'][0]["id"]
                    list_accessions.append(ID)
    return list_accessions


def is_html(output_file):
    """Assess that the output file is in html format."""
    assert output_file[-4:] == "html", print("Save the results in html file format!")


def Blast_sequence(seq, filename="temporary", file_type="html", keep_files=True):
    """seq: sequence to blast, filename: name of the output file"""

    p = query_Blast(seq)  # Submit the request to the BLAST site
    rid = extract_attribute(p, RID)  # Extract the request ID (rid) that will be used for next steps
    while True:
        Status = check_request_status(rid)
        if Status == 'READY':
            print("The status is ready")
            break
        else:
            print("Will wait and check in 60 seconds")
            time.sleep(60)  # Do not poll for any single RID more often than once a minute.
    s = get_results(rid, response_format=file_type)
    print("The sequence was successfully aligned!")
    if file_type == "json":
        filename = filename + ".zip"
        path = download_zip_json_file(s, filename)  # save the zip/json file
        list_matches = zip_reader(path)  # extract the hits IDs
    else:
        filename = filename + ".html"
        path = download_html_file(s, filename)  # save the zip/json file
        list_matches = html_reader(path)  # extract the hits IDs

    if keep_files == False:  # delete the result files
        os.remove(path)
    return list_matches


def Blast_orfs(OrfList, filename="temporary", file_type="html", keep_files=True):
    """Given a list of orfs it blasts them against pdb database and return a dictionary with the orf as key and the list of
    alignments as value"""

    n = 0
    dict_matches = {}
    for orf in OrfList:
        n += 1
        file_name = "orf" + str(n) + "_" + filename
        list_matches = Blast_sequence(orf, file_name, file_type, keep_files)
        dict_matches[orf] = list_matches
        time.sleep(10)  # Do not contact the server more often than once every 10 seconds.
    print("The job is complete for all the sequences!")
    return dict_matches
