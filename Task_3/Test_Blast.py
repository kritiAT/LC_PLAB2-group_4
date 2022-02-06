""" Test for the module Query_Blast"""

from Query_Blast import extract_attribute, check_request_status
from Statics_Query import *  # names imported: URL, PUT_Request, GET_Request, Program, Database, RID, GET_query_head, url_request_head,
                             # LOG_FILE_PATH, DATA_DIR, LOG_DIR    PROJECT_DIR   home_dir

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#self.cache_file_path = os.path.join(DATA_DIR, f"{symbol.lower()}.json")



class TestQueryBlast():
    
    def test_extract_attribute(self):
        
        toy_orf = "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLG" 
        
        # 1.Verify the validity of the aminoacidic sequence
        aa = {'C','D', 'S', 'Q', 'K','I', 'P', 'T', 'F', 'N', 'G', 'H', 'L', 'R', 'W', 'A', 'V', 'E', 'Y', 'M'}
        for letter in toy_orf:
            assert letter in aa, print ("The aminoacidic sequence contain invalid symbols: {}".format(letter))
            
        Query = "QUERY="+toy_orf+"&"
        PUT_query = URL + PUT_Request + Query + Program + Database 
        p = requests.put(PUT_query)
        
        # 2. Verify the status of the request:
        assert p.status_code==200, print ("The status of your request is not succesfull")
        
        rid = extract_attribute(p.text, RID)
        
        # 3. Verify that the RID was properly extracted:
        assert type(rid)==str, print ("The response ID was not extracted properly")

    def TestRequestStatus(self):
        
        # QUESTION: should I define here the name of the output file? 
        output_file= "test_file.html"
        
        
        GET_query = GET_query_head + rid 
        url_submit = url_request_head + rid 
        check_request_status(GET_query, url_submit, output_file)
        
        # 1. Assess that the file was created in the folder
        cache_file_path = os.path.join(DATA_DIR, output_file)
        assert cache_file_path.is_file()




