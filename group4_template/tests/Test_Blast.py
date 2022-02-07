""" Test for the module Query_Blast"""

from group4_template.group4.Blast import *
from group4_template.group4.utils import *

class TestBlast():
    
    def test_extract_attribute(self):
        
        toy_orf = "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLG" 
        
        # 1.Verify the validity of the aminoacidic sequence
        aa = {'C','D', 'S', 'Q', 'K','I', 'P', 'T', 'F', 'N', 'G', 'H', 'L', 'R', 'W', 'A', 'V', 'E', 'Y', 'M'}
        for letter in toy_orf:
            assert letter in aa, print("The aminoacidic sequence contain invalid symbols: {}".format(letter))
            
        query = "QUERY="+toy_orf+"&"
        PUT_query = URL + PUT_Request + query + Program + Database_PDB
        p = requests.put(PUT_query)

        # 2. Verify the status of the request:
        assert p.status_code==200, print ("The status of your request is not succesfull")

        rid = extract_attribute(p, rid_attr = RID)

        # 3. Verify that the RID was properly extracted as a string:
        assert type(rid)==str, print ("The response ID was not extracted properly, it is not a string object")

    def TestCheckRequestStatus(self):

        toy_orf = "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLG"
        query = "QUERY="+toy_orf+"&"
        PUT_query = URL + PUT_Request + query + Program + Database_PDB
        output_file= "test_file.html"
        p = requests.put(PUT_query)
        rid = extract_attribute(p.text, rid_attr = RID)
        url_submit = URL +'CMD=Get&ORMAT_OBJECT=SearchInfo&RID=' + rid
        query_status = extract_attribute(submit_request, "Status=")
        assert query_status=="READY" or query_status=="WAITING", print ("The query status is unknown")

    def TestGetResults(self):
        # 1. Assess that the file was created in the folder
        toy_orf = "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLG"
        filename ="test_file.html"
        Blast_sequence (toy_orf, filename = filename, file_type="html", keep_files = True)
        cache_file_path = os.path.join(DATA_DIR, "test_file.html")
        assert cache_file_path.is_file()




