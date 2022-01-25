# url endpoint
URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?"     

# Requests
PUT_Request = "CMD=put&"
GET_Request = "CMD=get&"

# PUT Parameters
Query = "QUERY="+toy_seq+"&"
Program = "PROGRAM=blastp&"
Database = "DATABASE=pdb"

# GET Parameters
RID = "RID = "
GET_query_head = URL+ GET_Request + "RID="
url_request_head = URL+'CMD=Get&FORMAT_OBJECT=SearchInfo&RID='

# Put query
PUT_query = URL + PUT_Request + Query + Program + Database

