import re
import requests
import time

def extract_attribute(put_response, rid_attr):
    
    """ will go line by line through a file, searching for mentions of the rid_attribute and extract/return the value 
    associated with it """
    
    for line in put_response.splitlines():
        if rid_attr in line:
            attribute_value = re.sub(rid_attr, "", line)
            attribute_value="".join(attribute_value.split())
            return attribute_value
        
def check_request_status(GET_query, url_submit):
    
    """ Check the status ID and when ready download the results as html file"""
    
    # Submit the request of the status:
    submit_request = requests.put(url_submit)
    query_status = extract_attribute(submit_request.text, "Status=")
    
    if query_status=='WAITING':
        time.sleep(30)                               # need to wait as required by the policy of Blast API 
        check_request_status(GET_query, url_submit)  # submit another request of status after 30 sec
        
    else:  # The status is READY       

        # Save the Submit result
        
        with open("results_trial.html", "w") as file_handle:
            s = requests.get(GET_query)
            file_handle.write(s.text)
        print ("The results have been saved in the current folder")
        
        return 