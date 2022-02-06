import click
from utils import *  # names imported: URL, PUT_Request, GET_Request, Program, Database_PDB, RID, GET_query_head, url_request_head
from Blast import Blast_orfs, check_request_status, Blast_sequence

@click.group()
def cli():
    pass

@click.command()
@click.argument('orf')  # NB This should be provided by Faiza code
@click.argument('filename')
def predict(orflist, filename):
    list_results = Blast_sequence (orf, filename)
    print (list_results)

cli.add_command(predict) 

if __name__ == '__main__':
    cli() 
