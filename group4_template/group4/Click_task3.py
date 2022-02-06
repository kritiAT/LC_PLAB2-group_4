import click
from utils import *  # names imported: URL, PUT_Request, GET_Request, Program, Database_PDB, RID, GET_query_head, url_request_head
from Blast import Blast_orfs, check_request_status, Blast_sequence
# ok
@click.group()
def cli():
    pass

@click.command()
@click.argument('orflist')  # NB This should be provided by Faiza code
@click.option('-k', '--keep_files', default=True, help="Decide if to store the file results or erase")
@click.option('-f', '--filename', default="temporary", help="Give a file name where the results will be saved.")
@click.option('-t', '--file_type', default="html", help="Decide in which format to save the results.Allowed format are html or zip")
def predict(orflist, filename, file_type, keep_files):
    list_results = Blast_sequence (orflist, filename, file_type, keep_files)
    print (list_results)

cli.add_command(predict) 

if __name__ == '__main__':
    cli() 
