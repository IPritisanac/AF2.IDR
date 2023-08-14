import os, sys
import numpy as np
import urllib.request

# Now read a list of UniProt IDs and automatically download & save the PDB files locally
def download_AFDB_list(file_in, database_name):
    ''' Function to automaticallys download PDB files from the AlphaFold Database
        Input:  an input file that contains a list of UniProt IDs in the 1st column '''

    cnt_success = 0             # initiate a counter to count how many PDB files were successfully downloaded
    cnt_total = 0               # initiate a counter to count the number of UniProt IDs
    fin = open(file_in, 'r')    # open the input file that contains UniProt IDs
    counter = 0
    for line in fin:
        counter += 1            # iterate the counter
        print(counter)          # print to see the number that the loop is on
        new = line.strip().split() # split the line in the file
        uniprot = new[0]        # extract the UniProt ID from the first column
        cnt_total +=1           # iterate the counter (total)
        try:
            url = 'https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v3.pdb' % uniprot       # create the URL based on the AFDB syntax + UniProt I
            urllib.request.urlretrieve(url, pdb_path+'/AF-%s-F1.pdb' % uniprot)             # use urllib to save the file to disk
            cnt_success+=1                                                                  # if successful, iterate the counter
        except urllib.error.HTTPError:                                                      # if unsuccesful, make a note of the UniProt ID
            print('!! error for %s' % uniprot)

    print('\n##### FINISHED %s! downloaded %i of %i UniProt entries' % (database_name, cnt_success, cnt_total))   # print stats


if __name__=='__main__':
    db_start_end_dir = sys.argv[1]        # path to the folder which contains Uniprot_start_end files for each database
    # Create a directory to store the output PDB files (if needed)
    cwd = os.getcwd()               # get path for current working directory
    for file in os.listdir(db_start_end_dir):
        curr_file_path = os.path.join(db_start_end_dir, file)
        out_dir = file[:-15]       # name of output directory
        pdb_path = cwd+'/'+out_dir      # full path for output directory
        if not os.path.exists(pdb_path):    # if the directory doesn't exist
            os.makedirs(pdb_path)           # create directory (else: pass)
        download_AFDB_list(curr_file_path, out_dir)
