import pandas as pd
import numpy as np
import os, sys
import csv

def get_pLDDT(scores, rnge, control):
    '''Function that takes as input the file with all the pLDDT scores of a current dataset (scores), the corresponding file containing the range of disordered regions
    (UniProt    start   end) and the control variable which is either 1 for the positive datasets and 0 for the negative dataset (CheZOD)
       The output of the function is the DB_final_pLDDT_socres.txt file for each positive dataset separately, which contains pLDDT scores of disordered regions from
    the current positive dataset concatenated with the pLDDT scores of the disordered regions in CheZOD dataset'''
    counter = 0
    final_pLDDT = np.zeros((0,4), int)                      #create an empty array for storing pLDDT scores only for DOR and DDR regions
    for uniprot in range(np.shape(rnge)[0]):    #iterate through every line in .txt file with Uniprot ID, start and end DOR/DDR region
        for ID in range(np.shape(scores)[0]):               #iterate through /.pLDDT_scores.txt file in search for the same uniprot name
            if ((rnge[uniprot][0] == scores[ID][0]) and (int(rnge[uniprot][1]) <= int(scores[ID][3])) and (int(rnge[uniprot][2]) >= int(scores[ID][3]))):
                counter +=1
                new_row = np.array([[scores[ID][0], scores[ID][3], control, scores[ID][1]]])                 #get me rows of AA residues that are in the DOR/DDR region
                final_pLDDT = np.append(final_pLDDT, new_row, axis = 0)                      #update my final array
            if ((rnge[uniprot][0] == scores[ID][0]) and (int(rnge[uniprot][2]) < int(scores[ID][3]))):
                counter = 0
                break
            if ((rnge[uniprot][0] != scores[ID][0]) and (counter !=0)):
                counter = 0
                break
    return final_pLDDT


if __name__ == '__main__':
    plddt_dir = sys.argv[1]     # name of the directory with all DB_pLDDT_scores (the pLDDT_scores directory)
    rnge_dir = sys.argv[2]      # name of the directory with all DB_Uniprot_start_end files (the DB_uniprot_start_end directory)
    cwd = os.getcwd()           # get current working directory name
    plddt_dir_path = os.path.join(cwd, plddt_dir)       # create a path to the directory with all DB_pLDDT_scores.txt files
    rnge_dir_path = os.path.join(cwd, rnge_dir)         # create a path to the directory with all DB_no_overlap.txt files
    out_dir = 'final_pLDDT_scores'
    plddt_out_path = os.path.join(cwd, out_dir)                # create a path for the output directory with final pLDDT scores
    if not os.path.exists(plddt_out_path):          # create the output directory if it doesn't already exist, else: pass
        os.makedirs(plddt_out_path)
    # first read the plddt scores of all datasets in separate arrays
    for file in os.listdir(plddt_dir_path):
        if file[:-17] == 'CheZOD':
            chezod_plddt_path = os.path.join(plddt_dir_path, file)
            negative_chezod_scores = np.loadtxt(chezod_plddt_path, skiprows = 1, dtype = 'str', delimiter = '\t')
        if file[:-17] == 'Anchor_flex_link':
            anchor_plddt_path = os.path.join(plddt_dir_path, file)
            negative_anchor_scores = np.loadtxt(anchor_plddt_path, skiprows = 1, dtype = 'str', delimiter = '\t')
        if file[:-17] == 'MFIB':
            MFIB_plddt_path = os.path.join(plddt_dir_path, file)
            positive_MFIB_scores = np.loadtxt(MFIB_plddt_path, skiprows = 1, dtype = 'str', delimiter = '\t')
    # read the ranges of IDRs in all dataset in separate arrays
    for f in os.listdir(rnge_dir_path):
        if f[:-15] == 'CheZOD':
            chezod_rnge_path = os.path.join(rnge_dir_path, f)
            negative_chezod_range = np.loadtxt(chezod_rnge_path, dtype = 'str', delimiter = ' ')
        if f[:-15] == 'Anchor_flex_link':
            anchor_rnge_path = os.path.join(rnge_dir_path, f)
            negative_anchor_range = np.loadtxt(anchor_rnge_path, dtype = 'str')
        if f[:-15] == 'MFIB':
            MFIB_rnge_path = os.path.join(rnge_dir_path, f)
            positive_MFIB_range = np.loadtxt(MFIB_rnge_path, dtype = 'str', delimiter = ' ')
    # extract only the pLDDT regions that are within the ranges of given IDRs in each database
    negative_CheZOD_plddt = get_pLDDT(negative_chezod_scores, negative_chezod_range, control = '0')
    negative_Anchor_plddt = get_pLDDT(negative_anchor_scores, negative_anchor_range, control = '0')
    positive_MFIB_plddt = get_pLDDT(positive_MFIB_scores, positive_MFIB_range, control = '1')

    final_plddt_scores_chezod = np.concatenate((positive_MFIB_plddt, negative_CheZOD_plddt), axis = 0)
    final_plddt_scores_anchor = np.concatenate((positive_MFIB_plddt, negative_Anchor_plddt), axis = 0)
    #write out a file containing Uniprot ID, residue number, expected class (0 for CheZOD and Anchor flex link database and 1 for MFIB database) and AlphaFold pLDDT
    header = ['UniProt_ID', 'Residue_number', 'true_pLDDT', 'AF_pLDDT']
    df = pd.DataFrame(final_plddt_scores_chezod, columns = header)
    df_name = os.path.join(cwd, out_dir, 'MFIB_CheZOD_final_pLDDT_scores.txt' )
    df.to_csv(df_name, header = True, index = None, sep = ' ', mode = 'w')
    df_anchor = pd.DataFrame(final_plddt_scores_anchor, columns = header)
    df_anchor_name = os.path.join(cwd, out_dir, 'MFIB_Anchor_flex_link_final_pLDDT_scores.txt' )
    df_anchor.to_csv(df_anchor_name, header = True, index = None, sep = ' ', mode = "w")
