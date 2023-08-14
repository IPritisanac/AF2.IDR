import pandas as pd
import numpy as np
import os, sys
import csv

def get_pLDDT(scores, rnge, control):
    '''Function that takes as input the file with all the pLDDT scores of a current dataset (scores), the corresponding file containing the range of disordered regions
    (UniProt    start   end) and the control variable whih is either 1 for the positive datasets and 0 for the negative dataset (CheZOD)
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
    # first read the plddt scores only of CheZOD dataset in the negative_dataest_scores array
    for file in os.listdir(plddt_dir_path):
        if file[:-17] == 'CheZOD':
            chezod_plddt_path = os.path.join(plddt_dir_path, file)
            negative_dataset_scores = np.loadtxt(chezod_plddt_path, skiprows = 1, dtype = 'str', delimiter = '\t')
    # read the ranges of IDRs in CheZOD dataset in the negative_dataset_range array
    for f in os.listdir(rnge_dir_path):
        if f[:-15] == 'CheZOD':
            chezod_rnge_path = os.path.join(rnge_dir_path, f)
            negative_dataset_range = np.loadtxt(chezod_rnge_path, dtype = 'str', delimiter = ' ')
    # extract only the pLDDT regions that are within the ranges of given IDRs in the CheZOD database
    negative_plddt = get_pLDDT(negative_dataset_scores, negative_dataset_range, control = '0')

    # now iterate over the directory with pLDDT scores of all remaining datasets (positive) and read the plddt scores of each positive dataset one by one
    # in the positive_dataset_scores_array, then find the corresponding DB_uniprot_start_end file to extract the pLDDT scores of denoted regions
    union_plddt_scores = np.zeros((0,4))
    union_plddt_scores = np.concatenate((union_plddt_scores, negative_plddt), axis = 0)     # create a union_plddt_scores array which will contain all positive datasets and a negative dataset
    for file in os.listdir(plddt_dir_path):
        curr_dataset = file[:-17]
        if curr_dataset != 'CheZOD':
            positive_plddt_path = os.path.join(plddt_dir_path, file)
            positive_dataset_scores = np.loadtxt(positive_plddt_path, skiprows = 1, dtype = 'str', delimiter = '\t')
            for f in os.listdir(rnge_dir_path):
                dataset = f[:-15]
                if curr_dataset == dataset:
                    positive_rnge_path = os.path.join(rnge_dir_path, f)
                    positive_dataset_range = np.loadtxt(positive_rnge_path, dtype = 'str', delimiter = ' ')
                    positive_plddt = get_pLDDT(positive_dataset_scores, positive_dataset_range, control = '1')
                    #concatanate the current positive dataset pLDDT scores with the CheZOD pLDDT scores
                    final_plddt_scores = np.concatenate((positive_plddt, negative_plddt), axis = 0)
                    union_plddt_scores = np.concatenate((union_plddt_scores, positive_plddt), axis = 0)
                    #write out a file containing Uniprot ID, residue number, expected class (0 for CheZOD database and 1 for all positive databases) and AlphaFold pLDDT
                    header = ['UniProt_ID', 'Residue_number', 'true_pLDDT', 'AF_pLDDT']
                    df = pd.DataFrame(final_plddt_scores, columns = header)
                    file_name = plddt_out_path + '\\' + dataset + '_CheZOD_final_pLDDT_scores.txt'
                    df.to_csv(file_name, header = True, index = None, sep = ' ', mode = 'w')
    union_df = pd.DataFrame(union_plddt_scores, columns = header)
    file_name = plddt_out_path + '\\' + 'Union_CheZOD_final_pLDDT_scores.txt'
    union_df.to_csv(file_name, header = True, index = None, sep = ' ', mode = 'w')
