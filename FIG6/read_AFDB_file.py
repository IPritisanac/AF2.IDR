import os, sys
import numpy as np
import matplotlib.pyplot as plt


def convert_AA(aa_list):
    """ convert 3-letter amino-acid codes into 1-letter codes"""
    aa_dict = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K',
               'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V',
               'TRP': 'W', 'TYR':'Y'}
    new_aa = ''
    for aa in aa_list:
        new_aa += aa_dict[aa]

    return new_aa


def read_AFDB(proteome):
    """ read pdb files in an AlphaFold DB Proteome directory, return per-residue pLDDT scores

    Input : str, name of directory that contains the PDB files of interest """

    model_cnt = 0
    uniprot_list = []
    pLDDT_list = []
    cwd = os.getcwd()
    os.chdir(cwd+'/'+proteome)
    res_cnt = 0
    long_cnt = 0
    long_list = []
    long_res_cnt = 0
    final_dict = {}
    for r, d, f in os.walk(os.getcwd()):
        for file in f:
            if file.endswith(".pdb"):                           # reading only .pdb files here, not .gz
                pLDDT_dict = {}                                 # empty dictionary to store per-residue pLDDT scores for each PDB file
                aa_seq = ''                                     # empty string to store the amino acid sequence for each PDB file
                res_list = []                                   # list to store residue numbers
                res_start = 0                                   # used to set the first residue number
                model_cnt += 1                                  # counter to iterate over the PDB files
                uniprot_ID = ''                                 # set an empty uniprot ID

                fragment_num = 1 #int(file.split('-')[2][1:])   # this is only needed for long proteins >2700aa; get the fragment number (F1 -> 1, F2 -> 2, ...)
                pdb_file = open(file, "rt")                     # read the pdb file
                first_res_flag = 0                              # identifies the first residue when == 1



                for line in pdb_file:                           # now loop over PDB file... the exceptions below are because sometimes I found that the AlphaFold2 people had different spacing in their files
                    new = line.strip().split()

                    # First get the UniProt ID
                    if len(new) > 1 and new[0] == 'DBREF':
                        uniprot_ID += new[6]                    # Get the UniProt ID

                        # Count all proteins with multiple PDB files (i.e., proteins > 2700 residues)
                        if fragment_num > 1:
                            long_cnt += 1
                            long_list.append(uniprot_ID)
                            long_res_cnt += int(new[4])

                        uniprot_list.append(uniprot_ID)         # Append UniProt ID to list
                        total_res = int(new[4])                 # Total number of residues
                        res_cnt += total_res                    # Update the total residue number count

                    # Then get the AA sequence
                    elif len(new) > 1 and new[0] == 'SEQRES':
                        aa_seq += convert_AA(new[4:])           # Extract the 3-letter amino acid code, convert to 1-letter, and append to string


                    # Then get per-residue pLDDT scores
                    elif len(new) >1 and new[0] == 'ATOM':
                        first_res_flag += 1                     # Counter that gets set to 1 such that the first residue number can be stored
                        if fragment_num == 1:
                            scalar = 0                          # Multiplicative factor that's set to 0 for the first PDB file associated with long proteins > 2700 residues
                        else:
                            scalar = fragment_num-1             # Multiplicative factor to correct residue numbers in long proteins > 2700 residues

                        # these elif len(new) == 11 and elif len(new) == 10 conditions are necessary to store all pLDDT scores
                        if len(new) == 12:
                            res_num = int(new[5]) + (199)*scalar + fragment_num-1    # Correct for long proteins with the 200*scalar
                            res_list.append(res_num)                # Store the residue number
                            pLDDT = float(new[10])                  # Store the pLDDT score
                            if first_res_flag == 1:
                                res_start = (res_num)           # Store the first residue in case first_res_flag is 1
                                if res_start == 0:
                                    res_start = 1
                        elif len(new) == 11:
                            try:
                                res_num = int(new[4].split('A')[1]) + int((199)*scalar) + fragment_num-1
                                res_list.append(res_num)
                                pLDDT = float(new[9]) # some bug in AFDB that fuses the y, z coordinates...
                                if first_res_flag == 1:
                                    res_start = (res_num)
                                    if res_start == 0:
                                        res_start = 1
                            except ValueError:
                                res_num = int(new[5]) + int((199)*scalar) + fragment_num-1
                                res_list.append(res_num)
                                pLDDT = float(new[9])
                                if first_res_flag == 1:
                                    res_start = (res_num)
                                    if res_start == 0:
                                        res_start = 1
                        elif len(new) == 10:
                            res_num = int(new[4].split('A')[1]) + int((199)*scalar) + fragment_num-1
                            res_list.append(res_num)
                            pLDDT = float(new[8]) # some bug in AFDB that fuses the y, z coordinates...
                            if first_res_flag == 1:
                                res_start = (res_num)
                                if res_start == 0:
                                    res_start = 1

                        # Store per-residue pLDDDT scores in a dictionary {residue number : pLDDT score}
                        if res_num in pLDDT_dict.keys():
                            pass
                        else:
                            pLDDT_dict[res_num] = pLDDT


                # Identify PDB files where the number of pLDDT scores does not match the amino-acid sequence length
                if len(aa_seq) != len(list(pLDDT_dict.values())):
                    print(uniprot_ID, len(aa_seq), len(aa_seq) - len(list(pLDDT_dict.values())))

                # Store values to a list
                pLDDT_list.append([uniprot_ID, list(pLDDT_dict.values()), aa_seq, res_start])   # UniProt ID, pLDDT score, amino acid type, amino acid number

                # Now correct for the proteins that were split into multiple fragments
                temp_dict = {}
                temp_dict[res_start] = [uniprot_ID,list(pLDDT_dict.values()), aa_seq]

                final_dict.setdefault(uniprot_ID,{})
                final_dict[uniprot_ID].setdefault(res_start,[uniprot_ID,list(pLDDT_dict.values()), aa_seq])

                pdb_file.close()    # close the PDB file

    duplicated_residues = long_res_cnt + len(list(set(long_list)))*int(1400)        # figure out how many residues would be double-counted if the long proteins >2700aa with multiple fragments are not corrected for
    fragment_count = long_cnt + len(list(set(long_list)))                           # figure out how many proteins have multiple fragments

    return final_dict, fragment_count, long_list, duplicated_residues, pLDDT_list


def fix_long_protein(long_prot_dict, long_protein_uniprot):
    """ Fix the issue associated with proteins > 2700 residues that are split across multiple PDB files"""
    # Now correct for the long proteins
    long_dict_fix = {}

    for uniprot in long_protein_uniprot:
        pLDDT_fix = []
        AA_fix = ''
        pLDDT_final = []
        for i,val in enumerate(sorted(long_prot_dict[uniprot].items(), key=lambda item: int(item[0]))):
            temp_pLDDT = val[1][1]
            temp_AA = val[1][2]
            if i == 0:
                pLDDT_fix.append(temp_pLDDT)
                AA_fix += temp_AA
            elif i > 0 and i < len(sorted(long_prot_dict[uniprot].items(), key=lambda item: int(item[0])))-2:
                pLDDT_fix.append(temp_pLDDT[-200:])
                AA_fix += temp_AA[-200:]
            elif i == len(sorted(long_prot_dict[uniprot].items(), key=lambda item: int(item[0])))-2:
                pLDDT_fix.append(temp_pLDDT[-200:])
                AA_fix += temp_AA[-200:]
                last_one = int(val[0])+len(val[1][2])-1
                print('last one = ',last_one)
            elif i == len(sorted(long_prot_dict[uniprot].items(), key=lambda item: int(item[0])))-1:
                delta = (int(val[0])+len(val[1][2])-1) - last_one
                pLDDT_fix.append(temp_pLDDT[-delta:])
                AA_fix += temp_AA[-delta:]

        final = []
        final.append([item for sublist in pLDDT_fix for item in sublist])
        final.append([AA_fix])
        final.append(np.arange(1, len(AA_fix)+1, 1))

        long_dict_fix[uniprot] = [final]

    return long_dict_fix


def update_pLDDT_dict(list_pLDDT, fixed_dict):
    """ Update a dictionary with {uniprot: [pLDDT, AA_seq]} with the fixed values for long proteins """

    # Convert list of lists (uniprot, pLDDT, AA seq) into dictionary
    final_dict = {}
    for i, j in enumerate(list_pLDDT):
        final_dict[j[0]] = j[1:]

    # Update the dictionary with the corrected values for long proteins
    for key, val in fixed_dict.items():
        for rkey, rval in final_dict.items():
            if key == rkey:
                final_dict[rkey] = val

    return final_dict



if __name__ == '__main__':
    in_prot = sys.argv[1]  # name of the folder where the .pdb files of current dataset are saved - the folder of each dataset needs to be called from the command line separately
    # First, read the AFDB
    curr_dir = os.getcwd()         # get the current directory where the read_AFDB_file.py script is located
    dat = read_AFDB(in_prot)        # now read all PDB files in said directory
    # Now, fix the data and correct for proteins that are split into fragments
    fixed = fix_long_protein(dat[0], list(set(dat[2])))
    final = update_pLDDT_dict(dat[4], fixed)
    # Change directory to the output folder
    if not os.path.exists(curr_dir+'/pLDDT_scores'):
        os.makedirs(curr_dir+'/pLDDT_scores')            # create an output directory
    os.chdir(curr_dir+'/pLDDT_scores')                   # change into the output directory
    output_dir = os.getcwd()                        # get output directory path
    # Create the outupt file
    fout = open((output_dir+'/'+'%s_pLDDT_scores.txt' % (in_prot)), 'w')    # open the output file
    fout.write('# UniProtID\tpLDDT_score\tAA_type\tResidue_number\n')       # write header for output file
    # Loop over the list of UniProt, per-residue pLDDT score, and amino acids to write out the information
    res_cnt = 0
    for i,j in final.items():
        uniprot = i
        print(uniprot)
        try:
            aa_sequence = j[1]
            pLDDT_scores = j[0]
        except IndexError:
            aa_sequence = j[0][1][0]
            pLDDT_scores = j[0][0]
        for k,score in enumerate(pLDDT_scores):
            fout.write('%s\t%s\t%s\t%.0f\n' % (uniprot, score, aa_sequence[k], k+1))
        res_cnt += k+1
    fout.close()
    # Now print out some stats to the screen
    print('\nTotal number of residues in the fragment-corrected AFDB: %.0f' % (res_cnt))
    print('\nTotal number of proteins in the fragment-corrected AFDB: %.0f' % len(dat[0].keys()))
    print('\nTotal number of fragments in the uncorrected AFDB: %.0f' % dat[1])
    print('\nTotal number of duplicated residues in the uncorrected AFDB: %.0f' % (dat[3]))
