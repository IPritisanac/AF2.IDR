import os, sys
import numpy as np
import matplotlib.pyplot as plt 
import gzip


def convert_AA(aa_list):
    aa_dict = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 
               'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V',
               'TRP': 'W', 'TYR':'Y'}
    
    new_aa = ''
    for aa in aa_list:
        new_aa += aa_dict[aa]
    
    return new_aa


cwd = os.getcwd()

def get_pLDDT(file):
    pLDDT_dict = {}                                 # empty dictionary to store per-residue pLDDT scores for each PDB file
    aa_seq = ''                                     # empty string to store the amino acid sequence for each PDB file
    res_list = []                                   # list to store residue numbers
    res_start = 0
    uniprot_ID = ''                                 # set an empty uniprot ID
    uniprot_list = []
    pdb_file = open(file, "r")                      # read the pdb file
    first_res_flag = 0
    res_cnt = 0
    for line in pdb_file:

        new = line.strip().split()

        # First get the UniProt ID
        if len(new) > 1 and new[0] == 'DBREF':
            uniprot_ID += new[6]                      # fix bug from new[7] --> new[6]
            uniprot_list.append(uniprot_ID)         # Append UniProt ID to list
            total_res = int(new[4])                 # Total number of residues  
            res_cnt += total_res                    # Update the total residue number count
        
        # Then get the AA sequence 
        elif len(new) > 1 and new[0] == 'SEQRES':
            aa_seq += convert_AA(new[4:])           # Extract the 3-letter code of amino acids and append to string

        
        # Then get per-residue pLDDT scores
        elif len(new) >1 and new[0] == 'ATOM':
            first_res_flag += 1                     # Counter that gets set to 1 such that the first residue number can be stored
            if len(new) == 12:
                res_num = int(new[5]) # + (199)*scalar + fragment_num-1    # Correct for long proteins with the 200*scalar
                res_list.append(res_num)                # Store the residue number
                pLDDT = float(new[10])                  # Store the pLDDT score
                if first_res_flag == 1:
                    res_start = (res_num)           # Store the first residue in case first_res_flag is 1
                    if res_start == 0:
                        res_start = 1
            elif len(new) == 11:
                try:
                    res_num = int(new[4].split('A')[1]) # + int((199)*scalar) + fragment_num-1 
                    res_list.append(res_num)
                    pLDDT = float(new[9]) # some bug with fusion of y, z coordinates...
                    if first_res_flag == 1:
                        res_start = (res_num)
                        if res_start == 0:
                            res_start = 1
                except ValueError:
                    res_num = int(new[5]) # + int((199)*scalar) + fragment_num-1 
                    res_list.append(res_num)
                    pLDDT = float(new[9])
                    if first_res_flag == 1:
                        res_start = (res_num)
                        if res_start == 0:
                            res_start = 1
            elif len(new) == 10:
                res_num = int(new[4].split('A')[1]) # + int((199)*scalar) + fragment_num-1 
                res_list.append(res_num)
                pLDDT = float(new[8]) # some bug with fusion of y, z coordinates...
                if first_res_flag == 1:
                    res_start = (res_num)
                    if res_start == 0:
                        res_start = 1
                        
            # Store per-residue pLDDDT scores in a dictionary {residue number : pLDDT score}
            if res_num in pLDDT_dict.keys():
                pass
            else:
                pLDDT_dict[res_num] = pLDDT # pLDDT_dict.get(res_num, []) + [pLDDT]
                
    return pLDDT_dict
                        
                        


files = []
for r, d, f in os.walk(os.getcwd()):
    for file in f:
        if "template" in file:
            pass
        else:
            if file.endswith(".pdb"):
                files.append(file)

dat = np.zeros(shape=(140,len(files)))
file_cnt = -1
labels = []
for in_file in files:
    text = in_file.split('_')
    if len(text) == 2:
        labels.append(text[0])
        labels.append(text[1].split('.pdb')[0])
    else:
        labels.append(text[0])
        labels.append(text[1])
    file_cnt += 1
    data = get_pLDDT(in_file)
    cnt = -1
    for key, value in data.items():
        cnt += 1
        print(cnt, key, value)
        dat[cnt,file_cnt] = value
        
def rmsd(in1, in2):
    return np.sqrt(np.sum((in1-in2)**2)/len(in1))
    
    
# Now fix the labels and sort such that AlphaFold-DB is first 
new_labels = [ '-'.join(x) for x in zip(labels[0::2], labels[1::2]) ]
sort_idx = []
cntr = -1
for text in new_labels:
    cntr += 1
    if text == 'AlphaFold-DB':
        break
length = np.arange(0,len(files),1)
label_array = []
label_cnt = -1
for val in length:
    label_cnt +=1
    if label_cnt != cntr:
        label_array.append(int(label_cnt))

label_array.insert(0,cntr)
final_labels = []
plt.plot([0,100], [0,100], 'k-', label='y = x')
for i, item in enumerate(label_array):
    #final_labels.append(new_labels[item])
    if i == 0:
        fix = item
        print('fix = ',fix)
    elif i>0:
        plt.plot(dat[:,fix], dat[:,item], 'o', mfc='w', label=new_labels[item]+', RMSD = %.1f' % (rmsd(dat[:,fix], dat[:,item])))



    


#plt.plot([30,100], [30,100], 'k-', label='y = x')
#plt.plot(dat[:,0], dat[:,1], 'o', mfc='w', mec='r', label='AlphaFold Colab, RMSD = %.1f' % (rmsd(dat[:,0], dat[:,1])))
#plt.plot(dat[:,0], dat[:,2], 'o', mfc='w', mec='b', label='ColabFold, RMSD = %.1f' % (rmsd(dat[:,0], dat[:,2])))
plt.xlabel('pLDDT score, AlphaFold DB')
plt.ylabel('pLDDT score, Colab')
plt.xlim(30,100)
plt.ylim(30,100)
plt.legend(loc='lower right')
#plt.show()
plt.savefig('actr_compare.pdf')