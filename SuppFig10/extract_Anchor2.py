import os, sys
import numpy as np
import matplotlib.pyplot as plt

def get_ranges(DB_ranges):
    # read DB file (UniProt start end)
    uniprot_dict = {}
    for line in open(DB_ranges, 'r'):
        new_line = line.strip().split()
        if new_line[0] != '#':
            uniprot = new_line[0]
            res_start = int(new_line[1])
            res_end = int(new_line[2])
            res_range = np.arange(res_start, res_end+1)
            if uniprot in uniprot_dict.keys():
                uniprot_dict[uniprot].extend(res_range)
            else:
                uniprot_dict[uniprot] = []
                uniprot_dict[uniprot].extend(res_range)
    # now remove duplicates
    clean_uniprot_dict = {}
    for key, val in uniprot_dict.items():
        keep = list(set(val))
        clean_uniprot_dict[key] = []
        clean_uniprot_dict[key].extend(keep)
    return clean_uniprot_dict

def extract_anchor_scores(DB_ranges_dict, all_anchor_scores, dataset_name, output_dir):
    # now read Anchor2 file and keep overlaps
    out_path = os.path.join(output_dir, dataset_name)
    out = open(out_path+'_anchor2.txt', 'w')
    out.write('%s\t%s\t%s\t%s\n'%('UniprotID', 'residue', 'expected_score', 'anchor2_score'))
    for line in open(all_anchor_scores, 'r'):
        new = line.strip().split()
        if new[0] != '#':
            uniprot = new[0]    # uniprot ID
            anch = float(new[2])    # Anchor2
            resn = int(new[4])      # Residue number
            for key, value in DB_ranges_dict.items():
                if key == uniprot:
                    for res in value:
                        if ((resn == int(res)) and (dataset_name != 'CheZOD') and (dataset_name != 'Anchor_flex_link')):
                            out.write('%s\t%i\t1\t%f\n' % (key, res, 100*anch))            #we multiply the anchor score by 100 so that we get a range of scores 0-100 isntead of 0-1
                        if ((resn == int(res)) and ((dataset_name == 'CheZOD') or (dataset_name == 'Anchor_flex_link'))):
                            out.write('%s\t%i\t0\t%f\n' % (key, res, 100*anch))
    out.close()


if __name__=='__main__':
    DB_ranges = sys.argv[1]     #the directory with Uniprot_start_end.txt files of each dataset
    all_anchor_scores = sys.argv[2] #the direcroty with all anchor scores (the output of Anchor2 run) of each dataset
    curr_dir = os.getcwd()
    DB_ranges_path = os.path.join(curr_dir, DB_ranges)
    all_anchor_scores_path = os.path.join(curr_dir, all_anchor_scores)
    out_dir = 'Anchor2_scores'
    out_path = os.path.join(curr_dir, out_dir)
    if not os.path.exists(out_path):          # create the output directory if it doesn't already exist, else: pass
        os.makedirs(out_path)
    for file in os.listdir(DB_ranges):
        curr_dataset = file[:-15]
        for scores_file in os.listdir(all_anchor_scores):
            curr_dataset_scores = scores_file[:-11]
            if curr_dataset == curr_dataset_scores:
                curr_dataset_DB_path = os.path.join(DB_ranges_path, file)
                uniprot_ranges = get_ranges(curr_dataset_DB_path)
                curr_dataset_anchor_scores_path = os.path.join(all_anchor_scores_path, scores_file)
                extract_anchor_scores(uniprot_ranges, curr_dataset_anchor_scores_path, curr_dataset, out_path)
