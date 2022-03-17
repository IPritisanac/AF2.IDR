import os, sys
import matplotlib.pyplot as plt 
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import time

# to solve the np.genfromtxt memory issue: https://stackoverflow.com/questions/8956832/python-out-of-memory-on-large-csv-file-numpy
# modified function to work on my files
def iter_loadtxt(filename, skiprows=0, dtype=float, flag='None'):
    def iter_func():
        with open(filename, 'r') as infile:
            for _ in range(skiprows):
                next(infile)
            for line in infile:
                line = line.strip().split()
                if line[0] == '#':
                    pass
                else:
                    yield dtype(line[1])

    data = np.fromiter(iter_func(), dtype=dtype)

    return data
    
def fractions(dat_array):
    """ compute the % of various pLDDT score thresholds """
    vconf = np.count_nonzero(dat_array >= 90)
    conf = np.count_nonzero(dat_array >= 70) - vconf
    low = np.count_nonzero(dat_array >= 50) - vconf - conf
    vlow = np.count_nonzero(dat_array < 50)
    
    tot = (vconf) + (conf) + (low) + (vlow)
    
    vconf_frac = (vconf)/tot 
    conf_frac = (conf)/tot 
    low_frac = (low)/tot 
    vlow_frac = (vlow)/tot
    
    return vconf_frac, conf_frac, low_frac, vlow_frac, tot
    
    
def make_AF_dict(pLDDT):
    af_dict = {}
    fin = open(pLDDT, 'r')
    for line in fin:
        new = line.strip().split()
        if new[0] == '#':
            pass
        else:
            uniprot = new[0].split('_')[0]
            resn = int(new[3])
            score = float(new[1])
            if uniprot in af_dict.keys():
                af_dict[uniprot].append(score)
            else:
                af_dict[uniprot] = [score]
    return af_dict
    
    
### problem with slicing pLDDT list 
def filter_AF(af_dict, cons_dict):
    pLDDT_list = []
    for key, value in af_dict.items():
        for other_key, other_value in cons_dict.items():
            if key == other_key:
                temp = [int(x)-1 for x in other_value[0]]
                filter = np.asarray(value)[temp]
                pLDDT_list.append(np.asarray(filter))
                
    a = np.asarray(pLDDT_list)
    b = np.concatenate( a, axis=0 )
    return b
    
def compare_AF(pLDDT, morf_dic_data):
    # Read per-residue pLDDT file and look for matches in UniProt ID
    pLDDT_list = []
    cnt = 0
    prot_list = []
    fin = open(pLDDT, 'r')
    for line in fin:
        new = line.strip().split()
        if new[0] == '#':
            pass
        else:
            uniprot = new[0].split('_')[0]
            resn = new[3]
            pLDDT_score = new[1]
            if uniprot in morf_dic_data.keys():
                prot_list.append(uniprot)
                for res in morf_dic_data[uniprot][0]:
                    if int(res) == int(resn):
                        cnt += 1
                        pLDDT_list.append(float(pLDDT_score))
    pLDDT_list = np.asarray(pLDDT_list)
    fin.close()

    print(cnt)
    return np.asarray(pLDDT_list[0])
    
    
#proteome = 'updated_UP000005640_9606_HUMAN_pLDDT_scores.txt'
#dat = iter_loadtxt(proteome)
#print(fractions(dat))

# (1) Load IDRs / PFAM-corrected IDRs
idr_file = 'out/UP000005640_9606_HUMAN_pLDDT_scores__DISORDERED.txt'
idrs = iter_loadtxt(idr_file)
print(fractions(idrs))


# (2) Load structured regions
structured_file = 'out/UP000005640_9606_HUMAN_pLDDT_scores__STRUCTURED.txt'
structured = iter_loadtxt(structured_file)
print('\nloaded structured\n')  
print(fractions(structured))


# (3) Create array for the entire proteome
all = np.zeros(shape=(len(idrs)+len(structured),))
all[0:len(idrs)] = idrs
all[len(idrs)::] = structured
print('\nloaded all\n')


# (4) Plot the histogram
idr_counts, idr_bins = np.histogram(idrs, bins=50)
plt.hist(all, bins=50, facecolor='None', edgecolor='black', align='right', histtype='step', linewidth=0.5, label='Proteome', density=False)
plt.hist(structured, bins=50, alpha=0.5, facecolor='darkblue', edgecolor='k', align='right', label='Ordered',density=False)
plt.hist(idrs, bins=50, alpha=0.5, facecolor='orange', edgecolor='k', align='right', label='Disordered',density=False)
plt.xticks(np.arange(0, 101, 10))
plt.xlim(0,100)
plt.xlabel('pLDDT score')
plt.ylabel('Count')
plt.legend(loc='upper left', fontsize=8)

# Make the /plot directory if it doesn't exist
if not os.path.exists(os.getcwd()+'/plot'):
    os.makedirs(os.getcwd()+'/plot')
output_dir = os.getcwd()+'/plot'
plt.savefig(output_dir+'/histogram_pLDDT_IDRs_vs_structured_vs_all.pdf')

# Plot the zoomed in histogram; update xlim and ylim as needed
plt.clf()
plt.hist(idrs, bins=50, alpha=0.5, facecolor='orange', align='right', label='Disordered',density=False)
plt.xlabel('pLDDT score')
plt.ylabel('Count')
plt.xticks(np.arange(0, 101, 10))
plt.xlim(70,100)
plt.ylim(0,40000)
plt.legend(loc='upper left')
plt.savefig(output_dir+'/histogram_pLDDT_IDRs_vs_structured_vs_all_zoom.pdf')

# Now write out the pLDDT stats to file 
fout = open(output_dir+'/stats.txt', 'w')
fout.write('# %_>=90 %_90<x>=70 %_70<x>=50 %_<50 total_residues\n')
fout.write('IDRs: %.3f %.3f %.3f %.3f %i\n' % fractions(idrs))
fout.write('ordered: %.3f %.3f %.3f %.3f %i\n' % fractions(structured))
fout.write('all: %.3f %.3f %.3f %.3f %i\n' % fractions(all))
fout.close()


sys.exit()