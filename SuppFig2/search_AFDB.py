import os, sys
import numpy as np
import matplotlib.pyplot as plt


pLDDT_dis = 'updated_UP000005640_9606_HUMAN_pLDDT_scores__DISORDERED.txt' #_THRESH_above_70.txt'
pLDDT = 'updated_UP000005640_9606_HUMAN_pLDDT_scores.txt'
morfs = 'morfs.txt'
disprot = 'human_disprot.txt'

morf_dict = {}
for line in open(morfs, 'r'):
    new = line.strip().split()
    uniprot = new[0]
    resn = int(new[1])
    if uniprot in morf_dict.keys():
        morf_dict[uniprot].append(resn)
    else:
        morf_dict[uniprot] = [resn]
        
all_morf_counter = 0
for key, value in morf_dict.items():
    all_morf_counter += len(value)
    
    
count = 0
disprot_dict = {}
for line in open(disprot, 'r'):
    new = line.strip().split()
    if len(new) > 1:
        count+=1
        #print(count)
        try:
            uniprot_temp = new[0]#.split('-')[0]
            if '-' in uniprot:
                uniprot = new[0].split('-')[0]
            else:
                uniprot = uniprot_temp
            start = int(new[1])
            end = int(new[2])
        except IndexError:
            print(line)
        if uniprot in disprot_dict.keys():
            disprot_dict[uniprot].append(np.arange(start,end+1,1))
        else:
            disprot_dict[uniprot] = [np.arange(start, end+1,1)]
            
            
for key, value in disprot_dict.items():
    temp_dict = {}
    temp_dict[key] = []
    for x in range(len(value)):
        for y in value[x]:
            if y in temp_dict.items():
                print('fail')
                sys.exit()
            else:
                temp_dict[key].append(y)
    
            
all_morf_counter = 0
for key, value in disprot_dict.items():
    for x in range(len(value)):
        all_morf_counter += len(value[x])
print(all_morf_counter)
        
#print(morf_dict)
#print(morf_dict['X6R8D5'])


def compare_AF(morf_dic_data, pLDDT_file = pLDDT):
    # Read per-residue pLDDT file and look for matches in UniProt ID
    pLDDT_list = []
    prot_list = []
    conf_list = []
    fin = open(pLDDT_file, 'r')
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
                #print('~~~ YES!!!! ---> ', uniprot)
                for x in range(len(morf_dic_data[uniprot])):
                #for res in morf_dic_data[uniprot][0]:
                    for res in morf_dic_data[uniprot][x]:
                        if int(res) == int(resn):
                            #print(uniprot, res, pLDDT_score)
                            pLDDT_list.append(float(pLDDT_score))
                            if float(pLDDT_score) >= 70:
                                conf_list.append(uniprot)
    pLDDT_list = np.asarray(pLDDT_list)
    #print('~~~~ --> ',list(set(prot_list)))
    print('\n\n')
    print('~~~~ --> ',len(list(set(prot_list))))
    print('conf_list = ', len(list(set(conf_list))))

    return pLDDT_list, list(set(prot_list))
    
    
def make_AF_dict(morf_dic_data, pLDDT_file = pLDDT):
    # Read per-residue pLDDT file and look for matches in UniProt ID
    pLDDT_list = []
    prot_list = []
    test = {}
    res_test = {}
    fin = open(pLDDT_file, 'r')
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
                #print('~~~ YES!!!! ---> ', uniprot)
                for res in morf_dic_data[uniprot][0]:
                    if int(res) == int(resn):
                        if uniprot not in test.keys():
                            test[uniprot] = [float(pLDDT_score)]
                            res_test[uniprot] = [int(resn)]
                        else:
                            test[uniprot].append(float(pLDDT_score))
                            res_test[uniprot].append(int(resn))
                        #print(uniprot, res, pLDDT_score)
                        pLDDT_list.append(float(pLDDT_score))
    pLDDT_list = np.asarray(pLDDT_list)
    #print('~~~~ --> ',list(set(prot_list)))
    print('\n\n')
    print('~~~~ --> ',len(list(set(prot_list))))

    return test, res_test, pLDDT_list
    
#new_dat =(make_AF_dict(morf_dict))
new_dat = make_AF_dict(disprot_dict)
new_dat_dis = make_AF_dict(disprot_dict, pLDDT_dis)

missing_res = {}
for key, value in new_dat[1].items():
    for key2, value2 in new_dat_dis[1].items():
        if key == key2:
            for res in value:
                if int(res) in value2:
                    pass
                else:
                    if key in missing_res.keys():
                        missing_res[key].append(int(res))
                    else:
                        missing_res[key] =[int(res)]
        else:
            pass
                
#for key, value in missing_res.items():
#    print(key,value)


#sys.exit()

counter = 0
for key, value in new_dat[0].items():
    #print(key, value,'\n')
    counter += len(value)
    
print(counter/all_morf_counter)
print(counter)
print(all_morf_counter)

tot = compare_AF(disprot_dict, pLDDT)
dis_tot = compare_AF(disprot_dict, pLDDT_dis)

print('number of proteins in set = ',len(tot[1]))
print(len(dis_tot[1]))

empty_list = []
for i in tot[1]:
    if i in dis_tot[1]:
        pass
    else:
        empty_list.append(i)
        
#print(empty_list)
print(len(empty_list))



print(len(tot[0]))

plt.hist(tot[0], bins=50)
plt.show()

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

print('\n here are fractions...\n')
print(fractions(tot[0]))
sys.exit()

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
        #iter_loadtxt.rowlength = len(line)

    data = np.fromiter(iter_func(), dtype=dtype)
    #data = data.reshape((-1, iter_loadtxt.rowlength))
    return data
    
idrs = iter_loadtxt('updated_UP000005640_9606_HUMAN_pLDDT_scores__DISORDERED.txt')


# (4) Plot
idr_counts, idr_bins = np.histogram(idrs, bins=50)
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax2.hist(idrs, bins=50, facecolor='None', edgecolor='orange', align='right', histtype='step', linewidth=0.5, label='SPOT-Disorder IDRs', density=False)
ax1.hist(tot[0], bins=idr_bins, alpha=0.5, facecolor='maroon', edgecolor='None', align='right', label='DisProt IDRs',density=False)
    
ax1.set_xlabel('pLDDT score')
ax1.set_ylabel('Count, DisProt IDRs', color='k')
ax2.set_ylabel('Count, SPOT-Disorder IDRs', color='orange')
ax1.set_xticks(np.arange(0, 101, 10))
ax1.set_xlim(0,100)
plt.legend(loc='upper left', fontsize=8)
plt.show()
