import os, sys
import numpy as np
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
import time


pLDDT = 'UP000005640_9606_HUMAN_pLDDT_scores.txt'

# Create dictionary with Kyte-Doolittle hydropathy values
kyte_doolittle = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
       'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
       'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
       'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }

# Normalize hydropathy values to be between 0-1 according to Uversky 2000 Proteins
norm_kyte_doolittle = {}       
values = []
for key, val in kyte_doolittle.items():
    values.append(val)
for key, val in kyte_doolittle.items():
    norm_kyte_doolittle[key] = (val-min(values)) / (max(values) - min(values))
    
    

def get_AF(in_pLDDT):
    # Read per-residue pLDDT file and look for matches in UniProt ID
    uniprot_seqs = {}
    fin = open(in_pLDDT, 'r')
    cnt = 0
    for line in fin:
        cnt += 1
        new = line.strip().split()
        if new[0] == '#':
            pass
        else:
            uniprot = new[0].split('_')[0]
            res_type = new[2]
            res_num = new[3]
            pLDDT_score = new[1]
            if uniprot in uniprot_seqs.keys():
                uniprot_seqs[uniprot] += str(res_type)
            else:
                uniprot_seqs[uniprot] = str(res_type)
                
    fin.close()
    
    return uniprot_seqs


    

low_disorder = get_AF('updated_UP000005640_9606_HUMAN_pLDDT_scores__DISORDERED_THRESH_below_50.txt')
high_disorder = get_AF('updated_UP000005640_9606_HUMAN_pLDDT_scores__DISORDERED_THRESH_above_90.txt')
structured = get_AF('updated_UP000005640_9606_HUMAN_pLDDT_scores__STRUCTURED.txt')

print(len(low_disorder))
print(len(high_disorder))
print(len(structured))



def hydropathy(in_seq, window=5, normed=True):
    """ calculate mean hydropathy using Kyte-Doolittle values"""
    #seq = ''.join(list(in_dict.values())[0])
    if normed:
        hydro_dic = norm_kyte_doolittle
    else:
        hydro_dic = kyte_doolittle
    hydro_array = np.asarray([hydro_dic[x] for x in in_seq])
    
    smoothed_array = []
    for i in range((len(hydro_array)-window)+1):
        smoothed_array.append( sum(hydro_array[i:i+window])/window )
    smoothed_array = np.asarray(smoothed_array)
    
    return np.mean(smoothed_array)
    
def charge(in_seq):
    #seq = ''.join(list(in_dict.values())[0])
    aa_cnt = Counter(in_seq)
    charge = 0
    charge -= aa_cnt['E']
    charge -= aa_cnt['D']
    charge += aa_cnt['K']
    charge += aa_cnt['R']
    
    return charge/float(len(in_seq))
    

def uversky(input_dat):
    """ Dictionary with UniProtID : [mean hydropathy, mean net charge at pH 7] """
    uversky_dict = {}
    for key, val in input_dat.items():
        uversky_dict[key] = [hydropathy(val,5), np.abs(charge(val))]
        
    return uversky_dict
        
low_uversky = (uversky(low_disorder))
high_uversky =uversky(high_disorder)
structured_uversky = uversky(structured)

def plot_uversky(input_dict):
    for key, val in input_dict.items():
        plt.plot(val[0], val[1], 'o')
    plt.xlabel('Mean hydrophhobicity')
    plt.ylabel('Mean net charge')
    plt.xlim(0,0.7)
    plt.ylim(-0.02,0.7)
    plt.plot(np.linspace(0,1,11), 2.785*np.linspace(0,1,11) - 1.151, 'k-')
    plt.show()
    #return fig
    

#plot_uversky(uversky)
def make_arrays(in_dict):
    hydro_vals = []
    charge_vals =[]
    for key, val in in_dict.items():
        hydro_vals.append(val[0])
        charge_vals.append(val[1])
    hydro_vals = np.asarray(hydro_vals)
    charge_vals = np.asarray(charge_vals)
    
    return hydro_vals, charge_vals

hydro_low, charge_low = make_arrays(low_uversky)
hydro_high, charge_high = make_arrays(high_uversky)
hydro_structured, charge_structured = make_arrays(structured_uversky)

print(len(hydro_low))
print(len(hydro_high))
print(len(hydro_structured))

#sys.exit()    
    
'''    
#print(np.histogram2d(hydro_vals, charge_vals, bins=10))
plt.hist2d(hydro_vals, charge_vals, bins=100, range=[ [np.nanmin(hydro_vals), np.nanmax(hydro_vals)], [np.nanmin(charge_vals), np.nanmax(charge_vals)] ] )
plt.show()
#plt.show()'''

#### PLOT scatter data with 1D histograms ####

def scatter_hist(x1,y1,x2,y2, bin_number=100, ranges=(0,1), labels=[r'IDR$_{low}$','Ordered'], colors=['lightsalmon', 'mediumblue'], density_flag=False):
    """ function to plot an x/y scatter with histograms projected onto x and y axes"""
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure(figsize=(6, 4))

    # set up plot boundaries/ticks/etc
    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', bottom=True, top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)

    # the scatter plot:
    ax_scatter.plot(x1, y1, 'o', color=colors[0], label=labels[0])
    ax_scatter.plot(x2, y2, 'o', color=colors[1], label=labels[1])

    # now determine nice limits by hand:
    ax_scatter.set_xlim(-0.05,1.1)
    ax_scatter.set_ylim(-0.05,1.1)

    bins = bin_number
    weights = [np.ones_like(z)/float(len(z)) for z in [x1, y1, x2, y2]]
    ax_histx.hist(x1, bins=bins, range=ranges, histtype='bar', color=colors[0], alpha=0.5, weights=weights[0]) # density=density_flag
    ax_histy.hist(y1, bins=bins, range=ranges, histtype='bar',orientation='horizontal', color=colors[0], alpha=0.5, weights=weights[1]) # density=density_flag
    ax_histx.hist(x2, bins=bins, range=ranges, histtype='bar', color=colors[1], alpha=0.5, weights=weights[2]) # density=density_flag
    ax_histy.hist(y2, bins=bins, range=ranges, histtype='bar',orientation='horizontal', color=colors[1], alpha=0.5, weights=weights[3]) # density=density_flag

    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    
    ax_scatter.set_xlabel('Mean hydropathy')
    ax_scatter.set_ylabel('Mean net charge')
    ax_histx.set_ylabel('Count')
    ax_histy.set_xlabel('Count')

    ax_scatter.legend(loc='upper right')
    
    #plt.tight_layout()
    
    #plt.show()
    return ax_scatter, ax_histx, ax_histy


# PLOT for IDR low vs. structured and IDR high vs. structured

x1,y1,x2,y2 = hydro_low, charge_low, hydro_structured, charge_structured
bin_number=100
ranges=(0,1)
labels=[r'IDR$_{low-pLDDT}$','Ordered']
colors=['lightsalmon', 'mediumblue']
density_flag=False

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a rectangular Figure
fig = plt.figure(figsize=(6, 4))

# set up plot boundaries/ticks/etc
ax_scatter = plt.axes(rect_scatter)
ax_scatter.tick_params(direction='in', bottom=True, top=True, right=True)
ax_histx = plt.axes(rect_histx)
ax_histx.tick_params(direction='in', labelbottom=False)
ax_histy = plt.axes(rect_histy)
ax_histy.tick_params(direction='in', labelleft=False)

# the scatter plot:
ax_scatter.plot(x1, y1, 'o', color=colors[0], label=labels[0])
ax_scatter.plot(x2, y2, 'o', color=colors[1], label=labels[1])

# now determine nice limits by hand:
ax_scatter.set_xlim(-0.05,1.1)
ax_scatter.set_ylim(-0.05,1.1)

bins = bin_number
weights = [np.ones_like(z)/float(len(z)) for z in [x1, y1, x2, y2]]
ax_histx.hist(x1, bins=bins, range=ranges, histtype='bar', color=colors[0], alpha=0.5, weights=weights[0]) # density=density_flag
ax_histy.hist(y1, bins=bins, range=ranges, histtype='bar',orientation='horizontal', color=colors[0], alpha=0.5, weights=weights[1]) # density=density_flag
ax_histx.hist(x2, bins=bins, range=ranges, histtype='bar', color=colors[1], alpha=0.5, weights=weights[2]) # density=density_flag
ax_histy.hist(y2, bins=bins, range=ranges, histtype='bar',orientation='horizontal', color=colors[1], alpha=0.5, weights=weights[3]) # density=density_flag

ax_histx.set_xlim(ax_scatter.get_xlim())
ax_histy.set_ylim(ax_scatter.get_ylim())

ax_scatter.set_xlabel('Mean hydropathy')
ax_scatter.set_ylabel('Mean net charge')
ax_histx.set_ylabel('Fraction')
ax_histy.set_xlabel('Fraction')

import string 
ax_histx.text(-0.1, 1.3, string.ascii_uppercase[0], transform=ax_scatter.transAxes, 
        size=14, weight='bold')

ax_scatter.legend(loc='upper right')
plt.show()
#plt.savefig('update_lowIDR_vs_structured_Uversky.png')


# PLOT for IDR high vs. structured and IDR high vs. structured
plt.clf()

x1,y1,x2,y2 = hydro_high, charge_high, hydro_structured, charge_structured
bin_number=100
ranges=(0,1)
labels=[r'IDR$_{high-pLDDT}$','Ordered']
colors=['orangered', 'mediumblue']
density_flag=False

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a rectangular Figure
fig = plt.figure(figsize=(6, 4))

# set up plot boundaries/ticks/etc
ax_scatter = plt.axes(rect_scatter)
ax_scatter.tick_params(direction='in', bottom=True, top=True, right=True)
ax_histx = plt.axes(rect_histx)
ax_histx.tick_params(direction='in', labelbottom=False)
ax_histy = plt.axes(rect_histy)
ax_histy.tick_params(direction='in', labelleft=False)

# the scatter plot:
ax_scatter.plot(x1, y1, 'o', color=colors[0], label=labels[0])
ax_scatter.plot(x2, y2, 'o', color=colors[1], label=labels[1])

# now determine nice limits by hand:
ax_scatter.set_xlim(-0.05,1.1)
ax_scatter.set_ylim(-0.05,1.1)

bins = bin_number
weights = [np.ones_like(z)/float(len(z)) for z in [x1, y1, x2, y2]]
ax_histx.hist(x1, bins=bins, range=ranges, histtype='bar', color=colors[0], alpha=0.5, weights=weights[0]) # density=density_flag
ax_histy.hist(y1, bins=bins, range=ranges, histtype='bar',orientation='horizontal', color=colors[0], alpha=0.5, weights=weights[1]) # density=density_flag
ax_histx.hist(x2, bins=bins, range=ranges, histtype='bar', color=colors[1], alpha=0.5, weights=weights[2]) # density=density_flag
ax_histy.hist(y2, bins=bins, range=ranges, histtype='bar',orientation='horizontal', color=colors[1], alpha=0.5, weights=weights[3]) # density=density_flag

ax_histx.set_xlim(ax_scatter.get_xlim())
ax_histy.set_ylim(ax_scatter.get_ylim())

ax_scatter.set_xlabel('Mean hydropathy')
ax_scatter.set_ylabel('Mean net charge')
ax_histx.set_ylabel('Fraction')
ax_histy.set_xlabel('Fraction')

import string 
ax_histx.text(-0.1, 1.3, string.ascii_uppercase[1], transform=ax_scatter.transAxes, 
        size=14, weight='bold')

ax_scatter.legend(loc='upper right')
plt.show()
#plt.savefig('updated_highIDR_vs_structured_Uversky.png')




sys.exit()