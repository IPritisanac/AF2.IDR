import os,sys
import numpy as np


"""
Given position-based conservation measure for every residue
    - take column 1 for all lines != #
    -
normalized: False
# align_column_number	score	column
0	0.222811	MMMMIMMLLVVVMMMMMMMMMVVMMMVMMMMVMVMVVMMLMMVMMVVMMMV-MMVMMMMMVMVLMVVTVVMMMMVVMMVMM-VMVMVLMMMMIVMMMMVM-K-SNMLTLPPLLLLV-NN-NNV-KKKKKKKKRFLFQ-FFLVVLLFLIVLKLV-QKLL-F-K-LFKKLKKWL-F-KLTKKKRNNKNKI-QR--KCS-K
    N files with if > 60 % of values > 30%
    count file

    get average positional conservation across the sequence
    output = SEQ ID AVERAGE_POS_CONSERVATION LOWEST_CONS HIGHEST_CONS
"""

#def modify_aln_dir(indir,inputd,refdict):
def read_cons_dir(indir):
    noutput=indir+"IDR_CONS_STATS.out.txt"
    fout=open(noutput,'w')
    fout.write("#SEQID\tAVERAGE_POS_CONS\tMAX_POS_CONS\tMIN_POS_CONS\n")
    all_files=0
    cons_files=0
    for f in os.listdir(indir):
        if not f.endswith("CONS.out.txt"): # go through CONS files only
            continue
        all_files+=1
        fpath=indir+os.sep+f
        fout.write("%s"%(f)) # record name of the file
        fin=open(fpath,'r')
        cons_array=[]
        for line in fin:
            if line.startswith("#"):
                continue
            splitted=line.split("\t")
            cons_array.append(float(splitted[1]))
        cnt=0
        for val in cons_array:
            if float(val)>0.3:
                cnt+=1
        if cnt>0.6*len(cons_array):
            cons_files+=1
        #print("BEFORE")
        #print(cons_array)
        cons_array=np.array(cons_array)
        fcons_array = np.where(cons_array == -1000, np.nan, cons_array)
        #print("AFTER")
        #print(fcons_array)

        average_cons=np.nanmean(fcons_array)
        max_cons=np.nanmax(fcons_array)
        min_cons=np.nanmin(fcons_array)
        #print(len(cons_array))
        nonzero_count=np.count_nonzero(~np.isnan(fcons_array))
        if nonzero_count<10:
            fout.write('\t%s\t%s\t%s\n'%('NA','NA','NA'))
        else:
            fout.write('\t%.2f\t%.2f\t%.2f\n'%(average_cons,max_cons,min_cons))
    print(cons_files,all_files,float(cons_files)/all_files)

if __name__=="__main__":
    inputdir=sys.argv[1] # directory with all alignment files
    read_cons_dir(inputdir)
