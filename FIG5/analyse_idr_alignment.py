import os,sys
from Bio import AlignIO
import copy
import math
import numpy as np

"""
Given sequence alignment of an IDR
    > calculate average alignment depth
    > N sequences in MSA
    > % MSA gaps
        - all gaps
        - only if >50% of sequences have a gap
"""

#def modify_aln_dir(indir,inputd,refdict):
def read_aln_dir(indir):
    fout=open(indir+"_ALN_STATS.out.txt",'w')
    fout.write("N_ALN_SEQ\tAVERAGE_ALN_DEPTH,NORM_AVERAGE_ALN_DEPTH,PERCENT_GAPS_ALL,PERCENT_GAPS_50,PERCENT_GAPS_80\n")
    skipped=0
    nocols=0
    for f in os.listdir(indir):
        if not f.endswith(".fa"): # skip all non fasta files
            continue
        fout.write("%s"%(f)) # record name of the file
        fpath=indir+os.sep+f
        unid=f.split("_")[0]
        #gname=refdict[unid] # get genename from UNIPROT ID
        align=read_in_alignment(fpath)
        if not align: # skip any "NoneType"
            skipped+=1
            continue
        nseq=len(align)
        fout.write("\t%s"%(nseq))
        ncols=align.get_alignment_length()
        column_depth=[]
        norm_column_depth=[]
        all_gaps=0
        gaps_50p=0
        gaps_80p=0
        for i in range(0,ncols):
            column=align[:,i]
            len_init=len(column)
            nogaps=column.replace("-","")
            len_nogaps=len(nogaps)
            normdepth=len_nogaps/float(len_init)
            norm_column_depth.append(normdepth)
            column_depth.append(len_nogaps)
            if len_nogaps<len_init: # count every column with any gap
                all_gaps+=1
            if float(len_nogaps)<0.5*len_init: # count only columns with gaps over >50% of sequences
                gaps_50p+=1
            if float(len_nogaps)<0.2*len_init: # count only columns with gaps over >50% of sequences
                gaps_80p+=1
        alndepth=np.array(column_depth)
        avdepth=np.mean(alndepth)

        normalndepth=np.array(norm_column_depth)
        normavdepth=np.mean(normalndepth)

        #print("\t%s\t%.2f\t%.2f\t%.2f\n"%(avdepth,all_gaps/ncols,gaps_50p/ncols,gaps_80p/ncols))
        #sys.exit()
        try:
            fout.write("\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\n"%(avdepth,normavdepth,all_gaps/ncols,gaps_50p/ncols,gaps_80p/ncols))
        except ZeroDivisionError:
            nocols+=1
            continue
    print("NO DATA FOR %s SKIPPED and %s with columns 0"%(skipped,nocols))
    fout.close()

def read_in_alignment(infile):
    try:
        alignment = AlignIO.read(open(infile), "fasta")
        #print("Alignment length %i" % alignment.get_alignment_length())
    except ValueError:
        alignment=None

    return alignment

if __name__=="__main__":
    inputdir=sys.argv[1] # directory with all alignment files
    read_aln_dir(inputdir)
