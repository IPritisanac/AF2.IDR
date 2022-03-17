import os,sys
from Bio import AlignIO
import copy
"""
Given IDR boundaries of the human sequence and a sequence alignment of the orthologues
    > cut-out IDRs from the orthologues
"""
def modify_aln_dir(indir,inputd,outputp):
    dbounds=read_in_diso(inputd)
    noutput=indir+"N_HOMOLOGS.out.txt"
    fn=open(noutput,'w')
    for f in os.listdir(indir):
        if f.endswith("_ALN.fa"):
            fpath=indir+os.sep+f
            unid=f.split("_")[0]
            align=read_in_alignment(fpath)
            if not align:
                continue
            try:
                for entry in dbounds[unid]:
                    startd=entry[0]
                    endd=entry[1]
                    foutpath=outputp+os.sep+f[:-3]+"_IDR_"+str(startd)+"_"+str(endd)+".fa"
                    fout=open(foutpath,'w')
                    unid,startd,endd=get_new_diso_bounds(align,startd,endd)
                    idralign=cut_diso_alignment(align,startd,endd)
                    fn.write("%s\t%s\n"%(f[:-3]+"_IDR_"+str(startd)+"_"+str(endd),len(idralign.keys())))
                    for key,value in idralign.items():
                        fout.write(">%s\n"%(key))
                        fout.write("%s\n"%(value))                    
            except KeyError:
                print("NO DISORDER BOUNDS FOR %s"%(unid))

def read_in_diso(infile):
    disobounds={}
    fin=open(infile,'r')
    for line in fin:
        stripped=line.strip()
        splitted=stripped.split('\t')
        uniprot=splitted[0].split(":")[0]
        disobounds.setdefault(uniprot,[]).append((int(splitted[1]),int(splitted[2])))
        
    return disobounds

def read_in_alignment(infile):
    try:
        alignment = AlignIO.read(open(infile), "fasta")
        print("Alignment length %i" % alignment.get_alignment_length())
    except ValueError:
        alignment=None
    
    return alignment
        
def get_new_diso_bounds(aln,dstart,dend):
    for record in aln:
        if "HUMAN" in record.id:
            uniprot=record.id.split("|")[-1]
            dstartnew,dendnew=enumerate_alignment(dstart,dend,record.seq)
            break
    
    return uniprot,dstartnew,dendnew

# cut alignment of orthologs to extract IDRs
def enumerate_alignment(diso_start,diso_end,humanseq):
    cnt_gap=0
    cnt_aa=0
    diso_start_new=copy.deepcopy(diso_start) # initiate 
    diso_end_new=copy.deepcopy(diso_end)
    #get new N for start
    for aa in list(humanseq):
        if aa!="-":
            cnt_aa+=1
        if aa=="-":
            cnt_gap+=1
        if cnt_aa==diso_start: # when at the start of disordered region -- add gap count
            diso_start_new+=cnt_gap
            break
    #get new N and new gap count for end
    cnt_gap=0
    cnt_aa=0   
    for aa in list(humanseq):
        if aa!="-":
            cnt_aa+=1
        if aa=="-":
            cnt_gap+=1
        if cnt_aa==diso_end:
            diso_end_new+=cnt_gap
            break
    #print(diso_start_new)
    #print(diso_end_new)
            
    return diso_start_new,diso_end_new

def cut_diso_alignment(aln,startd,endd):
    idrs={}
    for record in aln:
        start=startd-1
        end=endd
        
        idr=record.seq[start:end]
        
        idrs.setdefault(record.id,idr)
        
    return idrs
  
if __name__=="__main__":
    inputdir=sys.argv[1] # directory with all alignment files
    inputdiso=sys.argv[2]
    outputdir=sys.argv[3]
    try:
        os.makedirs(outputdir)
    except FileExistsError:
        print("DIR EXISTS!")
    modify_aln_dir(inputdir,inputdiso,outputdir)
    
