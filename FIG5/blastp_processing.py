import os,sys

"""
loop over all files
make two column files for each threshold:

1. uniprot Id gene name -- for IDRs that pass the
"""

fdir=sys.argv[1]
evalue_thresh=sys.argv[2]
fout1=open("60_COVERAGE_30_IDENTITY_ALL_EVALUES.txt","w")
fout2=open("60_COVERAGE_30_IDENTITY"+"_"+str(evalue_thresh),"w")

fout1.write("%s\t%s\t%s\t%s\n"%("IDRID","COVERAGE","IDENTITY","E-VALUE"))
fout2.write("%s\t%s\t%s\t%s\n"%("IDRID","COVERAGE","IDENTITY","E-VALUE"))

for f in os.listdir(fdir):
    idrinit=f.split("_SELE_HUMAN_SEQ.fa.blastp.out.txt")[0]
    idrbounds=idrinit.split("ENSEMBL_ORTHOLOGUES_ALN_IDR_")[-1]
    unid=f.split("_")[0]
    fin=open(fdir+os.sep+f,'r')
    evalues=[]
    coverages=[]
    identities=[]
    seqids=[]
    for line in fin:
        if line.startswith("#"):
            continue
        else:
            stripped=line.strip()
            splitted=stripped.split()
            #unid=splitted[0].split("|")[1]
            #genename=splitted[0].split("|")[-1].split(".diso:")[0]
            #disoid=splitted[0].split("|")[-1].split(".diso")[1]
            try:
                if float(splitted[2])>30 and float(splitted[-2])>60:
                    seqid=unid+"_"+idrbounds
                    seqids.append(seqid)
                    evalues.append(float(splitted[-4]))
                    coverages.append(float(splitted[-2]))
                    identities.append(float(splitted[2]))
            except IndexError:
		print("INDEX ERROR FOR LINE %s"%(line))
    if len(evalues)==0: # if nothing passed the filter
        continue

    else:
        min_value=min(evalues)
        min_indx = evalues.index(min_value) # output only one per seuqence with min_evalue
        fout1.write("%s\t%s\t%s\t%s\n"%(seqids[min_indx],coverages[min_indx],identities[min_indx],evalues[min_indx]))
        if float(evalues[min_indx])<float(evalue_thresh):
            fout2.write("%s\t%s\t%s\t%s\n"%(seqids[min_indx],coverages[min_indx],identities[min_indx],evalues[min_indx]))
fout1.close()
fout2.close()
