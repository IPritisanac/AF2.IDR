import os,sys

# Enter the file names for the per-residue pLDDT and disorder prediction files 
plddt_input = 'out/UP000005640_9606_HUMAN_pLDDT_scores.txt'       # per-residue pLDDT scores
diso_file = "UP000005640_9606_SPOTD.out.txt"                      # disorder prediction file 

# Check for correction of PFAM domains in the disorder prediction input file 
if "NO_PFAM" in diso_file:
    pfam_flag = 'PFAM-corrected'
else:
    pfam_flag = ''

# Store disorder predictions in a dictionary
diso_dict={}
fin2=open(diso_file,'r')
for line in fin2:
    stripped=line.strip()
    splitted=stripped.split()
    gene_name=splitted[0].split("|")[-1].split(".diso")[0]
    uniprot_name = splitted[0].split("|")[1]
    start=int(splitted[1])
    end=int(splitted[2]) 
    diso_dict.setdefault(uniprot_name,[]).append((start,end)) 
fin2.close()

# Open the pLDDT file for reading and create 2 new files that are separated into ordered/disordered regions 
fin=open(plddt_input,'r')
fout1=open((plddt_input[:-4]+"_%s_DISORDERED.txt" % pfam_flag),'w')
fout2=open((plddt_input[:-4]+"_%s_STRUCTURED.txt" % pfam_flag),'w')
skipped_list = []
skipped=0
for line in fin:
    if line.startswith("#"):
        continue
    stripped=line.strip()
    splitted=stripped.split()
    gene=splitted[0]            # this is UniProt ID
    
    try:
        covered=False
        for region in diso_dict[gene]:
            start=region[0]
            end=region[1]
            diso_range=range(start,end+1)
            if int(splitted[3]) in diso_range:
                fout1.write(line)
                covered=True
                break
        if not covered:
            fout2.write(line)
    except KeyError:
        skipped+=1
        skipped_list.append(gene)

print("SKIPPED %s"%(skipped))

# Write out a file with all of the regions that were skipped because no match was found 
fout3 = open((plddt_input[:-4]+"_%s_skipped.txt" % pfam_flag),'w')
filtered_skipped = []
for item in skipped_list:
    if item in filtered_skipped:
        pass
    else:
        filtered_skipped.append(item)
for skipped_gene in filtered_skipped:
    fout3.write('%s\n' % skipped_gene)
fout3.close()

