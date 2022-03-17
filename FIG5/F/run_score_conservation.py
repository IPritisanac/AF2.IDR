import os,sys

inputdir=sys.argv[1]

for infile in os.listdir(inputdir):
    fpath=inputdir+os.sep+infile

    os.system("python score_conservation_IDR.py -s %s -o %s %s"%('property_entropy',fpath[:-3]+"_CONS.out.txt",fpath))
    #outpath=inputdir+"_CONS"+os.sep+infile[:-3]+".txt"

    #print("python score_conservation_IP.py %s -o %s"%(fpath,outpath))
    #os.system("python score_conservation_IP.py -o %s %s"%(outpath,fpath))
