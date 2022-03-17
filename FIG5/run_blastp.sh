#!/bin/bash

dirs=`find . -maxdepth 1 -mindepth 1 -type d`
path_db=/neuhaus/iva/ncbi-blast-2.8.1+/bin/cullpdb_pc30.0_res0.0-3.5_noBrks_noDsdr_len40-10000_R1.0_Xray+Nmr+EM_d2021_12_07_chains7046.fasta

for dir in $dirs
do
  seqs=`find $dir -maxdepth 1 -type f`
  for s in $seqs
  do
    echo $s
    /neuhaus/iva/ncbi-blast-2.8.1+/bin/blastp -query $s -db $path_db -max_hsps 1 -evalue 1e-6 -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bi
tscore qcovs qcovhsp' -out $s.blastp.out.txt -num_threads 4
  done
done
