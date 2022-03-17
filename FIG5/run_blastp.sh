#!/bin/bash

dirs=`find . -maxdepth 1 -mindepth 1 -type d`
path_db=/neuhaus/iva/ncbi-blast-2.8.1+/bin/pdb_7046_under40_9858.fasta

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
