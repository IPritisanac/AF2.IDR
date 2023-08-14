1. Figure 8A ==> distribution of pLDDT scores of CheZOD and Anchor flexible linkers database
   Figure 8D ==> ROC curve of MFIB as positive dataset and either Chezod or Anchor flexible linkers as negative dataset
    -- first we need to extract the pLDDT scores of Anchor flexible linkers, Chezod and MFIB datasets
    -- to get the MFIB_CheZOD_final_pLDDT_scores.txt, and Anchor_flex_link_CheZOD_final_pLDDT_scores.txt we need to:

        i) download the .pdb files from AlphaFold database for all proteins based on their Uniprot IDs in the CheZOD, MFIB and Anchor flexible linkers datasets, using the download_AFDB.py script
          -- the download_AFDB.py takes as input the name of the folder which contains the .txt files of format "Uniprot start end" belonging to CheZOD, MFIB and the ANchor flexible linkers datasets
          -- the input files can be found in the DB_uniprot_start_end folder
          -- from the command line call: python download_AFDB.py ./path_to_DB_uniprot_start_end
          -- the script creates a separate folder for each dataset; the folders contain successfully downloaded .pdb files of proteins in both datasets
          -- additionally, as the script runs, it displays Uniprot IDs of proteins which don't exist in AFDB and were therefore not downloaded (the Uniprot IDs are displayed in the terminal)

      ii) extract just the pLDDT scores from all .pdb files in all dataset directories and create a single file containing all pLDDT scores from each dataset, using the read_AFDB_file.py script
          -- this script needs to be called from the command line for each dataset separately:
                  - python read_AFDB_file.py ./path_to_CheZOD
                  - python read_AFDB_file.py ./path_to_MFIB
                  - python read_AFDB_file.py ./path_to_Anchor_flex_link
          -- the input to this script are the .pdb files inside the current dataset folder (CheZOD/MFIB/Anchor flexible linkers)
          -- the output of the script is a file of format "UniprotID    pLDDT_score  AA_type   residue_number" - once we run this script on both datasets, we get a separate file for each of them:
                  - CheZOD_pLDDT_scores.txt, MFIB_pLDDT_scores.txt, Anchor_flex_link_pLDDT_scores.txt
          -- all the output files are saved in a "pLDDT_scores" directory that is created during the run
          -- since the .pdb files of proteins >2700 residues long need to be fixed, the script displays some additional information in the terminal:
                  - the number of residues in the fragment-corrected AFDB
                  - the number of proteins in the fragment-corrected AFDB
                  - the number of fragments in the uncorrected AFDB
                  - the number of duplicated residues in the uncorrected AFDB

      iii) extract only the pLDDT scores of the disordered regions denoted as "start end" in the files from the DB_uniprot_start_end folder, using the final_pLDDT_scores.py script
          -- the script needs to be called once from the command line with: python final_plddt_scores.py ./path_to_pLDDT_scores DB_uniprot_start_end
          -- during the run a "final_pLDDT_scores" directory is created where all the positiveDB_cheZOD_final_pLDDT_scores.txt are saved
          -- the format of the output files is: Uniprot_ID    Residue_number  true_pLDDT    AF_pLDDT", where true_pLDDT is '1' MFIB and '0' for CheZOD and Anchor_flex_link dataset, while AF_pLDDT
         is the pLDDT assigned by AlphaFold for each residue
         -- in the end we have 2 separate files: MFIB_CheZOD_final_pLDDT_scores.txt where MFIB is the positive dataset and CheZOD the negative dataset; and the MFIB_Anchor_flex_link_final_pLDDT_scores.txt where
         MFIB is also the positive dataset, but the negative dataset is Anchor flexible linkers dataset
          -- these files are the final output that we can than use for the ROC analysis (Fig 8D) and for the distribution plot (Supp Fig 8A)

2. For the distribution plot in Figure 8A:
    -- from the command line call once python plddt_distribution.py ./path_to_final_pLDDT_scores
    -- the input are the MFIB_CheZOD_final_pLDDT_scores.txt and the MFIB_Anchor_flex_link_final_pLDDT_scores.txt
    -- the output is the pLDDT distribution plot of CheZOD and Anchor flexible linkers datasets saved as the CheZOD_Anchor_flex_link_pLDDT_distribution.eps

3. For the ROC plot in Figure 8D:
    -- from the command line call once python plot_overal_roc_curve.py ./path_to_final_pLDDT_scores
    -- the output is saved as Overall_ROC.eps

4. To get the size distribution of IDRs in CheZOD and Anchor flexible linkers dataset:
    -- the size_distribution.py needs to be called once from the command line as python size_distribution.py ./path_to_DB_uniprot_start_end
    -- the plot of IDR size distribution of CheZOD and ANchor flexible linkers is saved as CheZOD_Anchor_flex_link_size_distribution.eps
