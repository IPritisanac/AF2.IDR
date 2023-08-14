1. To get the list of pLDDT scores of the disordered regions from the databases in positive dataset (MFIB, FuzDB, DIBS, DisProt, MoRF) and from the negative dataset (CheZOD):

  i) download the .pdb files from AlphaFold database for all proteins based on their Uniprot IDs in the positive and negative datasets, using the download_AFDB.py script
      -- the download_AFDB.py takes as input the name of the folder which contains the .txt files of format "Uniprot start end" belonging to each positive and the negative dataset
      -- the input files can be found in the DB_uniprot_start_end folder
      -- from the command line call: python download_AFDB.py ./path_to_DB_uniprot_start_end
      -- the script creates a separate folder for each dataset (positive and negative); the folders contain successfully downloaded .pdb files of proteins in positive and negative datasets
      -- additionally, as the script runs, it displays Uniprot IDs of proteins which don't exist in AFDB and were therefore not downloaded (the Uniprot IDs are displayed in the terminal)

  ii) extract just the pLDDT scores from all .pdb files in each of the dataset directories and create a single file containing all pLDDT scores, using the read_AFDB_file.py script
      -- this script needs to be called from the command line for each dataset separately:
              - python read_AFDB_file.py ./path_to_CheZOD
              - python read_AFDB_file.py ./path_to_MFIB
              - python read_AFDB_file.py ./path_to_MoRF
              - python read_AFDB_file.py ./path_to_DisProt
              - python read_AFDB_file.py ./path_to_DIBS
              - python read_AFDB_file.py ./path_to_FuzDB
      -- the input to this script are the .pdb files inside the current dataset folder (CheZOD/MFIB/MoRF/DisProt/DIBS/FuzDB)
      -- the output of the script is a file of format "UniprotID    pLDDT_score  AA_type   residue_number" - once we run this script on all datasets, we get a separate file for each of them:
              - CheZOD_pLDDT_scores.txt, MFIB_pLDDT_scores.txt, MoRF_pLDDT_scores.txt, DisProt_pLDDT_scores.txt, DIBS_pLDDT_scores.txt, FuzDB_pLDDT_scores.txt
      -- all the output files are saved in a "pLDDT_scores" directory that is created during the run
      -- since the .pdb files of proteins >2700 residues long need to be fixed, the script displays some additional information in the terminal:
              - the number of residues in the fragment-corrected AFDB
              - the number of proteins in the fragment-corrected AFDB
              - the number of fragments in the uncorrected AFDB
              - the number of duplicated residues in the uncorrected AFDB

  iii) extract only the pLDDT scores of the disordered regions denoted as "start end" in the files from the DB_uniprot_start_end folder, using the final_pLDDT_scores.py script
      -- the script needs to be called once from the command line with: python final_plddt_scores.py ./path_to_pLDDT_scores DB_uniprot_start_end
      -- during the run a "final_pLDDT_scores" directory is created where all the positiveDB_cheZOD_final_pLDDT_scores.txt are saved
      -- the format of the output files is: Uniprot_ID    Residue_number  true_pLDDT    AF_pLDDT", where true_pLDDT is '1' for all positive datasets and '0' for CheZOD dataset, while AF_pLDDT
     is the pLDDT assigned by AlphaFold for each residue
      -- a Union_CheZOD_final_pLDDT_scores.txt is also created, which contains all the positive datasets and a negative dataset
      -- these files are the final output that we can than use for the ROC analysis

2. To get the plot in Figure 7A and the values of ROC analysis, we need the final_pLDDT_scores directory with the positiveDB_CheZOD_final_pLDDT_scores.txt file of each positive dataset. All plots and output files can be found in ROC_analysis_output directory.

  i) With the plot_roc_curve.py we get the ROC curve and a confusion matrix of each positiveDB_CheZOD_final_pLDDT_scores.txt file separately saved as a .eps file. Additionally, this script outputs the values of ROC analysis in a positiveDB_ROC_output.txt. This table was used for Supplementary table 5.
    -- the script needs to be called from the command line once: python plot_roc_curve.py ./path_to_final_pLDDT_scores
  ii) With the plot_overall_curve.py we get the ROC curves of all positiveDB_CheZOD_final_pLDDT_scores.txt files in one plot, like the one in Fig7A.
    -- the script needs to be called once from the command line: python plot_overall_roc_curve.py ./path_to_final_pLDDT_scores
    -- the output is one plot of ROC curves from all datasets saved in an .eps file format
