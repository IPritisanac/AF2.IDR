To get the pLDDT distribution plots for each dataset, we use the distribution_plots.py script
  -- the script takes as input a directory where the final pLDDT .txt files are (final_pLDDT_scores)
  -- the script iterates over each file in the folder and plots the distribution of the pLDDT scores for each positive and the negative dataset
  -- the plots are saved as .eps files (db_distribution_plot.eps)
  -- the script is called from the command line once: python distribution_plots.py ./path_to_final_pLDDT_scores
