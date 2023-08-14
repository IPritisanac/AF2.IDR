import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import os,sys


plddt_scores = sys.argv[1]      # name of the folder with the final_plddt scores
cnt = 0     # a counter which we will use to plot the pLDDT distribution plot of CheZOD dataset only once (since this dataset is part of each positiveDB_CheZOD_final_pLDDT_scores.txt file)
curr_dir = os.getcwd()
for file in os.listdir(plddt_scores):
    dataset = file[:-30]        #name of the current dataset
    file_path = os.path.join(curr_dir, plddt_scores, file)
    print(file_path)
    curr_dataset = pd.read_csv(file_path, delimiter = ' ')
    x = curr_dataset.loc[curr_dataset['true_pLDDT'] == 1, 'AF_pLDDT']
    kwargs = dict(alpha=0.5, bins=70)
    plt.hist(x, **kwargs, color='r')
    title = dataset + '_distribution_plot'
    plt.gca().set(title = title, ylabel = 'Number of residues', xlabel = 'pLDDT scores')
    plt.xlim(0,100)
    plt.legend()
    plt.savefig(title + ".eps", format = 'eps')
    plt.show()
    cnt += 1
    if cnt == 1:
        x = curr_dataset.loc[curr_dataset['true_pLDDT'] == 0, 'AF_pLDDT']
        kwargs = dict(alpha=0.5, bins=70)
        plt.hist(x, **kwargs, color='r')
        title = 'CheZOD_distribution_plot'
        plt.gca().set(title = title, ylabel = 'Number of residues', xlabel = 'pLDDT scores')
        plt.xlim(0,100)
        plt.legend()
        plt.savefig(title + ".eps", format = 'eps')
        plt.show()
