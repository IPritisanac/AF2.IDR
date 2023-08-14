from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import os, sys



if __name__ == '__main__':
    plddt = sys.argv[1] # name of the directory where final_pLDDT scores are
    curr_dir = os.getcwd()
    kwargs = dict(alpha = 0.5, bins = 'auto')
    colors = ['mediumturquoise','lightsalmon']
    cnt = 0
    for file in os.listdir(plddt):
        database = file[5:-23]
        print(database)
        file_path = os.path.join(curr_dir, plddt, file)
        curr_data = pd.read_csv(file_path, sep = ' ')
        x1 = curr_data.loc[curr_data['true_pLDDT'] == 0, 'AF_pLDDT']
        plt.hist(x1, **kwargs, color=colors[cnt], label=database)
        cnt += 1
    plt.gca().set(title='Distribution of pLDDT scores', ylabel='Number of residues')
    plt.xlim(0,100)
    plt.legend()
    plt.savefig('CheZOD_Anchor_flex_link_pLDDT_distribution.eps', format = 'eps')
    plt.show()
