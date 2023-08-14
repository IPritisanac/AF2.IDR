import numpy as np
from matplotlib import pyplot as plt
import os,sys


if __name__ == '__main__':
    ranges = sys.argv[1]        # name of the directory with start-end residues (DB_Uniprot_start_end)
    curr_dir = os.getcwd()
    cnt = 0
    for file in os.listdir(ranges):
        if file[:-15] != 'MFIB':
            rng = []
            curr_dataset_path = os.path.join(curr_dir, ranges, file)
            curr_dataset = np.loadtxt(curr_dataset_path, dtype = 'str')
            for i in range(len(curr_dataset)):
                rng.append(int(curr_dataset[i][2])-int(curr_dataset[i][1]) + 1)
                if (int(curr_dataset[i][2])-int(curr_dataset[i][1]) + 1) > 250:
                    print(curr_dataset[i])
            colors = ['mediumturquoise', 'lightsalmon']
            print('mean %s = %.3f'%(file[:-15], np.mean(rng)))
            print('median %s = %.3f'%(file[:-15], np.median(rng)))
            if cnt == 0:
                kwargs = dict(alpha=0.7, bins=20)
                plt.hist(rng, **kwargs, color=colors[cnt], label = file[:-15])
            else:
                kwargs = dict(alpha=0.4, bins=50)
                plt.hist(rng, **kwargs, color=colors[cnt], label = file[:-15])
            cnt += 1
    plt.gca().set(title='Distribution of  IDR length', ylabel='Number of IDRs', xlabel = 'IDR length')
    plt.xlim(0,100)
    plt.legend()
    plt.savefig('CheZOD_Anchor_flex_link_size_distribution.pdf', format = 'pdf')
    plt.show()
