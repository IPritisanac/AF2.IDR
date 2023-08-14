import os,sys
import matplotlib.pyplot as plt
import numpy as np
import copy

def read_in_probs(infile):
    fin=open(infile,"r")
    name=infile[:-4]
    order=[]
    known_probs=[]
    given_probs=[]
    gene_order=[]
    cnt=0
    for line in fin:
        stripped=line.strip()
        splitted=stripped.split()
        cnt+=1
        if cnt==1:
            continue
        given_probs.append(float(splitted[3])/100.)
        known_probs.append(float(splitted[2]))
        gene=splitted[0].replace('"','')
        gene_order.append(gene+"_"+splitted[1])

    prob_array=np.array(given_probs,dtype=float)
    labels_array=np.array(known_probs,dtype=float)
    gene_res_order=np.array(gene_order)

    return name,prob_array,labels_array,gene_res_order


def binary_roc_curve(probs,test_labels,label,figpath,tag="_real_"):

    from sklearn.metrics import roc_curve
    from sklearn.metrics import roc_auc_score
    from sklearn.metrics import precision_score
    from sklearn.metrics import recall_score
    from sklearn.metrics import confusion_matrix

    fpr_all = []
    tpr_all = []
    tresholds_all = []
    gmeans_all = []
    ix_all = []
    for i in range(len(probs)):
        probabilities=probs[i]
        labels=test_labels[i]

        try:
            auc = roc_auc_score(labels, probabilities)
        except ValueError:
            print("VALUE ERROR FOR %s"%(label))
            print(labels)
            print(probabilities)
            return 0

        fpr, tpr, thresholds = roc_curve(labels, probabilities)
        fpr_all.append(fpr)
        tpr_all.append(tpr)
        gmeans = np.sqrt(tpr * (1-fpr))
        gmeans_all.append(gmeans)
        ix = np.argmax(gmeans)
        ix_all.append(ix)

    fig, axs = plt.subplots(1, 1)
    colors = ['lightsteelblue', 'lightcoral']
    labels = ['MFIB_Anchor_flex_link', 'MFIB_CheZOD']
    for i in range(len(probs)):
        axs.plot([0, 1], [0, 1], linestyle='--')
        fpr = fpr_all[i]
        ix = ix_all[i]
        tpr = tpr_all[i]
        axs.scatter(fpr[ix], tpr[ix], marker='o', s = 22, color='black')
        axs.set_xlabel("False positive rate")
        axs.set_ylabel("True positive rate")
        axs.plot(fpr_all[i], tpr_all[i], linestyle='-', color = colors[i], label = labels[i]) #  marker='.', markersize = 4.5, 

    fig.tight_layout()
    plt.legend()
    plt.savefig('Overall_ROC'+".eps",format='eps')
    plt.show()


if __name__=="__main__":
    prob_dir=sys.argv[1]        # directory with all _final_plddt.txt files
    ordg_all = []
    probp_all = []
    labs_all = []
    genes_all = []
    for f in os.listdir(prob_dir):
        if f.endswith(".txt"):
            file_probs=prob_dir+os.sep+f
            file_labels=""
            ordg,probp,labs,genes=read_in_probs(file_probs)
            ordg_all.append(ordg)
            probp_all.append(probp)
            labs_all.append(labs)
    binary_roc_curve(probp_all,labs_all,ordg_all,'/'.join(file_probs.split(os.sep)[:-1]))
