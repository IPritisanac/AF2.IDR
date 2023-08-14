import os,sys
import matplotlib.pyplot as plt
import numpy as np
import copy

## read in probabilities from file
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

def binary_roc_curve(probs,test_labels,label,figpath, out_file, output_dir, tag="_real_"):
    from sklearn.metrics import roc_curve
    from sklearn.metrics import roc_auc_score
    from sklearn.metrics import precision_score
    from sklearn.metrics import recall_score
    from sklearn.metrics import confusion_matrix
    # keep probabilities for the positive outcome only
    probabilities=probs
    labels=test_labels
    # calculate AUC
    try:
        auc = roc_auc_score(labels, probabilities)
    except ValueError:
        print("VALUE ERROR FOR %s"%(label))
        print(labels)
        print(probabilities)
        return 0
    # calculate roc curve
    fpr, tpr, thresholds = roc_curve(labels, probabilities)
    gmeans = np.sqrt(tpr * (1-fpr))
    ix = np.argmax(gmeans)
    opt_threshold=thresholds[ix]
    g_mean=gmeans[ix]
    precision_probs=[]
    for probab in probs:
        if probab>=opt_threshold:
            precision_probs.append(1)
        else:
            precision_probs.append(0)
    ppv=precision_score(labels,precision_probs)
    recall=recall_score(labels,precision_probs)
    npositives=np.sum(labels)

    # compute confusion matrix
    cm=confusion_matrix(labels,precision_probs)

    tn, fp, fn, tp = confusion_matrix(labels,precision_probs).ravel()
    specificity = tn / (tn+fp)
    accuracy = (tp+tn)/(tp+tn+fp+fn)

    label=label.split(os.sep)[-1]
    print('%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.1f' % (label,auc,specificity,ppv,recall,accuracy,opt_threshold,g_mean,npositives))
    out_file.write('%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.1f\n' % (label,auc,specificity,ppv,recall,accuracy,opt_threshold,g_mean,npositives))
    fig, axs = plt.subplots(1, 2)

    axs[0].plot([0, 1], [0, 1], linestyle='--')
    fig.suptitle('%s AUC: %.3f, RECALL: %.3f\nSPECIFICITY: %.3f, PRECISION: %.3f' % (label,auc,recall,specificity,ppv))
    axs[0].scatter(fpr[ix], tpr[ix], color='black', label='Best')
    axs[0].set_xlabel("False positive rate")
    axs[0].set_ylabel("True positive rate")
    axs[0].plot(fpr, tpr, marker='.', markersize = 1)

    """
    Print and plot the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    cmap=plt.cm.Blues
    normalize=True
    title=None
    from sklearn.metrics import confusion_matrix
    if not title:
        if normalize:
            title = 'Normalized CM'
        else:
            title = 'CM no normalization'
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    else:
        pass
    im = axs[1].imshow(cm, interpolation='nearest', cmap=cmap)
    axs[1].figure.colorbar(im, ax=axs[1])
    im.set_clim(0,1)
    axs[1].set(xticks=np.arange(cm.shape[1]),
            yticks=np.arange(cm.shape[0]),
            title=title,
            ylabel='True label',
            xlabel='Predicted label')
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 1.5
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            axs[1].text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    plt_file = os.path.join(output_dir, label)
    plt.savefig(plt_file+".eps",format='eps')
    plt.show()


if __name__=="__main__":
    prob_dir=sys.argv[1] # directory with all _final_plddt.txt files
    curr_dir = os.getcwd()
    print("GO_LABEL\tAUC\tSPECIFICITY\tPRECISION(PPV)\tRECALL\tACCURACY\tOPTIMAL THRESHOLD\tG-MEAN\tN_PROTEINS_TRAINING")
    if not os.path.exists(curr_dir+'/ROC_analysis_output'):
        os.makedirs(curr_dir+'/ROC_analysis_output')            # create an output directory
    os.chdir(curr_dir+'/ROC_analysis_output')                   # change into the output directory
    output_dir = os.getcwd()                        # get output directory path
    os.chdir(curr_dir)
    out_file = os.path.join(output_dir, 'ROC_output.txt')
    roc_out = open(out_file, "w")
    roc_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('Dataset','AUC','Specificity','PPV','Recall','Accuracy','Optimal_threshold','G_mean','N_IDRs'))
    for f in os.listdir(prob_dir):
        if f.endswith(".txt"):
            database_name = f[:-23]
            file_probs=prob_dir+os.sep+f
            file_labels=""
            ordg,probp,labs,genes=read_in_probs(file_probs)
            binary_roc_curve(probp,labs,ordg,'/'.join(file_probs.split(os.sep)[:-1]), roc_out, output_dir)
