import numpy as np
from matplotlib import pyplot as plt

#anchor = np.loadtxt(r'C:\Users\patri\Desktop\for_Jess_Aplhafold\Uniprot_start_end\ANCHOR2_TESTING_flexible_linkers.txt', skiprows = 1, dtype = 'str', delimiter = '\t')
#DIBS = np.loadtxt(r'C:\Users\patri\Desktop\for_Jess_Aplhafold\Uniprot_start_end\DIBS_start_end_final.txt', skiprows = 1, dtype = 'str', delimiter = '\t')
#DisProt = np.loadtxt(r'C:\Users\patri\Desktop\for_Jess_Aplhafold\NEW_DisProt.txt', skiprows = 1, dtype = 'str', delimiter = ' ')
#DisProt2 = np.loadtxt(r'C:\Users\patri\Desktop\for_Jess_Aplhafold\Uniprot_start_end\DisProt_start_end_final.txt', skiprows = 1, dtype ='str', delimiter = '\t')
#DisProt3 = np.loadtxt(r'C:\Users\patri\Desktop\for_Jess_Aplhafold\IUPRED_NEW_DisProt.txt', skiprows = 1, dtype = 'str', delimiter = ' ')
#FuzDB = np.loadtxt(r'C:\Users\patri\Desktop\for_Jess_Aplhafold\Uniprot_start_end\FuzDB_DOR.txt', skiprows = 1, dtype = 'str', delimiter = '\t')
MFIB = np.loadtxt(r'C:\Users\patri\github\real_for_Jess_Aplhafold\SUPP_FIG_8_A_B_D\DB_Uniprot_start_end\CheZOD_no_overlap.txt', dtype = 'str')

'''
#Percantage of <5 and <10 residues
less_5 = 0
less_10 = 0
all = 0
for i in range(len(anchor)):
    all += 1
    if ((int(anchor[i][2]) - int(anchor[i][1])) +1) <5:
        less_5 += 1
    if ((int(anchor[i][2]) - int(anchor[i][1]) +1 ) >= 5) and ((int(anchor[i][2]) - int(anchor[i][1]) +1 ) <= 10):
        less_10 += 1
perc_5 = (less_5/all)*100
perc_10 = (less_10/all)*100
print('Anchor')
print('percantage <5 = %.1f'%perc_5)
print('percantage >=5 and <=10 = %.1f'%perc_10)

less_5 = 0
less_10 = 0
all = 0
for i in range(len(DIBS)):
    all += 1
    if ((int(DIBS[i][2]) - int(DIBS[i][1])) +1) <5:
        less_5 += 1
    if ((int(DIBS[i][2]) - int(DIBS[i][1]) +1 ) >= 5) and ((int(DIBS[i][2]) - int(DIBS[i][1]) +1 ) <= 10):
        less_10 += 1
perc_5 = (less_5/all)*100
perc_10 = (less_10/all)*100
print('DIBS')
print('percantage <5 = %.1f'%perc_5)
print('percantage >=5 and <=10 = %.1f'%perc_10)

less_5 = 0
less_10 = 0
all = 0
for i in range(len(DisProt2)):
    all += 1
    if ((int(DisProt2[i][2]) - int(DisProt2[i][1])) +1) <5:
        less_5 += 1
    if ((int(DisProt2[i][2]) - int(DisProt2[i][1]) +1 ) >= 5) and ((int(DisProt2[i][2]) - int(DisProt2[i][1]) +1 ) <= 10):
        less_10 += 1
perc_5 = (less_5/all)*100
perc_10 = (less_10/all)*100
print('DisProt')
print('percantage <5 = %.1f'%perc_5)
print('percantage >=5 and <=10 = %.1f'%perc_10)

less_5 = 0
less_10 = 0
all = 0
for i in range(len(FuzDB)):
    all += 1
    if ((int(FuzDB[i][2]) - int(FuzDB[i][1])) +1) <5:
        less_5 += 1
    if ((int(FuzDB[i][2]) - int(FuzDB[i][1]) +1 ) >= 5) and ((int(FuzDB[i][2]) - int(FuzDB[i][1]) +1 ) <= 10):
        less_10 += 1
perc_5 = (less_5/all)*100
perc_10 = (less_10/all)*100
print('FuzDB')
print('percantage <5 = %.1f'%perc_5)
print('percantage >=5 and <=10 = %.1f'%perc_10)

less_5 = 0
less_10 = 0
all = 0
for i in range(len(MFIB)):
    all += 1
    if ((int(MFIB[i][2]) - int(MFIB[i][1])) +1) <5:
        less_5 += 1
    if ((int(MFIB[i][2]) - int(MFIB[i][1]) +1 ) >= 5) and ((int(MFIB[i][2]) - int(MFIB[i][1]) +1 ) <= 10):
        less_10 += 1
perc_5 = (less_5/all)*100
perc_10 = (less_10/all)*100
print('MFIB')
print('percantage <5 = %.1f'%perc_5)
print('percantage >=5 and <=10 = %.1f'%perc_10)





'''
x1 = []
for i in range(len(MFIB)):
    x1.append(int(MFIB[i][2])-int(MFIB[i][1]) + 1)

kwargs = dict(alpha=0.5, bins=100)
plt.hist(x1, **kwargs, color='y')
plt.gca().set(title='Distribution of IDR lengths in Anchor2', ylabel='Number of IDRs', xlabel = 'IDR length')
plt.xlim(0,350)
plt.legend()
plt.show()

'''x1 = []
for i in range(len(DIBS)):
    x1.append(int(DIBS[i][2])-int(DIBS[i][1]) + 1)

kwargs = dict(alpha=0.5, bins=200)
plt.hist(x1, **kwargs, color='y')
plt.gca().set(title='Distribution of IDR lengths in DIBS', ylabel='Number of IDRs', xlabel = 'IDR length')
plt.xlim(0,400)
#plt.legend()
plt.show()

x1 = []
x2 = []
for i in range(len(DisProt2)):
    x1.append(int(DisProt2[i][2])-int(DisProt2[i][1]) + 1)
for i in range(len(DisProt3)):
    x2.append(int(DisProt3[i][2])-int(DisProt3[i][1]) + 1)

kwargs = dict(alpha=0.5, bins=70)
plt.hist(x1, **kwargs, color='b', label = 'Old')
plt.hist(x2, **kwargs, color='r', label = 'IUPRED')
plt.gca().set(title='Distribution of IDR lengths in DisProt', ylabel='Number of IDRs', xlabel = 'IDR length')
plt.xlim(0,100)
plt.legend()
plt.show()


x1 = []
cnt = 0
for i in range(len(FuzDB)):
    x1.append(int(FuzDB[i][2])-int(FuzDB[i][1]) + 1)
    if ((int(FuzDB[i][2]) - int(FuzDB[i][1]) +1) < 6):
        cnt +=1

kwargs = dict(alpha=0.5, bins=50)
plt.hist(x1, **kwargs, color='y')
plt.gca().set(title='Distribution of IDR lengths in FuzDB', ylabel='Number of IDRs', xlabel = 'IDR length')
plt.xlim(0,10)
#plt.legend()
plt.show()
print(cnt)

x1 = []
for i in range(len(MFIB)):
    x1.append(int(MFIB[i][2])-int(MFIB[i][1]) + 1)

kwargs = dict(alpha=0.5, bins=200)
plt.hist(x1, **kwargs, color='y')
plt.gca().set(title='Distribution of IDR lengths in MFIB', ylabel='Number of IDRs', xlabel = 'IDR length')
plt.xlim(0,500)
#plt.legend()
plt.show()
'''
