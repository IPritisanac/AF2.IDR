import os, sys
import numpy as np
import matplotlib.pyplot as plt

# Load a SPARTA+ prediction file
fin = open('input/5bxv_B_pred.tab', 'r')


# extract the CS information from SPARTA+
cs_dict = {}
cs_list = []
for line in fin:
    new = line.strip().split()
    if len(new) < 2:
        pass
    else:
        if new[0] == 'REMARK':
            pass
        elif new[0] in ('DATA', 'VARS', 'FORMAT'):
            pass
        else:
            cs = new[2]
            if cs in ('N','C','CA','CB','HN'):
                cs_dict[cs] = cs_dict.get(cs, []) + [[int(new[0]), float(new[4]), float(new[8]), float(new[5])]] # res CS err RCshift
            elif cs in ('HA','HA2','HA3'):
                cs_dict['HA'] = cs_dict.get('HA', []) + [[int(new[0]), float(new[4]), float(new[8]), float(new[5])]] # res CS err RCshift

fin.close()


'''# Load 1H and 15N assigned shifts from 16543
bmrb_dict = {}
bmrb = open('bmrb_50253.txt', 'r')
for line in bmrb:
    new = line.strip().split()
    if new[0] == 'sequence':
        pass
    else:
        res = int(new[0])
        bmrb_dict[res] = bmrb_dict.get(res, []) + [float(new[2]), float(new[3])] # CS1 CS2
bmrb.close()  '''     


def sparta_seq(sparta_in, ID):
    ''' extract the amino-acid sequence from a SPARTA+ prediction file '''
    out_seq = ''
    fin = open(sparta_in, 'r')
    for line in fin:
        new = line.strip().split()
        if len(new) < 2:
            pass
        else:
            if new[0] == 'REMARK':
                pass
            elif new[0] == 'DATA' and new[1] == 'SEQUENCE':
                out_seq += ''.join(new[2:])
    fout = open(('%s.seq' % ID), 'w')
    fout.write('%s' % out_seq)
    fout.close()
     
    return out_seq
    

def sparta_to_ssp(sparta_dict, ID):
    ''' convert the chemical shifts in SPARTA+ to SSP format '''
    for key, value in sparta_dict.items():
        fout = open(('%s.%s' % (ID, str(key).lower())), 'w')
        data = np.asarray(value)
        for i,res in enumerate(data[:,0]):
            fout.write('%.0f %.3f\n' % (res, data[i,1]))#-data[i,2]))
        fout.close()

     
# Now run the example for 5bxv     
sparta_seq('input/5bxv_B_pred.tab', 'output/sparta_4ebp1_5bxv')        
sparta_to_ssp(cs_dict, 'output/sparta_4ebp1_5bxv')    


sys.exit()
    
# Load 13CO, 13CA, 13CB assigned chemical shifts from 6968
#carbon_dat = np.genfromtxt('output/output_50253.txt') # phospho 4eb-2
carbon_dat = np.genfromtxt('output/output_19905.txt')
print(carbon_dat)   
# AA res C CA CB N H HA



def bmrb_to_ssp(bmrb_in, id):
    cs_order = ['CO', 'CA', 'CB']
    for i in range(np.shape(bmrb_in)[1] - 5):
        cnt_bmrb = -1
        fout = open(('%s.%s' % (id, cs_order[i].lower())), 'w')
        for j, resval in enumerate(bmrb_in[:,1]):
            cnt_bmrb += 1
            if carbon_dat[j,i+2] > 0.00:
                print(carbon_dat[j,i+2])
                fout.write('%.0f %.3f\n' % (resval, carbon_dat[j,i+2]))
        fout.close()
                
bmrb_to_ssp(carbon_dat, '19905')
sys.exit()
    
    
    
                

fig, ax = plt.subplots(3,2, figsize=(9,9), sharex=True)
ax = ax.ravel()
cnt_plt = -1  
for key, value in cs_dict.items():
    if key  == 'HN':
        cnt_sparta = -1
        cnt_plt += 1
        data = np.asarray(value)
        bmrb_idx = []
        sparta_idx = []
        for res in data[:,0]:
            cnt_sparta += 1
            cnt_bmrb = -1
            for i, resval in enumerate(carbon_dat[:,1]):
                cnt_bmrb += 1
                if res == resval and carbon_dat[i,6] > 0.00:
                    bmrb_idx.append(cnt_bmrb)
                    sparta_idx.append(cnt_sparta)
        
        # SPARTA+
        print(data[:,1][sparta_idx])
        print(np.shape(data[:,1][sparta_idx]))
        
        # BMRB
        #print(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])
        #print(np.shape(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]))
        
        rms_HN = np.sqrt(np.mean((data[:,1][sparta_idx] - carbon_dat[:,6][bmrb_idx])**2))

        # Plot
        '''ax[0].errorbar(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx], data[:,1][sparta_idx], yerr=data[:,2][sparta_idx], fmt='o', mfc='white', elinewidth=0.5, label='RMSD = %.2f ppm' % rms_HN) 
        ax[0].plot([np.min(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])],[np.min(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])], 'k-', label='y=x', linewidth=1) 
        ax[0].legend(loc='upper left')
        ax[0].set_xlabel(r'Measured $^{1}$H (ppm)')
        ax[0].set_ylabel(r'Simulated $^{1}$H (ppm)')
        #plt.xlim(min(np.min(data[:,1][sparta_idx] - data[:,2][sparta_idx]), np.min(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])), max(np.max(data[:,1][sparta_idx] + data[:,2][sparta_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])))
        ax[0].set_ylim(min(np.min(data[:,1][sparta_idx] - data[:,2][sparta_idx]), np.min(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]))-0.1, max(np.max(data[:,1][sparta_idx] + data[:,2][sparta_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]))+0.1)         
        #ax[0].set_xlim(7.8,8.8)
        #ax[0].set_ylim(7.8,8.8)'''
        
        ax[0].bar(data[:,0][sparta_idx], data[:,1][sparta_idx] - data[:,3][sparta_idx], label='simulated')
        ax[0].bar(data[:,0][sparta_idx], carbon_dat[:,6][bmrb_idx] - data[:,3][sparta_idx], label='measured')
        #ax[0].set_xlabel('Residue number')
        ax[0].set_ylabel('secondary 1HN shift (ppm)')
        ax[0].legend(loc='upper left')
        
    elif key  == 'CA':
        cnt_sparta = -1
        cnt_plt += 1
        data = np.asarray(value)
        bmrb_idx = []
        sparta_idx = []
        for res in data[:,0]:
            cnt_sparta += 1
            cnt_bmrb = -1
            for i, resval in enumerate(carbon_dat[:,1]):
                cnt_bmrb += 1
                if res == resval and carbon_dat[i,3] > 0.00:
                    bmrb_idx.append(cnt_bmrb)
                    sparta_idx.append(cnt_sparta)

        
        # SPARTA+
        print(data[:,1][sparta_idx])
        print(np.shape(data[:,1][sparta_idx]))
        print(data[:,0][sparta_idx])
        
        # BMRB
        print(carbon_dat[:,2][bmrb_idx])
        print(carbon_dat[:,1][bmrb_idx])
        
        rms_C = np.sqrt(np.mean((data[:,1][sparta_idx] - carbon_dat[:,3][bmrb_idx])**2))

        # Plot
        #ax[2].errorbar(carbon_dat[:,3][bmrb_idx], data[:,1][sparta_idx], yerr=data[:,2][sparta_idx], fmt='o', mfc='white', elinewidth=0.5, label='RMSD = %.2f ppm' % rms_C) 
        #ax[2].plot([np.min(carbon_dat[:,3][bmrb_idx]), np.max(carbon_dat[:,3][bmrb_idx])], [np.min(carbon_dat[:,3][bmrb_idx]), np.max(carbon_dat[:,3][bmrb_idx])], 'k-', label='y=x', linewidth=1) 
        #ax[2].legend(loc='upper left')
        #ax[2].set_xlabel(r'Measured $^{13}$CO (ppm)')
        #ax[2].set_ylabel(r'Simulated $^{13}$CO (ppm)')
        
        ax[2].bar(data[:,0][sparta_idx], data[:,1][sparta_idx]-data[:,3][sparta_idx], label='simulated')
        ax[2].bar(data[:,0][sparta_idx], carbon_dat[:,3][bmrb_idx]-data[:,3][sparta_idx], label='measured')
        ax[2].set_xlabel('Residue number')
        ax[2].set_ylabel('secondary 13CA shift (ppm)')
        ax[2].legend(loc='upper left')
        
        
    elif key  == 'C':
        cnt_sparta = -1
        cnt_plt += 1
        data = np.asarray(value)
        bmrb_idx = []
        sparta_idx = []
        for res in data[:,0]:
            cnt_sparta += 1
            cnt_bmrb = -1
            for i, resval in enumerate(carbon_dat[:,1]):
                cnt_bmrb += 1
                if res == resval and carbon_dat[i,2] > 0.00:
                    bmrb_idx.append(cnt_bmrb)
                    sparta_idx.append(cnt_sparta)

        
        # SPARTA+
        print(data[:,1][sparta_idx])
        print(np.shape(data[:,1][sparta_idx]))
        print(data[:,0][sparta_idx])
        
        # BMRB
        print(carbon_dat[:,2][bmrb_idx])
        print(carbon_dat[:,1][bmrb_idx])
        
        rms_C = np.sqrt(np.mean((data[:,1][sparta_idx] - carbon_dat[:,2][bmrb_idx])**2))

        # Plot
        #ax[2].errorbar(carbon_dat[:,3][bmrb_idx], data[:,1][sparta_idx], yerr=data[:,2][sparta_idx], fmt='o', mfc='white', elinewidth=0.5, label='RMSD = %.2f ppm' % rms_C) 
        #ax[2].plot([np.min(carbon_dat[:,3][bmrb_idx]), np.max(carbon_dat[:,3][bmrb_idx])], [np.min(carbon_dat[:,3][bmrb_idx]), np.max(carbon_dat[:,3][bmrb_idx])], 'k-', label='y=x', linewidth=1) 
        #ax[2].legend(loc='upper left')
        #ax[2].set_xlabel(r'Measured $^{13}$CO (ppm)')
        #ax[2].set_ylabel(r'Simulated $^{13}$CO (ppm)')
        
        ax[3].bar(data[:,0][sparta_idx], data[:,1][sparta_idx]-data[:,3][sparta_idx], label='simulated')
        ax[3].bar(data[:,0][sparta_idx[1:]], carbon_dat[:,2][bmrb_idx[1:]]-data[:,3][sparta_idx[1:]], label='measured') # avoid Met1
        ax[3].set_xlabel('Residue number')
        ax[3].set_ylabel('secondary 13CO shift (ppm)')
        ax[3].legend(loc='upper left')

    elif key  == 'N':
        cnt_sparta = -1
        cnt_plt += 1
        data = np.asarray(value)
        bmrb_idx = []
        sparta_idx = []
        for res in data[:,0]:
            cnt_sparta += 1
            cnt_bmrb = -1
            for i, resval in enumerate(carbon_dat[:,1]):
                cnt_bmrb += 1
                if res == resval and carbon_dat[i,5] > 0.00:
                    bmrb_idx.append(cnt_bmrb)
                    sparta_idx.append(cnt_sparta)

        
        # SPARTA+
        print(data[:,1][sparta_idx])
        print(np.shape(data[:,1][sparta_idx]))
        
        # BMRB
        #print(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx])
        #print(np.shape(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx]))
        
        rms_N = np.sqrt(np.mean((data[:,1][sparta_idx] - carbon_dat[:,5][bmrb_idx])**2))

        # Plot
        '''ax[1].errorbar(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx], data[:,1][sparta_idx], yerr=data[:,2][sparta_idx], fmt='o', mfc='white', elinewidth=0.5, label='RMSD = %.2f ppm' % rms_N) 
        ax[1].plot([np.min(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx])],[np.min(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx])], 'k-', label='y=x', linewidth=1) 
        ax[1].legend(loc='upper left')
        ax[1].set_xlabel(r'Measured $^{15}$N (ppm)')
        ax[1].set_ylabel(r'Simulated $^{15}$N (ppm)')'''
        #plt.xlim(min(np.min(data[:,1][sparta_idx] - data[:,2][sparta_idx]), np.min(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])), max(np.max(data[:,1][sparta_idx] + data[:,2][sparta_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])))
        #ax[1].set_ylim(min(np.min(data[:,1][sparta_idx] - data[:,2][sparta_idx]), np.min(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]))-0.1, max(np.max(data[:,1][sparta_idx] + data[:,2][sparta_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]))+0.1)                
        
        ax[1].bar(data[:,0][sparta_idx], data[:,1][sparta_idx] - data[:,3][sparta_idx], label='simulated')
        ax[1].bar(data[:,0][sparta_idx], carbon_dat[:,5][bmrb_idx] - data[:,3][sparta_idx], label='measured')
        #ax[1].set_xlabel('Residue number')
        ax[1].set_ylabel('secondary 15N shift (ppm)')
        ax[1].legend(loc='upper left')
        
        
    elif key  == 'CB':
        cnt_sparta = -1
        cnt_plt += 1
        data = np.asarray(value)
        bmrb_idx = []
        sparta_idx = []
        for res in data[:,0]:
            cnt_sparta += 1
            cnt_bmrb = -1
            for i, resval in enumerate(carbon_dat[:,1]):
                cnt_bmrb += 1
                if res == resval and carbon_dat[i,4] > 0.00:
                    bmrb_idx.append(cnt_bmrb)
                    sparta_idx.append(cnt_sparta)

        
        # SPARTA+
        print(data[:,1][sparta_idx])
        print(np.shape(data[:,1][sparta_idx]))
        
        # BMRB
        #print(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx])
        #print(np.shape(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx]))
        
        rms_N = np.sqrt(np.mean((data[:,1][sparta_idx] - carbon_dat[:,4][bmrb_idx])**2))

        # Plot
        '''ax[1].errorbar(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx], data[:,1][sparta_idx], yerr=data[:,2][sparta_idx], fmt='o', mfc='white', elinewidth=0.5, label='RMSD = %.2f ppm' % rms_N) 
        ax[1].plot([np.min(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx])],[np.min(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx])], 'k-', label='y=x', linewidth=1) 
        ax[1].legend(loc='upper left')
        ax[1].set_xlabel(r'Measured $^{15}$N (ppm)')
        ax[1].set_ylabel(r'Simulated $^{15}$N (ppm)')'''
        #plt.xlim(min(np.min(data[:,1][sparta_idx] - data[:,2][sparta_idx]), np.min(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])), max(np.max(data[:,1][sparta_idx] + data[:,2][sparta_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])))
        #ax[1].set_ylim(min(np.min(data[:,1][sparta_idx] - data[:,2][sparta_idx]), np.min(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]))-0.1, max(np.max(data[:,1][sparta_idx] + data[:,2][sparta_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]))+0.1)                
        
        ax[4].bar(data[:,0][sparta_idx], data[:,1][sparta_idx] - data[:,3][sparta_idx], label='simulated')
        ax[4].bar(data[:,0][sparta_idx], carbon_dat[:,4][bmrb_idx] - data[:,3][sparta_idx], label='measured')
        #ax[1].set_xlabel('Residue number')
        ax[4].set_ylabel('secondary 13CB shift (ppm)')
        ax[4].legend(loc='upper left')
        
        
        
plt.tight_layout()        
plt.show()
#plt.savefig('aSyn_secondary_shifts.pdf')
#plt.clf()



sys.exit()
    
