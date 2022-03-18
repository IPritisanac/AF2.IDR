import os, sys
import numpy as np
import matplotlib.pyplot as plt


############### TRANSLATIONAL DIFFUSION ########################
    
    
# now simulate diffusion
exp = 5.71e-11 # m^2/s
exp_err = 0.2e-11 # m^2/s 
calc = 5.141e-11 # m^2/s

# D = (kb*T)/(6pi*nu*Rh)


def calc_diffusion(temp,Rh):
    kb = 1.38064852e-23 # m^2 kg s^-2 K^-1
    T = temp+273.15 # K
    Rh = Rh*1e-10 # m
    nu = 1136.6e-6 # uPa/s @ 15C Pa = kg m-1 s-2
    return (kb*T)/(6*np.pi*nu*Rh)
    
# diffusion of dioxane in 0.1 mM aSyn @ 15C, pH 6.4 = ~7.9e-10 m^2/s
# Rh of dioxane = 2.12 A --> (7.9e-10/5.71e-11)*2.12 = 29.3 Angstroms
# calculated diffusion from the equation above yields 6.36e-11 m^2/s diffusion constant --> 26.3 Angstrom Rh
    
print(calc_diffusion(15,29.2))

def calc_intensity(D):
    # Ij = I0*exp([-gamma^2*Gj^2*delta^2]*(Delta-delta/3-tau/2)*D)
    # max Grad = 66.8 G/cm --> % * maxGrad = G (in Tesla/m) --> 66.8 G/cm = 0.668 T/m
    # delays from 200 to 400 ms --> Delta
    # 2 x 1.5 ms for encode/decode grads --> delta = 3 ms
    I0 = 5000 # intensity, AU
    gammaH = 267.52218744e6 # gyromagnetic ratio, rad s-1 T-1
    delta = 3e-3 # total encode/decode grad duration, s
    tau = 0.2e-3 # gradient recovery delay, s
    Delta = 200e-3 + delta + 2*tau # delay time, s
    G = 0.668 # max gradient strenght, T m^-1
    
    grad_array = np.linspace(0,1,21)*G
    #print(grad_array)
    #print(-gammaH**2 * grad_array**2 * delta**2)
    #print((Delta-delta/3-tau/2))
    
    Ij = I0*np.exp((-gammaH**2 * grad_array**2 * delta**2) * (Delta-delta/3-tau/2)*D)
    
    return grad_array/G, Ij/I0
    
fig, ax = plt.subplots(1,2, figsize=(9,4))
ax = ax.ravel()
    
ax[0].plot(calc_intensity(exp)[0], calc_intensity(exp)[1], 'o-', mfc='white', label=r'measured')#, $D$ = 5.71 $\pm$ 0.02 $x$ 10$^{-11}$ m$^{2}$ $s^{-1}$ ')
ax[0].plot(calc_intensity(calc)[0], calc_intensity(calc)[1], '-', mfc='orange', mec='white', label='simulated')#, $D$ = 5.86 $x$ 10$^{-11}$ m$^{2}$ $s^{-1}$ ')
ax[0].legend(loc='upper right')
ax[0].set_xlabel(r'$G^{2}$/$G_{max}^{2}$', fontsize=12)
ax[0].set_ylabel(r'$I$/$I_{0}$', fontsize=12)

ax[1].errorbar(calc_intensity(exp)[0]**2, np.log(calc_intensity(exp)[1]), yerr = (exp_err/exp)*np.log(calc_intensity(exp)[1]), fmt='o-', mfc='white', label=r'measured')
ax[1].plot(calc_intensity(calc)[0]**2, np.log(calc_intensity(calc)[1]),'o-', mfc='orange', mec='white', label='simulated')
ax[1].set_ylabel('ln (I/I0)', fontsize=12)
ax[0].set_xlabel(r'$G^{2}$/$G_{max}^{2}$', fontsize=12)

plt.tight_layout(w_pad=1, h_pad=2) # pad=0.4, 
#plt.savefig('aSyn_translational_diffusion.pdf')
plt.show()

#sys.exit()
'''
ax[1].errorbar(calc_intensity(exp)[1], calc_intensity(exp)[1], xerr=0.02/5.71*calc_intensity(exp)[1], fmt='o-')
ax[1].plot([np.min(calc_intensity(exp)[1]), np.max(calc_intensity(exp)[1])],[np.min(calc_intensity(exp)[1]), np.max(calc_intensity(exp)[1])], 'k-', label='y=x')
ax[1].set_xlabel(r'Measured I/I0')
ax[1].set_ylabel(r'Simulated I/I0')
ax[1].legend(loc='upper left')

plt.tight_layout(w_pad=1, h_pad=2) # pad=0.4, 
#plt.savefig('aSyn_translational_diffusion.pdf')
plt.show()


sys.exit()'''

fin = open('pred.tab', 'r')

cs_dict = {}

# first get CS information from SPARTA+
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

print(cs_dict)


def sparta_seq(sparta_in, ID):
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
    
sparta_seq('pred.tab', 'sparta_aSyn_AF')


def sparta_to_ssp(sparta_dict, ID):
    for key, value in sparta_dict.items():
        fout = open(('%s.%s' % (ID, str(key).lower())), 'w')
        data = np.asarray(value)
        for i,res in enumerate(data[:,0]):
            fout.write('%.0f %.3f\n' % (res, data[i,1]-data[i,2]))
        fout.close()

        
        
sparta_to_ssp(cs_dict, 'sparta_aSyn_A')        
#sys.exit()


# Load 1H and 15N assigned shifts from 16543
bmrb_dict = {}
bmrb = open('16543_bmrb.tab', 'r')
for line in bmrb:
    new = line.strip().split()
    if new[0] == 'sequence':
        pass
    else:
        res = int(new[0])
        bmrb_dict[res] = bmrb_dict.get(res, []) + [float(new[2]), float(new[3])] # CS1 CS2
bmrb.close()        

    
# Load 13CO, 13CA, 13CB assigned chemical shifts from 6968
carbon_dat = np.genfromtxt('output_6968.txt')
print(carbon_dat)
print(np.shape(carbon_dat))

    
    
def bmrb_to_ssp(bmrb_in, id):
    cs_order = ['CA', 'CB'] # no co in 5744
    for i in range(np.shape(bmrb_in)[1] - 3):
        cnt_bmrb = -1
        fout = open(('%s.%s' % (id, cs_order[i].lower())), 'w')
        for j, resval in enumerate(bmrb_in[:,1]):
            cnt_bmrb += 1
            if carbon_dat[j,i+3] > 0.00:
                print(carbon_dat[j,i+3])
                fout.write('%.0f %.3f\n' % (resval, carbon_dat[j,i+3]))
        fout.close()
                
bmrb_to_ssp(carbon_dat, '5744')
#sys.exit()
                

fig, ax = plt.subplots(2,2, figsize=(9,9), sharex=True)
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
            for key2, value2 in bmrb_dict.items():
                cnt_bmrb += 1
                if res == key2:
                    bmrb_idx.append(cnt_bmrb)
                    sparta_idx.append(cnt_sparta)
        
        # SPARTA+
        print(data[:,1][sparta_idx])
        print(np.shape(data[:,1][sparta_idx]))
        
        # BMRB
        print(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])
        print(np.shape(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]))
        
        rms_HN = np.sqrt(np.mean((data[:,1][sparta_idx] - np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])**2))

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
        ax[0].bar(data[:,0][sparta_idx], np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx] - data[:,3][sparta_idx], label='measured')
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
        ax[2].bar(data[:,0][sparta_idx], carbon_dat[:,3][bmrb_idx]-0.15 - data[:,3][sparta_idx], label='measured')
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
            for key2, value2 in bmrb_dict.items():
                cnt_bmrb += 1
                if res == key2:
                    bmrb_idx.append(cnt_bmrb)
                    sparta_idx.append(cnt_sparta)

        
        # SPARTA+
        print(data[:,1][sparta_idx])
        print(np.shape(data[:,1][sparta_idx]))
        
        # BMRB
        print(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx])
        print(np.shape(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx]))
        
        rms_N = np.sqrt(np.mean((data[:,1][sparta_idx] - np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx])**2))

        # Plot
        '''ax[1].errorbar(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx], data[:,1][sparta_idx], yerr=data[:,2][sparta_idx], fmt='o', mfc='white', elinewidth=0.5, label='RMSD = %.2f ppm' % rms_N) 
        ax[1].plot([np.min(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx])],[np.min(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx])], 'k-', label='y=x', linewidth=1) 
        ax[1].legend(loc='upper left')
        ax[1].set_xlabel(r'Measured $^{15}$N (ppm)')
        ax[1].set_ylabel(r'Simulated $^{15}$N (ppm)')'''
        #plt.xlim(min(np.min(data[:,1][sparta_idx] - data[:,2][sparta_idx]), np.min(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])), max(np.max(data[:,1][sparta_idx] + data[:,2][sparta_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx])))
        #ax[1].set_ylim(min(np.min(data[:,1][sparta_idx] - data[:,2][sparta_idx]), np.min(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]))-0.1, max(np.max(data[:,1][sparta_idx] + data[:,2][sparta_idx]), np.max(np.asarray(list(bmrb_dict.values()))[:,0][bmrb_idx]))+0.1)                
        
        ax[1].bar(data[:,0][sparta_idx], data[:,1][sparta_idx] - data[:,3][sparta_idx], label='simulated')
        ax[1].bar(data[:,0][sparta_idx], np.asarray(list(bmrb_dict.values()))[:,1][bmrb_idx]-0.5 - data[:,3][sparta_idx], label='measured')
        #ax[1].set_xlabel('Residue number')
        ax[1].set_ylabel('secondary 15N shift (ppm)')
        ax[1].legend(loc='upper left')
        
        
        
plt.tight_layout()        
#plt.show()
plt.savefig('aSyn_secondary_shifts_-0.15CA.pdf')
#plt.clf()



#sys.exit()
    
    

