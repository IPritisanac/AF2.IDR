import numpy as np
import os, sys

class BMRB(object):
    """ Scripts to convert NMR-STAR files into Sparky peak lists and TALOS input files """

    def __init__(self, in_file):
        """
        Parameters
        ----------
        in_file	: str
            the input NMR-STAR file containing chemical shifts of interest
        """
        self.input = in_file  # initialize to use later for reading input

    def version(self):
        """
        Determines which version of the NMR-STAR file is used
        """
        cnt = -1
        with open(self.input, 'r+') as f:
            lines = f.readlines()
            for i in range(0, len(lines)):
                line = lines[i].strip().split()
                if len(line)>0 and line[0] == '_Entry.NMR_STAR_version':
                    self.version = line[1]

        f.close()

        return self.version

    def sequence(self):
        """
        Extracts the amino acid sequence of a protein from an NMR-STAR file

        Parameters
        ----------
        None (depends on in_file from __init__)

        Returns
        ----------
        self.aa_seq : array
            the amino acid sequence of the protein in an array format
        """

        aa_seq = []
        cnt = -1
        with open(self.input, 'r+') as f:
            lines = f.readlines()
            for i in range(0, len(lines)):
                line = lines[i].strip().split()
                if len(line)>0 and line[0] == '_Entity.Polymer_seq_one_letter_code':
                    aa_seq.append(lines[i+2].strip().split()[0])
                    x=1
                    while len(lines[i+2+x].strip().split()) > 0 and lines[i+2+x].strip().split()[0] != ';':
                        aa_seq.append(lines[i+2+x].strip().split()[0])
                        x+=1
                    self.aa_seq = ''.join(map(str, aa_seq))
                else:
                    pass

        f.close()

        return self.aa_seq


    def aa_to_number(self, aa_3):
        """
        Converts a three-letter amino acid code into a float

        Parameters
        ----------
        aa_3 : str
            the 3-letter amino acid code in string format

        Returns
        ----------
        the float associated with the inputted 3-letter amino acid code
        """
        self.aa_number_dict = {'ALA':1, 'CYS':2, 'ASP':3, 'GLU':4, 'PHE':5, 'GLY':6, 'HIS':7, 'ILE':8, 'LYS':9, 'LEU':10, 'MET':11, 'ASN':12, 'PRO':13, 'GLN':14, 'ARG':15, 'SER':16, 'THR':17, 'VAL':18, 'TRP':19, 'TYR':20}

        return self.aa_number_dict[aa_3]

    def number_to_aa(self, num):
        """
        Converts from a float to one-letter AA code

        Parameters
        ----------
        num : float or int
            the number or integer associated with the 3-letter amino acid code
        Returns
        ----------
        the string for the associated three-letter amino acid code
        """
        self.number_aa_dict = {1:'A', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I', 9:'K', 10:'L', 11:'M', 12:'N', 13:'P', 14:'Q', 15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y'}

        return self.number_aa_dict[num]

    def starting_residue(self):
        """
        Find the starting residue number from the BMRB NMR-STAR file
        """
        with open(self.input, 'r+') as f:
            lines = f.readlines()
            for i in range(0, len(lines)):
                line = lines[i].strip().split()
                if len(line) > 0 and line[0] == '_Entity_comp_index.Entity_ID':
                    try:
                        self.starting_residue = int(lines[i+2].strip().split()[1])
                    except ValueError:
                        self.starting_residue = int(lines[i+2].strip().split()[0])
                #if len(line) > 0 and line[0] == '_Entity_poly_seq.Entity_ID':
                    #self.starting_residue = int(lines[i+2].strip().split()[2])

        f.close()

        return self.starting_residue


    def available_shifts(self):
        """
        Find the which chemical shifts are available in the NMR-STAR file
        """
        self.cs_list = []
        with open(self.input, 'r+') as f:
            lines = f.readlines()
            for k in range(0, len(lines)):
                line = lines[k].strip().split()
                if len(line) > 0 and line[0] == '_Atom_chem_shift.Assigned_chem_shift_list_ID':
                    y=0
                    # (2) start while loop to iterate over the residue numbers using the counter "y"
                    while len(lines[k+2+y].strip().split()) > 1 and float(lines[k+2+y].strip().split()[5]) <= (len(self.aa_seq) + 2):
                        atom = str(lines[k+2+y].strip().split()[7])
                        if atom in self.cs_list:
                            pass
                        else:
                            self.cs_list.append(atom)
                        y+=1
        f.close()


        return self.cs_list


    def shifts(self, flag='default', residue_labels='Y'):
        """
        Extract the BMRB chemical shifts and store them in an array for later use
        Parameters
        ----------
        xyz : file, asdf
        cs : string, denote the chemical shifts to store in the array (e.g. 'H,N', 'H,N,CA', 'H,N,CA,CB', or 'all')
        """
        
        start_res = self.starting_residue()
        
        # (1) identify which set of chemical shifts to extract
        # if methyls are requested, create the necessary dictionary and reset the nuclei_flag parameter
        if flag.upper() in ('DEFAULT'):
            nuclei_flag = self.cs_list
        elif flag.upper() == 'CARBONS':
            nuclei_flag = ['CA', 'CB']      # no CO shift in SDS micelle structure
        else:
            print('flag variable in the shifts function is not recognized!')
            sys.exit()

        # (2) verify that the shifts in self.cs_list match those in the datafile
        print('\n(1) CHEMICAL SHIFT VALIDATION')
        for nucleus in nuclei_flag:
            if nucleus not in self.cs_list:
                print('\nNOTE the chemical shift {} is not in the BMRB file!'.format(nucleus))
                print('please modify your selected chemical shifts')
                print('exiting now..........................................')
                sys.exit()
            else:
                print('{} exists in BMRB file'.format(nucleus))
        print('\nPassed chemical shift validation\nContinuing......................')

        # (3) create an empty array to store chemical shifts and residue numbers
        # if methyls are desired, only store 1H and 13C shifts a la HMQC
        if flag.upper() == 'DEFAULT':
            self.shifts_array = np.zeros(shape=(len(self.aa_seq), len(nuclei_flag)+2))
        elif flag.upper() == 'CARBONS':
            self.shifts_array = np.zeros(shape=(len(self.aa_seq), 5))
        else:
            print('\nNOTE within the function "shifts", please enter flag=methyls or flag=MILV, or some other combination of methyl-bearing residue types')
            print('exiting now........')
            sys.exit()


        # (4) read the NMR-STAR file and extract the chemical shifts
        with open(self.input, 'r+') as f:
            lines = f.readlines()
            for k in range(0, len(lines)):
                line = lines[k].strip().split()
                if len(line) > 0 and line[0] == '_Atom_chem_shift.Assigned_chem_shift_list_ID':
                    y=0
                    methyl_res_cnt = -1
                    # start while loop to iterate over the residue numbers using the counter "y"
                    while len(lines[k+2+y].strip().split()) > 1 and float(lines[k+2+y].strip().split()[5]) <= len(self.aa_seq):
                        try:
                            resn = int(lines[k+2+y].strip().split()[18])
                        except ValueError:
                            resn = int(lines[k+2+y].strip().split()[5])
                        aa_type = str(lines[k+2+y].strip().split()[6])
                        atom = str(lines[k+2+y].strip().split()[7])
                        ppm = float(lines[k+2+y].strip().split()[10])
                        print(resn,atom,ppm)

                        # iterate over the chemical shifts and store selected values
                        if flag.upper() == 'DEFAULT':
                            for l, spin in enumerate(nuclei_flag):
                                #print(spin, l)
                                if atom == spin:
                                    self.shifts_array[resn-start_res, 0] = self.aa_to_number(aa_type) #self.shifts_array[resn-self.starting_residue, 0] = aa_type
                                    self.shifts_array[resn-start_res, 1] = resn
                                    self.shifts_array[resn-start_res, l+2] = ppm
                                #else:
                                #    self.shifts_array[resn-start_res, 0] = self.aa_to_number(aa_type) #self.shifts_array[resn-self.starting_residue, 0] = aa_type
                                #    self.shifts_array[resn-start_res, 1] = resn
                                
                        elif flag.upper() == 'CARBONS':
                            #for l, spin in enumerate(nuclei_flag):
                            if atom == 'C':
                                self.shifts_array[resn-start_res, 0] = self.aa_to_number(aa_type) #self.shifts_array[resn-self.starting_residue, 0] = aa_type
                                self.shifts_array[resn-start_res, 1] = resn
                                self.shifts_array[resn-start_res, 2] = ppm
                            elif atom == 'CA':
                                self.shifts_array[resn-start_res, 0] = self.aa_to_number(aa_type) #self.shifts_array[resn-self.starting_residue, 0] = aa_type
                                self.shifts_array[resn-start_res, 1] = resn
                                self.shifts_array[resn-start_res, 3] = ppm
                            elif atom == 'CB':
                                self.shifts_array[resn-start_res, 0] = self.aa_to_number(aa_type) #self.shifts_array[resn-self.starting_residue, 0] = aa_type
                                self.shifts_array[resn-start_res, 1] = resn
                                self.shifts_array[resn-start_res, 4] = ppm
                                
                        # (4) iterate the counter within the while loop to move to the next line
                        y+=1

        return self.shifts_array


    def sparky(self, cs='NH', w1='N', w2='H', output='sparky.txt'):
        """
        Converts the desired chemical shifts into a Sparky peak list

        Parameters
        ----------
        cs : str
            the desired chemical shifts for the Sparky peak list ("NH", "METHYL", "HACA", "HBCB", "NHCO")
        w1 : str
            denotes the 1st axis (w1) in the Sparky peak list ("N", "H", "CO", "CA", "C")
        w2 : str
            denotes the 2nd axis (w2) in the Sparky peak list ("H", "N", "CO", "CA", "C")
        output : str
            name of the outputted Sparky peak list (default is "sparky.txt")
        """

        # Figure out the indices of H and N within self.cs_list
        # Use this to slice out H and N chemical shift values from self.shifts_array
        try:
            H_loc = self.cs_list.index('H')
            N_loc = self.cs_list.index('N')
        except ValueError:
            H_loc = self.cs_list.index('HN')

        if cs.upper() in ('NH', 'HN', 'CN', 'NC', 'METHYL'):
            # (1) trim the chemical shift list to exlclude prolines and other non-assigned residues
            data_cnt = 0
            res_list = [] ; res_type = [] ; h_list = [] ; x_list = []
            for n in range(0, len(self.shifts_array)):
                if np.abs(self.shifts_array[n,2 + H_loc]) > 0:
                    data_cnt+=1
                    res_list.append(int(self.shifts_array[n,1]))   # append the residue number
                    if cs.upper() in ('NH', 'HN'):
                        res_type.append(self.number_to_aa(self.shifts_array[n,0]))  # append the residue type
                        h_list.append(self.shifts_array[n,2 + H_loc])  # append the 1HN chemical shift
                        x_list.append(self.shifts_array[n,2 + N_loc])  # append the 15N chemical shift
                else:
                    pass

            # reformat the arrays that store chemical shifts as numpy arrays
            h_list = np.asarray(h_list)  ;  x_list = np.asarray(x_list)
            final_res_list = []


            # combine the residue type and number --> skip if methyls are chosen, as this is implemented for them by default
            if cs.upper() not in ('METHYLS', 'METHYL'):
                for r,val in enumerate(res_list):
                    final_res_list.append(res_type[r]+str(int(val)))
            else:
                final_res_list = self.shifts_array[:,0]

            # (2) create numpy array from the trimmed list of chemical shifts and save to file
            if cs.upper() in ('NH', 'HN'):
                if w1 == 'N' and w2 == 'H':
                    self.trimmed_shifts = np.column_stack([final_res_list, x_list.astype(np.object), h_list.astype(np.object)])
                    np.savetxt(output, self.trimmed_shifts, fmt=['%s'+'N-H', '%.3f', '%.3f'], header='Assignment\tw1\tw2', comments='', delimiter='\t')
                elif w1 == 'H' and w2 == 'N':
                    self.trimmed_shifts = np.column_stack([final_res_list, h_list.astype(np.object), x_list.astype(np.object)])
                    np.savetxt(output, self.trimmed_shifts, fmt=['%s'+'H-N', '%.3f', '%.3f'], header='Assignment\tw1\tw2', comments='', delimiter='\t')
            else:
                print('\nNOTE within the function "sparky", please enter cs=HN')
                print('exiting now........')
                sys.exit()

        print('\n(2) The Sparky peak list for {} has been successfully created'.format(cs))


    def talos(self, file_out='talos.tab'):
        """
        Creates a TALOS-N input file based on the NMR-STAR file

        Returns
        ----------
        talos.tab : file
            a TALOS-formatted file containing the assigned chemical shifts
        """

        # (1) create a talos.tab file and write out the header
        fout = open('output/'+file_out, 'w')
        fout.write("%s\t%s\t%s\n" % ("DATA", "FIRST_RESID", self.starting_residue))
        fout.write("%s\t%s\t%s\n" % ("DATA", "SEQUENCE", self.aa_seq))
        fout.write("%s\t%s\t%s\t%s\t%s\n" % ("VARS", "RESID", "RESNAME",  "ATOMNAME", "SHIFT"))
        fout.write("%s\t%s\t%s\t%s\t%s\n" % ("FORMAT", "%4d", "%1s", "%4s", "%8.3f"))

        # (2) loop over the chemical shifts arrary (self.shifts_array) and the list of available nuclei (self.cs_list)
        # write out to the TALOS file
        for i, row in enumerate(self.shifts_array):
            for j in range(len(self.cs_list)):
                if row[0] != 0:  # check if the row has data (i.e. aa_to_num != 0)
                    if row[j+2] != 0:  # check if the chemical shift is available (e.g. for Pro HN, N)
                        fout.write("%.0f\t%s\t%s\t%.3f\n" % (int(row[1]), self.number_to_aa(row[0]), self.cs_list[j], float(row[j+2])))
        fout.close()

        
        
### TEST 

new = BMRB('bmrb_5744.txt')
print(new.sequence())
print(new.available_shifts())
#print(new.shifts('carbons'))
print('done')
np.savetxt('output_5744.txt', new.shifts('carbons'))