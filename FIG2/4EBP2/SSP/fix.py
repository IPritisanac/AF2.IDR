import os, sys
# to fix a 7-residue offset in BMRB ID 191141 (cloning scar)
fin = open('191141.cb', 'r')
fout = open('191141_fix.cb', 'w')
for line in fin:
    new = line.strip().split()
    res = float(new[0])
    res = res+7
    fout.write('%.0f\t%.3f\n' % (float(res), float(new[1])))
fin.close()
fout.close()

