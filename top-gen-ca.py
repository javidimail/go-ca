#!/usr/bin/python 

import sys
domain = sys.argv[1]

if domain =='ch1': pdbfile, n = ['ch1_mini_aligned_with_cl-H-1-93-ca.pdb',93]
elif domain == 'cl': pdbfile, n = ['1fh5-L-112-214-ca.pdb',103]
elif domain == 'ch2': pdbfile, n = ['1fh5-L-112-214-ca.pdb',100]
elif domain == 'ch3': pdbfile, n = ['1fh5-L-112-214-ca.pdb',101]

topfile = open('toppar/%s-ca.top'%domain, 'w')

topcomment = ['* Topology file for a single bead GO model system\n', '*\n','    28    6\n']

topfile.writelines(topcomment)

for i in xrange(1, n+1, 1): topfile.write('MASS    %d  A%d   100\n' %(i, i))

print >> topfile, "\nDECL -CA    ! '-'=preceding atom"
print >> topfile, "DECL +CA    ! '+'=next atom"
print >> topfile, "DECL #CA    ! '#'=second next atom\n"

print >> topfile, "AUTOGENERATE ANGLES DIHEDRAL\n"

for i in xrange(1, n, 1): topfile.write("\nRESI P%d 0.00\nGROUP\nATOM CA A%d 0.00\nBOND CA +CA\n" %(i, i))

#last residu topology
topfile.write("\nRESI P%d 0.00\nGROUP\nATOM CA A%d 0.00\n" %(n, n))

