#!/usr/bin/python
import numpy as np
import sys

domain = sys.argv[1]

#open parameter file to start writing

if domain =='ch1': pdbfile, n = ['ch1_mini_aligned_with_cl-H-1-93-ca.pdb',93]
elif domain == 'cl': pdbfile, n = ['1fh5-l-112-214-ca.pdb',103]
elif domain == 'ch2': pdbfile, n = ['ch2-aligned-cl-a-241-340-ca.pdb',100] 
elif domain == 'ch3': pdbfile, n = ['1cqk-ch3-aligned-cl-a-4-104-ca.pdb',101] #4 to 104
 
print pdbfile, domain, n

#coorfile = open('1fh5-H-123-215-ca.pdb', 'r')
coorfile = open(pdbfile, 'r')
paramfile = open('toppar/%s-ca.prm' %domain, "w")

# BONDS
paramcomment=['BOND !Bonded pptential\n','!K(CA-CA)=100*epsilon/a^2\n']
paramfile.writelines(paramcomment)
for i in xrange(1,n,1): paramfile.write("A%d A%d 100.0 3.8\n"%(i,i+1))

# ANGLES
paramcomment=['\nANGLE !Bond-angle pptential\n','!atom atom atom kappa theta0\n']
paramfile.writelines(paramcomment)
for i in xrange(1,n-1,1): paramfile.write("A%d   A%d   A%d   12.5   105.0\n"%(i, i+1, i+2))

# DIHEDRALS
paramcomment=['\nPHI   ! dihedral-angle potential\n' \
, '! dihedral_plot = kappa * [1 + col(n*phi-delta)]\n','!atom atom atom atom kappa n delta\n']
paramfile.writelines(paramcomment)
for i in xrange(1, n-2,1): 
    paramfile.write('A%d   A%d   A%d   A%d   0.90 3  -66.5\n'%(i, i+1, i+2, i+3))
    paramfile.write('A%d   A%d   A%d   A%d   2.27 2  -68.1\n'%(i, i+1, i+2, i+3))
    paramfile.write('A%d   A%d   A%d   A%d   2.91 1  -37.3\n'%(i, i+1, i+2, i+3))

# NONBONDED
paramcomment=['\nNBONDed CUTNB 16.0 CTOFNB 14.0 CTONNB 12.0 WMIN 1.5 E14Fac 1.0 EPS 1.0 -\n' \
, 'ATOM CDIE SWITCH VSWITCH VATOM BYGROUP ! non-bonded potential\n' \
, '!atom           e_min    r_min/2\n']
paramfile.writelines(paramcomment)
for i in xrange(1, n+1, 1): paramfile.write('A%d  0.0  -1e-12  22.71 \n'%i)

def ca_native(coorfile, paramfile):
  x=[]; y=[]; z=[]
  i=0; j=0
  paramcomment=['\nNBFIX\n' \
  ,'!modify certain interactions, here rmin is given instead of rmin/2\n' \
  ,'!atom atom    e_min   r_min\n']
  paramfile.writelines(paramcomment)

  for column in (raw.strip().split() for raw in coorfile):
    # Column[0-->11], an example: 
    # ATOM[0] 2508 CA THR H 123 8.891[6-x] -17.658[7-y]  13.330[8-z]  1.00 33.77  C
    if domain=='ch1':
       x.append(float(column[5])); y.append(float(column[6])); z.append(float(column[7]))
    else:
       x.append(float(column[6])); y.append(float(column[7])); z.append(float(column[8]))
 
   
  for i in xrange(1, len(x)+1,1):
      for j in xrange(i+3, len(x)+1, 1):
          p1 = np.array((x[i-1], y[i-1], z[i-1]))
          p2 = np.array((x[j-1], y[j-1], z[j-1]))
          dist = np.linalg.norm(p1-p2)
          if dist <= 8 : paramfile.write('A%d  A%d  -1.25  %f\n'%(i, j, dist))
          if dist > 8 : paramfile.write('A%d  A%d  -0.125  %f\n'%(i, j, dist)) 
  print len(x)

ca_native(coorfile, paramfile)
