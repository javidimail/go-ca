#!/bin/bash

#Thsi short code cuts CAs from pdb file and put them in another file! That's it!
# please follow this format: ./pdb-ca-cutter.sh pdbfile chainID lowerlimit upperlimit

pdbfile=$1
chainID=$2
lowerLimit=$3
upperLimit=$4

pdbout=${pdbfile//.pdb}-${chainID,,}-$lowerLimit-$upperLimit-ca.pdb

# FOR STANDART PDB FILE
grep -e "ATOM.*CA" $pdbfile | awk -v id=$chainID -v l=$lowerLimit -v u=$upperLimit '$5 == id && $6>=l && $6<=u' > $pdbout

# FOR CHARMM PDB FILE
#grep -e "ATOM.*CA" $pdbfile | awk -v id=$chainID -v l=$lowerLimit -v u=$upperLimit '$5>=l && $5<=u' > $pdbout

# NOTE
# ADD "TER" to the last line of file! coz charmm os not able to terminate PDB file without that!!
