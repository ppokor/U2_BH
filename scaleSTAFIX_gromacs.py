"""
this is script to stafix-scale Gromacs topology
requirements:
	python 3 (remove input/output type definitions for is_to_be_scaled function to use this also with python 2)
    gromacs topology containing [ defaults ], [ atomtypes ], [ atoms ] and [ bonds ] sections. 
        [ bonds ], [ pairs ], [ angles ] and [ dihedrals ] had to be defined rather than [ dihedraltypes ] etc., see Gromacs manual for further details
usage:
	python scaleSTAFIX_gromacs.py parm_file scaling_factor [mask]
		[mask] is to be given in format resid1,resid2,resid3,... (e.g. 1,2,3,4). If no mask is specified, all RNA residues will be scaled"
flow:
	read original topology, get definitions of all present atomtypes from there and write a temporary topology that will have atomtypes of atoms that are to be scaled changed by appending "Y" to their atomtypename
    pick atomtypes that are to be scaled
    write new topology where [ nonbond_params ] section is added after [ atomtypes ] and contains definition of LJ params for scaled atomtypes (those with "Y" appended)
to generate Gromacs topology from Amber top with parmed run in shell:
    parmed -e <<EOF
    import parmed as pmd
    amber = pmd.load_file('amber.parm7', 'amber.rst7')
    amber.save('gromacs.top')
    amber.save('gromacs.gro')
    EOF
"""

import sys
import os
from math import sqrt
from typing import List

#definition of RNA residue names and of atom names for atoms that are NOT to be scaled
residue_namesA="A","A3","A5","RA","RA3","RA5"
residue_namesG="G","G3","G5","RG","RG3","RG5"
residue_namesC="C","C3","C5","RC","RC3","RC5"
residue_namesU="U","U3","U5","RU","RU3","RU5"
atom_namesSugar="O2\'","HO2\'","HO3\'","HO5\'"
atom_namesA="N3","N6","H61","H62"
atom_namesG="N1","H1","N2","H21","H22","N3"
atom_namesC="O2","N4","H41","H42"
atom_namesU="O2","N3","H3"

#
present_atom_types=[] #that are to be scaled, will be renamed by appending "Y" in the new topology
atomtypesdef=[] #deffinition of atomtypes read from original topology
topfileTMP="stafixTMP.top"

#takes array of nonbond_params lines, scales the epsilon and returns string to write to topology file
def scale_nonbond_params(scalefact: float, input_arr: List[str]) -> str:    
    data=input_arr
    outstr="\n[ nonbond_params ]\n; i    j func      sigma      epsilon\n"
    for i in range(0,len(data)):
        for j in range(0,i+1):
            datai=data[i].split()
            dataj=data[j].split()
            sigma=0.5*(float(datai[5])+float(dataj[5]))
            eps=sqrt(float(datai[6])*float(dataj[6]))
            eps=eps*scalefact
            outstr+=" %3s   %2s   1   %10.8f    %10.8f"%(datai[0],dataj[0],sigma,eps)+"\n"
    return outstr+"\n"

#This decides whether to scale the atom of a given name (atname) in a residue with name resname or not and returns True or False
def is_to_be_scaled(resname: str,atname: str) -> bool:
    scaling=False
    if resname in residue_namesA and atname not in atom_namesA+atom_namesSugar:
        scaling=True
    elif resname in residue_namesG and atname not in atom_namesG+atom_namesSugar:
        scaling=True
    elif resname in residue_namesC and atname not in atom_namesC+atom_namesSugar:
        scaling=True
    elif resname in residue_namesU and atname not in atom_namesU+atom_namesSugar:
        scaling=True
    if ("3" in resname and atname=="O3\'") or ("5" in resname and atname=="O5\'"): #dont scale terminal hydroxyl groups
        scaling=False
    return scaling

# check if all arguments were given, if topology file exists, if scaling factor is a valid number, if mask is not provided in 1-5 format. Exit if problem
intro="\nthis is script to STAFIX-scale parmed-generated gromacs topology\nusage:\n\tpython scaleSTAFIX_gromacs.py topology_file.top scaling_factor [mask]\n\n \
\t\t\n *.top file is to be generated via parmed\
\t\t[mask] is to be given in format resid1,resid2,resid3,... (e.g. 1,2,3,4). If no mask is specified, all RNA residues will be scaled\n" 
try:
	inputTOP=sys.argv[1]
	if not os.path.exists(inputTOP):
		print("\n!invalid input file!")
		raise
	scalefact=float(sys.argv[2])
	try:
		mask=sys.argv[3].split(",")
	except:
		mask="*"
	if "-" in mask:
		raise
except:
	print(intro)
	sys.exit()

inputTOPfilemane,inputTOPextension=os.path.splitext(inputTOP)
topfile_new=inputTOPfilemane+"STAFIX"+str(scalefact)+inputTOPextension
if mask=="*":
        mask_print="all RNA residues"
        mask = [str(i) for i in range(10000)]
else:
    mask_print=",".join(mask)


#read input topology file, make a new one (topfileTMP) where atoms that are to be scaled will be assigned new types and get list of them and load atomtype data
with open(inputTOP, "r") as procfile: 
    defaultsread=False
    atomread=False
    atomtyperead=False
    with open(topfileTMP,"w") as tmpfile:
        for pline in procfile:
            if defaultsread and ";" not in pline:
                combrule=pline.split()[1]
                if combrule!="2":
                    # here if comb-rule is not 2 we can't use sigma-epsilon read/write as implemented in this script. see gromacs manual Parameter files -> Non-bonded parameters for more details
                    print("\nplease use topology with combination rule 2 ([ defaults ] section)\n")
                    sys.exit()
                else:  
                    tmpfile.write(pline)   
                    defaultsread=False
            elif "[ defaults ]" in pline:
                defaultsread=True
                tmpfile.write(pline)
            elif atomtyperead and pline[0]=="[":
                atomtyperead=False
                tmpfile.write(pline)
            elif atomtyperead:
                if pline[0]==";":
                    tmpfile.write(pline)
                else:
                    try:
                        atomtype=pline.split()[0]
                        atomtypesdef.append(atomtype+"Y"+pline.split(atomtype)[1][1:])  #load definitions of all present atomtypes, later sort out those that are not to be scaled
                        tmpfile.write(pline)
                    except:
                        tmpfile.write(pline)
            elif "[ atomtypes ]" in pline:
                atomtyperead=True
                tmpfile.write(pline)
            elif atomread and (pline[0]!=";"):
                if "[ bonds ]" in pline:
                    atomread=False
                    tmpfile.write(pline)
                else:
                    try:
                        #    1         HO      1     U5   HO5'      1   0.429500     3.0240   ; qtot 0.4295
                        data=pline.split()
                        atomtype=data[1]
                        residue_number=data[2]
                        residue=data[3]
                        atom=data[4]
                        if residue_number in mask and is_to_be_scaled(residue,atom):
                            if atomtype not in present_atom_types:	#add new atom type to list of those atomstypes, that are to be scaled
                                present_atom_types.append(atomtype)
                            tmpfile.write(pline.split(atomtype, 1)[0][:-1]+atomtype+"Y"+pline.split(atomtype, 1)[1])
                        else: 
                            tmpfile.write(pline)
                    except:
                        tmpfile.write(pline)
            elif "[ atoms ]" in pline:
                atomread=True
                tmpfile.write(pline)
            else:
                tmpfile.write(pline)

#remove from atomtypesdef list those atom types, that belong to atoms that are not to be scaled
atomtypesdef2=[]
for atomtype in atomtypesdef:
    if atomtype.split()[0][:-1] in present_atom_types:
        atomtypesdef2.append(atomtype)
atomtypesdef=atomtypesdef2

#write new file with definitions of added atomtypes and nonbond_params for their pairs
with open(topfileTMP,"r") as tmpfile:
    with open(topfile_new,"w") as outputfile:
        outputfile.write(";\tModified topology file with stafix scaling factor "+str(scalefact)+" applied on residues: "+mask_print+"\n")
        atomread=False
        for tline in tmpfile:
            if atomread:    
                if tline and tline[0]!="[":
                    outputfile.write(tline)
                elif tline[0]=="[":
                    atomread=False
                    outputfile.write("; STAFIX atom types\n")
                    for scaled_type in atomtypesdef:
                        outputfile.write(scaled_type)
                    outputfile.write(scale_nonbond_params(scalefact, atomtypesdef)) 
                    outputfile.write(tline)
                else:
                    atomread=False
                    outputfile.write("; STAFIX atom types\n")
                    for scaled_type in atomtypesdef:
                        outputfile.write(scaled_type)
                    outputfile.write(scale_nonbond_params(scalefact, atomtypesdef)) 
            elif "[ atomtypes ]" in tline:
                atomread=True
                outputfile.write(tline)
            else:
                outputfile.write(tline)

os.remove(topfileTMP)
print(("\nnew file generated: " +topfile_new+"\nscaled residues: "+mask_print+"\n"))
