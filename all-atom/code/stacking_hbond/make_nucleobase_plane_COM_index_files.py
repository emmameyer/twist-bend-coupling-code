"""
    Author: Rachel Bricker
    Year:   2025
"""

"""
    Creates index files 'nucleobase_vec_atoms.ndx' and 'nucleobase_COM_atoms.ndx' which are used as input for GROMACS utility 'traj'.
    
    THIS CODE REQUIRES:
        * ENERGY MINIMIZATION TO BE PERFORMED. HENCE, IT USES THE .gro FILE OUTPUTTED BY THE
          MDRUN THAT CARRIED OUT ENERGY MINIMIZATION
"""

"""
   usage: python3 make_nucleobase_plane_COM_index_files.py
      1. total number of nucleotides
      2. path to .gro file outputted by the mdrun that carried out energy minimization
      3. output directory
   
   example: python3 make_nucleobase_plane_COM_index_files.py \
            42 \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/dsDNA1/em.gro \
            /mnt/c/Users/brick/Documents/alkyl_chain_stuff/GROMACS_files/AMBER/dsDNA1/
"""

import sys
import numpy as np

# command line input
max_residue_id = int(sys.argv[1])
path           = str(sys.argv[2])
output_dir     = str(sys.argv[3])

# add trailing forward slash to directory path if necessary
output_dir_split = output_dir.split("/")
if output_dir_split[-1] != "":
    output_dir = output_dir + "/"

################################################################################################
#
# FUNCTIONS
#
################################################################################################

def make_index_files(path, max_residue_id):
    nucleobase_atoms = ["C2", "C4", "C5", "C6", "C7", "C8", "C5M", "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"]
    atoms_COM        = {}
    thymine          = ["O2", "O4"]
    adenine          = ["N6", "C8"]
    guanine          = ["O6", "C8"]
    cytosine         = ["O2", "N4"]
    atoms_vec        = {}

    with open(path, "r") as f:
        residue_id   = 1
        residue_name = None
        for line in f:
            line_split = line.split()
            if len(line_split) > 1:
                if line_split[0] != residue_name and residue_name != None:
                    residue_id += 1
                    if residue_id >= max_residue_id+1:
                        break
                if line_split[0][0:len(str(residue_id))] == str(residue_id):
                    residue_name = line_split[0]
                    
                    if line_split[1] in nucleobase_atoms:
                        atoms_COM.setdefault(residue_name, []).append(line_split[2])
                    
                    if ( ('A' in residue_name and line_split[1] in adenine)  or
                         ('T' in residue_name and line_split[1] in thymine)  or
                         ('G' in residue_name and line_split[1] in guanine)  or
                         ('C' in residue_name and line_split[1] in cytosine) ):
                        atoms_vec.setdefault(residue_name, []).append(line_split[2])

    file_name = "nucleobase_COM_atoms.ndx"
    with open(output_dir + file_name, "w+") as output:
        resi_ID = 1
        for residue in atoms_COM:
            output.write("[ RESI_" + str(resi_ID) + " ]\n")
            for atom in atoms_COM[residue]:
                if atom != atoms_COM[residue][-1]:
                    output.write(atom + " ")
                else:
                    output.write(atom + "\n")
            resi_ID += 1

    file_name = "nucleobase_vec_atoms.ndx"
    with open(output_dir + file_name, "w+") as output:
        resi_ID = 1
        for residue in atoms_vec:
            output.write("[ RESI_" + str(resi_ID) + " ]\n")
            for atom in atoms_vec[residue]:
                if atom != atoms_vec[residue][-1]:
                    output.write(atom + " ")
                else:
                    output.write(atom + "\n")
            resi_ID += 1

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    make_index_files(path, max_residue_id)

if __name__ == "__main__": 
    main()
