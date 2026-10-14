#Helper program that attempts to find all complementary pairs in a system by checking multiple configurations. If not provided by input, searches all configurations for complementary pairs, or, if max_num_confs is provided, will check only the requisite number of snapshots.
#If flagged to search only for fully hybridized configurations, instead of providing the mode of determined complementary pairs across all checked configurations, it will exclusively search for configurations where all nucleotides are bonded. Be wary doing so without searching many configurations (as well as using this program in general!).

#Requires older version of oxdna_analysis_tools (before RyeReader). Newer version modifies output_bonds.py such that it is unusable for this purpose.
#Because of this, you **must** hard-code path to older oxDNA_analysis_tools src folder (as done below with `sys.path.insert` for my case) an old version if want to minimize path conflicts with newest installed versions of oxDNA analysis tools.
#If still conflicts, try using a virtual environment without oxDNA installed and hard-code to the src folder of the old version of OAT.

#For usage, run "python3 get_comp_list.py --help".

import sys
sys.path.insert(0,"/mnt/c/Users/thosw/Downloads/oxdna_analysis_tools-master/src/")
import argparse
from oxDNA_analysis_tools.UTILS.base import *
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2,get_input_parameter
from oxDNA_analysis_tools.output_bonds import output_bonds
from typing import TextIO, List, Tuple

#Below contains code thanks to Daniel Hollis. 
#Builds list of complementary pairs for a given configuration.
#assumes read_H_bonds_output_bonds already ran to map output of output_bonds to system object
def find_comp_pairs(system, N_nuc: int) -> List[Tuple[int,int]]:
    comp_pairs: List[Tuple[int,int]] = [] 
    for ID in range(N_nuc):
        if system._nucleotides[ID].interactions: #if ID H-bonded
            comp_ID = system._nucleotides[ID].interactions[0] #record complement
            if (comp_ID, ID) not in comp_pairs: #lazy way to not double-count pairs
                comp_pairs.append((ID, comp_ID)) #record pair

    return comp_pairs

#Below contains code thanks to Daniel Hollis. 
def find_unpaired_IDs(comp_pairs_list: List[Tuple[int,int]], N_nuc: int) -> List[int]:
    #this searching algorithm could be made faster
    unpaired_IDs: List[int] = []
    for i in range(N_nuc):
        pair_found = False
        for j in range(0, i):
            if((i, j) in comp_pairs_list or (j, i) in comp_pairs_list):
                pair_found = True
                break
        if(not(pair_found)):
            for k in range(i+1, N_nuc):
                if((i, k) in comp_pairs_list or (k, i) in comp_pairs_list):
                    pair_found = True
                    break

        if(not(pair_found)):
            unpaired_IDs.append(i)
    
    return unpaired_IDs

#Below contains code thanks to Daniel Hollis. 
#Attempts to read in configurations up until `max_time` in the trajectory, and find all complementary pairs in each configuration.
#In the end, it also does some processing to determine what the most common (mode) complementary pairs list was across all configurations, and also produces some diagnostics (says if all nucleotides are bonded, and, if not, which aren't bonded).
#Important note: Alters state of reader -- must refresh reader after use as far as I can tell (couldn't find an easy way not to alter state) else will start pointer from time step after finds a hybridized ring. 
def get_mode_comp_list(input_file: str, traj_file: str, top_file: str, max_time: int, check_hybridized: bool) -> List[Tuple[int,int]]:
    reader = LorenzoReader2(traj_file, top_file)
    system = reader._get_system()

    N_nuc: int = system.get_N_Nucleotides()
    
    #For all configurations up to `max_time`, find complementary pairs in the configuration and record for comparison later.
    all_comp_pairs_lists: List[List[Tuple[int,int]]] = [] #an array of lists of complementary pairs, which we will take the mode of later to be the best guess for intended complementary pairs
    while system and system._time <= max_time: #the first check may be redundant, but keeping just in case
        system.map_nucleotides_to_strands()
        out = output_bonds(input_file, system)
        system.read_H_bonds_output_bonds(out)

        #Find this system's complementary pairs.
        comp_pairs: List[Tuple[int,int]] = find_comp_pairs(system, N_nuc)

        #Append this configuration's complementary pairs and keep for later, unless flagged to check only for hybridized ring in which case an additional check is performed.
        if check_hybridized:
            if not find_unpaired_IDs(comp_pairs, N_nuc): #if conf is not hybridized, continue
                all_comp_pairs_lists.append(comp_pairs)
        else:
            all_comp_pairs_lists.append(comp_pairs)
        
        system = reader._get_system()

    #Check if found any systems with complementary pairs (should definitely have, but just in case run into problems).
    if len(all_comp_pairs_lists) == 0:
        raise Exception("Unable to find a configuration where any nucleotides were bonded. Check if input files are all okay (and, if set `-c`, make sure it makes sense for this system).")
    
    #Determine the most frequently recorded complementary pairs lists.
    #thanks to: https://stackoverflow.com/questions/10797819/finding-the-mode-of-a-list (simple one-liner solution, rather than iterating the list and keeping counts)
    mode_comp_pairs = max(all_comp_pairs_lists, key=all_comp_pairs_lists.count)

    if len(all_comp_pairs_lists) == all_comp_pairs_lists.count(mode_comp_pairs):
        print("All complementary pairs lists across all configurations searched were consistent.")
    else:
        print(f"{len(all_comp_pairs_lists) - all_comp_pairs_lists.count(mode_comp_pairs)} configurations had complementary pairs lists that differed from the mode (this is likely if search a lot of configurations due to fraying).")

    unpaired_IDs: List[int] = find_unpaired_IDs(comp_pairs, N_nuc)
    if not unpaired_IDs:
        print("Provided complementary pairs list contains all nucleotide IDs (i.e., is hybridized).")
    else:
        print(f"Provided complementary pairs list does not contain all nucleotide IDs (i.e., it is not hybridized). The following IDs were unpaired: {unpaired_IDs}")
    
    return mode_comp_pairs
        
#Grabs the total number of steps of simulation from input file.
#This is a *silly* function due to an oversight in how oxDNA's get_input_parameter works.
#When call get_input_parameter(input_file,"steps"), it grabs newtonian_steps value at the very least if
#steps key precedes newtonian_steps key in input file. Kludged by modifying get_input_parameter from base.py.
def get_num_steps(input_file):
    fin = open(input_file)
    value = ''
    for line in fin:
        line = line.lstrip()
        if not line.startswith('#'):
            if "steps" in line:
                if "newtonian_steps" not in line:
                    value = line.split('=')[1].replace(' ','').replace('\n','')
    
    fin.close()
    if value == '':
        print("ERROR: Key {} not found in input file {}".format(parameter, input_file))
    return value

#Parse CLI input parameters.
prog = "get_comp_list.py"
parser = argparse.ArgumentParser(prog = prog, description='''Attempts to produce a space-separated list of complementary pairs by searching a trajectory. 
                                 This is intended for small structures (duplexes or minicircles). The program will create a list of all complementary pairs in each read-in configuration, 
                                 and output the mode list containing the most frequently recorded list of complementary pairs. If you only want a hybridized ring, then you can also set a flag below to instead ignore the mode if it is not hybridized and attempt to override with a hybridized ring. 
                                 It will also inform you if all nucleotides in the system are H-bonded in this outputted list. If any nucleotides are not H-bonded in this list, it will inform you which were not paired.''')
parser.add_argument("input_file", type = str, nargs = 1, help = "OxDNA input file path to be read in.")
parser.add_argument("-m", "--max_num_confs", type = int, nargs = 1, dest = "max_num_confs", help = '''optional: Max number of configurations to search in trajectory. If not provided, will search all configurations of the trajectory. 
                    Else, will only search up to the number provided (assuming does not exceed number in trajectory). Typically 80-90 percent of configurations for nondefective minicircles are fully hybridized (i.e., all nts are H-bonded).
                    Ideally should search at minimum ~100 for safety.''')
parser.add_argument("-c", "--check_hybridized", dest = "check_hybridized", action='store_true', default = False, help = "optional: If flagged, will *only* search for configurations that appear to be hybridized. This is likely to be the mode with the flag unset, but in cases where trajectory may involve lots of fraying, this may be useful as a sanity check.")
parser.add_argument("output_file", type = str, nargs = 1, help = "Output file containing complementary pairs. Every two entries comprises a pair of nucleotide indices.")

#Process CLI inputs with some minor diagnostics/prep.
args = parser.parse_args()
input_file: str = args.input_file[0]
steps: int = int(float(get_num_steps(input_file)))
print_conf_interval: int = int(float(get_input_parameter(input_file, "print_conf_interval")))
num_confs: int = int(steps//print_conf_interval)
max_num_confs: int
if args.max_num_confs[0]:
    max_num_confs = int(args.max_num_confs[0])
    if max_num_confs > num_confs:
        max_num_confs = num_confs
        print(f"Provided `max_num_confs` is greater than number of snapshots in simulation. Searching all {num_confs} configurations for complementary pairs.")
    else:
        print(f"Searching provided {max_num_confs} configurations for complementary pairs.")
else:
    max_num_confs = num_confs
    print(f"Searching all {num_confs} configurations for complementary pairs")
check_hybridized: bool = False
if args.check_hybridized:
    check_hybridized = True
    print("Checking only for configurations that appear to be hybridized (i.e., all nts are bonded).")

output_file_txt: TextIO = open(args.output_file[0],"w")
max_time: int = max_num_confs*print_conf_interval

traj_file: str = get_input_parameter(input_file, "trajectory_file")
top_file: str = get_input_parameter(input_file, "topology")

#Main program call for determining the native pairs.
comp_pairs_list: List[Tuple[int,int]] = get_mode_comp_list(input_file, traj_file, top_file, max_time, check_hybridized)

for i in range(len(comp_pairs_list)):
    if i != len(comp_pairs_list) - 1:
        output_file_txt.write(f"{comp_pairs_list[i][0]} {comp_pairs_list[i][1]} ")
    else:
        output_file_txt.write(f"{comp_pairs_list[i][0]} {comp_pairs_list[i][1]}")

output_file_txt.close()

print(f"Wrote space-separated list of pairs to {args.output_file[0]}.")
print("**Please be cautious here as it is possible that the pairs did not detected correctly. Must check oxview to be absolutely sure.")
print("**If fraying is common for your system, and do not set `-c` (for example, if because there are structural defects such as bulges that prohibit this), be aware that the mode may still be missing certain complementary pairs. Again, be wary.")