#Helper program to kink_finder_Harrison_method.py that attempts to automatically give an input to the program (a list of all the complementary pairs near the nick site in a space-separated list).
#This has a dual purpose as an input to oat align.py function for focusing camera on nick site, but you may need to up max_shift to focus on a larger region of the ring.
#However, must be careful with this as there is a chance ring will not hybridize as intended. Must check output with oxview. 

#You must provide the indices of what I call nuc_1 and nuc_2:
#nuc_1: lower ID of the two sites opposite of the nick usually (second ID from 5' to 3' in the case of wrapping)
#nuc_2: higher ID of the two sites opposite of the nick usually; usually nuc_1 + 1 (first ID from 5' to 3' in the case of wrapping)

#Searches all configurations for hybridized rings, or, if max_num_confs provided, the provided number of snapshots.

#Requires older version of oxdna_analysis_tools (before RyeReader). Newer version modifies output_bonds.py such that it is unusable for this purpose.
#Because of this, you **must** hard-code path to older oxDNA_analysis_tools src folder (as done below with sys.path.append for my case) an old version if want to minimize path conflicts with newest installed versions of oxDNA analysis tools.
#If still conflicts, try using a virtual environment without oxDNA installed and hard-code to the src folder of the old version of OAT.

#Usage: python3 [input file] [program output (space-separated list of indices)] [nuc_1] [nuc_2] [optional: max_num_confs for redundant search (will search max_num_confs for hybridized rings to make sure comp_pairs_list is consistent for each hybridized ring in case of errors)]

import sys
sys.path.insert(0,"/mnt/c/Users/thosw/Downloads/oxdna_analysis_tools-master/src/")
from oxDNA_analysis_tools.UTILS.base import *
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2,get_input_parameter
from oxDNA_analysis_tools.output_bonds import output_bonds

#PROGRAM PARAMETERS:
#Toggle the max range about nick site to check for kinking (max_shift = 2 checks 6-base window)
max_shift = 2

#Below contains code thanks to Daniel Hollis.
#Very slightly modified is_hybridized function of Daniel's to check if just the nucleotides near the nick site are hybridized. Returns true if current configuration is hybridized near nick site.
def is_nick_hybridized(rel_IDs, system):
    N_nuc = system.get_N_Nucleotides()
    
    hybridized = True

    #this searching algorithm could be made faster
    for i in rel_IDs:
        pair_found = False
        for j in range(0, i):
            if(is_paired(system, i, j)):
                pair_found = True
                break
        if(not(pair_found)):
            for k in range(i+1, N_nuc):
                if(is_paired(system, i, k)):
                    pair_found = True
                    break

        if(not(pair_found)):
            hybridized = False
            break
    return(hybridized)

#Below contains code thanks to Daniel Hollis.
#returns true if nucleotides with ids ID_1 and ID_2 are a base pair.
#assumes read_H_bonds_output_bonds already ran to map output of output_bonds to system object
def is_paired(system, ID_1, ID_2):
    bonded = False
    if ID_2 in system._nucleotides[ID_1].interactions:
        bonded = True

    return (bonded)

#Below contains code thanks to Daniel Hollis. 
#Finds the most commonly reported native pairs for nucleotides near the nick site as determined by get_rel_pairs_list. Searches configurations up to max_time (determined by optional runtime argument max_num_confs, else searches all currently).
#Returns a list of all of the necessary native pairs near the nick site. Each entry in outputted list is a complementary pair [nuc, comp].
#Important note: Alters state of reader -- must refresh reader after use as far as I can tell (couldn't find an easy way not to alter state) else will start pointer from time step after finds a hybridized ring. 
def get_comp_list(input_file, traj_file, top_file, rel_pairs, max_time):
    # find and save a configuration in the trajectory which has no broken base pairs. this will be used to determine
    # which nucleotides are complementary to which.
    temp_reader = LorenzoReader2(traj_file, top_file)
    temp_system = temp_reader._get_system()
    
    rel_IDs = []
    for pair in rel_pairs:
        if pair[0] not in rel_IDs:
            rel_IDs.append(pair[0])
        if pair[1] not in rel_IDs:
            rel_IDs.append(pair[1])
    
    #Attempt to find all hybridized rings up to max_time step in simulation. 
    all_comp_pairs_lists = [] #an array of lists of complementary pairs, which we will take the mode of later to be the best guess for actual hybridized ring
    while temp_system and temp_system._time <= max_time: #the first check may be redundant, but keeping just in case
        temp_system.map_nucleotides_to_strands()
        out = output_bonds(input_file, temp_system)
        temp_system.read_H_bonds_output_bonds(out)
        if(is_nick_hybridized(rel_IDs,temp_system)): #determine if this time step has the nick sites fully hybridized
            #Build list of complementary pairs for this time step.
            comp_pairs = [] 
            for rel_ID in rel_IDs:
                comp_ID = temp_system._nucleotides[rel_ID].interactions[0]
                comp_pairs.append([rel_ID, comp_ID])
            all_comp_pairs_lists.append(comp_pairs)
        temp_system = reader._get_system()
    if len(all_comp_pairs_lists) == 0:
        raise Exception("Could not find a configuration where the nucleotides were all bonded.")
    
    #thanks to: https://stackoverflow.com/questions/10797819/finding-the-mode-of-a-list (simple one-liner solution, rather than iterating the list and keeping counts)
    mode_comp_pairs = max(all_comp_pairs_lists, key=all_comp_pairs_lists.count)
    print("Found", len(all_comp_pairs_lists), "hybridized rings.")
    if len(all_comp_pairs_lists) == all_comp_pairs_lists.count(mode_comp_pairs):
        print("All comp pairs lists were consistent.")
    else:
        print(len(all_comp_pairs_lists) - all_comp_pairs_lists.count(mode_comp_pairs), "configurations had comp pairs lists that differed from the mode.")
    return mode_comp_pairs
        
#Returns a list of consecutive bases whose angles between their base-normal vectors are checked for kinking in main program
#Should account wrapping around any of the indices to the corresponding max/min index of the current strand
#Builds the list in such a way that the pairs closest to the nick site are earlier in the list.
def get_rel_pairs_list(nuc_1, nuc_2, max_shift, system):
    N_nuc = system.get_N_Nucleotides()
    
    #Determine internal strand ID for nuc_1 (and nuc_2, assuming user provided correctly indices on same strand)
    for strand in system._strands:
        for nucleotide in strand._nucleotides:
            if nuc_1 == nucleotide.index:
                strand_index = strand.index
                
    #For simplicity, I'm constructing my own list of all indices on the current strand.
    ID_list = []
    for nucleotide in system._strands[strand_index]._nucleotides:
        ID_list.append(nucleotide.index)
    
    min_nuc_index = min(ID_list)
    max_nuc_index = max(ID_list) 
    
    rel_pairs = [[nuc_1,nuc_2]] #nuc_1 and nuc_2 assumed as first pair and build rest of list from here
    
    for shift in range(0,max_shift):
        nuc_a_pos_shift = nuc_2 + shift #first nucleotide in positive direction
        nuc_b_pos_shift = nuc_2 + shift + 1 #and second nucleotide in positive direction to be checked pair-wise with first for stacking

        if nuc_a_pos_shift > max_nuc_index:
            diff = nuc_a_pos_shift - max_nuc_index 
            nuc_a_pos_shift = min_nuc_index + diff - 1
        if nuc_b_pos_shift > max_nuc_index:
            diff = nuc_b_pos_shift - max_nuc_index
            nuc_b_pos_shift = min_nuc_index + diff - 1

        rel_pairs.append([nuc_a_pos_shift,nuc_b_pos_shift])

        nuc_a_neg_shift = nuc_1 - shift - 1
        nuc_b_neg_shift = nuc_1 - shift

        if nuc_a_neg_shift < min_nuc_index:
            diff = min_nuc_index - nuc_a_neg_shift
            nuc_a_neg_shift = max_nuc_index - diff + 1
        if nuc_b_neg_shift < min_nuc_index:
            diff = min_nuc_index - nuc_b_neg_shift
            nuc_b_neg_shift = max_nuc_index - diff + 1

        rel_pairs.append([nuc_a_neg_shift,nuc_b_neg_shift])
    
    return rel_pairs

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

#Handle input arguments. 
input_file = sys.argv[1]
output_file = open(sys.argv[2],"w")
nuc_1 = int(sys.argv[3]) #lower ID of the two sites opposite of the nick usually (second ID from 5' to 3' in the case of wrapping)
nuc_2 = int(sys.argv[4]) #higher ID of the two sites opposite of the nick usually; usually nuc_1 + 1 (first ID from 5' to 3' in the case of wrapping)
steps = int(float(get_num_steps(input_file)))
print_conf_interval = int(float(get_input_parameter(input_file, "print_conf_interval")))
num_confs = int(steps//print_conf_interval)
print("Simulation has", num_confs, "snapshots.")
try: 
    max_num_confs = int(sys.argv[5])
    if max_num_confs > num_confs:
        max_num_confs = num_confs
        print("Provided max_num_confs is greater than number of snapshots in simulation. Searching all", num_confs, "configurations for hybridized rings.")
    else:
        print("Searching", max_num_confs, "configurations for hybridized rings.")
except:
    max_num_confs = num_confs #if unprovided, assume want to search all configurations currently
max_time = max_num_confs*print_conf_interval

traj_file = get_input_parameter(input_file, "trajectory_file")
top_file = get_input_parameter(input_file, "topology")

reader = LorenzoReader2(traj_file, top_file)
system = reader._get_system()

#Main program calls for determining the native pairs near the nick site.
rel_pairs_list = get_rel_pairs_list(nuc_1, nuc_2, max_shift, system)
print("Relevant pairs of consecutive indices on in-tact strand:",rel_pairs_list)
comp_pairs_list = get_comp_list(input_file, traj_file, top_file, rel_pairs_list, max_time)
print("Reported complementary pairs near the nick site:",comp_pairs_list)

for i in range(len(comp_pairs_list)):
    if i != len(comp_pairs_list) - 1:
        output_file.write(str(comp_pairs_list[i][0]) + " " + str(comp_pairs_list[i][1]) + " ")
    else:
        output_file.write(str(comp_pairs_list[i][0]) + " " + str(comp_pairs_list[i][1]))

output_file.close()

print("Wrote space-separated list of indices to", sys.argv[2] + ".")
print("**If you want to be certain these are correct, I highly suggest checking two things before use in kink finder:")
print("1) The first and third index of the list are nuc_1 and nuc_2 as you provided.")
print("nuc_1: lower ID of the two sites opposite of the nick usually (second ID from 5' to 3' in the case of wrapping)")
print("nuc_2: higher ID of the two sites opposite of the nick usually; usually nuc_1 + 1 (first ID from 5' to 3' in the case of wrapping)")
print("2) To be sure, each consecutive pair of indices is near the nick site and complementary as you'd expect.")
print("Issues are unlikely but, depending on the issue, may not be obvious in output of kink finder.")