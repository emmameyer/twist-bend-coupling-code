#Program that reads in the output of get_comp_list.py or get_nick_comp_list.py, and determines if any complementary pairs in this list are frayed (melted) *energetically* (|V_HB| < .1(max|V_HB|) I think). Alternatively, you can manually write out a text file with a list of space-separated complementary bases all on one line that you want to check.
#It is slow, but works on trajectories that have already been run.

#Requires older version of oxdna_analysis_tools (before RyeReader). Newer version modifies output_bonds.py such that it is unusable for this purpose.
#Because of this, you **must** hard-code path to older oxDNA_analysis_tools src folder (as done below with sys.path.append for my case) an old version if want to minimize path conflicts with newest installed versions of oxDNA analysis tools.
#If still conflicts, try using a virtual environment without newest version of oxDNA installed and/pr hard-code path to the src folder of the old version of OAT.

#Usage: python3 [simulation input file] [list of space-separated pairs to check for melting (auto-generating attempted by get_comp_list.py or get_nick_comp_list.py)] [main output: text file with a list of frayed time steps and offending pairs] 

import sys
#sys.path.insert(0, "<your path to the older version of analysis tools, if necessary>")
from oxDNA_analysis_tools.UTILS.base import *
from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2,get_input_parameter
from oxDNA_analysis_tools.output_bonds import output_bonds

#Handle input arguments. 
input_file = sys.argv[1]
comp_list_txt = open(sys.argv[2])
output_file_txt = open(sys.argv[3], "w")

#Parse comp_list_txt input. 
comp_pairs_1D_str = comp_list_txt.readline().split() #raw read-in; reads in as strings and as a 1D list
comp_pairs_list = []
for i in range(int(len(comp_pairs_1D_str)/2)):
    comp_pairs_list.append([int(comp_pairs_1D_str[2*i]),int(comp_pairs_1D_str[2*i + 1])])

#Below contains code thanks to Daniel Hollis.
#Modified is_hybridized function of Daniel's to check if any of the provided list of complementary pairs is non-interacting.
#This version in particular returns the pairs found to be non-interacting instead of a boolean flag. If the list contains any elements, will evaluate truthy.
def get_broken_pairs(comp_pairs, system):
    broken_pairs = []
    for pair in comp_pairs:
        if not is_paired(system, pair[0], pair[1]):
            broken_pairs.append(pair)
    return broken_pairs

#Below contains code thanks to Daniel Hollis.
#returns true if nucleotides with ids ID_1 and ID_2 are a base pair.
#assumes read_H_bonds_output_bonds already ran to map output of output_bonds to system object
def is_paired(system, ID_1, ID_2):
    bonded = False
    if ID_2 in system._nucleotides[ID_1].interactions:
        bonded = True

    return (bonded)

#Read in trajectory and topology file names from keys in input file.
traj_file = get_input_parameter(input_file, "trajectory_file")
top_file = get_input_parameter(input_file, "topology")

#Fully set up initial configuration.
reader = LorenzoReader2(traj_file, top_file)
system = reader._get_system()
system.map_nucleotides_to_strands()
out = output_bonds(input_file, system)
system.read_H_bonds_output_bonds(out)
system.inbox()

while system:
    #Determine if any broken pairs in comp_pairs_list for current configuration.
    broken_pairs = get_broken_pairs(comp_pairs_list, system)
    if broken_pairs:
        output_file_txt.write(str(int(system._time)) + " " + "FRAYED " + str(broken_pairs) + "\n")
    else:
        output_file_txt.write(str(int(system._time)) + " " + "STACKED" + "\n")
        
    #Set up for next iteration (read in next configuration and determine H-interactions).
    system = reader._get_system()
    if system:
        system.map_nucleotides_to_strands()
        out = output_bonds(input_file, system)
        system.read_H_bonds_output_bonds(out)
        system.inbox()
    
remove_chars = len(os.linesep)
try:
    output_file_txt.truncate(output_file_txt.tell() - remove_chars)
    print("Time steps where provided base pairs were not all hybridized were printed to", sys.argv[3] + ".")
except:
    print("Failed to print any frayed times... best guess empty trajectory or indices file.")
output_file_txt.close()
