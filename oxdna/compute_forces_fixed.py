#Reads simulation external force file and determines the pulling force on a base pair with nucleotides trapped to one coordinate in space (i.e., "frozen in" bases if stiffness is sufficiently strong); i.e., type = 'trap' with rate = '0'.
#Intended to determine the pulling force on just a single base pair, but written so that can handle multiple. The forces on base pairs must be provided sequentially in external force file (i.e., first trap with type = trap is pair_0 nt_0, next type = trap is pair_0 nt_1, then third type = trap pair_1 nt_0, followed by type = trap pair_1 nt_1, etc.) 
#This is done innately in my prep_pull_forces.py program to my knowledge.

#Usage: python3 [input_file for oxDNA simulation] [main output file with forces printed on every time step] [OPTIONAL: times output of kink_finder_Harrison_method_new.py, energetic kink_finder.py, or fray_finder.py (for pruning to check stacked configurations only!)]

import sys
from oxDNA_analysis_tools.UTILS.data_structures import *
from oxDNA_analysis_tools.UTILS.RyeReader import *
from oxDNA_analysis_tools.external_force_utils.forces_modified import *
from oxDNA_analysis_tools.external_force_utils.force_reader_modified import *
import numpy as np
import os

input_file = sys.argv[1]
#comp_list_txt = open(sys.argv[2])
output_file_txt = open(sys.argv[2], "w")
try:
    times_out = sys.argv[3]
except:
    times_out = "" #this is falsy so can be used to check if times file was provided

#Used to parse initial conditions file for the "frozen in" base initial coordinates (rather than repeating the read-in of initial configuration done in prep_slippery_forces.py).
#Should correctly parse the external force file for only the type = "trap" bases.
def read_forces_from_file(file_dir):
    forces = read_force_file(file_dir)
    ID_list = []
    stiff_list = []
    trap_coords_list = []
    for force in forces:
        print(force.get("type"))
        if force.get("type") == "trap":
            ID_list.append(int(force.get("particle")))
            stiff_list.append(float(force.get("stiff")))
            trap_coords_list.append([float(coord) for coord in force.get("pos0")])
    return ID_list, stiff_list, trap_coords_list

#Builds a list of all stacked configuration time steps from times output of kink_finder_Harrison_method.py or energetic kink_finder.py.
#Splits line in text file on whitespace, and if string "STACKED" appears, will add the time to the list if there is no string "FRAYED" also in that list.
#Assumes the first entry of the list is the time step of the configuration.
#(Why check both? My structural kink finder has a chance though rare to print both STACKED and FRAYED. It is not likely, but because configurations can fall in weird gray areas
# where the ends look stacked, but the ends are frayed, this is a possibility.)
def get_stacked_times(file_dir):
    times_list = []

    file = open(file_dir)
    while True:
        curLine = file.readline()
        if not curLine:
            break
        curLineSplit = curLine.split()
        if "STACKED" in curLineSplit:
            if "FRAYED" not in curLineSplit:
                times_list.append(int(float(curLineSplit[0])))
    file.close()

    return times_list
    
#Grabs the total number of steps of simulation from input file.
#This is a silly function due to an oversight (afaik, anyway) in how oxDNA's get_input_parameter works.
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

#Determines the normal of the helix currently by drawing a line from the center of helix of the previous base pair from 5' to 3' end to the center of the helix for the current base pair.
#Assumes for now that provided ID is on the strand containing indices from 0 to <N_bp - 1>. Should only affect the sign of the helix normal direction in this case currently.
#This function only really intended for ends of a duplex. If not, would probably average the normal direction of ID, ID + 1 and ID, ID - 1.
#It's unclear also to me yet if this is a better method than just averaging the base-normal vectors of the two normals of ID, comp_ID. 
#def get_helix_normal(conf,system,top_info,ID,COM_coords,comp_pairs_list):
    # #Find the strand index for ID
    # for strand in system.strands:
        # for monomer in strand.monomers:
            # if ID == monomer.id:
                # strand_index = strand.id - 1
    
    # if strand_index != 0:
        # print("Possible error: get_helix_normal was fed an index that was not on the strand from 0 to <N_bp -1>; unsure if will affect sign of torque.")
                
    # #For simplicity, I'm constructing my own list of all indices on the current strand.
    # rel_ID_list = []
    # for monomer in system.strands[strand_index].monomers:
        # rel_ID_list.append(monomer.id)
        
    # #Let ID_2 = ID + 1, unless ID = max_nuc_index on in-tact strand
    # ID_2 = ID + 1
    # if ID_2 not in rel_ID_list:
        # ID_2 = ID - 1
    
    # #Find complement to ID_2:
    # for pair in comp_pairs_list:
        # if ID_2 == pair[0]:
            # comp_ID_2 = pair[1]
        # elif ID_2 == pair[1]:
            # comp_ID_2 = pair[0]
    
    # #Find center of helix of ID_2.
    # COM_coords_2 = np.average([conf.positions[ID_2],conf.positions[comp_ID_2]], axis = 0) #COM of base pair (i.e., center of helix)
    
    # #Finally, determine normal to helix for this end of strand.
    # if ID_2 > ID:
        # helix_normal = np.subtract(COM_coords,COM_coords_2)
    # else: 
        # helix_normal = np.subtract(COM_coords_2,COM_coords)
    
    # #Normalize vector here.
    #helix_normal_dir = helix_normal/np.linalg.norm(helix_normal)
    
    #return helix_normal_dir

#Define helix normal to be the average of the normals of the two bases of the pair. 
#This is e_3 vector on page 11 of: https://www.semanticscholar.org/reader/43a1c055243ab732c37d2cd48bbd1ebf962aadec (supplement to DNA elasticity from coarse-grained simulations: the effect of groove asymmetry; I learned this after writing this)
def get_helix_normal(conf,pair):
    rel_a3s = [-1* conf.a3s[pair[0]], conf.a3s[pair[1]]] #as these are on two different strands, will face opposite directions, so must take the negative of one of the base-normal vectors
    avg_a3 = np.average(rel_a3s, axis = 0)
    norm_avg_a3 = avg_a3 / np.linalg.norm(avg_a3)
    
    return norm_avg_a3

#Main function that computes the forces acting on a complementary pair of indices that each have an external force of type = 'trap' (currently also assumes rate = 0). 
#These traps have spatial trap position at trap_pair_coords and each trap having stiffnesses stiff_pair. 
#This force vector is then projected onto the local tangent of the helix (this component I believe is representative of the local pulling force on the nucleotides along the helix normal).
def compute_forces_fixed(conf, pair, trap_pair_coords, stiff_pair):
    COM_coords = np.average([conf.positions[pair[0]],conf.positions[pair[1]]], axis = 0) #COM of base pair (i.e., center of helix)
    if pair[0] < pair[1]:
        #helix_normal = get_helix_normal(conf,system,top_info,pair[0],COM_coords,comp_pairs_list)
        helix_normal = get_helix_normal(conf,pair)
    else:
        #helix_normal = get_helix_normal(conf,system,top_info,pair[1],COM_coords,comp_pairs_list)
        pair.reverse() #reverse pair order
        helix_normal = get_helix_normal(conf,pair)
    proj_Fs = [] #force projected onto helix normal axis
    for ID, trap_coords, stiff in zip(pair, trap_pair_coords, stiff_pair):
        cur_coords = conf.positions[ID]
        dis = np.subtract(trap_coords, cur_coords) #displacement from trap position (vector from trap to current position; feel force along dis direction :))
        F = stiff * dis #afaik from documentation... it gives reasonable values for forces
        proj_Fs.append(np.dot(F,helix_normal))
    
    return proj_Fs
            
#Create TopInfo and TrajInfo objects that are used to read in trajectory.
traj_file = get_input_parameter(input_file, "trajectory_file")
top_file = get_input_parameter(input_file, "topology")
top_info, traj_info = describe(top_file, traj_file) 

#Read in the topology (creates a system with same hierarchy of System -> Strands -> Monomer in same way as old analysis code)
#system,monomers = strand_describe(top_file) #Creates System and list of Monomers; see data_structures.py for instance variables and believe monomers list is indexed exactly the same as oxview indexes nucleotides

#Find number of configurations in current simulation (assumes a linear time scale).
#'steps' holds the total number of time steps of the simulation (the 'steps' key in input file), and the print_conf_interval is directly grabbed from input file as well.
steps = int(float(get_num_steps(input_file)))
print_conf_interval = int(float(get_input_parameter(input_file,"print_conf_interval")))
num_confs = steps/print_conf_interval

#Read in from the input file the name of the external forces file to be read in.
external_forces_file = get_input_parameter(input_file, "external_forces_file")

#Read in IDs of base pairs and their initial trap coordinates/stiffnesses.
ID_list, stiff_list, trap_coords_list = read_forces_from_file(external_forces_file)

#Parse comp_list_txt input. 
#comp_pairs_1D_str = comp_list_txt.readline().split() #raw read-in; reads in as strings and as a 1D list
#comp_pairs_list = []
#for i in range(int(len(comp_pairs_1D_str)/2)):
#    comp_pairs_list.append([int(comp_pairs_1D_str[2*i]),int(comp_pairs_1D_str[2*i + 1])])
    
#Create a more formal list of the base pairs and the coordinates of each of the traps of the two bases in the base pair.
#In other words, rel_comp_pairs_list holds pairs of form [ID, comp_ID] and trap_pair_coords_list holds the corresponding trap coordinates of [ID_trap_coords, comp_ID_trap_coords]
#This structure feels needlessly confusing but I can't think of a better way just yet.
rel_comp_pairs_list = []
trap_pair_coords_list = []
stiff_pair_list = []
for i in range(int(len(ID_list)/2)):
    rel_comp_pairs_list.append([ID_list[2*i],ID_list[2*i + 1]])
    trap_pair_coords_list.append([trap_coords_list[2*i],trap_coords_list[2*i + 1]])
    stiff_pair_list.append([stiff_list[2*i],stiff_list[2*i + 1]])

print("Complementary pairs, corresponding trap coordinates, and spring stiffnesses as provided in external forces file:")
print(rel_comp_pairs_list, "\n" + str(trap_pair_coords_list), "\n" + str(stiff_pair_list))

#Build list of times that are stacked to analyze.
if times_out:
    stacked_times = get_stacked_times(times_out)
    if not stacked_times:
        raise Exception("A times file was provided but no stacked configurations found. Check times file.")
    print("Reading in", len(stacked_times), "stacked configurations.")
else:
    stacked_times = [] #this will count as falsy
    print("No times file provided. Reading in all", int(num_confs), "configurations.")

increment = 10000

#This loop structure could be made more readable... read in an increment number of confs at a time, and determine the phase of each conf in this increment. Then, loop to next increment. Don't read all configurations at once in case of memory issues (**may not actually be an issue? another bug was giving segfaults/memory issues but not sure if read-ins of large # of configs can be problematic still).
for i in range(int(num_confs//increment)):
    #Read in an increment number of configurations
    cur_confs = get_confs(top_info, traj_info, i*increment, increment)
    
    #Prune this batch of configurations for just the stacked configurations. If no stacked_times provided, assume all configurations are to be read in as stacked (vacuously true! :)).
    stacked_confs = []
    if stacked_times:
        for conf in cur_confs:
            if conf.time in stacked_times:
                stacked_confs.append(conf)
    else:
        stacked_confs = cur_confs
    
    #Loop over every conf in this batch:
    for conf in stacked_confs:
        #print("Time step:", conf.time)
        #print("ID || cur_coords trap_coords COM_coords dis r F np.cross(r,F)")
        
        #DO NOT INBOX HERE. Was confused for a long while why the coordinates do not match that of trajectory file. It is because I inboxed here previously.
        #You do *not* want to inbox here, as the external forces should maintain the ring relatively stationary at center of box, and otherwise screws up the analysis as want the trap coordinates to remain stationary (and otherwise does not if inboxing shifts the coordinate system).
        
        #Attempt to compute the force on the bases for each conf.
        for pair, trap_pair_coords, stiff_pair in zip(rel_comp_pairs_list, trap_pair_coords_list, stiff_pair_list):
            proj_Fs = compute_forces_fixed(conf, pair, trap_pair_coords, stiff_pair)
            output_file_txt.write(str(conf.time) + " " + str(pair[0]) + "," + str(pair[1]) + " " + str(proj_Fs[0]) + " " + str(proj_Fs[1]) + " " + str(proj_Fs[0]/2 + proj_Fs[1]/2) + "\n")
                
remove_chars = len(os.linesep)
try:
    output_file_txt.truncate(output_file_txt.tell() - remove_chars) #removes trailing newline in text file
    print("Main output with torques on each base printed to", sys.argv[3] + ".")
except: #if above fails, nothing was printed (as no trailing newline)
    print("No configurations were found in trajectory file. Make sure provided correct directory (else, unsure why this branch would get hit).")
output_file_txt.close()
