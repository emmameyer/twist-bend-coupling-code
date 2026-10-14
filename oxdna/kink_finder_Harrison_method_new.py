#Newest version of kink finder inspired by Harrison et al's method outlined in "Identifying physical causes of apparent enhanced
#cyclization of short DNA molecules with a coarse-grained model."
#This program uses the same principal of computing angles between base-normal vectors of consecutive pairs of bases in a 6-base window about the nick site on the in-tact strand. The paper may describe this better than I can.
#What's different: We change the max angle that counts as kinked from 90 -> 40. It also attempts to distinguish kinking at the nick site from kinking away from the nick site by computing the max angle computed between any pair in this window and comparing to the angle between the two nucleotides directly across the nick.
#It incorporates a distance check to determine fraying events.
#It reports fractions for each phase and gives an error bar by dividing the simulation into intervals and taking a STDEV of the fractions in each interval.
#On user end, must provide input file (with trajectory and topology files in corresponding directories) and provide the two nucleotides directly across the nick site as an input (nuc_1 and nuc_2). See usage at bottom of header.
#Outputs a file with the phase at every time step, and two separate files for the max angle and nuc_1 nuc_2 angle at every time step.

#Less detailed description of how it works (TL;DR):
#Determines phase of every printed time step and in almost all cases (literally ~99% of a typical simulation) are determined via the following criteria:
#STACKED if max_angle_deg < kink_thres_angle_deg (small max angle)
#KINKED at nick site if max_angle_deg >= kink_thres_angle_deg and nuc_1_nuc_2_angle_deg > kinked_away_thres_angle_deg (both large max angle and large nuc_1 nuc_2 angle) 
#KINKED_AWAY from nick site if max_angle_deg >= kink_thres_angle_deg and nuc_1_nuc_2_angle_deg <= kinked_away_thres_angle_deg (large max angle but small nuc_1 nuc_2 angle)
#and fraying events are just done via a distance check (if any two native pairs close to the nick site are far enough apart, frayed).

#More detailed description:
#Currently sets the maximum threshold angle for minicircles to be >=40 degrees to be kinked.
#A structure is considered to be kinked away from the nick if the angle between nuc_1 and nuc_2 (the two sites immediately opposite of the nick)
#is *small* (currently <=30 degrees) while the max angle is larger than the threshold angle.
#There is also a small edge case where the ends are parallel but disjointed like a broken bone; there is a distance check currently implemented that checks the nucleotides directly next to the nick to see if the backbone spheres are sufficiently far apart.
#If it is deemed kinked, it tries to differentiate kinked from kinked away by looking at difference between nuc_1 and nuc_2 angle from the max angle (generally, nuc_1 and nuc_2 angle is the same as max angle when kinked at the nick site even in this edge case but admittedly not perfect).
#A structure is considered to be stacked if the max angle is <40 degrees and passes the distance check (i.e., backbones are <=2 sim units apart).
#These angles were determined via analyzing histograms and many configurations to see what the distribution looked like and what angles reflected stacking vs. kinking (and if it was kinked at nick site or not).
#Distances chosen to best reflect what I saw viewing many configurations as well.

#A structure is considered to also be frayed if any native pairs near the nick site are distanced >=.7 sim units currently.
#It should be noted that kinking away from the nick always seems to inherently be associated with a frayed motif. As such, you can subtract the total frayed fraction from the kinked away fraction
#to estimate the number of frays that do not occur from this kinking away motif. Kinking away fraying events tend to be longer lasting than these sort of transient frays from what I've seen.

#WORD OF CAUTION: Wrapping of indices should be properly handled in all cases if provide nuc_1 and nuc_2 for a minicircle with no
#gap/bulge/handle defects as the program was originally designed for. When these defects are present, I'm not 100% positive that there aren't cases where it will improperly handle indices and I don't know if there
#is a good way to handle these except hard-coding the relevant pairs list instead of calling get_angles_list_deg or more sophisticated way to determine them than currently implemented.
#In most cases it'll still grab them correctly for these, but be wary and check the program's printed pairs in oxview to make sure it is checking pairs of bases that you want it to.

#Usage: python3 [comma-separated list of input files (input1,input2,etc.)] [main output ("times file" -- prints phase at each time step) directory] [nuc_1 nuc_2 angles output directory] [max angles output directory] [list of nucleotides near the nick site (use program get_nick_comp_list.py to try to generate automatically... but must check these indices manually to be sure)]

from oxDNA_analysis_tools.UTILS.data_structures import *
from oxDNA_analysis_tools.UTILS.RyeReader import *
import sys
import numpy as np
import os
from statistics import stdev

#Command line inputs.
input_files_unsplit = sys.argv[1]
input_files = input_files_unsplit.split(",")
times_out_txt = open(sys.argv[2], "w")
angles_out_txt = open(sys.argv[3], "w")
max_angles_out_txt = open(sys.argv[4], "w")
comp_list_txt = open(sys.argv[5])

#PROGRAM PARAMETERS:
#Toggle the max range about nick site to check for kinking (max_shift = 2 checks 6-base window)
#Number of intervals sets how many intervals to divide simulation into for determining error bar (if choose 4 intervals, divides simulation into quarters and determines fractions in each quarter to estimate an STDEV error).
#Toggle kink_thres_angle_deg to change the *minimum* angle threshold that counts as kinked (kink_thres_angle_deg = 90 is what is used in paper!)
#Toggle kinked_away_thres_angle_deg to change the *maximum* angle threshold between nuc_1 and nuc_2 that counts as kinked away (this is new; see program description at top)
#Toggle frayed_distance_thres to change the minimum distance (in sim units) between native pairs near nick that counts as frayed.
#Toggle stacked_distance_thres to change the minimum distance (in sim units) between the two ends closest to nick site to count as stacked.
#Toggle max_angle_diff_thres for a niche case: If determine a false positive on being stacked, check difference between max angle and nuc_1 nuc_2 angle. This sets maximum difference that counts as kinked at nick site (vs. kinked away). Generally, the difference is very often zero if kinked at the nick site but not always.
max_shift = 2
num_intervals = 4
kink_thres_angle_deg = 40
kinked_away_thres_angle_deg = 30
frayed_distance_thres = .7
stacked_distance_thres = 2
max_angle_diff_thres = 12

#Parse comp_list_txt input. 
comp_pairs_1D_str = comp_list_txt.readline().split() #raw read-in; reads in as strings and as a 1D list
comp_pairs_list = []
for i in range(int(len(comp_pairs_1D_str)/2)):
    comp_pairs_list.append([int(comp_pairs_1D_str[2*i]),int(comp_pairs_1D_str[2*i + 1])])
nuc_1 = comp_pairs_list[0][0] #lower ID of the two sites opposite of the nick usually (second ID from 5' to 3' in the case of wrapping)
comp_1 = comp_pairs_list[0][1]
nuc_2 = comp_pairs_list[1][0] #higher ID of the two sites opposite of the nick usually; usually nuc_1 + 1 (first ID from 5' to 3' in the case of wrapping)
comp_2 = comp_pairs_list[1][1]

#Returns a list of consecutive bases whose angles between their base-normal vectors are checked for kinking in main program
#Should account wrapping around any of the indices to the corresponding max/min index of the current strand
#Builds the list in such a way that the pairs closest to the nick site are earlier in the list.
def get_rel_pairs_list(nuc_1, nuc_2, max_shift, system):
    #Determine internal strand ID for nuc_1 (and nuc_2, *assuming* user provided correctly indices on same strand)
    for strand in system.strands:
        for monomer in strand.monomers:
            if nuc_1 == monomer.id:
                strand_index = strand.id - 1
                
    #For simplicity, I'm constructing my own list of all indices on the current strand.
    ID_list = []
    for monomer in system.strands[strand_index].monomers:
        ID_list.append(monomer.id)
    
    #Determine max/min index of strand for handling wrapping of indices when increment/decrement past min/max
    min_nuc_index = min(ID_list)
    max_nuc_index = max(ID_list) 
    
    rel_pairs = [[nuc_1,nuc_2]] #nuc_1 and nuc_2 assumed as first pair and build rest of list from here
    
    for shift in range(0,max_shift):
        nuc_a_pos_shift = nuc_2 + shift #first nucleotide in positive direction
        nuc_b_pos_shift = nuc_2 + shift + 1 #and second nucleotide in positive direction to be checked pair-wise with first for stacking

        if nuc_a_pos_shift > max_nuc_index:
            diff = nuc_a_pos_shift - max_nuc_index #amount past max index
            nuc_a_pos_shift = min_nuc_index + diff - 1 #newly wrapped index
        if nuc_b_pos_shift > max_nuc_index:
            diff = nuc_b_pos_shift - max_nuc_index
            nuc_b_pos_shift = min_nuc_index + diff - 1

        rel_pairs.append([nuc_a_pos_shift,nuc_b_pos_shift])

        nuc_a_neg_shift = nuc_1 - shift - 1
        nuc_b_neg_shift = nuc_1 - shift

        if nuc_a_neg_shift < min_nuc_index:
            diff = min_nuc_index - nuc_a_neg_shift #amount past min index
            nuc_a_neg_shift = max_nuc_index - diff + 1 #newly wrapped index
        if nuc_b_neg_shift < min_nuc_index:
            diff = min_nuc_index - nuc_b_neg_shift
            nuc_b_neg_shift = max_nuc_index - diff + 1

        rel_pairs.append([nuc_a_neg_shift,nuc_b_neg_shift])
    
    return rel_pairs

#Computes the dot product between base-normal vectors of the provided nucleotide IDs and returns the angle in degrees.
#Should note that a3 (base-normal vector in oxDNA) is normalized already but has round-off error when printed to trajectory.
def get_angle_deg(nuc_1, nuc_2, conf):
    rad_to_deg_conv = 180/np.pi
    dot_prod = np.clip(conf.a3s[nuc_1].dot(conf.a3s[nuc_2]),-1.0,1.0) #np.clip handles round-off errors due to binary math/trajectory rounding
    angle_rad = np.arccos(dot_prod)
    angle_deg = angle_rad*rad_to_deg_conv
    
    return angle_deg

#Returns the list of all the base-normal angles between all of the relevant pairs from get_rel_pairs_list. 
def get_angles_list_deg(rel_pairs, conf):
    angles_list_deg = []
    
    for pair in rel_pairs:
        angles_list_deg.append(get_angle_deg(pair[0],pair[1],conf))
    
    return angles_list_deg

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

#see https://dna.physics.ox.ac.uk/index.php/Documentation#Geometry_of_the_Model for below functions
def get_backbone_vector(nuc,conf,interaction_type):
    cm_pos = conf.positions[nuc]
    _a1 = conf.a1s[nuc]
    if interaction_type == "DNA2":
        _a3 = conf.a3s[nuc]
        _a2 = np.cross(_a3,_a1)
        return cm_pos - .34*_a1 + .3408*_a2
    else:
        return cm_pos - .4*_a1   

#same coordinates for both oxDNA1 and oxDNA2
def get_H_bond_site(nuc, conf):
    cm_pos =  conf.positions[nuc]
    _a1 = conf.a1s[nuc]
    return cm_pos + .4*_a1

for input_file in input_files: #simple loop over all trajectory files if needed
    #Create TopInfo and TrajInfo objects that are used to read in trajectory.
    traj_file = get_input_parameter(input_file, "trajectory_file")
    top_file = get_input_parameter(input_file, "topology")
    top_info, traj_info = describe(top_file, traj_file) 
    
    #Read in the topology (creates a system with same hierarchy of System -> Strands -> Monomer in same way as old analysis code)
    system,monomers = strand_describe(top_file) #Creates System and list of Monomers; see data_structures.py for instance variables and believe monomers list is indexed exactly the same as oxview indexes nucleotides
    
    #Find number of configurations in current simulation (assumes a linear time scale).
    #'steps' holds the total number of time steps of the simulation (the 'steps' key in input file), and the print_conf_interval is directly grabbed from input file as well.
    steps = int(float(get_num_steps(input_file)))
    print_conf_interval = int(float(get_input_parameter(input_file,"print_conf_interval")))
    num_confs = steps/print_conf_interval
    
    #########################################################################################################################
    #NEXT SEVERAL LINES ARE MOSTLY PRINTED OUT DIAGNOSTICS -- useful to read to understand some of the internals of program #
    #########################################################################################################################
    #I choose to print out some information about what is checked. Some of the math described here is important for the actual analysis.
    N_nuc = top_info.nbases
    
    #Determine internal strand ID for nuc_1 (and nuc_2, assuming user provided correctly indices on same strand)
    for strand in system.strands:
        for monomer in strand.monomers:
            if nuc_1 == monomer.id:
                strand_index = strand.id - 1 #for some reason, strand indices now do not count up from 0!
    
    #For simplicity, I'm constructing my own list of all indices on the current strand.
    ID_list = []
    for monomer in system.strands[strand_index].monomers:
        ID_list.append(monomer.id)
    
    min_nuc_index = min(ID_list)
    max_nuc_index = max(ID_list) 
    
    print("For this version of the code, a configuration is deemed \"KINKED\" if max angle between any two consecutive base-normal vectors in window is >=" + str(kink_thres_angle_deg) + "deg.\nIt also prints if it is \"KINKED_AWAY\" from the nick site (same max angle condition but with angle between nucleotides directly across the nick having angle <=" + str(kinked_away_thres_angle_deg) + "deg).\nIt is \"STACKED\" otherwise.\nHowever, there is a small edge case where the ends are nearly parallel but disjointed and kinked. The program attempts to remedy this by doing a distance check of the backbone spheres of the two nt closest to the nick site.\nFinally, any time step is determined to be frayed if any native pair near the nick site have H-bond sites sufficiently far apart.")
    print("The header of the program also describes this more.")
    
    print("Number of nucleotides:",N_nuc)
    print("Current strand ID:", strand_index)
    print("Min nucleotide index on strand:", min_nuc_index)
    print("Max nucleotide index on strand:", max_nuc_index)
    
    rel_pairs_list = get_rel_pairs_list(nuc_1, nuc_2, max_shift, system) #contains all of the relevant pairs whose base-normal angles are checked for kinking
    print("Below are the checked pairs.")
    print("Original IDs:", rel_pairs_list[0][0],rel_pairs_list[0][1])
    print("Positive shift pairs of nucleotide IDs:")
    for shift in range(0,max_shift):
        print(rel_pairs_list[2*shift + 1][0],rel_pairs_list[2*shift + 1][1])

    print("Negative shift pairs of nucleotide IDs:")
    for shift in range(0,max_shift):
        print(rel_pairs_list[2*shift + 2][0],rel_pairs_list[2*shift + 2][1])
    
    print("Raw comp_pairs_list:",comp_pairs_list)
    
    #Determine if oxDNA1 or oxDNA2 and print diagnostics.
    interaction_type = get_input_parameter(input_file, "interaction_type")
    if interaction_type == "DNA2":
        print("Auto-detected oxDNA2 from input file.")
    elif interaction_type == "DNA":
        print("Auto-detected oxDNA1 from input file.")
    else:
        print("Could not determine from input file interaction_type. Defaulting to oxDNA1.")
        
    #####################################################################################################################################################################################################################
    #MAIN PROGRAM. Loop over all configurations and determine which phase the each configuration is in (criteria described in header), and output the max angle and nuc_1 nuc_2 angles to separate text files currently.#     
    #####################################################################################################################################################################################################################
    
    #Set-up for statistics. 
    #rel_times is used to break the simulation into num_intervals intervals. The program determines the fractions in each of these intervals
    #and computes a STDEV in the end to propagate an error bar for each of the total simulation fractions.
    rel_times = [i * int(((steps/num_intervals)//print_conf_interval) * print_conf_interval) for i in range(1,num_intervals)] #algebra here is an attempt to make times stored here actual printed time steps
    rel_times.append(steps) #last rel_time is last time step of simulation
   
    kinked_away_count = [0] * num_intervals
    kinked_count = [0] * num_intervals
    frayed_count = [0] * num_intervals
    configs_count = [0] * num_intervals
    
    increment = 10000

    #This loop structure could be made more readable... read in an increment number of confs at a time, and determine the phase of each conf in this increment. Then, loop to next increment. Don't read all configurations at once in case of memory issues (**may not actually be an issue? another bug was giving segfaults/memory issues but not sure if read-ins of large # of configs can be problematic still).
    for i in range(int(num_confs//increment)):
        #Read in an increment number of configurations
        cur_confs = get_confs(top_info, traj_info, i*increment, increment)
        
        #Loop over every conf in this batch:
        for conf in cur_confs:
            #print("Time step:", conf.time)
            
            #Inbox every configuration (so that PBCs are not an issue for distance checks)
            conf = inbox(conf)
            
            #Determine which interval of the simulation we are in for statistics.
            for j in range(len(rel_times)):
                if conf.time <= rel_times[j]:
                    index = j
                    break
            
            #Compute the angles between base-normals of nuc_1 and nuc_2 and max of all pairs checked.
            angles_list_deg = get_angles_list_deg(rel_pairs_list, conf)
            max_angle_deg = max(angles_list_deg)
            nuc_1_nuc_2_angle_deg = angles_list_deg[0]
            
            #Check for normal kinking vs kinking away from the nick:
            if max_angle_deg >= kink_thres_angle_deg:
                if nuc_1_nuc_2_angle_deg <= kinked_away_thres_angle_deg:
                    times_out_txt.write(str(conf.time) + " KINKED_AWAY") #mark when kinked away from nick
                    kinked_away_count[index] += 1
                else:
                    times_out_txt.write(str(conf.time) + " KINKED") #write the current time of the configuration when kinked
                    kinked_count[index] += 1
            
            else: #edge case on what otherwise would be deemed stacked: possible to have kinking events where both ends point parallel are disjointed
                #Check if there is a false positive on stacking via checking distances of nucleotides directly at nick site:
                backbone_distance = np.linalg.norm(get_backbone_vector(comp_1, conf, interaction_type) - get_backbone_vector(comp_2, conf, interaction_type))
                
                #If past stacked distance threshold, do a small check to see if think is kinked at nick site or kinked away.
                #Often, when kinked at nick site the max angle is equal to nuc_1 and nuc_2 angle or very small.
                #When kinked away, it tends to be larger. As with the others, there is a degree of arbitrariness in what you set the threshold
                #but I've found that this works for many many configurations.
                #In any case, this is less surefire than the first check for kinked/kinked away... but it is better than saying stacked when it isn't for our purposes!
                if backbone_distance >= stacked_distance_thres: #see if ends of nick site are far apart
                    if max_angle_deg - nuc_1_nuc_2_angle_deg <= max_angle_diff_thres: 
                        times_out_txt.write(str(conf.time) + " KINKED2 " + str(round(max_angle_deg - nuc_1_nuc_2_angle_deg,2)))
                        kinked_count[index] += 1
                    else: #assume kinked_away
                        times_out_txt.write(str(conf.time) + " KINKED_AWAY2 " + str(round(max_angle_deg - nuc_1_nuc_2_angle_deg,2)))
                        kinked_away_count[index] += 1
                else: #stacked if ends of nick site are sufficiently close
                    times_out_txt.write(str(conf.time) + " STACKED")
                
                            
            #Check in addition for any fraying events near the nick site. Just checks if H-bond sites are sufficiently far apart.
            frayed = False
            for pair in comp_pairs_list:            
                base_distance = np.linalg.norm(get_H_bond_site(pair[0], conf) - get_H_bond_site(pair[1],conf))
                if base_distance >= frayed_distance_thres:
                    frayed = True
                    break
            if frayed:    
                times_out_txt.write(" || FRAYED " + str(round(base_distance,2)) + "\n")
                frayed_count[index] += 1
            else:
                #times_out_txt.write(" || NOT FRAYED " + str(round(base_distance,2)) + "\n")
                times_out_txt.write("\n") #setup for next line of text file          
                    
            #Write out both of the pertinent angles checked:
            angles_out_txt.write(str(conf.time) + " " + str(nuc_1_nuc_2_angle_deg) + "\n")
            max_angles_out_txt.write(str(conf.time) + " " + str(max_angle_deg) + "\n")
            
            #Add to configuration count (this turns out to count correctly even though count when !system at end; we didn't count initial configuration)
            configs_count[index] += 1

remove_chars = len(os.linesep)
try:
    times_out_txt.truncate(times_out_txt.tell() - remove_chars) #removes trailing newline in text file
    print("Main output with phase classifications at each time step was written to", sys.argv[2] + ".")
except: #if above fails, nothing was printed (as no trailing newline)
    print("No configurations were found in trajectory file. Make sure provided correct directory (else, unsure why this branch would get hit).")
times_out_txt.close()

remove_chars = len(os.linesep)
try:
    angles_out_txt.truncate(angles_out_txt.tell() - remove_chars)
    print("Simulation time and angle between normal vectors of immediate two bases on in-tact strand printed to", sys.argv[3] + ".")
except:
    print("Failed to properly print base-normal angles... likely due to previous problem.")
angles_out_txt.close()

remove_chars = len(os.linesep)
try:
    max_angles_out_txt.truncate(max_angles_out_txt.tell() - remove_chars)
    print("Simulation time and max angle between checked bases on in-tact strand printed to", sys.argv[4] + ".")
except:
    print("Failed to properly print base-normal angles... likely due to previous problem.")
max_angles_out_txt.close()

print("Raw relevant times:", rel_times)
print("Raw kinked_fractions:", [i / j for i,j in zip(kinked_count,configs_count)])
print("Raw kinked_away_fractions:", [i / j for i,j in zip(kinked_away_count,configs_count)])
print("Raw stacked_fractions:", [1- i/k - j/k for i,j,k in zip(kinked_count,kinked_away_count,configs_count)])
print("Raw frayed_fractions:", [i / j for i,j in zip(frayed_count,configs_count)])
print("Raw tran_fray_fractions:", [i/k - j/k for i,j,k in zip(frayed_count, kinked_away_count,configs_count)])
print("Raw configs_count:", configs_count)

print()
print("STATISTICS:")

#Print out main statistics. 
#Results for total fractions are computed over entirety of simulation, and error determined via STDEV of fractions across ~even divisions of the simulation (if this is unclear, read the last print statement in case it's any clearer).
print("Estimate of fraction of KINKED configs:", sum(kinked_count)/sum(configs_count), "±", stdev([i / j for i,j in zip(kinked_count,configs_count)]))
print("Estimate of fraction of KINKED_AWAY configs:", sum(kinked_away_count)/sum(configs_count), "±", stdev([i / j for i,j in zip(kinked_away_count,configs_count)]))
print("Estimate of fraction of STACKED configs:", 1 - sum(kinked_count)/sum(configs_count) - sum(kinked_away_count)/sum(configs_count), "±", stdev([1 - i/k - j/k for i,j,k in zip(kinked_count, kinked_away_count, configs_count)]))
print("Estimate of fraction of configs with a FRAYED event (independent of above classifications):", sum(frayed_count)/sum(configs_count), "±", stdev([i / j for i,j in zip(frayed_count,configs_count)]))
print("Estimate of fraction of configs FRAYED but **NOT** KINKED_AWAY:", sum(frayed_count)/sum(configs_count) - sum(kinked_away_count)/sum(configs_count), "±", stdev([i/k - j/k for i,j,k in zip(frayed_count, kinked_away_count,configs_count)]))
print("Total number of configs checked:", sum(configs_count))
print("To determine reported error, simulation was split into", num_intervals, "approx. equal intervals and used the STDEV of the corresponding fractions across each of these intervals.")