#!/usr/bin/env python3

########################################################################################################################
# TITLE:       kink_finder_precise.py
# AUTHOR:      Daniel Hollis
# CREATED:     August 11, 2023
# LAST EDITED: July 17, 2024
# DESCRIPTION: Takes as inputs the trajectory and topology file of a DNA minicircle simulation, as well as nick location
#              (if applicable). Prints to the console a list of all the configurations in which the energetic criterion
#              for a kink is met, along with the locations on the circle of these kinking events, identified by two
#              nucleotides, one on either side of the region of bond disruption. If the kink is at the nick, this is
#              reported in the console. Writes final statistics to the stat file. Assumes the circle is entirely of
#              dsDNA with no abnormalities besides possibly one nick. Assumes nucleotide ID numbers are arranged such
#              that one strand comprises 0-->strand_length-1 and the other comprises strand_length-->total_#_bases-1.
#
#              Approximately, this energetic criterion corresponds to two instances of broken stacking on either strand,
#              either adjacent to each other or separated by a stretch of broken hydrogen bonds. This specific
#              combination of disrupted bonds compromises the structural integrity of the circle and allows for the
#              possibility of sharp bending. Consult the code or the author for a more detailed definition.
#
#              Note: To avoid ambiguity in the case of large kinked regions (>= strand_length/2), boundaries are given
#              in a specific order. If you look at the strand containing nucleotides 0-->strand_length-1, such that the
#              ID numbers run in a clockwise direction, the kinked region will extend between the first and second
#              boundaries in a clockwise direction. A kink with one boundary represents the rare case where the circle
#              is kinked, but has only one intact base pair.
#
#              A suggested improvement to this code would be adding a physical criterion, similar to the one found in
#              Harrison, Romano et al. 2019, in addition to the energetic criterion implemented here. This would check
#              for sharp bending in the region between the boundaries and report a kink only if the bending exceeded
#              some threshold angle.
#
#              As of 7/17/24, the program should also be applicable to a linear duplex without defects. The program
#              will auto-detect if the system appears to be linear or circular by checking if any circular strands
#              are present in the system, and adjust the logic slightly to be appropriate for linear duplexes.
#
#              Thanks to Thomas Swope for implementing the stacking check, and Rachel Bricker for patching some index
#              errors.
#
########################################################################################################################

#This program relies on the old version of oxDNA analysis tools (OAT) which is deprecated since circa May 2022.
#Point to your copy of this version (or older) below.
#import sys
#sys.path.insert(0, "<your path to the aforementioned version of oxDNA analysis tools>")
import argparse
#readers_modified is a version of readers from the OAT UTILS, modified to be compatible with base_modified--see below.
from oxDNA_analysis_tools.UTILS.readers_modified import LorenzoReader2, get_input_parameter
from oxDNA_analysis_tools.UTILS import geom
#Below is a version of base.py from the OAT UTILS folder which, among other things, has been modified to support
#The reading of stacking interactions.
from oxDNA_analysis_tools.UTILS import base_modified
from oxDNA_analysis_tools.output_bonds import output_bonds

# run system checks
from oxDNA_analysis_tools.config import check_dependencies
check_dependencies(["python", "numpy"])

#returns true if nuc_1 and nuc_2 are stacked. Assumes nuc_1 and nuc_2 are adjacent and connected by backbone.
#assumes read_stack_bonds_output_bonds already ran to map output of output_bonds to system object
def is_stacked(system, ID_1, ID_2):
    stacked = False

    # for i in system._nucleotides[ID_1].interactions:
    #     print(i)
    if ID_2 in system._nucleotides[ID_1].stack_interactions:
        stacked = True

    return (stacked)

#returns true if nuc_1 and nuc_2 are coaxially stacked. Assumes nuc_1 and nuc_2 are adjacent but not connected by backbone.
#assumes read_cxstack_bonds_output_bonds already ran to map output of output_bonds to system object.
def is_cxstacked(system, ID_1, ID_2):
    cxstacked = False

    # for i in system._nucleotides[ID_1].interactions:
    #     print(i)
    if ID_2 in system._nucleotides[ID_1].cxstack_interactions:
        cxstacked = True

    return (cxstacked)
    
#def is_cr_stacked(system, ID_1, ID_2):
#    crstacked = False
#
#    if ID_2 in system._nucleotides[ID_1].crstack_interactions:
#        crstacked = True
#
#    return (crstacked)

#returns true if nuc_1 and nuc_2 are stacked, either normally or coaxially.
#assumes read_stack_bonds_output_bonds and read_cx stack_bonds_output_bonds already ran to map output of output_bonds
#to system object
def is_adj_stacked(system, ID_1, ID_2):
    stacked = False

    if(is_stacked(system, ID_1, ID_2)):
        stacked = True
    elif(is_cxstacked(system, ID_1, ID_2)):
        stacked = True

    return (stacked)

#returns true if nucleotides with ids ID_1 and ID_2 are a base pair.
#assumes read_H_bonds_output_bonds already ran to map output of output_bonds to system object
def is_paired(system, ID_1, ID_2):
    bonded = False
    if ID_2 in system._nucleotides[ID_1].H_interactions:
        bonded = True

    return (bonded)

#returns true if there are no broken base pairs in the DNA ring system
def is_hybridized(system, num_nuc):
    hybridized = True

    #this searching algorithm could be made faster
    for i in range(num_nuc):
        pair_found = False
        for j in range(0, i):
            if(is_paired(system, i, j)):
                pair_found = True
                break
        if(not(pair_found)):
            for k in range(i+1, num_nuc):
                if(is_paired(system, i, k)):
                    pair_found = True
                    break

        if(not(pair_found)):
            hybridized = False
            break
    return(hybridized)

# returns 2d array of nucleotide IDs representing two complementary strands (coined ref_ring)
# assumes read_H_bonds_output_bonds already ran to map output of output_bonds to system object
# assumes system has no broken base pairs
def get_base_pairs(ref_sys, num_nuc, strand_length):
    ref_ring = [[], []]
    for i in range(strand_length):
        ref_ring[0].append(i)
        for j in range(strand_length, num_nuc):
            if(is_paired(system, i, j)):
                ref_ring[1].append(j)
                break
    return ref_ring

#Returns None if no kinks are present in system; otherwise, returns pairs of boundaries for each kinked region in the
# circle (going from ID 0-->strand_length).
#These boundary nucleotides will always be on strand 0, and will come in the form of 2x1 arrays bundled into a 2D array.
#Assumes read_H_bonds_output_bonds, read_stack_bonds_output_bonds, read_cxstack_bonds_output_bonds already ran to map
#output of output_bonds to system object.
def find_kinks(system, ref_ring, num_nuc, strand_length, is_circular):
    all_bounds = []
    #We need to start our search with a hybridized pair.
    start_nuc = first_base_pair(system, ref_ring, strand_length)
    #Intact hydrogen bond which is currently the subject of focus.
    current_nuc = start_nuc
	#Find succeeding intact hydrogen bond.
	#It's possible this will return current_nuc if there's only one hybridized base pair.
    next_nuc = next_base_pair(system, ref_ring, strand_length, current_nuc)
    while(next_nuc != start_nuc):
        #Determine boundaries of kink in the region between current_nuc and next_nuc. If no kink, bounds will equal None.
        bounds = is_kink(system, ref_ring, num_nuc, strand_length, current_nuc, next_nuc)
        if(bounds != None):
            all_bounds.append(bounds)
        current_nuc = next_nuc
        next_nuc = next_base_pair(system, ref_ring, strand_length, current_nuc)
    if(is_circular):
        bounds = is_kink(system, ref_ring, num_nuc, strand_length, current_nuc, next_nuc)
        if(bounds != None):
            all_bounds.append(bounds)
    
    if(len(all_bounds) == 0):
        return None
    else:
        return(all_bounds)
    
#Returns the nucleotide which is located on strand 0 and is a part of the first hybridized base pair in the circle when
# going in the direction of 0-->strand_length. Used by find_kinks()
def first_base_pair(system, ref_ring, strand_length):
    for i in range(strand_length):
        if(is_paired(system, i, ref_ring[1][i])):
            first_nuc = i
            break
    else: #if no break
        raise Exception("No intact hydrogen bonds at t = %s. Ambiguous case." % str(system._time))
    return(first_nuc)

#Finds the hybridized base pair which succeeds the base pair containing nuc. Goes from ID 0-->strand_length on strand 0.
#Nuc and next_nuc on strand 0. Will return nuc argument if there's only one hybridized pair. Used by find_kinks().
def next_base_pair(system, ref_ring, strand_length, nuc):
    for i in range(nuc+1, strand_length):
        if(is_paired(system, i, ref_ring[1][i])):
            next_nuc = i
            return(next_nuc) #if found, return value and exit function

    for j in range(0, nuc):
        if(is_paired(system, j, ref_ring[1][j])):
            next_nuc = j
            return(next_nuc)

    next_nuc = nuc #special case
    return(next_nuc)

#If there is a kink located between the two hybridized base pairs, represented by their nucleotides which are located
#on strand 0, this function returns the boundaries of this kinked region. Works in the direction of 0-->strand_length.
#Boundaries are located on strand 0. Returns None if no kink is found. In case where nuc1=nuc2, region will span the
#entire circle. Used by find_kinks().
#
#This function accomplishes two things at once: (1.) It tests whether both strands have at least one instance of broken
#stacking between base pair containing nuc1 and base pair containing nuc2. This is a condition for a kink. And (2.) it
#determines bounds for the *potential* kinked region, to return in the event that the region is indeed kinked. If the
#region is not kinked, these bounds will be irrelevant and disregarded. These two tasks are combined here to increase
#efficiency, since they both involve the reading of the same group of stacking interactions.
def is_kink(system, ref_ring, num_nuc, strand_length, nuc1, nuc2):
    bound1 = -1 #bounds will remain at -1 unless potential bounds are located, i.e. broken stacking is identified
    bound2 = -1
    strand0broken = False
    strand1broken = False
    both_broken = False

    #begin by searching for first bound while simultaneously checking if stacking condition is met.
    if(nuc2 > nuc1):
        endsearch = nuc2-1 #final nucleotide to check (code will check stacking between endsearch and endsearch+1)
    else: #if nuc2 <= nuc1
        endsearch = strand_length-2

    for i in range(nuc1, endsearch+1):
        if(not is_adj_stacked(system, i, i+1)):
            strand0broken = True
        if(not is_adj_stacked(system, ref_ring[1][i], ref_ring[1][i+1])):
            strand1broken = True
        if(strand0broken or strand1broken):
            #bound1 corresponds to the first case of broken stacking when moving from nuc1-->nuc2, so if bound1 is still
            #at its default value and we find an example of broken stacking, set i to bound1.
            if(bound1 == -1):
                bound1 = i
            #If both strands have at least one case of broken stacking, the region is kinked and bound1 has already been
            #identified. It is pointless to keep going.
            if(strand0broken and strand1broken):
                both_broken = True
                break

    #if nuc2 <= nuc1, test region between 0 and nuc2 which was left out when we searched nuc1-->strand_length-2
    while(nuc2 <= nuc1 and not both_broken): #this while loop functions as a "breakable if statement"
        #first, let's test the edge case
        if(not is_adj_stacked(system, strand_length-1, 0)):
            strand0broken = True
        if(not is_adj_stacked(system, ref_ring[1][strand_length-1], ref_ring[1][0])):
            strand1broken = True
        if(strand0broken or strand1broken):
            if(bound1 == -1):
                bound1 = strand_length-1
            if(strand0broken and strand1broken):
                both_broken = True
                break

        #now we'll test the aforementioned region
        for j in range(0, nuc2):
            if(not is_adj_stacked(system, j, j+1)):
                strand0broken = True
            if(not is_adj_stacked(system, ref_ring[1][j], ref_ring[1][j+1])):
                strand1broken = True
            if(strand0broken or strand1broken):
                if(bound1 == -1):
                    bound1 = j
                if(strand0broken and strand1broken):
                    both_broken = True
                    break

        break #so the while loop cannot repeat

    #At this point we're still not sure whether both_broken is True or False
    if(not both_broken):
        return None #stacking condition is not met, region is not kinked; bounds are irrelevant.
    #else:
    #stacking condition is met and region is kinked; now we must find second bound.
    #bound2 is located by going in the opposite direction, from nuc2-->nuc1.
    if(nuc2 > nuc1):
        endsearch2 = nuc1 + 1  #final nucleotide to check. endsearch2 to avoid confusion with endsearch.
    else:  #if nuc2 <= nuc1
        endsearch2 = 1

    for k in range(nuc2, endsearch2-1, -1):
        if(not is_adj_stacked(system, k, k-1) or not is_adj_stacked(system, ref_ring[1][k], ref_ring[1][k-1])):
            bound2 = k
            break

    #if nuc2 <= nuc1 and bound2 is not found, test region which was left out when we searched between nuc2 and 1
    while(nuc2 <= nuc1 and bound2 == -1): #this while loop serves as a "breakable if statement"
        #first, test edge case
        if(not is_adj_stacked(system, 0, strand_length-1) or not is_adj_stacked(system, ref_ring[1][0], ref_ring[1][strand_length-1])):
            bound2 = 0
            break

        #now, test the aforementioned region
        for l in range(strand_length-1, nuc1, -1):
            if(not is_adj_stacked(system, l, l-1) or not is_adj_stacked(system, ref_ring[1][l], ref_ring[1][l-1])):
                bound2 = l
                break

        break #so the while loop cannot repeat

    #bound2 must have been found by now, since we know there is at least one instance of broken stacking somewhere.
    if(bound2 == -1):
        raise Exception("bound2 not found; there's a bug in the program!")

    #bound1==bound2 is possible if there's only one hydrogen bond in the circle. This is unlikely, but be aware.
    return[bound1, bound2]

#checks for all strands if all IDs are stacked/coaxially stacked with closest neighbor ID_2 on strand
def is_all_intact(system, ref_ring, is_circular):
    flag = "" #this flag will determine if it is intact; if it is empty, no broken stacking was found -> intact (and empty string will evaluate falsy)
              #if it is not intact, it will say ADJ if broken stacking between two adjacent pairs or CR if broken cross stacking (and nonempty string will evaluate truthy)
              #these will be concatenated if both are broken and printed to stdout along with UNSTACKED
    
    adj_intact = True
    cr_intact = True
              
    broken_pairs = [] #list holding any broken pairs found
    for strand_ID in range(len(ref_ring)):
        ID_list = ref_ring[strand_ID] #list of all IDs on current strand
            
        #Check for any broken adjacent stacking on either strand.
        for ID, ID_2 in zip(ID_list, ID_list[1:]): #check every nucleotide on current strand and it's adjacent pair for any broken stacking
            if not is_adj_stacked(system, ID, ID_2):
                broken_pairs.append([ID,ID_2])
                adj_intact = False #if any adjacent nucleotides have broken stacking, not intact
        
        if is_circular: #if circular, check also the wrapped indices
            ID = ID_list[-1]
            ID_2 = ID_list[0]
            if not is_adj_stacked(system, ID, ID_2):
                broken_pairs.append([ID,ID_2])
                adj_intact = False #if any adjacent nucleotides have broken stacking, not intact

    # #Check for any broken cross-stacking interactions. Do this only on one strand as otherwise double-check pairs (i.e., [82,160] and [160,82] will both be reported if check both strands)
    # #Note: This was written before the program was adjusted to handle linear DNA, assuming it is functional at all. Rework/recheck if want to implement cross-stacking analysis.
    # comp_ID_list = ref_ring[0] #ID_list from last iteration of loop is second strand index; comp_ID_list will be the list of complementary indices to ID_list and is parallel to ID_list

    # for i in range(len(ID_list)):
        # #Search the neighbor to comp_ID in either direction to see if cross-stacking is broken
        # ID = ID_list[i]
        # comp_ID = comp_ID_list[i]
        # upper_comp_ID = comp_ID + 1
        # if upper_comp_ID > max(comp_ID_list):
            # upper_comp_ID = min(comp_ID_list)
        # if not is_cr_stacked(system,ID,upper_comp_ID):
            # broken_pairs.append([ID,upper_comp_ID])
            # cr_intact = False
        # lower_comp_ID = comp_ID - 1
        # if lower_comp_ID < min(comp_ID_list):
            # lower_comp_ID = max(comp_ID_list)
        # if not is_cr_stacked(system,ID,lower_comp_ID):
            # broken_pairs.append([ID,lower_comp_ID])
            # cr_intact = False
    
    #Concatenate flag to contain the type of broken stacking.
    if not adj_intact:
        flag = flag + "ADJ"
        
    if not cr_intact:
        flag = flag + "CR"
    
    return flag, broken_pairs

#Returns True if kinked region is considered "at the nick."
#Assumes the bounds are listed in the order specified in the program description.
#Assumes circle is nicked between nick1 and nick2, and that these are given in the order specified in the help prompt.
#Assumes all the nucleotide IDs are arranged as specified in the description.
def kink_at_nick(system, ref_ring, strand_length, bounds, nick1, nick2):
    bound1 = bounds[0]
    bound2 = bounds[1]
    nick_in_bounds = False

    #If nick is on strand 1, set pair1 and pair2 to complements of nick2 and nick1, respectively. This simplifies our
    #calculations. This may be counter-intuitive, but remember that the IDs of strand 1 run in the direction opposite
    #of strand 2. We just need to know that pair1 is first and pair2 is second when moving clockwise on strand 0.
    if(nick1 >= strand_length): #nick is on strand 1
        pairsfound = 0 #when this equals 2, break loop below
        for i in range(strand_length):
            if(ref_ring[1][i] == nick1):
                pair2 = i
                pairsfound += 1
            if(ref_ring[1][i] == nick2):
                pair1 = i
                pairsfound += 1
            if(pairsfound == 2):
                break
    else: #nick is on strand 0
        pair1 = nick1
        pair2 = nick2

    #Now, let's test!
    #By the way, this could be re-written in a more elegant way. This is just how I worked through this problem.
    if(bound2 > bound1):
        #if pair1 and pair2 (and by extension, nick1 and nick2) are sequentially numbered, not straddling the edge
        if(abs(pair1 - pair2) == 1):
            if(pair1 >= bound1 and pair2 <= bound2):
                nick_in_bounds = True
    elif(bound2 < bound1):
        if(pair1 >= bound1 and pair2 <= bound2):
            nick_in_bounds = True
    else: #if bound2 == bound1, an extremely rare case where a circle has only one intact base pair
        nick_in_bounds = True

    #Below are the two conditions for the kink to be considered "at the nick"
    if(nick_in_bounds):
        if(not is_cxstacked(system, nick1, nick2)):
            return True
    return False #catch-all

########################################################################################################################
#MAIN

# handle arguments entered in the command line
parser = argparse.ArgumentParser(description='''Finds timestamps and locations of kinks in the simulation of a nicked dsDNA circle, non-defective dsDNA circle, or a nondefective linear dsDNA duplex.
                                 For a nicked dsDNA circle, if `nick1` and `nick2` arguments are provided, 
                                 the statistics file will be more detailed differentiating when kinks/unstacking events are at the nick vs. away from the nick.''')
parser.add_argument('input_file', type=str, nargs=1, help="The input file used to run the simulation")
parser.add_argument('traj_file', type=str, nargs=1, help="The trajectory of the simulation")
parser.add_argument('stat_file', type=str, nargs=1, help="The name of the file to which to write the final statistics (i.e., the fractions of stacking events, fractions of kinking events, etc.)")
parser.add_argument('data_file', type=str, nargs=1, help="The name of the file to which to write the data (i.e., the readout from the main console containing the verdicts of either STACKED, KINKED, or UNSTACKED for each configuration)")
parser.add_argument('nick1', type=str, nargs='?', help="(For nicked circles; optional) The ID of the nucleotide which is adjacent to the nick and is first when moving clockwise from 0-->strand_length-1")
parser.add_argument('nick2', type=str, nargs='?', help="(For nicked circles; optional) The ID of the nucleotide which is adjacent to the nick and is second when moving clockwise from 0-->strand_length-1")
parser.add_argument('indices_file', type=str, nargs='?', help='''(optional) A text file with a list of space-separated nucleotide IDs which form base pairs (otherwise, complementary base pairs will attempt 
                    to be automatically created by the program -- in certain cases, a configuration with all nucleotides unmelted will not be found and must be provided manually).''')
args = parser.parse_args()

# parse command line arguments
input_file   = args.input_file[0]
traj_file    = args.traj_file[0]
stat_file    = args.stat_file[0]
reader_times = open(args.data_file[0], "w")
top_file     = get_input_parameter(input_file, "topology")

#test nick location arguments
if(args.nick1 is None and args.nick2 is None):
    print("Nick-adjacent nucleotides not provided. Assuming linear duplex or fully ligated circle.")
    nicked = False
elif((args.nick1).isdigit() and (args.nick1).isdigit()):
    print("Nick-adjacent nucleotides provided. Kinks at the nick will be reported to console.")
    nicked = True
else:
    raise Exception("Nick-adjacent nucleotide IDs must be positive integer values!")
if nicked:
    nick1 = int(args.nick1)
    nick2 = int(args.nick2)

#read in a list of complementary pairs if provided
if(args.indices_file is not None):
    indices_txt   = open(args.indices_file, "r")
else:
    indices_file  = ""

reader = LorenzoReader2(traj_file, top_file)
system = reader._get_system()
num_nuc = system.get_N_Nucleotides()
strand_length = int(num_nuc / 2)
is_circular = False #False if neither strand is circular, True if at least one strand is circular
for strand in system._strands:
    if strand._circular:
        is_circular = True
        print("Found a circular strand in the system. Treating the system as circular.")
if not is_circular:
    print("Found no circular strands in the system. Treating the system as linear.")

# find and save a configuration in the trajectory which has no broken base pairs. This will be used to determine
# which nucleotides are complementary to which, assuming this has not been user-provided.
if(args.indices_file is None):
    while(system != False):
        system.map_nucleotides_to_strands()
        out = output_bonds(input_file, system)
        system.read_H_bonds_output_bonds(out)
        system.inbox()
        if(is_hybridized(system, num_nuc)):
            ref_sys = system
            break
        system = reader._get_system()

    if(system == False):
        raise Exception("ERROR: For this program to work, the trajectory must contain at least one configuration which has no broken base pairs, or complementary indices must be user-provided as an input to the program.")

if(args.indices_file is None):
    ref_ring = get_base_pairs(ref_sys, num_nuc, strand_length)
    # reset reader
    reader = LorenzoReader2(traj_file, top_file)
    system = reader._get_system()

else:
    ref_ring = [[],[]] #two parallel lists representing each strand; parallel entries are complementary pairs
    #parse space-separated list of indices representing complementary pairs (this structure is messy but was chosen so they were usable files also for align.py in oat)
    indices_1D_str = indices_txt.readline().split() #raw read-in; reads in as strings and as a 1D list
    for i in range(int(len(indices_1D_str)/2)):
        ref_ring[0].append(int(indices_1D_str[2*i]))
        ref_ring[1].append(int(indices_1D_str[2*i + 1]))

kink_total_count = 0
kink_nick_count = 0 #number of confs with kinks at nick ONLY
kink_lig_count = 0 #number of confs with kinks away from nick ONLY
kink_both_count = 0 #number of confs with kinks at both the nick and away from the nick
stack_count = 0
unstack_count = 0
unstack_nick_count = 0 #number of confs with broken stacking across nick ONLY
unstack_away_count = 0 #number of confs with broken stacking away from nick ONLY
unstack_both_count = 0 #confs with broken stacking both across and away from the nick
melt_count = 0 #number of confs with all but 1 broken H-bond
total_count = 0

# check every configuration for kinking
while (system != False):
    # run output_bonds, map output to system object
    system.map_nucleotides_to_strands()
    out = output_bonds(input_file, system)
    system.read_stack_bonds_output_bonds(out)
    system.read_cxstack_bonds_output_bonds(out)
    #system.read_crstack_bonds_output_bonds(out)
    system.read_H_bonds_output_bonds(out)
    system.inbox()

    all_bounds_ref = find_kinks(system, ref_ring, num_nuc, strand_length, is_circular)
    if(all_bounds_ref is not None):
        #Consolidate kinks whose boundaries coincide; this criterion is subjective and depends on how you define a kink.
        #You can think of this code as combining kinked regions which share a hydrogen bond boundary. This could be
        #modified to only combine kinks which are separated by one hydrogen bond, two, three, etc.
        #Assumes that bounds in all_bounds are listed consecutively.
        num_bounds = len(all_bounds_ref)
        all_bounds = []
        if(num_bounds > 1):
            skip_iteration = False
            for i in range(num_bounds-1):
                if not skip_iteration:
                    # if right boundary of kink at i coincides with left boundary of kink at i+1, consolidate.
                    if(all_bounds_ref[i][1] == all_bounds_ref[i+1][0]):
                        all_bounds_ref[i][1] = all_bounds_ref[i+1][1]
                        skip_iteration = True
                    all_bounds.append(all_bounds_ref[i])
                else:
                    skip_iteration = False

                # get last element during last iteration if necessary
                if (i == num_bounds-2) and (not skip_iteration):
                    all_bounds.append(all_bounds_ref[i+1])

            #check edge case
            if(all_bounds[-1][1] == all_bounds[0][0]): 
                all_bounds[-1][1] = all_bounds[0][1]
                all_bounds.pop(0)
        else:
            while all_bounds_ref:
                all_bounds.append(all_bounds_ref.pop(0))

        #The below if statement handles a rare case where find_kinks() returns two identical bounds. This may happen
        #if there is only one intact hydrogen bond in the circle.
        is_near_melt = False
        if(num_bounds == 1 and all_bounds[0][0] == all_bounds[0][1]):
            is_near_melt = True

        #Separate kinks into three categories: nick, ligated, and extremely rare case where duplex is mostly melted. Create output strings for each and write to console.
        nick_str    = "NICK"
        ligated_str = "LIG"
        melt_str =    "MELT" 
        for i in range(len(all_bounds)):
            bounds = all_bounds[i]
            if is_near_melt:
                melt_str += " " + str(bounds[0])
            elif(nicked and kink_at_nick(system, ref_ring, strand_length, bounds, nick1, nick2)):
                nick_str += " " + str(bounds)
            else:
                ligated_str += " " + str(bounds) 

        if(melt_str != "MELT"): #catastrophic case where `melt_str != "MELT"`
            melt_count +=1
            print(str(int(system._time)) + " KINKED || " + melt_str)
            reader_times.write(str(int(system._time)) + " KINKED || " + melt_str + "\n")
        elif(nick_str == "NICK"): #if no kinks were reported at the nick and not catastrophic case (see previous check)
            kink_lig_count += 1
            print(str(int(system._time)) + " KINKED || " + ligated_str)
            reader_times.write(str(int(system._time)) + " KINKED || " + ligated_str + "\n")
        elif(ligated_str == "LIG"): #if no kinks were reported anywhere else except the nick and not catastrophic case (see first check)
            kink_nick_count += 1
            print(str(int(system._time)) + " KINKED || " + nick_str)
            reader_times.write(str(int(system._time)) + " KINKED || " + nick_str + "\n")
        else: #if kinks were reported both at the nick and not
            kink_both_count += 1
            print(str(int(system._time)) + " KINKED || " + nick_str + " " + ligated_str)
            reader_times.write(str(int(system._time)) + " KINKED || " + nick_str + " " + ligated_str + "\n")
        
        kink_total_count += 1
    else:
        #Not kinked, but now we will check for any broken stacking in the ring.
        #This includes normal and coaxial stacking, but not cross stacking, anywhere on the ring (both strands)
        
        flag, broken_pairs = is_all_intact(system, ref_ring, is_circular)
        if not flag: #this is counter-intuitive, but if not flag (empty string) means fully intact -> STACKED
            reader_times.write(str(int(system._time)) + " STACKED\n")
            print(str(int(system._time)) + " STACKED")
            stack_count +=1
            
        else:
            unstacked_str = ""
            for i in range(len(broken_pairs)):
                if(i != 0):
                    unstacked_str += " "
                unstacked_str += str(broken_pairs[i])
            
            unstack_count += 1
            if nicked:
                if([nick1, nick2] in broken_pairs or [nick2, nick1] in broken_pairs):
                    if(len(broken_pairs) == 1):
                        unstack_nick_count += 1
                        unstacked_str += " ( NICK )" # tells you the broken stacking interaction is at the nick
                    else:
                        unstack_both_count += 1
                        # this means that the nick-adjacent nucleotides are unstacked AND some other adjacent nucleotides
                        # which are not both paired to the nick-adjacent nucleotides are unstacked
                        unstacked_str += " ( NICK AND AWAY )"
                else:
                    unstack_away_count += 1
                    unstacked_str += " ( AWAY )" # tells you some other adjacent nucleotides which are not both paired to the nick-adjacent nucleotides are unstacked
            reader_times.write(str(int(system._time)) + " UNSTACKED || " + str(flag) + " " + str(unstacked_str) + "\n")
            print(str(int(system._time)) + " UNSTACKED ||", flag, unstacked_str)
    total_count += 1
    system = reader._get_system()
    # below code can be used to truncate analysis for debugging
    # if(system._time > 10000000):
    #     system = False

with open(stat_file, 'w') as reader_out:
    reader_out.write("Total: " + str(total_count) + "\n")
    reader_out.write("\n")
    reader_out.write("Kinked (total): " + str(kink_total_count) + "\n")
    reader_out.write("Percent: " + str((kink_total_count/total_count)*100) + "\n")
    if nicked:
        reader_out.write("Kinked (nick only): " + str(kink_nick_count) + "\n")
        reader_out.write("Percent: " + str((kink_nick_count/total_count)*100) + "\n")
        reader_out.write("Kinked (away from nick only): " + str(kink_lig_count) + "\n")
        reader_out.write("Percent: " + str((kink_lig_count/total_count)*100) + "\n")
        reader_out.write("Kinked (nick and away from nick): " + str(kink_both_count) + "\n")
        reader_out.write("Percent: " + str((kink_both_count/total_count)*100) + "\n")
    reader_out.write("\n")
    reader_out.write("Stacked: " + str(stack_count) + "\n")
    reader_out.write("Percent: " + str((stack_count/total_count)*100) + "\n")
    reader_out.write("\n")
    reader_out.write("Unstacked (total): " + str(unstack_count) + "\n")
    reader_out.write("Percent: " + str((unstack_count/total_count)*100) + "\n")
    if melt_count != 0: #extremely rare, so will only print this if actually occurs
        reader_out.write("\n")
        reader_out.write("Near-melted (total): " + str(melt_count) + "\n")
        reader_out.write("Percent: " + str((melt_count/total_count)*100) + "\n")
    if nicked:
        reader_out.write("Unstacked (across nick only): " + str(unstack_nick_count) + "\n")
        reader_out.write("Percent: " + str((unstack_nick_count/total_count)*100) + "\n")
        reader_out.write("Unstacked (away only): " + str(unstack_away_count) + "\n")
        reader_out.write("Percent: " + str((unstack_away_count/total_count)*100) + "\n")
        reader_out.write("Unstacked (across nick and away): " + str(unstack_both_count) + "\n")
        reader_out.write("Percent: " + str((unstack_both_count/total_count)*100) + "\n")
