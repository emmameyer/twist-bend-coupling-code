#Compute the bend and twist angles between consecutive nucleotides following the arithmetic in supplement to https://pubs.aip.org/aip/jcp/article-abstract/146/21/214902/974005/DNA-elasticity-from-coarse-grained-simulations-The?redirectedFrom=fulltext (supplement to DNA elasticity from coarse-grained simulations: the effect of groove asymmetry).
#The bend and twist angles are contained in the vector in the paper defined as captial `Theta`. `Theta[1]` and `Theta[2]` representing bending degrees of freedom, and `Theta[3]` represent the twist degree of freedom.
#Using input parameter `num_prune`, you can opt to exclude a certain number of bases on each end of the duplex from the analysis and statistics.
#If choose to prune bases, regardless if the strand is determined to be linear or circular, the angle between the first index and last index will be excluded (i.e., won't wrap the indices). This is surely what you want.
#If don't prune bases, will still exclude the angle between the first index and last index for completely **linear** DNA. This is also likely what you want.
#Finally, if don't prune bases, and at least one strand is **circular**, the angle between the first and last index *will* be included in the analysis.
#This should be currently written so that all pruning is done before analyzing (i.e., will not burn time computing the angles at the pruned bases first before discarding).

#Credit for much of this goes to Daniel Hollis.

#For usage, run "python3 get_bend_and_twist_angles.py --help".

#import sys
#sys.path.insert(0, "<your path to the latest version of oxDNA analysis tools, if necessary>")
from oxDNA_analysis_tools.UTILS.data_structures import *
from oxDNA_analysis_tools.UTILS.RyeReader import *
import numpy as np
from typing import List, Tuple, TextIO
import numpy.typing as npt
import argparse

#Parse CLI input parameters.
prog = "get_bend_and_twist_angles.py"
parser = argparse.ArgumentParser(prog = prog, description='''Calculates the vector Theta as described in the supplement to article here: https://pubs.aip.org/aip/jcp/article-abstract/146/21/214902/974005/DNA-elasticity-from-coarse-grained-simulations-The?redirectedFrom=fulltext. 
                                 The first and second component of Theta represent bending angles, while the third component represents the twist angle (in the supplement, represented in units of radians). 
                                 Returned are text files with the averages of these angles across all dinucleotide pairs in the simulation, along with the raw angles for each nucleotide (all written in degrees).''')
parser.add_argument("input_file", type = str, nargs = 1, help = "The input file of the oxDNA simulation.")
parser.add_argument("comp_list_dir", type = str, nargs = 1, help = '''The indices file from get_comp_list.py (contains a space-separated list of complementary indices). 
                    Necessary for defining triad vectors in the theory.''')
parser.add_argument("avg_angles_dir", type = str, nargs = 1, help = "Output text file giving the average bend and twist angles across all dinucleotide steps of a given configuration line-by-line.")
parser.add_argument("all_angles_dir", type = str, nargs = 1, help = '''Output text file giving all of the raw measured angles across each dinucleotide step of every configuration. In principle can be used to reconstruct the average bend and twist angles. 
                    The first entry on any given line is the first model index of the dinucleotide step, and every set of three values thereafter are the values of Theta_1, Theta_2, and Theta_3 repeating for each configuration.''')
parser.add_argument("-n", "--num_prune", type = int, nargs = 1, dest = "num_prune", help = '''optional: Number of pruned bases to exclude on either end of the duplex for measuring the bend and twist angles. 
                    For example, if `num_prune == 3`, will exclude computation of the twist angles involving the three base pairs on the first end and the three base pairs on the second end of the duplex (for nicked circles, will attempt to exclude six base pairs centered by the nick from the analysis). 
                    Primarily this is to exclude base pairs that are prone to fraying (ends of linear duplex), or for excluding the twist angle across a kink (possibly ill-defined). 
                    Defaults to 0.''')
parser.add_argument("-c", "--check_kinked", dest = "check_kinked", action='store_true', default = False, help = "optional: If `times_dirs` is provided, flag this to check kinked configurations only. If unflagged, the program will search for unkinked configurations only.")
parser.add_argument("-t", "--times_dirs", type = str, nargs = 1, dest = "times_dirs", help = '''optional: Comma-separated list of files from kink detectors or fray finder for filtering configurations out on the basis of being kinked or unkinked. If not included, all configurations in the trajectory will be checked. 
                    If multiple are included, an intersection of the unkinked/kinked lines included in all files will be performed to determine the final list (i.e., if provided a fray finder and kink detector file, will exclude all frayed AND kinked confs from analysis).''')

args = parser.parse_args()
input_file: str = args.input_file[0]
comp_list_dir: str = args.comp_list_dir[0]
all_angles_txt: TextIO = open(args.all_angles_dir[0], "w")
avg_angles_txt: TextIO = open(args.avg_angles_dir[0], "w")
num_prune: int; times_dirs: List[str]; check_unkinked: bool
if args.times_dirs:
    times_dirs = args.times_dirs[0].split(",")
else:
    times_dirs = [] #counts as falsy for future checks
if args.check_kinked:
    check_unkinked = not args.check_kinked
else:
    check_unkinked = True
if args.num_prune:
    num_prune = args.num_prune[0]
else:
    num_prune = 0

#Parse `comp_list_txt` to find all complementary pairs on the ring.
def parse_comp_list_txt(comp_list_dir: str) -> List[Tuple[int,int]]:
    comp_list_txt: TextIO = open(comp_list_dir, "r")
    comp_sites_str: List[str] = comp_list_txt.readline().split()
    comp_pairs_list: List[Tuple[int,int]] = []
    for i in range(int(len(comp_sites_str))//2):
        comp_pairs_list.append((int(comp_sites_str[2*i]),int(comp_sites_str[2*i+1])))

    return comp_pairs_list

#Grabs the total number of steps of simulation from input file.
#This is a silly function due to an oversight (afaik, anyway) in how oxDNA's get_input_parameter works.
#When call get_input_parameter(input_file,"steps"), it grabs newtonian_steps value at the very least if
#steps key precedes newtonian_steps key in input file. Kludged by modifying get_input_parameter from base.py.
def get_num_steps(input_file: str) -> str:
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

#Reads in times files, and either builds a list of relevant ("rel") unkinked or kinked time steps depending on provided boolean flag `check_unkinked`.
#The times files currently are intended to come from programs `new_kink_finder_Harrison_method.py` or `kink_finder_precise.py`, whose structure is
#space-separated with the first column being the time step and following that keywords ("STACKED", "KINKED", "FRAYED", etc.) to identify the type of deformation.
def get_rel_times(file_dir: str, check_unkinked: bool = True) -> List[int]:
    times_list: List[int] = []

    file: TextIO = open(file_dir)
    while True:
        curLine: str = file.readline()
        if not curLine:
            break
        curLineSplit: List[str] = curLine.split()
        if check_unkinked:
            if "STACKED" in curLineSplit or "UNSTACKED" in curLineSplit:
                if "FRAYED" not in curLineSplit:
                    times_list.append(int(float(curLineSplit[0])))
        else: #check for kinked configurations (but exclude fraying as twist angle becomes more ill-defined)
            if "STACKED" not in curLineSplit and "UNSTACKED" not in curLineSplit:
                if "FRAYED" not in curLineSplit:
                    times_list.append(int(float(curLineSplit[0])))
    file.close()

    return times_list

#Returns a list of nucleotides on one strand to be analyzed for the bend and twist angles, attempting to strip 'num_prune' IDs off each end of the strand from the full `ID_list`.
#Essentially removes num_prune indices on both ends of the sorted `ID_list` containing `nuc`
def get_pruned_ID_list(system, top_file: str, num_prune: int) -> List[int]:
    #Determine nuc in the special case of nicked minicircle, otherwise just set 0 as irrelevant
    nuc: int = determine_nuc(top_file)
    print(f"Determined \'nuc\': {nuc}")
    
    #Determine internal strand ID for nuc.
    for strand in system.strands:
        for monomer in strand.monomers:
            if nuc == monomer.id:
                strand_index: int = strand.id - 1 #-1 as `strand.id` apparently isn't 0-indexed

    #Construct a list of all indices on the strand containing `nuc`.
    ID_list: List[int] = []
    for monomer in system.strands[strand_index].monomers:
        ID_list.append(monomer.id)

    #Attempt to prune the ends. Literally removes all indices from min_nuc_index (inclusive) to min_nuc_index + num_prune (exclusive), and max_nuc_index (inclusive) to max_nuc_index - num_prune (exclusive).
    ID_list = sorted(ID_list) #sort the list in order in case for whatever reason oxDNA would read them in out-of-order
    pruned_ID_list: List[int] = ID_list[num_prune:(len(ID_list) - num_prune)] #slice the list to exclude num_prune indices on either end of the list

    return pruned_ID_list

#Get the list of all complementary nucleotide IDs to those in provided `ID_list` using a full list of complementary pairs `comp_pairs_list`. This returns a parallel list to `ID_list`.
def get_comp_ID_list(ID_list: List[int], comp_pairs_list: List[List[int]]) -> List[int]:
    comp_ID_list: List[int] = []
    for ID in ID_list:
        for pair in comp_pairs_list:
            if ID == pair[0]:
                comp_ID_list.append(pair[1])
                break
            elif ID == pair[1]:
                comp_ID_list.append(pair[0])
                break
    return comp_ID_list
    
#Naive way of determining if duplex contains at least one circular strand by reading in the topology file (I couldn't find a method for this in new oxDNA libraries, but surely missed it).
#If there is an entry of '-1' for 5' or 3' neighbors, implies topologically not connected to another nucleotide in that direction.
#So if there are two strands (which is what this program is intended for), it is only a fully linear duplex if there are four '-1's in the last two columns, else at least one strand is circular.
def is_circular(top_file: str) -> bool:
    file: TextIO = open(top_file)
    counts: int = 0 #increment each time find a '-1' in split list
    file.readline() #ignore header line
    while True:
        curLine: str = file.readline()
        if not curLine:
            break
        curLineSplit: List[str] = curLine.split()
        if curLineSplit[2] == "-1":
            counts += 1
        if curLineSplit[3] == "-1":
            counts += 1
    file.close()
    if counts == 4:
        return False
    else:
        return True

#Determines 'nuc', which is only relevant when pruning bases. This is necessary for one particular case.
#For linear DNA and fully ligated circular DNA, it does not matter which 'nuc' you choose, though pruning on a fully ligated ring makes little sense for our current purposes.
#For nicked minicircles in particular, you must choose a 'nuc' to be on the nicked strand if you want it to prune the ends closest to the nick site. This is what we want. 
def determine_nuc(top_file: str) -> int:
    file: TextIO = open(top_file)
    nuc: int = 0
    file.readline() #ignore header line
    while True:
        curLine: str = file.readline()
        if not curLine:
            break
        curLineSplit: List[str] = curLine.split()
        if curLineSplit[2] == "-1" or curLineSplit[3] == "-1":
            return nuc #if it is a linear duplex, this will return nuc = 0, else will return nuc = 0 or nuc = N_bp/2 if nicked minicircle I believe
        nuc += 1
    file.close()
    return 0 #if it is a fully ligated minicircle

#Determine the model base-normal vectors of the provided `pair` of indices of the given `conf`, and calculate the normalized triad vector e_3 described in supplement to paper in header (essentially an average of the two base-normal vectors).
def get_e_3(conf, pair: List[int]) -> npt.NDArray:
    rel_a3s: List[npt.NDArray] = [conf.a3s[pair[0]], -1* conf.a3s[pair[1]]] #as these are on two different strands, the base-normal vectors face opposite directions, so must take the negative of one of the two vectors before averaging
    avg_a3: npt.NDArray = np.average(rel_a3s, axis = 0)
    e_3: npt.NDArray = avg_a3 / np.linalg.norm(avg_a3)
    
    return e_3

#Determine the model COM locations of the provided `pair` of indices of the given `conf`, and calculate the normalized e_2 vector in the Skoruppa supplement.
#y_hat is the vector stemming from the COM of `pair[0]` to COM of `pair[1]`, which is then slightly modulated using Gram-Schmidt to get orthogonal vector e_2 to e_3
#This is defined according to Triad II definition (used in Skoruppa main paper; claim is most robust)
def get_e_2(conf, pair: List[int], e_3: npt.NDArray) -> npt.NDArray:
    rel_COMs: List[npt.NDArray] = [conf.positions[pair[1]],conf.positions[pair[0]]]
    dis_COMs: npt.NDArray = rel_COMs[0]-rel_COMs[1]
    y_hat: npt.NDArray = dis_COMs/np.linalg.norm(dis_COMs) #this is y-hat as defined in paper
    
    e_2: npt.NDArray = (y_hat - np.dot(y_hat,e_3)*e_3)/np.linalg.norm(y_hat - np.dot(y_hat,e_3)*e_3)
    
    return e_2
    
#Determine Theta as defined in the Skoruppa supplement, which essentially involves linear algebra of T-matrices with form [e_1, e_2, e_3] (triad vectors form the columns of these matrices) for each base pair.
#All formulas are given in that supplement and reproduced here.
def get_Theta_vec(T_i: npt.NDArray, T_i_1: npt.NDArray) -> List[float]:
    R: npt.NDArray = np.matmul(T_i.transpose(), T_i_1) #R is the rotation matrix that transforms i -> i + 1 space iirc
    Theta_mag: float = np.arccos(np.trace(R)/2. - 1./2.) #this gives scalar Theta as described in supplement
    coeff: float = Theta_mag/(2.*np.sin(Theta_mag)) #leading scalar in formula for vector Theta
    vec: List[float] = [R[2][1] - R[1][2], R[0][2] - R[2][0], R[1][0] - R[0][1]] #vector in formula for Theta 
    Theta: List[float] = [coeff*entry for entry in vec] #vector Theta, which represent the bend (1st and 2nd components) and twist (third component) angles
    
    return Theta

#Driving analytic function that computes the bend and twist angles for all unpruned dinucleotide steps in the duplex for a provided `conf`. 
#Returns parallel lists `total_angles` and `avg_angles`, which each contain three entries: one for the total/avg Omega_1, one for the total/avg Omega_2, and one for the total/avg Omega_3 across all IDs in `pruned_ID_list` for the given configuration.
#Also returns `ID_Thetas`, which contains `len(pruned_ID_list)` lists containing the Theta vector for each ID. Essentially, this is a list of all measured angles for this configuration structured as a list of `Theta`s.
def get_bend_and_twist_angles(conf, pruned_ID_list: List[int], pruned_comp_ID_list: List[int], circular: bool) -> Tuple[npt.NDArray, List[float], List[List[float]]]:
    #Initialize main construct for holding all measured angles in the analysis.
    #It will form a parallel array to `pruned_ID_list` (i.e., `all_Theta_degs[3]` holds the Theta vector for ID `pruned_ID_list[3]`). 
    all_Theta_degs: List[List[float]] = []

    #Loop over every index on in-tact strand (every index in `pruned_ID_list`)
    for i, (ID, comp_ID) in enumerate(zip(pruned_ID_list, pruned_comp_ID_list)):
        #Determine ID_2, the neighboring base to ID (generally `ID_2 == ID + 1`, except if circular may wrap around to beginning of strand)
        ID_2: int; comp_ID_2: int
        if i < len(pruned_ID_list) - 1:
            ID_2 = pruned_ID_list[i+1]
            comp_ID_2 = pruned_comp_ID_list[i+1]
        elif circular: #wrap indices
            ID_2 = pruned_ID_list[0]
            comp_ID_2 = pruned_comp_ID_list[0]
        else: #noncircular and done
            break
        
        #Compute triad vectors e_i for both `ID` and `ID_2`. These form the bases represented as matrices.
        e_3_ID: npt.NDArray = get_e_3(conf, [ID, comp_ID])
        e_2_ID: npt.NDArray = get_e_2(conf, [ID, comp_ID], e_3_ID)
        e_1_ID: npt.NDArray = np.cross(e_2_ID, e_3_ID)
        basis_ID: List[npt.NDArray] = [e_1_ID, e_2_ID, e_3_ID]
        
        e_3_ID_2: npt.NDArray = get_e_3(conf, [ID_2, comp_ID_2])
        e_2_ID_2: npt.NDArray = get_e_2(conf, [ID_2, comp_ID_2], e_3_ID_2)
        e_1_ID_2: npt.NDArray = np.cross(e_2_ID_2, e_3_ID_2)
        basis_ID_2: List[npt.NDArray] = [e_1_ID_2, e_2_ID_2, e_3_ID_2]
        
        #Compute bend and twist angles in radians, and convert to degrees.
        T_i: npt.NDArray = np.array(basis_ID).transpose() #this creates a matrix of form [e_1, e_2, e_3]; transpose is needed to form column vectors out of the row-like list
        T_i_1: npt.NDArray = np.array(basis_ID_2).transpose() #this also creates a matrix of same form for the neighboring nucleotide
        Theta_vec: List[float] = get_Theta_vec(T_i, T_i_1) #compute Theta for given dinucleotide step
        rad2deg: float = 180 / np.pi
        Theta_vec_deg: List[float] = [Theta_comp*rad2deg for Theta_comp in Theta_vec] #convert Theta entries to degrees
        
        #Store measured bend and twist angle for current ID in this configuration.
        all_Theta_degs.append(Theta_vec_deg)
    
    #Sum all of the like angles across all IDs (i.e., contains the total Theta_1, total Theta_2, and total Theta_3 across all confs).
    total_angles: npt.NDArray = np.sum(all_Theta_degs, axis = 0)
    
    #Divide each entry of `total_angles` by number of measurements to get the `average_angles` (parallel to `total_angles`)
    avg_angles: List[float]
    if circular:
        avg_angles = [total_angle/len(pruned_ID_list) for total_angle in total_angles]
    else:
        avg_angles = [total_angle/(len(pruned_ID_list) - 1) for total_angle in total_angles]

    return total_angles, avg_angles, all_Theta_degs

#MAIN PROGRAM.
#Determine the values of Theta_1, Theta_2, and Theta_3 between consecutive nts for every configuration (or only unkinked/kinked configs if times files are provided and according to `check_unkinked` flag) excluding dinucleotide steps from the analysis as indicated and write it out.
#Particularly, all measured angles will be written out to `all_angles_txt`, while averages and totals will be written to `avg_angles_txt`. 
def main():
    #Create TopInfo and TrajInfo.
    traj_file: str = get_input_parameter(input_file, "trajectory_file")
    top_file: str = get_input_parameter(input_file, "topology")
    top_info, traj_info = describe(top_file, traj_file)

    #Read in the topology (creates a system with same hierarchy of System -> Strands -> Monomer in same way as old analysis code)
    system, monomers = strand_describe(top_file) #Creates System and list of Monomers; see data_structures.py for instance variables and believe monomers list is indexed exactly the same as oxview indexes nucleotides

    #Fatal error check.
    if len(system.strands) != 2:
        raise Exception("System does not contain exactly two strands. Program currently requires two strands as it is written. Aborting...")

    #Read topology file to see if either of the strands is circular. If it is, and num_prune == 0, will consider the angle between the wrapped indices.
    circular: bool = is_circular(top_file)

    #Finally, read in complementary pairs from file.
    comp_pairs_list: List[Tuple[int,int]] = parse_comp_list_txt(comp_list_dir)

    #Print some diagnostics and perform some checks. If `num_prune == 0` and there is at least one circular strand, consider the angles between the wrapped indices; else, override and discard the ends.
    if num_prune != 0:
        print(f"Pruning {num_prune} bases on each end of the duplex.")
    else: 
        print("Not pruning any bases on each end of duplex.")

    if circular:
        if num_prune == 0:
            print("Flagged duplex as circular, and chose not to prune ends, so program will include the angles between the wrapped indices on each end of the duplex.")
        else:
            print("Duplex auto-detected as having a circular strand, yet chose to prune ends, so will force analysis to ignore the angle between the wrapped indices of the pruned lists.")
            circular = False
    else:
        print("Flagged duplex as linear. This means the program will *not* check the angle between the ending indices of the duplex.")

    if times_dirs:
        type_readin_str: str = "unkinked" if check_unkinked else "kinked"
        print(f"Including only {type_readin_str} configurations in the analysis, with the exception that they will still be excluded if they are also indicated as frayed{' (in this case, this tends to mean the strong majority of `KINKED_AWAY` configurations)' if not check_unkinked else ''}.")

    #Build list of list of IDs to be considered after pruning. These will be the indices whose angles will be analyzed.
    pruned_ID_list: List[int] = get_pruned_ID_list(system, top_file, num_prune)
    pruned_comp_ID_list: List[int] = get_comp_ID_list(pruned_ID_list, comp_pairs_list)

    #Find number of configurations in current simulation (assumes a linear time scale).
    #'steps' holds the total number of time steps of the simulation (the 'steps' key in input file), and the print_conf_interval is directly grabbed from input file as well.
    steps: int = int(float(get_num_steps(input_file)))
    print_conf_interval: int = int(float(get_input_parameter(input_file,"print_conf_interval")))
    num_confs: int = int(steps/print_conf_interval)

    #Build list of the relevant times to analyze. If times files are provided, either check for unkinked confs or kinked confs depending on set `check_unkinked` flag, and perform an intersection across all files if multiple are provided. 
    #If no times file provided, read in all configurations.
    rel_times: List[int]
    if times_dirs: #read files if provided
        cur_times: List[int]
        for i, times_dir in enumerate(times_dirs):
            cur_times: List[int] = get_rel_times(times_dir, check_unkinked)
            print(f"Found {len(cur_times)} {type_readin_str} configurations in file: {times_dir}.")
            if i == 0:
                rel_times = cur_times
            else:
                rel_times = list(set(rel_times) & set(cur_times)) #intersection of current list with all unkinked/unfrayed confs already found in other files
        if not rel_times: #intersection is empty -> no times found meeting `type` criterion
            raise Exception(f"Times files were provided but no {type_readin_str} configurations found in the intersection of all provided files. Check times files.")
        else: #found unkinked/kinked configurations in file
            print(f"Analyzing {len(rel_times)} {type_readin_str} configurations.")
            
    else: #no times file provided
        rel_times = [] #this will count as falsy
        print(f"No times file provided. Reading in all {num_confs} configurations.")

    #Minor diagnostics and setup for main loop. Prepare all angles file for writing.
    print(f"Final list of nucleotide indices considered (post-pruning if applicable):\n{pruned_ID_list}")
    print(f"Determined complementary list of nucleotide indices to above list:\n{pruned_comp_ID_list}")
    increment: int = 10000
    print_ID_list: List[int] = []
    if circular:
        print_ID_list = pruned_ID_list
    else:
        print_ID_list = pruned_ID_list[0:len(pruned_ID_list) - 1] #slice the last entry, as do not wrap when num_prune != 0 (and otherwise prints an empty line)
    all_angles_lines: List[str] = [] #will hold all lines to be placed in all_angles_txt
    for ID in print_ID_list:
        all_angles_lines.append(f"{ID}\n")

    #Read in an `increment` number of configurations at a time, analyze them, and move onto next batch.
    for i in range(int(num_confs//increment)):
        #Read in configurations and, if `times_dirs` provided, only analyze the corresponding configurations.
        cur_confs: List = get_confs(top_info, traj_info, i*increment, increment)
        rel_confs: List = []
        if rel_times:
            for conf in cur_confs:
                if conf.time in rel_times:
                    rel_confs.append(conf)
        else:
            rel_confs = cur_confs
    
        #Analyze each "relevant" conf and write out/prepare to write out to files:
        for conf in rel_confs:
            conf = inbox(conf) #inbox the configuration so that PBCs do not affect vector analysis
            total_angles: List[float]; avg_angles: List[float]; all_Theta_degs: List[List[float]] 
            total_angles, avg_angles, all_Theta_degs = get_bend_and_twist_angles(conf, pruned_ID_list, pruned_comp_ID_list, circular) #main fn: compute all bend/twist angles and store in `all_Theta_degs`, also compute totals/averages of each component of Theta across all dinuc steps
            avg_angles_txt.write(f"{conf.time} {avg_angles[0]} {avg_angles[1]} {avg_angles[2]} {total_angles[0]} {total_angles[1]} {total_angles[2]}\n")

            #Prepare lines for all angles file output.
            for i, Theta_vec_deg in enumerate(all_Theta_degs):
                all_angles_lines[i] = "".join([all_angles_lines[i].strip(), f" {Theta_vec_deg[0]} {Theta_vec_deg[1]} {Theta_vec_deg[2]}\n"])
    
    #Write out all prepared lines for all angles file.
    all_angles_txt.writelines(all_angles_lines)

    print(f"Wrote all measured angles to {args.all_angles_dir[0]} and averages/totals for each configuration to {args.avg_angles_dir[0]}.")

if __name__ == "__main__":
    main()
