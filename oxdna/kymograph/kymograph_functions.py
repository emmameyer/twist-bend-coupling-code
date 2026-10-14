# Function file for plotting kymograph
# usage: from kymograph_functions import *
# Author: Rachel Bricker
# Acknowledgments: Thank you Tommy Swope for catching/fixing mistakes and updating the code

from oxDNA_analysis_tools.UTILS.RyeReader import *
import matplotlib.pyplot as plt
from math import log10, floor

def get_colors(file) -> list[tuple[int]]:
    """
        Returns a list which contains the "color" of each configuration

        ---------------------------------
        | Color   | Meaning             |
        ---------------------------------
        | Red     | Stacked at nick     | 
        | Blue    | Kinked at nick      |
        ---------------------------------

        For a configuration to be considered "Kinked at nick", the duplex must be
        kinked at the nick.
 
        The configuration is considered "Stacked at nick" if the duplex is not kinked at
        the nick.

        Parameters:
            file       (text file)        : precise_kink_finder.py times output file

        Returns:
            color_list (list[tuple[int]]) : a list of RGB values corresponding to each
                                            configuration in ascending order
    """

    stacked_color   = (1, 0, 0)    # red for configurations stacked at nick
    unstacked_color = (0, 0, 1)    # dark blue for configurations which are not stacked at the nick

    color_list = []

    line = file.readline()

    while line:
        line_split = line.split()

        # On July 13, 2024, I removed a section of code that would never be executed
        # (which was brought to my attention by Tommy Swope... thank you); however,
        # I also accidentally removed necessary lines of code, too! D:
        # Tommy fixed this on January 24, 2025 by adding:
        # `or "UNSTACKED" in line_split`
        # to the following if statement... thank you again! :D

        if "STACKED" in line_split or "UNSTACKED" in line_split:    # ring is unkinked
            color_list.append(stacked_color)
        else:    # the ring is kinked somewhere
            if "NICK" not in line_split:    # check if the kink is at the nick site
                color_list.append(stacked_color)
            else:
                color_list.append(unstacked_color)

        line = file.readline()

    file.close()

    return color_list

def get_avg_angles(file) -> list[float]:
    """
        Returns a list containing the average twist angle of each configuration.
        
        Parameters:
            file       (text file)   : get_bend_and_twist_angles.py average bend and twist angle output file.
                                       Preferably, at least the two base pairs closest to the
                                       nick site should be pruned.

        Returns:
            avg_angles (list[float]) : a list of the average twist angle corresponding
                                       to each configuration in ascending order
    """

    avg_angles = []

    line = file.readline()
    while line:
        line_split = line.split()

        # UPDATE MADE BY TOMMY (THANK YOU):
        # Each line of the file follows the following format:
        # (configuration index) (average Theta_1) (average Theta_2) (average Theta_3) (total Theta_1) (total Theta_2) (total Theta_3)
        # Therefore, grab the fourth element (representing average twist angle).
        avg_angles.append(float(line_split[3]))

        line = file.readline()

    file.close()

    return avg_angles

def get_base_pairs_per_turn_list(avg_angles_list : list[float]) -> list[float]:
    """
        Returns a list containing the number of base pairs per turn of each configuration given
        a list containing the average twist angle of each cnfiguration: avg_angles_list

        Parameters:
            avg_angles_list (list[float]) : list obtained from get_avg_angles()
    """

    return [(360.0 / avg_angle) for avg_angle in avg_angles_list]

def get_n_confs(colors : list[str], avg_angles : list[float]) -> int:
    """
        Retrieves the number of configurations from the lists "colors" and "avg_angles". Both
        lists are used as a safety check to make sure they both contain the same number of
        elements.
    """

    # make sure the user inputted corresponding files
    if (len(avg_angles) == len(colors)):
        n_confs = len(colors)    # the number of configurations
    else:
        print("Error: The number of configurations in the average twist angles file (" + str(len(avg_angles)) +
              ") does not match the number of configurations recorded in the kink detectection times file ("
              + str(len(colors)) + ")")
        exit()

    return n_confs

def get_time_list(print_time_interval : float, n_confs : int) -> list[float]:
    """
        Creates a time list given the amount of time elapsed between snapshots (print_interval)
        and the number of configurations (n_confs)
    """

    time = []    # unit depends on print_interval
    for i in range(n_confs):
        time.append(i*print_time_interval)

    return time

def find_where_to_slice_list(time : list[float], lower_limit : float, upper_limit : float) -> tuple[int, int]:
    """
        Finds the indexes of the elements which are equal or else closest to lower_limit and upper_limit.

        Parameters:
            time        (list[float]) : time list
            lower_limit (float)       : the lower bound of the range (must be in the same units as the
                                        elements in time list)
            upper_limit (float)       : the upper bound of the range (must be in the same units as the
                                        elements in time list)

        Returns:
            index_1     (int)         : the index of the element which is either equal to lower_limit
                                        or else the element which is the closest to lower_limit and less
                                        than lower_limit
            index_2     (int)         : the index of the element which is either equal to upper_limit
                                        or else the element which is the closest to upper_limit and
                                        greater than upper_limit
    """

    index_1 = -1
    index_2 = -1
    for i in range(len(time)):
        if (time[i] <= lower_limit):
            index_1 = i     # this will update until time[i] is no longer less than or equal to lower_bound
        elif (time[i] >= upper_limit):
            index_2 = i
            break

    return index_1, index_2

def set_size(w, h, ax=None):
    """
        Sets explicitly the size of the axes. This function was obtained from
        https://stackoverflow.com/questions/44970010/axes-class-set-explicitly-size-width-height-of-axes-in-given-units

        Parameters: 
            w (float) : desired width of horizontal axis (in inches)
            h (float) : desired height of vertical axis (in inches)
    """

    # get axis
    if not ax: ax=plt.gca()

    l    = ax.figure.subplotpars.left
    r    = ax.figure.subplotpars.right
    t    = ax.figure.subplotpars.top
    b    = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)

    ax.figure.set_size_inches(figw, figh)

def round_sig(num, sig):
    """
        Rounds a number (num) to a specific number of significant digits (sig). Function obtained by indgar
        on Stack Overflow:
        https://stackoverflow.com/questions/3410976/how-to-round-a-number-to-significant-figures-in-python
    """

    return round(num, sig - int(floor(log10(abs(num)))) - 1)

def get_print_conf_interval_and_dt(input_file : str) -> tuple[int, float]:
    """
        Parameters:
            input_file          (str)   : path to the input file

        Returns:
            print_conf_interval (int)   : the number of configurations between snapshots
            dt                  (float) : time step of the integration in ps
    """

    # get print interval
    print_conf_interval = int(float(get_input_parameter(input_file, "print_conf_interval")))

    dt = float(get_input_parameter(input_file, "dt"))    # unit: time simulation units

    # 1 simulation time unit = 3.03 ps
    dt = round_sig(dt*3.03, 2)    # unit: picoseconds

    return print_conf_interval, dt

def compute_print_time_interval(print_conf_interval : int, dt : float) -> float:
    """
        Calculate the time elasped (in ns) between snapshots

        Parameters:
            print_conf_interval (int)   : the number of configurations between snapshots
                                          (obtained from get_print_conf_interval_and_dt())
            dt                  (float) : time step of the integration in ps (obtained from
                                          get_print_conf_interval_and_dt())

        Returns:
            print_time_interval (float) : the time elasped (in ns) between snapshots
    """

    print_time_interval = (dt*print_conf_interval) / 1000    # unit: nanoseconds

    return print_time_interval

def convert_conf_index_to_time(conf_index : int, print_time_interval : float, print_conf_interval : int) -> float:
    """
        Computes the amount of time elasped (measured in the same units as print_time_interval)
        up to the specified configuration (conf_index).

        Parameters:
            conf_index          (int)   : the index of the configuration
            print_time_interval (float) : the time elasped between snapshots (obtained from
                                          compute_print_time_interval())
            print_conf_interval (int)   : the number of configurations between snapshots
                                          (obtained from get_print_conf_interval_and_dt())
    """

    return (conf_index*print_time_interval) / (print_conf_interval)
