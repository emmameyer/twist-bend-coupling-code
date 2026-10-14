# Function file for plotting graphs
# usage: from functions_for_plots import *
# Author: Rachel Bricker

import matplotlib.pyplot as plt
import statistics
import matplotlib.font_manager as font_manager

def set_rcParameters():
    """
        Sets rcParameters
        
        Returns:
            font_leg : for setting the font for legends
    """
    
    font_size    = 6
    font_leg     = font_manager.FontProperties(family="Arial", style='normal', size=font_size)
    widths       = 1
    line_width   = widths+(widths*0.5)
    rcParameters = {'font.size'        : font_size,                                  # font size
                    'font.family'      : "sans-serif",                               # specify the font
                    'font.sans-serif'  : "Arial",
                    'axes.linewidth'   : widths,                                     # set width of axes
                    'xtick.major.width': widths,                                     # set tick widths
                    'xtick.minor.width': widths,
                    'ytick.major.width': widths,
                    'ytick.minor.width': widths,
                    'grid.color'       : 'lightgray',                                # set grid color
                    'grid.linewidth'   : widths,                                     # set grid width
                    'grid.linestyle'   : '--',                                       # dashed grid
                    'lines.linewidth'  : line_width,                                 # line width
                    'patch.linewidth'  : widths,                                     # legend line width
                    'axes.axisbelow'   : True,                                       # plot grid behind data
                    'axes.titlesize'   : font_size,                                  # font size of plot title
                    'axes.prop_cycle'  : plt.cycler("color", plt.cm.Set1.colors)}    # color cycle
    plt.rcParams.update(rcParameters)
    
    return font_leg

def read_xvg_file(file):
    """
        Gets data from a .xvg file from GROMACS. Function is based off of code from:
        https://github.com/kulasinski/python-for-gromacs/blob/master/xvgplot.py

        Parameters:
            file (xvg file)          : MD simulation analysis script from GROMACS

        Returns:
            data (list[list[float]]) : list which stores each line of the input file as a list which
                                       stores each column of that line as an element
    """

    data = []
    for line in file:
        if line[0] in ['#','@']: # ignore title, labels, and comments in file
            continue
        s = line.split()
        data.append([float(x) for x in s])
    
    return(data)

def moving_average(data, max_window_size):
    """
        The moving average is calculated with a varying window size to:
            (i)  force the initial data point to be plotted (the initial configuration of all DNA
                 models are generated from the same structure, meaning that all DNAs initially
                 have a similar value, e.g. radius of gyration value, in the original series...
                 we want this to be shown in our plots); and
            (ii) keep the size of the output series the same as the original series.
        
        In the first iteration, the window size = 1.
        In each following iteration, the window size will increase by two (i.e. increase the number of
        elements we take from both the left and the right of the i-th element by one) until the window
        size = max_window_size.
        When we approach the rightmost boundary of the original data series and the window size can no
        longer equal max_window_size, we decrease the window size by one (i.e. decrease the number of
        elements we take from the right by one) each iteration.

        Parameters:
            data            (list[int or float]) : data points to compute moving average from
            max_window_size (int)                : the maximum size of the window, i.e., the 
                                                   maximum number of data points used to calculate
                                                   the moving average

        Returns:
            moving_averages (list[float])        : series of averages of data (same size as the input
                                                   series)
    """

    moving_averages = []
    
    # number of data points needed to the LEFT of the i-th element to have a window size of
    # max_window_size
    n_data_left = max_window_size//2
    
    # number of data points needed to the RIGHT of the i-th element to have a window size of
    # max_window_size
    n_data_right = max_window_size-1-max_window_size//2
    
    # list of the number of data points to the LEFT that you put into the window when there
    # aren't enough to make a window of size max_window_size
    varying_n_data_left = list(range(n_data_left+1))
    
    # list of the number of data points to the RIGHT that you put into the window when there
    # aren't enough to make a window of size max_window_size
    varying_n_data_right = list(reversed(range(n_data_right+1)))
    
    # allows us to iterate varying_n_data_right
    offset = None

    for i in range(len(data)):
        if i < n_data_left: 
            # if we are near the leftmost boundary and there aren't enough data points to the left and right
            # of the i-th element to make the window size = max_window_size, use a window size < max_window_size
            
            # get window: the number of data points we take from both the left and right of the i-th element vary
            window = data[i-varying_n_data_left[i] : i+varying_n_data_right[-1-i]+1]
            
            # compute window average
            window_avg = sum(window) / len(window)

            moving_averages.append(window_avg)
        elif i >= (len(data)-n_data_right):
            # if we are near the rightmost boundary and there aren't enough data points to the left and right
            # of the i-th element to make the window size = max_window_size, use a window size < max_window_size

            if offset == None:
                offset = i
            index = i-offset    # allows us to iterate varying_n_data_right

            # get window: only the number of data points we take from the right of the i-th element varies
            window = data[i-n_data_left : i+varying_n_data_right[index]+1]
            
            # compute window average
            window_avg = sum(window) / len(window)

            moving_averages.append(window_avg)
        else:
            # else we have enough data points to create a window of size max_window_size

            # get window: the window size = max_window_size and the number of data points we take from both
            # the left and right of the i-th element does NOT vary
            window = data[i-n_data_left : i+n_data_right+1]

            # compute window average
            window_avg = sum(window) / len(window)

            moving_averages.append(window_avg)

    return moving_averages
    
def plot_data(x, y, x_label, y_label, title, legend, file_name, fig_width, fig_height):
    # set rcParams
    font_leg = set_rcParameters()

    # set figure dimensions
    fig, ax = plt.subplots(1, figsize=(fig_width, fig_height))

    # plot data
    for i in range(len(y)):
        plt.plot(x, moving_average(moving_average(y[i], 500), 100), label=legend[i])

    # position legend to the left
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop=font_leg) 

    # set x-axis label
    plt.xlabel(x_label)
        
    # set y-axis label
    plt.ylabel(y_label)

    # set title
    if title != None:
        plt.title(title)

    # show grid
    plt.grid()

    plt.tight_layout()

    # save figure
    plt.savefig(file_name, bbox_inches="tight", dpi=600)
