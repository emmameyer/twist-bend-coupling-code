import numpy as np

# input file outputted from process_twist.py
input_file = "twist_internal_avg.xvg"


print(f"Running block averaging on {input_file}...")


try:
    # load the data, skipping comment lines
    data = np.loadtxt(input_file, comments=["#"])
    
    # column 0 is time (ps), column 1 is twist (deg)
    time = data[:, 0]
    twist = data[:, 1]
    
    # filter out any NaN values (from masked array operations)
    valid_mask = ~np.isnan(twist) # this creates a boolean mask where twist is not NaN, so True where valid and False where NaN
    time = time[valid_mask] # this defines new time arrays that only have valid data
    twist = twist[valid_mask] # this defines new twist arrays that only have valid data
    
    print(f"Loaded {len(time)} valid timepoints from {time[0]:.0f} to {time[-1]:.0f} ps")
    
    # define the time boundaries for each 50 ns (50,000 ps) block: for block averaging
    boundaries = [0, 50000, 100000, 150000, 200000]
    block_averages = [] # empty list to store block averages
    block_std_devs = [] # empty list to store block standard deviations
    
    print("\n--- Block Analysis Results ---")
    print("Block | Time Range (ns) |  Average (deg)  | Std Dev (deg) | N points")
    print("-----------------------------------------------------------------------")
    
    # this is saying for each block average defined by the boundaries above, for the length of the boundaries -1 (so if there are 5 boundaries, there are 4 blocks)
    for i in range(len(boundaries) - 1):
        start_time = boundaries[i] #the start_time is that boundary time
        end_time = boundaries[i+1] # the end_time is the next boundary time
        
        # create a boolean mask to select data within the time range
        if i == 0:  # include start time for the first block
            mask = (time >= start_time) & (time <= end_time) # boolean array, True where time is within the block range and false where its not
        else:
            mask = (time > start_time) & (time <= end_time)
        
        # select the twist values for this block
        block_data = twist[mask] # apply the mask to get twist values (from the twist array) within this time block and store in block_data
        
        if len(block_data) > 0:  # only calculate if there is data in this block, amount of data points greater than 0
            avg = np.mean(block_data) # avg of the block is found from the mean of the data
            std = np.std(block_data, ddof=1)  # sample standard deviation with a degree of freedom of 1
            block_averages.append(avg) # append the average to the list of block averages
            block_std_devs.append(std) # append the std dev to the list of block std devs
            
            range_ns = f"{start_time/1000:.0f} - {end_time/1000:.0f}" # time range in ns for printing, this just converts ps to ns
            print(f"  {i+1}   | {range_ns:<15} | {avg:^15.3f} | {std:^15.3f} | {len(block_data):^8d}") # print block results in the terminal
        else:
            range_ns = f"{start_time/1000:.0f} - {end_time/1000:.0f}" # time range in ns for printing
            print(f"  {i+1}   | {range_ns:<15} | --- NO DATA --- | --- NO DATA --- | 0") # print if no data in block
    
    # final calculations
    # final average is the mean of the block averages
    final_avg = np.mean(block_averages) # simply averages the block averages
    
    # standard error of the mean from block averages
    final_std_err = np.std(block_averages, ddof=1)
    

    # prints stuff in the terminal
    print("-----------------------------------------------------------------------")
    print("\n--- Final Results ---")
    print(f"Overall Average (Mean of blocks): {final_avg:.3f} deg")
    print(f"Standard Error (Std Dev of blocks): {final_std_err:.3f} deg")
    print(f"Helical Repeat: {360/final_avg:.2f} bp/turn")
    print(f"\nReport as: {final_avg:.3f} ± {final_std_err:.3f} deg")
    
    # Save the block averages for plotting
    np.savetxt("block_averages_for_plotting.dat", block_averages, fmt="%.3f")
    
except Exception as e:
    print(f"An error occurred: {e}")
