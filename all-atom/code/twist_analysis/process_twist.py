import numpy as np

# Define input and output file names
input_file = "Twist__twist.xvg" 
output_file = "twist_internal_avg.xvg"


print(f"Processing {input_file}...")

# Loading data
try:
    data = np.loadtxt(input_file, comments=["@", "#"]) # skip comment / header lines
    time_column = data[:, 0] # define the first column as time
    
    # extract columns 1-42 (the 42 valid base pair steps), this is hard coded. change it per system. (in python indexing, if u want to go until column 42, you must add one, so 43)
    all_twist_steps = data[:, 1:43]
    
    # exclude 2 terminal steps from each end (keeps steps 3-40, indices 2-39) (we only want internal base pair steps, as the terminals flip out. this will change for nonlinear DNA)
    internal_twist_steps = all_twist_steps[:, 2:-2]
    
    # mask invalid values (this defines values that are invalid, like those above 100 degrees, equal to 999.0 (my data had a bunch of these at the end for some reason), or negative values)
    masked_data = np.ma.masked_where(
        (internal_twist_steps > 100) | (internal_twist_steps == 999.0) | (internal_twist_steps < 0), # so basically, these values are bad, "mask them" inside this defined data set internal_twist_steps
        internal_twist_steps
    )
    
    # calculate mean per frame, now our data is named masked_data, since we are using masked arrays to hide bad data
    avg_internal_twist_per_frame = np.ma.mean(masked_data, axis=1) # axis=1 means we want to average across columns (i.e., across the base pair steps) for each row (time frame), if its set to 0, it would average across rows (i.e., across time frames) for each column (base pair step)
    
    # save output
    output_data = np.column_stack((time_column, avg_internal_twist_per_frame)) # create 2D array with time and average twist
    header_text = (f"Time(ps)   Avg_Internal_Twist(deg)\n"
                   f"Excluding 2 terminal steps from each end (steps 1-2 and 41-42)") # define header
    np.savetxt(output_file, output_data, header=header_text, fmt="%10.3f  %10.4f") # save to file with formatting
    

    # this is just printed output in the terminal to help you know what happened
    print(f"✓ Saved to {output_file}")
    print(f"  Internal steps: {internal_twist_steps.shape[1]}")
    print(f"  Invalid values filtered: {np.sum(masked_data.mask)} ({100*np.sum(masked_data.mask)/masked_data.size:.2f}%)")
    print(f"  Overall average: {np.mean(avg_internal_twist_per_frame):.3f}°")
    print(f"  Helical repeat: {360/np.mean(avg_internal_twist_per_frame):.2f} bp/turn")

# this will catch errors like file not found, or data loading issues and it will print in the terminal    
except Exception as e:
    print(f"Error: {e}")
