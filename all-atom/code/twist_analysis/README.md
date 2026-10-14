# README file for twist/helical repeat analysis
This file lays out the general plan for extracting twist angles from DNA md simulations. Most of these scripts are pretty hard coded so they must be edited for each scenario and are designed for linear DNA, however, it can be used for circular, just keep in mind the twist step between the ligated/nick site (like between 1 and 80 for an 80 bp circle).

1. Extract data from the simulation using do_x3dna: do_x3dna -f *.xtc -s *.gro -lbpsm -noavg -name _twist

2. run [process_twist.py](process_twist.py) to extract the twist angles from the outputted file (the file used here has twist values for each time step for each base pair step, this script averages the twist values so that there is just one averaged twist angle for each time step)
## Then you can go two different ways. I have block average -> plot_twist, and also straight from process_twist to plot_h.py
3. [plot helical repeat over the simulation](../plotters/plot_h.py) using the outputted file from process_twist.py. this is the preferred route.
## Or go the block average route.
3. run [block_average.py](block_average.py): this block averages the twists over the simulation and outputs averages for four segments of the simulation, then finds the average of that.

4. run [plot_twist.py](plot_twist.py) (very much so hardcoded, probably not ideal for analysis of one scenario, it is a comparison plot, block averages)
