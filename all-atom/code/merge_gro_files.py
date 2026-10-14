import MDAnalysis as mda
from MDAnalysis import Merge

# Load the two systems
ref = mda.Universe("dna_processed.gro")
mob = mda.Universe("histone_processed.gro")

# Merge the two universes
merged = Merge(ref.atoms, mob.atoms)

# Write to a new GRO file
merged.atoms.write("merged_system.gro")

print("Merged system written to merged_system.gro")


# i also have used the following code too:

#import MDAnalysis as mda
# Load the two gromacs files
#u1 = mda.Universe("/Users/emma/Documents/kent/minicircles/BSC1/82bp_Tw8_deltaTw0/strand1_processed.gro")
#u2 = mda.Universe("/Users/emma/Documents/kent/minicircles/BSC1/82bp_Tw8_deltaTw0/strand2_processed.gro")
# Combine atoms from both universes
#combined = mda.Merge(u1.atoms, u2.atoms)
# Save to the same directory
#output_path = "/Users/emma/Documents/kent/minicircles/BSC1/82bp_Tw8_deltaTw0/merged.gro"
#combined.atoms.write(output_path)
