Input files used for the OL21 forcefield with CUFIX.
1. [ions.mdp](ions.mdp): simple ions mdp file
2. [em.mdp](em.mdp): energy minimization, simple, has xyc pbc and h-bond constraints
3. [nvt.mdp](nvt.mdp): 100 ps npt, v-rescale is used
4. [npt.mdp](npt.mdp): 100 ps npt, c-rescale and v-rescale are used
5. [md.mdp](md.mdp): 50 ns production run, 0.002 dt
6. [pulling_rotation_example.mdp](pulling_rotation_example.mdp): example production run mdp file used for DNA rotation and enforced h-bond bonds.
