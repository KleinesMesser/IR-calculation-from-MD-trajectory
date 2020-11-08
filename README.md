# IR-calculation-from-MD-trajectory
This project contains codes I currently use to calculate IR spectrum from MD simulations.
I and my coworkers developed the first frequency maps and coupling models within the nucleobase carbonyl stretch region. These models form a complete framework enabling one to predict the vibrational Hamiltonian from coordinate files generated from MD simulations.

To understand the theoretical background, please refer to the following papers:
1. Frequency map development and IR calculation: Y. Jiang and L. Wang, "Development of vibrational frequency maps for nucleobases," J. Phys. Chem. B, vol. 123, no. 27, pp. 5791-5804, 2019.
2. Coupling model development: Y. Jiang and L. Wang, "Modeling the vibrational couplings of nucleobases," J. Chem. Phys., vol. 152, no. 8, p. 084114, 2020.

In this project, I include all the codes needed to calculate IR spectrum from MD simulations. This code can be applied to systems containing multiple identical oligonucleotide strands with primary bases (AGCTU). Here I use an RNA oligomer as an example. Changes are needed if one wants to consider DNA oligomers because of the atom number difference between ribose and deoxyribose. The codes were written in fortran. The source codes (.f90) as well as the compiled files (.out) are both included. Besides, some necessary parameter files are included too.

Example input and output files (.dat) of this RNA oligomer are provided. The direct input coordinate file should be an xyz file. I also provide the corresponding pdb file so that one can know more detailed information concerning each atom. The charge file (1_Charge.dat) contains partial charges in AMBER internal unit and they are converted to electron charges accordingly in the electric field calculation (1_EF_traj.f90). The box parmeters (1_box.dat) are in Angstrom. These example files only contain 10 snaptshots. To acquire a reasonable final IR spectrum, much more snapshots are needed. For example, as indicate in these codes, I need 200,000 snapshots with a time step of 10 fs to get an IR spectrum. Finally, I provide an example script file to execute those calculation steps (3_script_IR).
