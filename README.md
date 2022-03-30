# EFDA-of-Persistence-Landscapes

This repository is meant to accompany the paper 'Robust Persistence Homology Through Amplitude-Phase Separation in Persistence Landscapes'.

It contains code and data to reproduce the results for simulated data presented in the paper. 

There are 3 scripts and 2 directories, which are summarized as follows:

Scripts:
	main_replicatio1.m - Used to reproduce the simulations 1 and 2 in the paper, 
			     and replicate figures 3 and 5. 
	main_replication2.m - Used to reproduce the simulations 1 and 2 with noise in the paper, 
			     and replicate figures 4 and 6.
	main_replication3.m - Used to reproduce the circles, spirals, and torus examples in the paper,
			     and replicate figures 7, 8 and 9.

	All scripts call code MATLAB functions the 'code' directory and .mat files from the 'data' directory.
	Before running any examples, the section entitled 'Load code and mex c file for alignment' in each script must be ran. 

Directories:

	code: this directory contains MATLAB functions that perform alignment of persistence landscapes. There are 2 main 
	     functions ran in the scripts: 
		
		align_landscapes.m: input of this function are a time index set, t, and an array of landscapes, land. The
				   function outputs, a mean SRVF, reparameterizations, and aligned SRVFs. 
		plotLandscape.m: input of this function are a time index set, t, and a matrix corresponding to a discretized
			      landscape, landscape, and z axis limits, zl. 
			      Calling the function produces a plot of a landscape, e.g., Figure 5(a).
		The remaining functions are called within align_landscapes.m. 

	data: this directory contains .mat files that contain the data used in analyses. The correspondence between example
	     and .mat file is as follows.
		
		Figure 3: main_sim_1.mat
		Figure 5: main_sim_2.mat
		Figure 4: main_sim_1_noise.mat
		Figure 6: main_sim_2_noise.mat 
		Figure 7: main_circles.mat 
		Figure 8: main_spirals.mat 
		Figure 9: main_torus.mat 


