# 5th-proj

5th year project.
Monte Carlo code.
		--------------------------------------
		3D SCATTERED LIGHT CODE: POINT SOURCES
		--------------------------------------

This directory contains the source codes for the 3D Cartesian grid. The code can be set up to simulate emission from an isotropic point source or other sources such as a gaussian beam or pencil beam. The density structure is set to be a uniform density sphere discretised into the 3D grid. The code outputs the 3D density grid as an unformatted file. At the end of the simulation the code outputs to the screen the average number of scatterings per Monte Carlo photon packet and outputs diffuse reflectance, transmission and depth resolved fluence.

The FORTRAN files are:

	  density.f
	  gridset.f
        iarray.f
        mcpolar.f
        ran2.f
        stokes.f
        sourceph.f
        tauint2.f
        fresnel.f
Include files are:

	grid.txt  
	photon.txt  

Input parameters are in:

	input.params

The file that compiles the code and creates the executable file 'mcgrid' is:

	Makefile


