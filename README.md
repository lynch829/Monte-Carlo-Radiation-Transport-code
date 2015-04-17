# Monte-Carlo-Radiation-Transport-code
            --------------------------------------
            3D SCATTERED LIGHT CODE
            --------------------------------------

This Git contains the source codes for the 3D Cartesian grid. The code can be set up to simulate emission from an isotropic point source, gaussian beam, pencil beam, and uniform irradiation. The density structure is set to be a uniform density discretised into the 3D grid. The code outputs the 3D density grid as an unformatted file. Various features are present in the code, such as: Frsenel reflections, photon weight and russian roulette system[L.Wang et al. 1995], support for semi infinte media and rough surfaces(currently in development).

At the end of the simulation the code outputs to the screen the average number of scatterings per Monte Carlo photon packet,total diffuse reflectance and total transmission. Other ouputs are the depth resolved fluence, Diffuse reflectance as afunction of radius, total trasmission as function of radius and various slices through the media of fluence.

The FORTRAN files are:

        density.f
        gridset.f
        iarray.f
        mcpolar.f
        ran2.f
        stokes.f
        search.f
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


Wang, Lihong, Steven L. Jacques, and Liqiong Zheng. "MCML—Monte Carlo modeling of light transport in multi-layered tissues." Computer methods and programs in biomedicine 47.2 (1995): 131-146.
