#*************** Makefile created by Mike Wolff ****************************

#******************************** G77/Linux Fortran ************************
FC     =       gfortran-4.8  -fbounds-check
#EXTRA_OPT =     -mpentium -malign-double -fforce-mem -fforce-addr \
#                -ffast-math -funroll-all-loops
# May want to experiment by adding the extra optimization flags to get
## better runtime. But then again, maybe not.
#FFLAGS  =       -Ofast $(EXTRA_OPT)
FFLAGS  =       -Ofast
LDFLAGS = 
time_it         = get_cpu_sun

#****************************************************************************


OBJSB =     density.o \
            gridset.o \
            iarray.o \
            mcpolar.o \
            ran2.o \
            stokes.o \
            sourceph.o \
            tauint2.o \
            search.o \
            fresnel.o

mcgrid:	$(OBJSB)
		$(FC) $(OBJSB) $(LDFLAGS) -o mcgrid

clean:;		/bin/rm -f *.o

