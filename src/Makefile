SHELL = /bin/bash
MAIN = FLEXDUST_Iceland

FC       = gfortran
INCPATH  = /xnilu_wrk/flex_wrk/bin64/grib_api/include
LIBPATH1 = /xnilu_wrk/flex_wrk/bin64/grib_api/lib
LIBPATH2 =   /usr/lib/x86_64-linux-gnu/
FFLAGS   =    -O3 -Wunused -fbounds-check -m64 -ffree-form -mcmodel=medium -fconvert=little-endian -frecord-marker=4 -I$(INCPATH)
LDFLAGS  = $(FFLAGS) -L$(LIBPATH2) -L$(LIBPATH1) -lgrib_api_f90 -lgrib_api -lm -ljasper


MODOBJS = \
par_mod.o \
com_mod.o \
cmapf_mod.o \
conv_mod.o \
dust_mod.o

OBJECTS = \
FLEXDUST.o \
subroutinesSoil.o \
readInput.o \
juldate.o \
caldate.o \
shift_field.o \
shift_field_0.o \
pbl_profile.o \
psim.o \
psih.o \
ew.o \
getfields.o \
readwind_nests.o \
calcpar.o \
calcpar_nests.o \
verttransform.o \
verttransform_nests.o \
richardson.o \
scalev.o \
getvdep.o \
getvdep_nests.o \
obukhov.o \
calcpv.o \
calcpv_nests.o \
getrb.o \
getrc.o \
raerod.o \
partdep.o \
qvsat.o \
readavailable.o \
readwind.o \
gridcheck.o \
dustEmission_HSO.o \
getSizeDistribution.o \
writeRELEASEfile.o \
subroutines.o \
dustEmission_MB95.o \
dustEmission_MB95_tau.o\
dustEmission_Kok14.o \
gridcheck_nests.o \
calcThreshold.o \
calcThreshold_tau.o



$(MAIN): $(MODOBJS) $(OBJECTS)
	$(FC) *.o -o $(MAIN) $(LDFLAGS)

$(OBJECTS): $(MODOBJS)

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm *.o *.mod

cleanall:
	rm *.o *.mod $(MAIN)
