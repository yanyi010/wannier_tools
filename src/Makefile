OBJ =  module.o sparse.o wt_aux.o math_lib.o mgmres.o symmetry.o \
		  readHmnR.o inverse.o proteus.o \
       eigen.o ham_qlayer2qlayer.o psi.o unfolding.o rand.o \
  		 ham_slab.o ham_bulk.o ek_slab.o ek_bulk_polar.o ek_bulk.o \
       readinput.o fermisurface.o surfgreen.o surfstat.o \
  		 mat_mul.o ham_ribbon.o ek_ribbon.o \
       fermiarc.o berrycurvature.o \
  		 wanniercenter.o dos.o  orbital_momenta.o \
  		 landau_level_sparse.o landau_level.o lanczos_sparse.o \
		 berry.o wanniercenter_adaptive.o \
  		 effective_mass.o findnodes.o \
		 sigma_OHE.o sigma.o Boltz_transport_anomalous.o \
		 2D_TSC.o \
	main.o

# compiler
F90  = mpiifort -fpp -DMPI -DINTELMKL -fpe3

INCLUDE = -I${MKLROOT}/include
WFLAG = -nogen-interface
OFLAG = -O3 -static-intel
FFLAG = $(OFLAG) $(WFLAG)
LFLAG = $(OFLAG)

# ARPACK LIBRARY
ARPACK=/data/home/yanyi2024/0-tools/0-Basic/ARPACK/libarpack_MAC.a

# blas and lapack libraries
# static linking
#LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
        ${MKLROOT}/lib/intel64/libmkl_sequential.a \
        ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl \
		  ${ARPACK}

# dynamic linking
LIBS = ${ARPACK} -L/${MKLROOT}/lib/intel64 -lmkl_core -lmkl_sequential -lmkl_intel_lp64 -lpthread
 
main : $(OBJ)
	$(F90) $(LFLAG) $(OBJ) -o wt-dev.x $(LIBS)
	cp -f wt-dev.x ../bin

.SUFFIXES: .o .f90

.f90.o :
	$(F90) $(FFLAG) $(INCLUDE) -c $*.f90

clean :
	rm -f *.o *.mod *~ wt-dev.x

