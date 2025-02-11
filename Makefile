# ~~~~~~~~~~~~~~~~~
# PFLARE - Steven Dargaville
# Makefile for PFLARE
#
# Must have defined PETSC_DIR and PETSC_ARCH before calling
# PETSc must be at least version 3.15
# ~~~~~~~~~~~~~~~~~
 
# Compiler options - need both a fortran and c compiler
# with appropriate mpi wrappings

# If you want to override these just pass in the values when you call this make
# e.g., make CC=cc, FC=ftn
export FC := mpif90
# This needs to be exported into the sub-makefile for Cython to see the MPI headers
export CC := mpicc
# The routine to launch an executable with mpi (sometimes mpiexec, sometimes mpirun, etc)
# exported into the sub-makefiles for tests
export MPIEXEC := mpiexec

# Compile flags - we deliberately don't export these here
# If we want to compile PFLARE with omp, we can set the CFLAGS and FFLAGS variables
# externally to have -fopenmp
# When we compile the tests, we need to add -gomp to the LDFLAGS and not have
# -fopenmp in the compile flags, as:
# 1) They don't explicitly have any omp in them
# 2) On cray machines, linking happens when we compile/link the tests
#    and we want don't want PFLARE to be linked to the threaded libraries
#    This is because we have to set the OMP_NUM_THREADS/=1 (and sort out the pinning)
#    to have the cray runtime correctly run in oversubscribed mode 
#    but that would then trigger the threaded blas/lapack
#    The only omp we want is internal to PFLARE
CFLAGS := ${CFLAGS} -O3 -fPIC
FFLAGS := ${FFLAGS} -O3 -fPIC

SHARED_FLAG := -shared
FORTMOD     := -J

# ~~~~~~~~~~~~~~~~~~~~~~~~
# Compiler specific flags
# If on a Cray machine, just add whichever flags are necessary to FFLAGS before compiling
# There is no easy way to work out what compiler is used given the ftn wrapper
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Intel
ifneq ($(filter ifx mpiifx,$(FC)),)
FORTMOD     := -module
FFLAGS      := ${FFLAGS} -fpscomp logicals
endif
# GNU
ifneq ($(filter gfortran mpif90,$(FC)),)
FFLAGS      := ${FFLAGS} -ffixed-line-length-none -ffree-line-length-none
endif
# LLVM 
ifneq ($(filter flang amdflang,$(FC)),)
endif	
# ~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~

INCLUDEDIR  := include
SRCDIR      := src
OBJDIR      := obj
# This needs to be exported into the sub-makefile for Cython
export LIBDIR	:= $(CURDIR)/lib

# Include directories
INCLUDE		:= -I$(INCLUDEDIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include

# Output the library
OUT = $(LIBDIR)/libpflare.so

# Define file extensions
.SUFFIXES: .F90 .o .mod .c

# All the files required by PFLARE
OBJS := $(OBJDIR)/NonBusyWait.o \
		  $(OBJDIR)/Binary_Tree.o \
		  $(OBJDIR)/TSQR.o \
		  $(OBJDIR)/Gmres_Poly_Data_Type.o \
		  $(OBJDIR)/AIR_Data_Type.o \
		  $(OBJDIR)/Matshell_Data_Type.o \
		  $(OBJDIR)/Matshell.o \
		  $(OBJDIR)/Sorting.o \
		  $(OBJDIR)/C_PETSc_Interfaces.o \
		  $(OBJDIR)/PCPFLAREINV_Interfaces.o \
		  $(OBJDIR)/PCAIR_Data_Type.o \
		  $(OBJDIR)/PETSc_Helper.o \
		  $(OBJDIR)/Gmres_Poly.o \
		  $(OBJDIR)/Gmres_Poly_Newton.o \
		  $(OBJDIR)/AIR_MG_Stats.o \
		  $(OBJDIR)/SAI_Z.o \
		  $(OBJDIR)/Constrain_Z_or_W.o \
		  $(OBJDIR)/PMISR_DDC.o \
		  $(OBJDIR)/Aggregation.o \
		  $(OBJDIR)/CF_Splitting.o \
		  $(OBJDIR)/Repartition.o \
		  $(OBJDIR)/Timers.o \
		  $(OBJDIR)/Weighted_Jacobi.o \
		  $(OBJDIR)/Neumann_Poly.o \
		  $(OBJDIR)/Approx_Inverse_Setup.o \
		  $(OBJDIR)/AIR_MG_Setup.o \
		  $(OBJDIR)/PCAIR_Shell.o \
		  $(OBJDIR)/PCAIR_Interfaces.o \
		  $(OBJDIR)/PFLARE.o \
		  $(OBJDIR)/C_PETSc_Routines.o \
		  $(OBJDIR)/C_Fortran_Bindings.o \
		  $(OBJDIR)/PCAIR_C_Fortran_Bindings.o \
		  $(OBJDIR)/PCAIR.o \
		  $(OBJDIR)/PCPFLAREINV.o		  

# Build rules
all: $(OBJDIR) $(OUT)

# Create the build directories
$(OBJDIR):
	mkdir -p $(OBJDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(INCLUDEDIR)

# Fortran
# Place the .o and .mod files in the $(OBJDIR) directory
$(OBJDIR)/%.o: $(SRCDIR)/%.F90
	$(FC) $(FFLAGS) -c $(INCLUDE) $^ -o $@ $(FORTMOD) $(INCLUDEDIR)

# C files
# Place the .o files in the $(OBJDIR) directory
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $(INCLUDE) $^ -o $@

# Make our shared library
$(OUT): $(OBJS) 
	$(FC) $(SHARED_FLAG) -fPIC $^ -o $(OUT)

# Build the tests
build_tests: $(OUT)
	$(MAKE) -C tests

# Build and run all the tests
tests: $(OUT)
	$(MAKE) -C tests
	$(MAKE) -C tests run_tests_load
	$(MAKE) -C tests run_tests_no_load

# Build and run only the tests that don't load 
# the example matrix (which uses 32 bit ints)
tests_no_load: $(OUT)
	$(MAKE) -C tests
	$(MAKE) -C tests run_tests_no_load	

# Build the Python module with Cython
.PHONY: python
python: $(OUT)
	$(MAKE) -C python python

# Run the python tests
tests_python: $(OUT)
	$(MAKE) -C python run_tests

# Cleanup
clean:
	$(RM) -r $(OBJDIR)
	$(RM) -r $(LIBDIR)
	$(RM) $(INCLUDEDIR)/*.mod
	$(RM) $(SRCDIR)/*.mod
	$(MAKE) -C tests clean
	$(MAKE) -C python clean

.PHONY: tests
