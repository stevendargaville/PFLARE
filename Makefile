# ~~~~~~~~~~~~~~~~~
# PFLARE - Steven Dargaville
# Makefile for PFLARE
#
# Must have defined PETSC_DIR and PETSC_ARCH before calling
# PETSc must be at least version 3.15
# Copied from $PETSC_DIR/share/petsc/Makefile.basic.user
# This uses the compilers and flags defined in the PETSc configuration
# ~~~~~~~~~~~~~~~~~

# Get the flags we have on input
CFLAGS_INPUT := $(CFLAGS)
FFLAGS_INPUT := $(FFLAGS)
CPPFLAGS_INPUT := $(CPPFLAGS)
FPPFLAGS_INPUT := $(FPPFLAGS)

# Read in the petsc compile/linking variables and makefile rules
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Directories we want
INCLUDEDIR  := include
SRCDIR      := src
# This needs to be exported into the sub-makefile for Cython
export LIBDIR := $(CURDIR)/lib

# Include directories - include top level directory in case compilers output modules there
INCLUDE := -I$(CURDIR) -I$(INCLUDEDIR)

# Output the library
OUT = $(LIBDIR)/libpflare.so

# All the files required by PFLARE
OBJS := $(SRCDIR)/NonBusyWait.o \
		  $(SRCDIR)/Binary_Tree.o \
		  $(SRCDIR)/TSQR.o \
		  $(SRCDIR)/Gmres_Poly_Data_Type.o \
		  $(SRCDIR)/AIR_Data_Type.o \
		  $(SRCDIR)/Matshell_Data_Type.o \
		  $(SRCDIR)/Matshell.o \
		  $(SRCDIR)/Sorting.o \
		  $(SRCDIR)/C_PETSc_Interfaces.o \
		  $(SRCDIR)/PCPFLAREINV_Interfaces.o \
		  $(SRCDIR)/PCAIR_Data_Type.o \
		  $(SRCDIR)/PETSc_Helper.o \
		  $(SRCDIR)/Gmres_Poly.o \
		  $(SRCDIR)/Gmres_Poly_Newton.o \
		  $(SRCDIR)/AIR_MG_Stats.o \
		  $(SRCDIR)/SAI_Z.o \
		  $(SRCDIR)/Constrain_Z_or_W.o \
		  $(SRCDIR)/PMISR_DDC.o \
		  $(SRCDIR)/Aggregation.o \
		  $(SRCDIR)/CF_Splitting.o \
		  $(SRCDIR)/Repartition.o \
		  $(SRCDIR)/Timers.o \
		  $(SRCDIR)/Weighted_Jacobi.o \
		  $(SRCDIR)/Neumann_Poly.o \
		  $(SRCDIR)/Approx_Inverse_Setup.o \
		  $(SRCDIR)/AIR_MG_Setup.o \
		  $(SRCDIR)/PCAIR_Shell.o \
		  $(SRCDIR)/PCAIR_Interfaces.o \
		  $(SRCDIR)/PFLARE.o \
		  $(SRCDIR)/C_PETSc_Routines.o \
		  $(SRCDIR)/C_Fortran_Bindings.o \
		  $(SRCDIR)/PCAIR_C_Fortran_Bindings.o \
		  $(SRCDIR)/PCAIR.o \
		  $(SRCDIR)/PCPFLAREINV.o	

# Define a variable containing all the tests
export TEST_TARGETS = ex12f ex6f ex6f_getcoeffs ex6 adv_1d adv_diff_2d ex6_cf_splitting

# Add the pflare include files
PETSC_FC_INCLUDES += $(INCLUDE)
PETSC_CC_INCLUDES += $(INCLUDE)

# Include any additional flags we input
CFLAGS += $(CFLAGS_INPUT)
FFLAGS += $(FFLAGS_INPUT)
CPPFLAGS += $(CPPFLAGS_INPUT)
FPPFLAGS += $(FPPFLAGS_INPUT)

# ~~~~~~~~~~~~~~~~~~~~~~~~
# Has petsc has been configured with 64 bit integers
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the petscconf.h
PETSC_HEADER_FILE := $(PETSC_DIR)/$(PETSC_ARCH)/include/petscconf.h
CONTENTS := $(file < $(PETSC_HEADER_FILE))
PETSC_USE_64BIT_INDICES := 0
ifneq (,$(findstring PETSC_USE_64BIT_INDICES 1,$(CONTENTS)))
PETSC_USE_64BIT_INDICES := 1
endif  		  
		  	
all: $(OUT)

# Create our directory structure and build the shared library
$(OUT): $(OBJS)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(INCLUDEDIR)
	$(LINK.F) -shared -o $(OUT) $(OBJS) $(LDLIBS)

# Build the tests
build_tests: $(OUT)
	@for t in $(TEST_TARGETS); do \
		$(MAKE) -C tests $$t; \
	done

# Build and run all the tests
# Only run the tests that load the 32 bit test matrix in /tests/data
# if PETSC has been configured without 64 bit integers
.PHONY: tests
tests: $(OUT)
	$(MAKE) build_tests
ifeq ($(PETSC_USE_64BIT_INDICES),0)
	$(MAKE) -C tests run_tests_load
endif	
	$(MAKE) -C tests run_tests_no_load

# Build the Python module with Cython
.PHONY: python
python: $(OUT)
	$(MAKE) -C python python

# Run the python tests
tests_python: $(OUT)
	$(MAKE) -C python run_tests

# Cleanup
clean::
	$(RM) -r $(LIBDIR)
	$(RM) $(INCLUDEDIR)/*.mod
	$(RM) $(SRCDIR)/*.mod
	$(RM) $(SRCDIR)/*.o
	$(RM) $(CURDIR)/*.mod
	$(MAKE) -C tests clean
	$(MAKE) -C python clean
