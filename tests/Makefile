# ~~~~~~~~~~~~~~~~~
# PFLARE - Steven Dargaville
# Makefile for tests
# ~~~~~~~~~~~~~~~~~

# Have to include the modules in ../include
INCLUDEDIR  := ../include
INCLUDE		:= -I$(INCLUDEDIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include

# These are the same in the top level makefile, but we redefine them here rather than use an export
# in the top level makefile
# That way we can compile PFLARE with openmp
# But to build the tests we then reset the CFLAGS and FFLAGS to empty and include 
# include -lgomp in the LDFLAGS
# That way we don't link to threaded libraries on cray machines (as that happens automatically
# if the flags include -fopenmp)
CFLAGS := ${CFLAGS} -O3 -fPIC
FFLAGS := ${FFLAGS} -O3 -fPIC

# ~~~~~~~~~~~~~~~~~~~~~~~~
# Compiler specific flags
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Intel
ifneq ($(filter ifx mpiifx,$(FC)),)
FORTMOD     := -module
FFLAGS      := ${FFLAGS} -real-size 64 -fpscomp logicals
endif
# GNU
ifneq ($(filter gfortran mpif90,$(FC)),)
FFLAGS      := ${FFLAGS} -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-none -ffree-line-length-none
endif
# LLVM 
ifneq ($(filter flang amdflang,$(FC)),)
endif	
# ~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~

# Get any library flags we have on input
LDFLAGS_INPUT := $(LDFLAGS)
# Include pflare in the linking
LDFLAGS 		:= -L$(LIBDIR) -lpflare -Wl,-rpath,$(LIBDIR)
# Include petsc
LDFLAGS     := $(LDFLAGS) -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib
# Have to link to blas, lapack
# Don't have to do this on cray machines as the compiler wrappers include links
ifeq ($(filter ftn,$(FC)),)
LDFLAGS     := $(LDFLAGS) -lblas -llapack -lm
endif
# Given we have C files in our library need to add gfortran if GNU compiler
ifneq ($(filter gfortran mpif90,$(FC)),)
LDFLAGS     := $(LDFLAGS) -lgfortran
endif
# Have to link to mpi
LDFLAGS     := $(LDFLAGS) -lmpi
# Include any input flags we had - these go last as order is important
LDFLAGS     := $(LDFLAGS) $(LDFLAGS_INPUT)


# Define file extensions
.SUFFIXES: .F90 .o .mod .c

# Test examples
OBJS := ex12f.o ex6f.o ex6f_getcoeffs.o ex6.o ex6_cf_splitting.o ex86.o adv_diff_2d.o

# Build rules
all: $(OBJS) 

# Fortran
%.o: %.F90
	$(FC) $(FFLAGS) $(INCLUDE) $^ $(LDFLAGS) $(LDLIBS) -o $@

# C
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(LDFLAGS) $(LDLIBS) -o $@

run_tests_load:
#
	@echo ""
	@echo "Test AIRG with GMRES polynomials for hyperbolic streaming problem"
	./ex12f.o -f data/mat_stream_2364
	@echo "Test AIRG with GMRES polynomials for hyperbolic streaming problem, matrix-free smoothing"
	./ex12f.o -f data/mat_stream_2364 -pc_air_matrix_free_polys
#
	@echo ""
	@echo "Test AIRG with GMRES polynomials for hyperbolic streaming problem in C"
	./ex6.o -f data/mat_stream_2364
	@echo "Test AIRG with GMRES polynomials for hyperbolic streaming problem, matrix-free smoothing in C in parallel"
	$(MPIEXEC) -n 2 ./ex6.o -f data/mat_stream_2364 -pc_air_matrix_free_polys	
#
	@echo ""
	@echo "Test lAIR with GMRES polynomial smoothing for hyperbolic streaming problem"
	./ex12f.o -f data/mat_stream_2364 -pc_air_z_type lair
	@echo "Test lAIR with strong R tolerance with GMRES polynomial smoothing for hyperbolic streaming problem"
	./ex12f.o -f data/mat_stream_2364 -pc_air_z_type lair	-pc_air_strong_r_threshold 0.01
# 
	@echo ""
	@echo "Test single level GMRES polynomial preconditioning for hyperbolic streaming problem in C"
	./ex6.o -f data/mat_stream_2364 -pc_type pflareinv -pc_pflareinv_type power	
	@echo "Test single level GMRES polynomial preconditioning for hyperbolic streaming problem in C in parallel"
	$(MPIEXEC) -n 2 ./ex6.o -f data/mat_stream_2364	-pc_type pflareinv -pc_pflareinv_type power
# 
	@echo ""
	@echo "Test single level GMRES polynomial preconditioning with the Newton basis matrix-free for hyperbolic streaming problem in C"
	./ex6.o -f data/mat_stream_2364 -pc_type pflareinv -pc_pflareinv_type newton -pc_pflareinv_matrix_free	
	@echo "Test single level GMRES polynomial preconditioning with the Newton basis matrix-free for hyperbolic streaming problem in C in parallel"
	$(MPIEXEC) -n 2 ./ex6.o -f data/mat_stream_2364	-pc_type pflareinv -pc_pflareinv_type newton -pc_pflareinv_matrix_free
#
	@echo ""	 
	@echo "Test AIRG as an exact solver, truncating hierarchy and using high order mf GMRES poly in the\
	 Arnoldi basis as a coarse solver for hyperbolic streaming"
	./ex12f.o -f data/mat_stream_2364 -pc_air_strong_threshold 0.0 -pc_air_a_drop 0.0 -pc_air_r_drop 0.0 \
	 -pc_air_inverse_type jacobi -mg_coarse_ksp_type richardson -mg_coarse_ksp_max_it 5 -ksp_type richardson -ksp_norm_type unpreconditioned \
	 -pc_air_max_levels 30 -pc_air_coarsest_poly_order 18 \
	 -pc_air_coarsest_matrix_free_polys -pc_air_coarsest_inverse_type arnoldi
# 
	@echo ""		 
	@echo "Test AIRG as an exact solver, heavily truncating hierarchy and using high order mf GMRES poly in the\
	 Newton basis as a coarse solver for hyperbolic streaming"
	./ex12f.o -f data/mat_stream_2364 -pc_air_strong_threshold 0.0 -pc_air_a_drop 0.0 -pc_air_r_drop 0.0 \
	 -pc_air_inverse_type jacobi -ksp_type richardson -ksp_norm_type unpreconditioned \
	 -pc_air_max_levels 10 -pc_air_coarsest_poly_order 60 \
	 -pc_air_coarsest_matrix_free_polys -pc_air_coarsest_inverse_type newton -pc_air_max_luby_steps 3
	@echo "Test AIRG as an exact solver, heavily truncating hierarchy and using high order mf GMRES poly in the\
	 Newton basis as a coarse solver for hyperbolic streaming in parallel"
	$(MPIEXEC) -n 2 ./ex6.o -f data/mat_stream_2364 -pc_air_strong_threshold 0.0 -pc_air_a_drop 0.0 -pc_air_r_drop 0.0 \
	 -pc_air_inverse_type jacobi -ksp_type richardson -ksp_norm_type unpreconditioned \
	 -pc_air_max_levels 10 -pc_air_coarsest_poly_order 60 \
	 -pc_air_coarsest_matrix_free_polys -pc_air_coarsest_inverse_type newton -pc_air_max_luby_steps 3	
#
	@echo ""
	@echo "Test AIRG with GMRES polynomials for hyperbolic streaming problem with coefficients calculated on subcomms"
	$(MPIEXEC) -n 2 ./ex12f.o -f data/mat_stream_2364 -pc_air_subcomm -pc_air_inverse_type arnoldi -pc_air_coarsest_subcomm \
	 -pc_air_coarsest_inverse_type arnoldi
# 
	@echo ""
	@echo "Test PMISR DDC CF splitting in C"
	./ex6_cf_splitting.o  -f data/mat_stream_2364
	@echo "Test PMISR DDC CF splitting in C in parallel"
	$(MPIEXEC) -n 2 ./ex6_cf_splitting.o  -f data/mat_stream_2364   

run_tests_no_load:
#
	@echo ""
	@echo "Test AIRG with GMRES polynomials for 2D finite difference stencil"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air
	@echo "Test AIRG with GMRES polynomials Arnoldi basis for 2D finite difference stencil"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air -pc_air_inverse_type arnoldi -pc_air_coarsest_inverse_type arnoldi
	@echo "Test AIRG with GMRES polynomials for 2D finite difference stencil in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air
	@echo "Test AIRG with GMRES polynomials Arnoldi basis for 2D finite difference stencil in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air -pc_air_inverse_type arnoldi -pc_air_coarsest_inverse_type arnoldi
	@echo ""
	@echo "Test lAIR with GMRES polynomial smoothing for 2D finite difference stencil in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air -pc_air_z_type lair
# 
	@echo ""
	@echo "Test single level GMRES polynomial preconditioning for 2D finite difference stencil"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type pflareinv -pc_pflareinv_type power
	@echo "Test single level GMRES polynomial preconditioning for 2D finite difference stencil in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type pflareinv -pc_pflareinv_type power
#
	@echo ""
	@echo "Test AIRG with Neumann polynomials for 2D finite difference stencil"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air -pc_air_inverse_type neumann
	@echo "Test AIRG with Neumann polynomials for 2D finite difference stencil in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air  -pc_air_inverse_type neumann
	@echo "Test AIRG with Neumann polynomials for 2D finite difference stencil, matrix-free smoothing"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air -pc_air_inverse_type neumann -pc_air_matrix_free_polys
	@echo "Test AIRG with Neumann polynomials for 2D finite difference stencil, matrix-free smoothing in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air  -pc_air_inverse_type neumann -pc_air_matrix_free_polys
#
	@echo ""
	@echo "Test AIRG with SAIs for 2D finite difference stencil"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air -pc_air_inverse_type sai
	@echo "Test AIRG with SAI for 2D finite difference stencil in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air  -pc_air_inverse_type sai
# 
	@echo ""
	@echo "Test AIRG with ISAIs for 2D finite difference stencil"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air -pc_air_inverse_type isai
	@echo "Test AIRG with ISAIs for 2D finite difference stencil in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air  -pc_air_inverse_type isai
# 
	@echo ""
	@echo "Test AIRG with Weighted Jacobi for 2D finite difference stencil"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air -pc_air_inverse_type wjacobi
	@echo "Test AIRG with Weighted Jacobi for 2D finite difference stencil in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air  -pc_air_inverse_type wjacobi
# 
	@echo ""
	@echo "Test AIRG with Unweighted Jacobi for 2D finite difference stencil"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air -pc_air_inverse_type jacobi
	@echo "Test AIRG with Unweighted Jacobi for 2D finite difference stencil in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air  -pc_air_inverse_type jacobi
# 
	@echo ""
	@echo "Test AIRG as an exact solver for 2D finite difference stencil"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air \
	 -pc_air_strong_threshold 0.0 -pc_air_a_drop 0.0 -pc_air_r_drop 0.0 -pc_air_inverse_type jacobi \
	 -mg_coarse_ksp_type richardson -mg_coarse_ksp_max_it 10 -ksp_type richardson -ksp_norm_type unpreconditioned
	@echo "Test AIRG as an exact solver for 2D finite difference stencil in parallel"
	$(MPIEXEC) -n 2 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 8 -da_grid_y 8 -pc_type air \
	 -pc_air_strong_threshold 0.0 -pc_air_a_drop 0.0 -pc_air_r_drop 0.0 -pc_air_inverse_type jacobi \
	 -mg_coarse_ksp_type richardson -mg_coarse_ksp_max_it 10 -ksp_type richardson -ksp_norm_type unpreconditioned
# 
	@echo ""
	@echo "Test AIRG with GMRES polynomials with PC reused with no sparsity change"
	./ex6f.o -m 10 -n 10
	@echo "Test AIRG with GMRES polynomials with PC regenerated with no sparsity change"
	./ex6f.o -m 10 -n 10 -regen -pc_air_reuse_sparsity
	@echo "Test AIRG with GMRES polynomials with PC regenerated with no sparsity change and strong R threshold"
	./ex6f.o -m 10 -n 10 -regen -pc_air_reuse_sparsity -pc_air_strong_r_threshold 0.01	
	@echo "Test lAIR with GMRES polynomials with PC regenerated with no sparsity change"
	./ex6f.o -m 10 -n 10 -regen -pc_air_z_type lair -pc_air_reuse_sparsity
	@echo "Test lAIR SAI with GMRES polynomials with PC regenerated with no sparsity change"
	./ex6f.o -m 10 -n 10 -regen -pc_air_z_type lair_sai -pc_air_reuse_sparsity	
# 
	@echo ""
	@echo "Test AIRG with GMRES polynomials with PC regenerated with no sparsity change and polynomial coeffs stored"
	./ex6f_getcoeffs.o -m 10 -n 10
	@echo "Test AIRG with GMRES polynomials with PC regenerated with no sparsity change and polynomial coeffs stored and lumping"
	./ex6f_getcoeffs.o -m 10 -n 10 -pc_air_a_lump
# 
	@echo ""
	@echo "Test single level GMRES polynomials with PC reused with no sparsity change"
	./ex6f.o -m 10 -n 10 -pc_type pflareinv
	@echo "Test single level GMRES polynomials with PC regenerated with no sparsity change"
	./ex6f.o -m 10 -n 10 -pc_type pflareinv -regen	
# 
	@echo ""
	@echo "Test solving isotropic diffusion with fast coarsening and near-nullspace"
	./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 50 -da_grid_y 50 -pc_type air -pc_air_z_type lair \
	 -pc_air_one_c_smooth -pc_air_symmetric -pc_air_constrain_z \
	 -pc_air_cf_splitting_type agg -pc_air_a_drop 1e-5 -pc_air_a_lump\
	 -pc_air_r_drop 0 -ksp_rtol 1e-10 -ksp_pc_side right	
	@echo "Test solving isotropic diffusion with fast coarsening and near-nullspace in parallel"
	$(MPIEXEC) -n 4 ./adv_diff_2d.o -u 0 -v 0 -alpha 1.0 -da_grid_x 50 -da_grid_y 50 -pc_type air -pc_air_z_type lair \
	 -pc_air_one_c_smooth -pc_air_symmetric -pc_air_constrain_z \
	 -pc_air_cf_splitting_type agg -pc_air_a_drop 1e-5 -pc_air_a_lump\
	 -pc_air_r_drop 0 -ksp_rtol 1e-10 -ksp_pc_side right		
# 
	@echo ""
	@echo "Test AIRG on steady 1D advection"
	./ex86.o -n 1000 -ksp_rtol 1e-10 -ksp_atol 1e-50 -ksp_pc_side right \
	 -pc_air_coarsest_inverse_type newton -pc_air_coarsest_poly_order 10 -pc_air_coarsest_matrix_free_polys
	@echo "Test AIRG on steady 1D advection in parallel"
	$(MPIEXEC) -n 2 ./ex86.o -n 1000 -ksp_rtol 1e-10 -ksp_atol 1e-50 -ksp_pc_side right \
	 -pc_air_coarsest_inverse_type newton -pc_air_coarsest_poly_order 10 -pc_air_coarsest_matrix_free_polys
# 
	@echo ""
	@echo "Test lAIR on steady 2D structured advection"
	./adv_diff_2d.o -da_grid_x 50 -da_grid_y 50 -pc_type air -ksp_pc_side right \
	 -ksp_rtol 1e-10 -ksp_atol 1e-50 -pc_air_a_lump -pc_air_a_drop 1e-4 -pc_air_one_c_smooth \
	 -pc_air_z_type lair -pc_air_inverse_type wjacobi
# 
	@echo ""
	@echo "Test AIRG on steady 2D structured advection with 0th order sparsity GMRES poly C smooth and faster coarsening"
	./adv_diff_2d.o -da_grid_x 50 -da_grid_y 50 -pc_type air -ksp_pc_side right \
	 -ksp_rtol 1e-10 -ksp_atol 1e-50 -pc_air_a_lump -pc_air_a_drop 1e-4 -pc_air_one_c_smooth \
	 -pc_air_c_inverse_sparsity_order 0 -pc_air_strong_threshold 0.99 -pc_air_ddc_fraction 0.0
# 
	@echo ""
	@echo "Test high order GMRES polynomials in small linear system"
	./adv_diff_2d.o -da_grid_x 5 -da_grid_y 5 -pc_type pflareinv -pc_pflareinv_type newton \
	 -pc_pflareinv_matrix_free -pc_pflareinv_order 16
	@echo "Test high order GMRES polynomials in small linear system - slightly bigger"
	./adv_diff_2d.o -da_grid_x 10 -da_grid_y 10 -pc_type pflareinv -pc_pflareinv_type newton \
	 -pc_pflareinv_matrix_free -pc_pflareinv_order 50
# 
	@echo ""
	@echo "Test auto truncation with matrix-free Newton polynomials as coarse grid solver" 
	./adv_diff_2d.o -da_grid_x 10 -da_grid_y 10 -ksp_type richardson -pc_type air \
	 -pc_air_coarsest_inverse_type newton -pc_air_coarsest_poly_order 10 -pc_air_coarsest_matrix_free_polys \
	 -pc_air_auto_truncate_start_level 1 -pc_air_auto_truncate_tol 1e-2

clean:
	$(RM) *.o *.mod *.out