<img align="right" img src="PFLARE_logo.png" width="300" height="300" />

# PFLARE library
### Author: Steven Dargaville

This library contains methods which can be used to solve linear systems in parallel with PETSc, with interfaces in C/Fortran/Python.
   
PFLARE can scalably solve:
1) Hyperbolic problems implicitly without Gauss-Seidel methods, such as advection equations, streaming operators from Boltzmann transport applications, multigrid in time discretisations, etc. This includes time dependent or independent equations, with structured or unstructured grids, with lower triangular structure or without.
3) Other asymmetric problems such as heavily anisotropic Poisson/diffusion problems.
4) Symmetric Poisson/diffusion problems, including those with varying coefficients.

## New methods in PFLARE

1) **PMISR DDC**: A new parallel CF splitting algorithm. When used on a matrix $\mathbf{A}$ this returns a set of "fine" and "coarse" points. PMISR DDC is similar to a PMIS CF splitting, but the resulting fine-fine submatrix $\mathbf{A}_ \textrm{ff}$ is more diagonally dominant than $\mathbf{A}$. $\mathbf{A}_ \textrm{ff}$ can also be made strongly diagonally dominant if desired.

2) **PCPFLAREINV**: A new PETSc PC type, containing methods for computing approximate inverses, most of which can be applied as assembled matrices or matrix-free. PCPFLAREINV can be used with the command line argument ``-pc_type pflareinv``, with several different PFLAREINV types available with ``-pc_pflareinv_type``:

   | Command line type  | Flag | Description |
   | ------------- | -- | ------------- |
   | power  |  PFLAREINV_POWER  | GMRES polynomial with the power basis  |
   | arnoldi  |  PFLAREINV_ARNOLDI  | GMRES polynomial with the Arnoldi basis  |
   | newton  |  PFLAREINV_NEWTON  | GMRES polynomial with the Newton basis  |
   | neumann  |  PFLAREINV_NEUMANN  | Neumann polynomial  |
   | sai  |  PFLAREINV_SAI  | Sparse approximate inverse  |
   | isai  |  PFLAREINV_ISAI  | Incomplete sparse approximate inverse (equivalent to a one-level RAS)  |
   | wjacobi  |  PFLAREINV_WJACOBI  | Weighted Jacobi  |
   | jacobi  |  PFLAREINV_JACOBI  | Jacobi  |

3) **PCAIR**: A new PETSc PC type, containing different reduction multigrids. PCAIR can be used with the command line argument ``-pc_type air``. Several different CF splittings are available, including the PMISR DDC method above, along with several different grid-transfer operators and smoothers. There are several features used to improve the parallel performance of PCAIR:

   - The number of active MPI ranks on lower levels can be reduced where necessary. If this is used then:
     - Repartitioning with graph partitioners can be applied.
     - OpenMP can be used in the polynomial inverse assembly (i.e., AIRG or nAIR) to reduce setup time (without requiring support for non-busy waits in the MPI library).
     - Calculation of polynomial coefficients can be done on subcommunicators.
   - The PCPFLAREINV methods above can be used as parallel coarse grid solvers, allowing heavy truncation of the multigrid hierarchy.
   - The sparsity of the multigrid hierarchy (and hence the CF splitting, repartitioning and symbolic matrix-matrix products) can be reused during setup. 

   The combination of ``-pc_air_z_type`` and ``-pc_air_inverse_type`` (defined by the PCPFLAREINV types above) defines several different reduction multigrids:

   | ``-pc_air_z_type``  | ``-pc_air_inverse_type`` | Description |
   | ------------- | -- | ------------- |
   | product  |  power, arnoldi or newton  | AIRG  |
   | product  |  neumann  | nAIR with Neumann smoothing  |
   | product  |  sai  | SAI reduction multigrid  |
   | product  |  isai  | ISAI reduction multigrid  |
   | product  |  wjacobi or jacobi  | Distance 0 reduction multigrid  |
   | lair  |  wjacobi or jacobi  | lAIR  |
   | lair_sai  |  wjacobi or jacobi  | SAI version of lAIR  |

   Different combinations of these types can also be used, e.g., ``-pc_air_z_type lair -pc_air_inverse_type power`` uses a lAIR grid transfer operator and GMRES polynomial smoothing with the power basis.

## Building PFLARE

This library depends on MPI, BLAS, LAPACK (>= 3.4) and PETSc (3.14 to 3.22) configured with a graph partitioner (e.g., ParMETIS). Please compile PETSc directly from the source code, as PFLARE requires access to some of the PETSc types only available in the source.

1) Set `PETSC_DIR` and `PETSC_ARCH` environmental variables.
2) Call ``make`` in the top level directory (you may need to modify the Makefile).
3) Call ``make python`` in the top level directory to build the Python module.

Then if desired:

4) Call ``make tests`` in the top level directory to check the build worked with some simple Fortran and C tests.
5) Call ``make tests_python`` in the top level directory to check the Python build worked with some simple Python tests.

PFLARE is compatible with 64-bit integers if PETSc has been configured with 64-bit integers. Note however that some of the tests use the PETSc Mat ``/data/mat_stream_2364``, which was output in PETSc format with 32-bit integers. Hence tests that load this Mat will crash if PETSc has been configured with 64-bit integers.

A Dockerfile is also provided which builds all the dependencies, compiles the library and runs all the tests. To run this Docker image, from the top level directory use:

     docker build -t "pflare" .
     docker run pflare

## Linking to PFLARE

1) For Fortran/C, link the library `libpflare.so` to your application; it is output to `lib/`. For Fortran, you must include `pflare.h` in `include/finclude/`, for C you must include `pflare.h` in `include/`.
2) For Python, ensure the full path to `python/` is in your `PYTHONPATH` environmental variable along with the path for PETSc. Your `LD_LIBRARY_PATH` must include the `lib/` directory (along with paths for PETSc, BLAS and LAPACK).

Using the components of PFLARE in an existing PETSc code is very simple. For C/Fortran, the user must call a single function which registers the new PC types with PETSc, while in Python this is handled by the import statement:

in Fortran:

     #include "finclude/pflare.h"
     ...
     call PCRegister_PFLARE()
     ...

or in C:

     #include "pflare.h"
     ...
     PCRegister_PFLARE();
     ...

or in Python with petsc4py:

     import pflare
     ...

## Using PFLARE

There are several different ways to use the methods in the PFLARE library.

### Using the new PC types:

Once the new PC types have been registered, they can then be used like native PETSc types, either by writing code to set the PETSc type/options, or through command line arguments. A few examples include:

#### 1) Using PCAIR with default options, namely AIRG with parameters tuned for a time independent advection equation on 2D unstructured triangles:

in Fortran:

     call KSPGetPC(ksp, pc, ierr)
     call PCSetType(pc, PCAIR, ierr)

     ...[e.g., KSPSolve somewhere here]

or in C:

     ierr = KSPGetPC(ksp, &pc);
     ierr = PCSetType(pc, PCAIR);

     ...[e.g., KSPSolve somewhere here]

or in Python with petsc4py:

     pc = ksp.getPC()
     pc.setType("air")

     ...[e.g., KSPSolve somewhere here]
     
or via the command line: ``-pc_type air``. 

#### 2) Using PCAIR to apply lAIR with FFC Weighted Jacobi smoothing:

in Fortran:

     call KSPGetPC(ksp, pc, ierr)
     call PCSetType(pc, PCAIR, ierr)

     call PCAIRSetZType(pc, AIR_Z_LAIR, ierr)
     call PCAIRSetInverseType(pc, PFLAREINV_WJACOBI, ierr)
     call PCAIRSetOneCSmooth(pc, PETSC_TRUE, ierr)

     ...[e.g., KSPSolve somewhere here]

or in C:

     ierr = KSPGetPC(ksp, &pc);
     ierr = PCSetType(pc, PCAIR);

     ierr = PCAIRSetZType(pc, AIR_Z_LAIR);
     ierr = PCAIRSetInverseType(pc, PFLAREINV_WJACOBI);
     ierr = PCAIRSetOneCSmooth(pc, PETSC_TRUE);

     ...[e.g., KSPSolve somewhere here]

or in Python with petsc4py:

     pc = ksp.getPC()
     pc.setType("air")

     petsc_options = PETSc.Options()
     petsc_options['pc_air_z_type'] = 'lair'
     petsc_options['pc_air_inverse_type'] = 'wjacobi'
     petsc_options['pc_air_one_c_smooth'] = ''     

     ...[e.g., KSPSolve somewhere here]
     
or via the command line: ``-pc_type air -pc_air_z_type lair -pc_air_inverse_type wjacobi -pc_air_one_c_smooth``.

#### 3) Using PCFLAREINV to apply a 20th order GMRES polynomial in the Newton basis matrix free:

in Fortran:

     call KSPGetPC(ksp, pc, ierr)
     call PCSetType(pc, PCPFLAREINV, ierr)

     call PCPFLAREINVSetOrder(pc, 20, ierr)
     call PCPFLAREINVSetType(pc, PFLAREINV_NEWTON, ierr)
     call PCPFLAREINVSetMatrixFree(pc, PETSC_TRUE, ierr)

     ...[e.g., KSPSolve somewhere here]

or in C:

     ierr = KSPGetPC(ksp, &pc);
     ierr = PCSetType(pc, PCPFLAREINV);

     ierr = PCPFLAREINVSetOrder(pc, 20);
     ierr = PCPFLAREINVSetType(pc, PFLAREINV_NEWTON);
     ierr = PCPFLAREINVSetMatrixFree(pc, PETSC_TRUE);

     ...[e.g., KSPSolve somewhere here]

or in Python with petsc4py:

     pc = ksp.getPC()
     pc.setType("pflareinv")

     petsc_options = PETSc.Options()
     petsc_options['pc_pflareinv_type'] = 'newton'
     petsc_options['pc_pflareinv_order'] = '20'
     petsc_options['pc_pflareinv_matrix_free'] = ''

     ...[e.g., KSPSolve somewhere here]
     
or via the command line: ``-pc_type pflareinv -pc_pflareinv_type newton -pc_pflareinv_order 20 -pc_pflareinv_matrix_free``.

## Reuse with PCAIR:

If solving multiple linear systems, there are several ways to reuse components of the PCAIR preconditioner to reduce setup times:

#### 1) Setup once 

In PETSc the ``KSPSetReusePreconditioner`` flag can be set to ensure the preconditioner is only setup once, even if the entries or sparsity of the underlying matrix in the KSP is changed. This can be useful in time stepping problems where the linear system (or Jacobian for non-linear problems) does not change, or if the changes are small. 

#### 2) Reuse sparsity during setup

When solving a linear system where the matrix has the same sparsity pattern as in a previous solve, PCAIR can can reuse the CF splitting, repartitioning, symbolic matrix-matrix products and resulting sparsity throughout the hierarchy during the setup. This takes more memory (typically 2-5x the storage complexity) but significantly reduces the time required in subsequent setups (10-20x). If we have a PETSc matrix $\mathbf{A}$:

in Fortran:

     call KSPGetPC(ksp, pc, ierr)
     call PCSetType(pc, PCAIR, ierr)
     
     ! Before the first setup, we tell it to 
     ! store any data we need for reuse
     call PCAIRSetReuseSparsity(pc, PETSC_TRUE, ierr)

     ! First solve - the PCAIR will be setup
     call KSPSolve(ksp, b, x, ierr)
     
     ...[Modify entries in A but keep the same sparsity]

     ! Second solve - the PCAIR will be setup
     ! but reusing the same sparsity
     call KSPSolve(ksp, b, x, ierr)     

or in C:

     ierr = KSPGetPC(ksp, &pc);
     ierr = PCSetType(pc, PCAIR);

     // Before the first setup, we tell it to 
     // store any data we need for reuse
     ierr = PCAIRSetReuseSparsity(pc, PETSC_TRUE);

     // First solve - the PCAIR will be setup
     ierr = KSPSolve(ksp, b, x)
     
     ...[Modify entries in A but keep the same sparsity]

     // Second solve - the PCAIR will be setup 
     // but reusing the same sparsity
     ierr = KSPSolve(ksp, b, x)

or in Python with petsc4py:

     pc = ksp.getPC()
     pc.setType("air")

     petsc_options = PETSc.Options()
     petsc_options['pc_air_reuse_sparsity'] = ''

     # First solve - the PCAIR will be setup
     ksp.solve(b,x)
     
     ...[Modify entries in A but keep the same sparsity]

     # Second solve - the PCAIR will be setup 
     # but reusing the same sparsity
     ksp.solve(b,x)
     
or via the command line: ``-pc_type air -pc_air_reuse_sparsity``.

#### 3) Reuse sparsity and GMRES polynomial coefficients during setup with AIRG

When using AIRG in PCAIR (which is the default), in addition to reusing the sparsity in the setup, the GMRES polynomial coefficients can also be reused. This is useful if the matrix is very similar to that used in a previous solve and saves parallel reductions in the setup. Note that the GMRES polynomial coefficients are very sensitive to changes in the matrix, so this may not always give good results. To enable this we simply set another reuse flag before (or after) the first solve:

in Fortran:

     ! Before the first setup, we tell it to 
     ! store any data we need for reuse
     call PCAIRSetReuseSparsity(pc, PETSC_TRUE, ierr)

     ! Tell it to store the gmres polynomial coefficients
     call PCAIRSetReusePolyCoeffs(pc, PETSC_TRUE, ierr)        
    
     ! First solve - the PCAIR will be setup
     call KSPSolve(ksp, b, x, ierr)  
     
     ...[Modify entries in A but keep the same sparsity]

     ! Second solve - the PCAIR will be setup
     ! but reusing sparsity & polynomial coefficients
     call KSPSolve(ksp, b, x, ierr)     

or in C:

     // Before the first setup, we tell it to 
     // store any data we need for reuse
     ierr = PCAIRSetReuseSparsity(pc, PETSC_TRUE);

     ! Tell it to store the gmres polynomial coefficients
     ierr = PCAIRSetReusePolyCoeffs(pc, PETSC_TRUE);

     // First solve - the PCAIR will be setup
     ierr = KSPSolve(ksp, b, x)
     
     ...[Modify entries in A but keep the same sparsity]

     // Second solve - the PCAIR will be setup
     // but reusing sparsity & polynomial coefficients
     ierr = KSPSolve(ksp, b, x)

or in Python with petsc4py:

     pc = ksp.getPC()
     pc.setType("air")

     petsc_options = PETSc.Options()
     petsc_options['pc_air_reuse_sparsity'] = ''
     petsc_options['pc_air_reuse_poly_coeffs'] = ''

     # First solve - the PCAIR will be setup
     ksp.solve(b,x)
     
     ...[Modify entries in A but keep the same sparsity]

     # Second solve - the PCAIR will be setup 
     # but reusing sparsity & polynomial coefficients
     ksp.solve(b,x)

or via the command line: ``-pc_type air -pc_air_reuse_sparsity -pc_air_reuse_poly_coeffs``.

#### 4) Reuse when solving the same linear system with AIRG

Often an outer loop (e.g., eigenvalue) will require solving a series of different linear systems one after another, and then returning to the first linear system to start the loop again. If the linear systems are close enough to each other to get good performance from reusing the sparsity, but different enough from each other that reusing the GMRES polynomials coefficients from a different linear system gives poor performance, PCAIR allows the GMRES polynomial coefficients to be saved externally and then restored before a solve. 

This means we can cheaply generate the exact same preconditioner when periodically solving the same linear system. An example of this is given in ``tests/ex6f_getcoeffs.F90``, where we solve a linear system, store the resulting GMRES polynomial coefficients, reuse the sparsity to solve a different linear system, then reuse the sparsity and restore the GMRES polynomial coefficients to solve the first linear system again. We should note this type of reuse is not yet available in Python. 

## OpenMP with PCAIR

If processor agglomeration has been enabled (it is by default) and:
1) ``-pc_air_matrix_free_polys`` has not been set (default).
2) ``-pc_air_inverse_type`` is one of: power (default), arnoldi, newton or neumann.
3) ``-pc_air_inverse_sparsity_order=1`` (default)

then OpenMP can be used to reduce the time required to assemble the fixed-sparsity polynomial inverses. The idle MPI ranks that have been left with 0 unknowns after repartitioning will be used to thread. 

The number of threads will automatically increase on lower grids with processor agglomeration. For example, if the number of active MPI ranks is halved each time processor agglomeration is triggered (``-pc_air_processor_agglom_factor 2``), then 2 threads will be used after the first processor agglomeration, with 4 threads used after the second, etc. 

This relies on a repartitioning where the MPI ranks that are made inactive are on the same node (and NUMA region), i.e., each node has fewer active MPI ranks after repartitioning (typically called an "interleaved" partitioning), rather than some nodes being full and others empty. This will depend on how the ranks are numbered by the MPI library. Good performance is also dependent on appropriate pinning of MPI ranks and OpenMP threads to hardware cores. 

Given that the time required to assemble the approximate polynomial inverses is typically small, using OpenMP often has little impact on the overall setup time. If however reuse has been enabled with PCAIR (see above), then the approximate polynomial inverse time is often the largest component of the (reused) setup time and OpenMP can help.

To build PFLARE with OpenMP, add the appropriate OpenMP compiler flag before calling make, e.g., with GNU

     export CFLAGS="-fopenmp"
     export FFLAGS="-fopenmp" 
     make

and then before running a problem with PFLARE, set the ``OMP_NUM_THREADS`` environmental variable to be the maximum number of threads to use; typically this would be the number of hardware cores per NUMA region.

It is recommended that PFLARE be linked with unthreaded BLAS/LAPACK libraries, as there is often little gain from using threaded libraries with PFLARE and that ensures the ``OMP_NUM_THREADS`` environmental variable is used purely to control the OpenMP described above. 

On Cray machines, adding the OpenMP compiler flag typically triggers linking to threaded libraries (e.g., BLAS/LAPACK). To work around this, it is recommended that PFLARE is built as a shared library (which is the default) with the compiler flag enabled, then any code which uses PFLARE is compiled without the compiler flag but with the OpenMP libraries included in the ``LDFLAGS``. This will ensure PFLARE can use OpenMP, but the unthreaded libraries are used to link. 

For example, on ARCHER2 (a HPE Cray EX machine) to build PFLARE and the tests with GNU:

     # Cray MPI compiler wrappers
     export CC=cc
     export FC=ftn

     # Build PFLARE
     export CFLAGS="-fopenmp"
     export FFLAGS="-fopenmp" 
     make CC=${CC} FC=${FC}

     # Build the tests
     export CFLAGS=""
     export FFLAGS="" 
     export LDFLAGS="-lgomp"
     make build_tests CC=${CC} FC=${FC} 

Running a test with OpenMP then requires setting the ``OMP_NUM_THREADS`` variable, ensuring the MPI ranks and OpenMP threads are correctly pinned and telling the queueing system that the problem is oversubscribed. For example, the (partial) slurm script to run ex62f on one node is:

      #SBATCH --time=0:01:0
      #SBATCH --nodes=1
      #SBATCH --ntasks-per-node=128
      #SBATCH --overcommit
      #SBATCH --partition=standard
      #SBATCH --qos=standard

      # Give it a maximum of 4 threads
      export OMP_NUM_THREADS=4

      [ Set the hex mask which describes the pinning ]

      srun --oversubscribe --distribution=block:block --cpu-bind=mask_cpu:$mask ex62f.o


## CF splittings

The CF splittings in PFLARE are used within PCAIR to form the multigrid hierarchy. There are several different CF splittings available with ``-pc_air_cf_splitting_type``:

   | Command line type  | Flag | Description |
   | ------------- | -- | ------------- |
   | pmisr_ddc  |  CF_PMISR_DDC  | Two-pass splitting giving diagonally dominant $\mathbf{A}_\textrm{ff}$ |
   | pmis  |  CF_PMIS  | PMIS method with symmetrized strength matrix |
   | pmis_dist2  |  CF_PMIS_DIST2  | Distance 2 PMIS method with strength matrix formed by S'S + S and then symmetrized |
   | agg  |  CF_AGG  | Aggregation method with root-nodes as C points. In parallel this is processor local aggregation  |
   | pmis_agg  |  CF_PMIS_AGG  | PMIS method with symmetrized strength matrix on boundary nodes, then processor local aggregation.  |   

The CF splittings can also be called independently from PCAIR. The CF splittings are returned in two PETSc IS's representing the coarse and fine points. For example, to compute a PMISR DDC CF splitting of a PETSc matrix $\mathbf{A}$:

in Fortran:

     IS :: is_fine, is_coarse
     ! Threshold for a strong connection
     PetscReal :: strong_threshold = 0.5
     ! Second pass cleanup
     PetscReal :: ddc_fraction = 0.1
     ! As many steps as needed
     PetscInt :: max_luby_steps = -1
     ! PMISR DDC
     integer :: algorithm = CF_PMISR_DDC
     ! Is the matrix symmetric?
     logical :: symmetric = .FALSE.
     ...
     call compute_cf_splitting(A, &
           symmetric, &
           strong_threshold, max_luby_steps, &
           algorithm, &
           ddc_fraction, &
           is_fine, is_coarse) 

or in C (please note the slightly modified name in C):

     IS is_fine, is_coarse;
     // Threshold for a strong connection
     PetscReal strong_threshold = 0.5;
     // Second pass cleanup
     PetscReal ddc_fraction = 0.1;
     // As many steps as needed
     PetscInt max_luby_steps = -1;
     // PMISR DDC
     int algorithm = CF_PMISR_DDC;
     // Is the matrix symmetric?
     int symmetric = 0;

     compute_cf_splitting_c(&A, \
         symmetric, \
         strong_threshold, max_luby_steps, \
         algorithm, \
         ddc_fraction, \
         &is_fine, &is_coarse);

or in Python with petsc4py:

     # Threshold for a strong connection
     strong_threshold = 0.5;
     # Second pass cleanup
     ddc_fraction = 0.1;
     # As many steps as needed
     max_luby_steps = -1;
     # PMISR DDC
     algorithm = pflare.CF_PMISR_DDC;
     # Is the matrix symmetric?
     symmetric = False;

     [is_fine, is_coarse] = pflare.pflare_defs.compute_cf_splitting(A, \
           symmetric, \
           strong_threshold, max_luby_steps, \
           algorithm, \
           ddc_fraction)
     
## More examples

For more ways to use the library please see the Fortran/C examples and the Makefile in `tests/`, along with the Python examples in `python/`. All the available options and default values are listed in ``src/PCAIR.c`` and ``src/PCPFLAREINV.c``.

## References

Please see these references for more details:

[1] S. Dargaville, R. P. Smedley-Stevenson, P. N. Smith, C. C. Pain, AIR multigrid with GMRES polynomials (AIRG) and additive preconditioners for Boltzmann transport, Journal of Computational Physics 518 (2024) 113342.  
[2] S. Dargaville, R. P. Smedley-Stevenson, P. N. Smith, C. C. Pain, Coarsening and parallelism with reduction multigrids for hyperbolic Boltzmann transport, arXiv preprint arXiv:2408.08262, (2024).  
[3] T. A. Manteuffel, S. Münzenmaier, J. Ruge, B. Southworth, Nonsymmetric Reduction-Based Algebraic Multigrid, SIAM Journal on Scientific Computing 41 (2019) S242–S268.  
[4] T. A. Manteuffel, J. Ruge, B. S. Southworth, Nonsymmetric algebraic multigrid based on local approximate ideal restriction (lAIR), SIAM Journal on Scientific Computing 40 (2018) A4105–A4130.  
[5] A. Ali, J. J. Brannick, K. Kahl, O. A. Krzysik, J. B. Schroder, B. S. Southworth, Constrained local approximate ideal restriction for advection-diffusion problems, SIAM Journal on Scientific Computing (2024) S96–S122.  
[6] T. Zaman, N. Nytko, A. Taghibakhshi, S. MacLachlan, L. Olson, M. West, Generalizing reduction-based algebraic multigrid, Numerical Linear
Algebra with Applications 31 (3) (2024) e2543.
