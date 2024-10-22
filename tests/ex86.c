static char help[] = "Solves a one-dimensional steady upwind advection system with KSP.\n\n";

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h    - base PETSc routines   petscvec.h - vectors
     petscmat.h    - matrices              petscpc.h  - preconditioners
     petscis.h     - index sets
     petscviewer.h - viewers
*/
#include <petscksp.h>
#include "pflare.h"

int main(int argc, char **args)
{
  Vec         x, b, work_vec; /* approx solution, RHS, work vector */
  Mat         A;              /* linear system matrix */
  KSP         ksp;            /* linear solver context */
  PC          pc;             /* preconditioner context */
  PetscInt    i, j, n = 10, col[2], its;
  PetscScalar work_scalar, value[2];
  PetscRandom r;
  KSPConvergedReason reason;

  PetscFunctionBeginUser;
  PetscInitialize(&argc, &args, (char*)0, help);

  PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL);

  // Register the pflare types
  PCRegister_PFLARE();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed.
  */
  VecCreate(PETSC_COMM_WORLD, &x);
  PetscObjectSetName((PetscObject)x, "Solution");
  VecSetSizes(x, PETSC_DECIDE, n);
  VecSetFromOptions(x);
  VecDuplicate(x, &b);
  VecDuplicate(x, &work_vec);

  /*
     Create matrix.  When using MatCreate(), the matrix format can
     be specified at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetFromOptions(A);
  MatSetUp(A);

  /*
     Assemble matrix
  */
  value[0] = -1.0;
  value[1] = 1.0;
  for (i = 1; i < n; i++) {
    col[0] = i - 1;
    col[1] = i;
    MatSetValues(A, 1, &i, 2, col, value, INSERT_VALUES);
  }
  i           = 0;
  j           = 0;
  work_scalar = 1;
  MatSetValues(A, 1, &i, 1, &j, &work_scalar, INSERT_VALUES);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  /*
     Create x, b
  */
  PetscRandomCreate(PETSC_COMM_WORLD, &r);
  VecSetRandom(x, r);
  VecSetRandom(work_vec, r);
  MatMult(A, work_vec, b);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  KSPCreate(PETSC_COMM_WORLD, &ksp);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the matrix that defines the preconditioner.
  */
  KSPSetOperators(ksp, A, A);

  /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following four statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions();
  */
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCAIR);

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  KSPSetFromOptions(ksp);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  KSPSolve(ksp, b, x);

  KSPGetIterationNumber(ksp,&its);
  KSPGetConvergedReason(ksp,&reason);
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 17)
      PetscCheck(reason > 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Didn't converge");
#else
      if (reason < 0) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Didn't converge");
#endif   
  printf("iterations %d \n", its);

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  VecDestroy(&x);
  VecDestroy(&b);
  VecDestroy(&work_vec);
  MatDestroy(&A);
  KSPDestroy(&ksp);
  PetscRandomDestroy(&r);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_view).
  */
  PetscFinalize();
  return 0;
}
