#if !defined (PETSC_LEGACY_H)
#define PETSC_LEGACY_H

! This header file should be included in fortran modules that use petsc
! directly after the implicit none. It provides legacy support for building
! petsc versions older than the latest released petsc. Where names
! have changed this #defines the old name as its newer equivalent, so that new
! names can be used in the code everywhere. Where interfaces have changed we
! still need #ifdef PETSC_VERSION>... in the main code
#include "petscversion.h"
#include "petsc/finclude/petsc.h"

! Defined in PETSC 3.22 and above
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<22)
#define PETSC_NULL_INTEGER_ARRAY PETSC_NULL_INTEGER
#define PETSC_NULL_SCALAR_ARRAY PETSC_NULL_SCALAR
#define PETSC_NULL_REAL_ARRAY PETSC_NULL_REAL
#endif

! PETSc 3.22 changed how to test for null things
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<22)
#define PetscMatIsNull(mat) mat==PETSC_NULL_MAT
#define PetscISIsNull(is) is==PETSC_NULL_IS
#define PetscNullspaceIsNull(nullspace) nullspace==PETSC_NULL_MATNULLSPACE
#else
#define PetscMatIsNull(mat) PetscObjectIsNull(mat)
#define PetscISIsNull(is) PetscObjectIsNull(is)
#define PetscNullspaceIsNull(nullspace) PetscObjectIsNull(nullspace)
#endif

#endif