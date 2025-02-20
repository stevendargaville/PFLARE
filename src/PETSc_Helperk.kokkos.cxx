#include <petscvec_kokkos.hpp>
#include <petscksp.h>
#include <iostream>

using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;
using DefaultMemorySpace    = Kokkos::DefaultExecutionSpace::memory_space;
using PetscScalarKokkosView = Kokkos::View<PetscScalar *, DefaultMemorySpace>;

// 
PETSC_INTERN void generate_identity_is_kokkos(Mat *input_mat, IS *indices, Mat *output_mat)
{
   MPI_Comm MPI_COMM_MATRIX;
   PetscInt local_indices_size, local_rows, local_cols, global_rows, global_cols;
   PetscInt global_row_start, global_row_end_plus_one;
   MatType mat_type;
   const PetscInt *is_pointer;
   PetscInt *oor, *ooc;   

   // Get the comm
   PetscObjectGetComm((PetscObject)*input_mat, &MPI_COMM_MATRIX);

   // Get the local sizes
   ISGetLocalSize(*indices, &local_indices_size);

   MatGetLocalSize(*input_mat, &local_rows, &local_cols);
   MatGetSize(*input_mat, &global_rows, &global_cols);
   // This returns the global index of the local portion of the matrix
   MatGetOwnershipRange(*input_mat, &global_row_start, &global_row_end_plus_one);

   MatCreate(MPI_COMM_MATRIX, output_mat);
   MatSetSizes(*output_mat, local_rows, local_cols, \
                       global_rows, global_cols);
   // Match the output type
   MatGetType(*input_mat, &mat_type);
   MatSetType(*output_mat, mat_type);
   MatSetUp(*output_mat);

   // Don't set any off processor entries so no need for a reduction when assembling
   MatSetOption(*output_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
   MatSetOption(*output_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);   

   // Get the indices we need
   ISGetIndices(*indices, &is_pointer);

   // Row and column indices for COO assembly
   // MatSetPreallocationCOO could modify the values so have to copy the values of 
   // is_pointer
   PetscMalloc2(local_indices_size, &oor, local_indices_size, &ooc);
   for (int i = 0; i < local_indices_size; i++)
   {
      oor[i] = is_pointer[i];
      ooc[i] = is_pointer[i];
   }

   // Set the diagonal
   MatSetPreallocationCOO(*output_mat, local_indices_size, oor, ooc);
   // Can delete oor and ooc now
   PetscFree2(oor, ooc);   
   ISRestoreIndices(*indices, &is_pointer);

   {
      // coo_v is stored on the device
      PetscScalarKokkosView coo_v("coo_v", local_indices_size);
      // Set coo_v to one
      Kokkos::deep_copy(coo_v, 1.0);   
      // This should all happen on the gpu
      MatSetValuesCOO(*output_mat, coo_v.data(), INSERT_VALUES);    
   }

   return;
}