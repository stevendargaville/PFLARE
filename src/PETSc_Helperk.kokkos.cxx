#include <petscvec_kokkos.hpp>
#include <petscksp.h>
#include <iostream>
#include <../src/mat/impls/aij/seq/aij.h>
#include <../src/mat/impls/aij/mpi/mpiaij.h>
#include <../src/mat/impls/aij/seq/kokkos/aijkok.hpp>

using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;
using DefaultMemorySpace    = Kokkos::DefaultExecutionSpace::memory_space;
using PetscScalarKokkosView = Kokkos::View<PetscScalar *, DefaultMemorySpace>;

// Generate identity but with kokkos - trying to keep everything on the device where possible
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

//------------------------------------------------------------------------------------------------------------------------

// Drop according to a tolerance but with kokkos - trying to keep everything on the device where possible
PETSC_INTERN void remove_small_from_sparse_kokkos(Mat *input_mat, PetscReal tol, Mat *output_mat, \
                  int relative_max_row_tolerance_int, int lump_int, int allow_drop_diagonal_int)
{
   MPI_Comm MPI_COMM_MATRIX;
   PetscInt local_rows, local_cols, global_rows, global_cols;
   PetscInt global_row_start, global_row_end_plus_one;
   PetscInt global_col_start, global_col_end_plus_one;
   MatType mat_type;
   PetscReal *rel_row_tol;
   PetscInt max_nnzs, max_nnzs_total, ncols, ncols_seq_local, ncols_seq_nonlocal;
   PetscInt nnzs_local, nnzs_nonlocal;
   const PetscScalar *local_data, *nonlocal_data;

   MatGetType(*input_mat, &mat_type);
   bool mpi = strcmp(mat_type, MATMPIAIJKOKKOS) == 0;

   // ~~~~~~~~~~~
   // This gets a little messy, as kokkos seems to be *very very* slow 
   // when we call matgetrow
   // So we try and avoid that in this routine
   // When we have a kokkos matrix, the ->A and ->B are of type AIJMat
   // as that is the copy on the host that is kept in sync with the device mat
   // ~~~~~~~~~~~
   Mat_MPIAIJ *mat_mpi;
   Mat_SeqAIJ *mat_seq_local, *mat_seq_nonlocal;
   Mat mat_local, mat_nonlocal;
   if (mpi)
   {
      mat_mpi = (Mat_MPIAIJ *)(*input_mat)->data;
      mat_seq_local = (Mat_SeqAIJ *)mat_mpi->A->data;
      mat_seq_nonlocal = (Mat_SeqAIJ *)mat_mpi->B->data;
      mat_local = mat_mpi->A;
      mat_nonlocal = mat_mpi->B;
   }
   else
   {
      mat_seq_local = (Mat_SeqAIJ *)(*input_mat)->data;
      mat_local = *input_mat;
   }

   // Get the comm
   PetscObjectGetComm((PetscObject)*input_mat, &MPI_COMM_MATRIX);
   MatGetLocalSize(*input_mat, &local_rows, &local_cols);
   MatGetSize(*input_mat, &global_rows, &global_cols);
   // This returns the global index of the local portion of the matrix
   MatGetOwnershipRange(*input_mat, &global_row_start, &global_row_end_plus_one);
   MatGetOwnershipRangeColumn(*input_mat, &global_col_start, &global_col_end_plus_one);

   // Get the number of nnzs
   max_nnzs = 0;
   max_nnzs_total = 0;
   nnzs_local = 0;
   nnzs_nonlocal = 0;
   for (int i = 0; i < local_rows; i++)   
   {
      ncols = mat_seq_local->i[i + 1] - mat_seq_local->i[i];
      nnzs_local += ncols;
      if (mpi)
      {
         ncols += mat_seq_nonlocal->i[i + 1] - mat_seq_nonlocal->i[i];
         nnzs_nonlocal += mat_seq_nonlocal->i[i + 1] - mat_seq_nonlocal->i[i];
      }

      if (ncols > max_nnzs) max_nnzs = ncols;
   }   
   max_nnzs_total = nnzs_local + nnzs_nonlocal;

   // We know we never have more to do than the original nnzs
   PetscInt *row_indices, *col_indices;
   PetscMalloc2(max_nnzs_total, &row_indices, max_nnzs_total, &col_indices);
   // By default drop everything
   for (int i = 0; i < max_nnzs_total; i++)
   {
      row_indices[i] = -1;
   }

   MatCreate(MPI_COMM_MATRIX, output_mat);
   MatSetSizes(*output_mat, local_rows, local_cols, global_rows, global_cols);
   // Match the output type
   MatSetType(*output_mat, mat_type);
   MatSetUp(*output_mat);

   // Don't set any off processor entries so no need for a reduction when assembling
   MatSetOption(*output_mat, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
   MatSetOption(*output_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
   MatSetOption(*output_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);     
  
   // Gets a pointer to the vals data
   MatSeqAIJGetArrayRead(mat_local, &local_data);
   if (mpi) MatSeqAIJGetArrayRead(mat_nonlocal, &nonlocal_data);

   // Relative dropping tolerance per row
   PetscMalloc1(local_rows, &rel_row_tol);   

   // Let's calculate the relative row tolerances on the device
   {
      // Let's create device memory with a dual view
      auto &exec = PetscGetKokkosExecutionSpace();
      MatScalarKokkosViewHost rel_row_tol_h(rel_row_tol, local_rows);
      auto rel_row_tol_d = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, exec, rel_row_tol_h);
      auto rel_row_tol_dual = MatScalarKokkosDualView(rel_row_tol_d, rel_row_tol_h);

      // Set device data to tol
      Kokkos::deep_copy(rel_row_tol_d, tol);   

      // Get the i, j and vals device pointers 
      const PetscInt *device_local_i, *device_local_j, *device_nonlocal_i, *device_nonlocal_j;
      PetscMemType mtype;
      PetscScalar *device_local_vals, *device_nonlocal_vals;
      MatSeqAIJGetCSRAndMemType(mat_local, &device_local_i, &device_local_j, &device_local_vals, &mtype);   
      if (mpi) MatSeqAIJGetCSRAndMemType(mat_nonlocal, &device_nonlocal_i, &device_nonlocal_j, &device_nonlocal_vals, &mtype);  
      
      if (relative_max_row_tolerance_int) 
      {       
         Kokkos::parallel_for(
            Kokkos::RangePolicy<>(0, local_rows), KOKKOS_LAMBDA(int i) {

               PetscInt ncols_local = device_local_i[i + 1] - device_local_i[i];
               PetscScalar max_val = -1.0;

               for (int j = 0; j < ncols_local; j++)
               {
                  if (abs(device_local_vals[device_local_i[i] + j]) > max_val) max_val = abs(device_local_vals[device_local_i[i] + j]);
               }

               if (mpi)
               {
                  PetscInt ncols_nonlocal = device_nonlocal_i[i + 1] - device_nonlocal_i[i];               

                  for (int j = 0; j < ncols_nonlocal; j++)
                  {
                     if (abs(device_nonlocal_vals[device_nonlocal_i[i] + j]) > max_val) max_val = abs(device_nonlocal_vals[device_nonlocal_i[i] + j]);
                  }  
               }              

               rel_row_tol_d(i) *= max_val;
            });
      }

      // We've modified the data on the device
      rel_row_tol_dual.modify_device();
      // Let's copy the data back to the host
      rel_row_tol_dual.sync_host(exec);

      // Now go and fill the new matrix
      // These loops just set the row and col indices to not be -1
      // if we are including it in the matrix
      // Loop over global row indices
      PetscInt counter = 0;

      // Go and do the local part
      for (int i = 0; i < local_rows; i++)
      {
         ncols_seq_local = mat_seq_local->i[i + 1] - mat_seq_local->i[i];

         // Do the local part
         for (int j = 0; j < ncols_seq_local; j++)
         {
            // Set the row/col to be included (ie not -1) 
            // if it is bigger than the tolerance
            if (abs(local_data[mat_seq_local->i[i] + j]) >= rel_row_tol[i])
            {
               row_indices[counter + j] = i + global_row_start;
               // Careful here to use global_col_start in case we are rectangular
               col_indices[counter + j] = mat_seq_local->j[mat_seq_local->i[i] + j] + global_col_start;
            }
            // If the entry is small and we are lumping, then add it to the diagonal
            // or if this is the diagonal and it's small but we are not dropping it 
            else if (lump_int || \
                  (!allow_drop_diagonal_int && \
                     mat_seq_local->j[mat_seq_local->i[i] + j] + global_col_start == i + global_row_start))
            {
               row_indices[counter + j] = i + global_row_start;
               col_indices[counter + j] = i + global_row_start;
            }
         }     
         counter = counter + ncols_seq_local;
      }
      MatSeqAIJRestoreArrayRead(mat_local, &local_data);

      // Do the non-local part
      if (mpi)
      {
         for (int i = 0; i < local_rows; i++)
         {
            ncols_seq_nonlocal = mat_seq_nonlocal->i[i + 1] - mat_seq_nonlocal->i[i];

            for (int j = 0; j < ncols_seq_nonlocal; j++)
            {
               // Set the row/col to be included (ie not -1) 
               // if it is bigger than the tolerance
               if (abs(nonlocal_data[mat_seq_nonlocal->i[i] + j]) >= rel_row_tol[i])
               {
                  row_indices[counter + j] = i + global_row_start;
                  // garray is the colmap
                  col_indices[counter + j] = mat_mpi->garray[mat_seq_nonlocal->j[mat_seq_nonlocal->i[i] + j]];
               }                 
               // Be careful to use the colmap to get the off-diagonal global column index
               else if (lump_int || \
                     (!allow_drop_diagonal_int && \
                        mat_mpi->garray[mat_seq_nonlocal->j[mat_seq_nonlocal->i[i] + j]] \
                           == i + global_row_start))
               {
                  row_indices[counter + j] = i + global_row_start;
                  col_indices[counter + j] = i + global_row_start;
               }     
            }
            counter = counter + ncols_seq_nonlocal;
         } 
         MatSeqAIJRestoreArrayRead(mat_nonlocal, &nonlocal_data);         
      }

      // Now set the sparsity
      MatSetPreallocationCOO(*output_mat, counter, row_indices, col_indices);
      PetscFree2(row_indices, col_indices);   
   }
   PetscFree(rel_row_tol);

   // For the values, we can just directly get access to the device pointers
   {
      const PetscInt *i, *j;
      PetscMemType mtype;
      PetscScalar *device_vals;
      MatSeqAIJGetCSRAndMemType(mat_local, &i, &j, &device_vals, &mtype);

      // Have to copy the local and then the nonlocal data into one array on the device
      if (mpi)
      {         
         // coo_v is stored on the device
         PetscScalarKokkosView coo_v("coo_v", max_nnzs_total);    

         // Copy the local values - could also do this with views, subviews
         // and calling MatSeqAIJGetKokkosView, but we know our device_vals 
         // and coo_v share the same execution space
         Kokkos::parallel_for(
            Kokkos::RangePolicy<>(0, nnzs_local), KOKKOS_LAMBDA(int i) {
               coo_v(i) = device_vals[i];
            });    

         MatSeqAIJGetCSRAndMemType(mat_nonlocal, &i, &j, &device_vals, &mtype);
         // Copy the nonlocal values
         Kokkos::parallel_for(
            Kokkos::RangePolicy<>(0, nnzs_nonlocal), KOKKOS_LAMBDA(int i) {
               coo_v(nnzs_local + i) = device_vals[i];
            });           

         // This should all happen on the gpu
         MatSetValuesCOO(*output_mat, coo_v.data(), INSERT_VALUES);             
      }
      // Don't have to do a copy if we have no non-local data
      else
      {
         // This should all happen on the gpu
         MatSetValuesCOO(*output_mat, device_vals, INSERT_VALUES);           
      }
   }  

   return;
}