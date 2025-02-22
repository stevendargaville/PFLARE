#include <petscvec_kokkos.hpp>
#include <petscksp.h>
#include <iostream>
#include <../src/mat/impls/aij/seq/aij.h>
#include <../src/mat/impls/aij/mpi/mpiaij.h>
#include <../src/mat/impls/aij/seq/kokkos/aijkok.hpp>

using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;
using DefaultMemorySpace    = Kokkos::DefaultExecutionSpace::memory_space;
using PetscScalarKokkosView = Kokkos::View<PetscScalar *, DefaultMemorySpace>;
using PetscIntKokkosView     = Kokkos::View<PetscInt *, DefaultMemorySpace>;
using PetscIntKokkosViewHost = Kokkos::View<PetscInt *, Kokkos::HostSpace>;
using PetscIntKokkosDualView = Kokkos::DualView<PetscInt *>;

// Horrid copy of this code and MatSetSeqAIJKokkosWithCSRMatrix_mine as they're declared PETSC_INTERN
// so I can't get at them
PetscErrorCode MatSeqAIJSetPreallocation_SeqAIJ_mine(Mat B, PetscInt nz, const PetscInt *nnz)
{
  Mat_SeqAIJ *b              = (Mat_SeqAIJ *)B->data;
  PetscBool   skipallocation = PETSC_FALSE, realalloc = PETSC_FALSE;
  PetscInt    i;

  PetscFunctionBegin;
  if (B->hash_active) {
    B->ops[0] = b->cops;
    PetscCall(PetscHMapIJVDestroy(&b->ht));
    PetscCall(PetscFree(b->dnz));
    B->hash_active = PETSC_FALSE;
  }
  if (nz >= 0 || nnz) realalloc = PETSC_TRUE;
  if (nz == MAT_SKIP_ALLOCATION) {
    skipallocation = PETSC_TRUE;
    nz             = 0;
  }
  PetscCall(PetscLayoutSetUp(B->rmap));
  PetscCall(PetscLayoutSetUp(B->cmap));

  if (nz == PETSC_DEFAULT || nz == PETSC_DECIDE) nz = 5;
  PetscCheck(nz >= 0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "nz cannot be less than 0: value %" PetscInt_FMT, nz);
  if (nnz) {
    for (i = 0; i < B->rmap->n; i++) {
      PetscCheck(nnz[i] >= 0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "nnz cannot be less than 0: local row %" PetscInt_FMT " value %" PetscInt_FMT, i, nnz[i]);
      PetscCheck(nnz[i] <= B->cmap->n, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "nnz cannot be greater than row length: local row %" PetscInt_FMT " value %" PetscInt_FMT " rowlength %" PetscInt_FMT, i, nnz[i], B->cmap->n);
    }
  }

  B->preallocated = PETSC_TRUE;
  if (!skipallocation) {
    if (!b->imax) { PetscCall(PetscMalloc1(B->rmap->n, &b->imax)); }
    if (!b->ilen) {
      /* b->ilen will count nonzeros in each row so far. */
      PetscCall(PetscCalloc1(B->rmap->n, &b->ilen));
    } else {
      PetscCall(PetscMemzero(b->ilen, B->rmap->n * sizeof(PetscInt)));
    }
    if (!b->ipre) PetscCall(PetscMalloc1(B->rmap->n, &b->ipre));
    if (!nnz) {
      if (nz == PETSC_DEFAULT || nz == PETSC_DECIDE) nz = 10;
      else if (nz < 0) nz = 1;
      nz = PetscMin(nz, B->cmap->n);
      for (i = 0; i < B->rmap->n; i++) b->imax[i] = nz;
      PetscCall(PetscIntMultError(nz, B->rmap->n, &nz));
    } else {
      PetscInt64 nz64 = 0;
      for (i = 0; i < B->rmap->n; i++) {
        b->imax[i] = nnz[i];
        nz64 += nnz[i];
      }
      PetscCall(PetscIntCast(nz64, &nz));
    }

    /* allocate the matrix space */
    PetscCall(MatSeqXAIJFreeAIJ(B, &b->a, &b->j, &b->i));
    PetscCall(PetscShmgetAllocateArray(nz, sizeof(PetscInt), (void **)&b->j));
    PetscCall(PetscShmgetAllocateArray(B->rmap->n + 1, sizeof(PetscInt), (void **)&b->i));
    b->free_ij = PETSC_TRUE;
    if (B->structure_only) {
      b->free_a = PETSC_FALSE;
    } else {
      PetscCall(PetscShmgetAllocateArray(nz, sizeof(PetscScalar), (void **)&b->a));
      b->free_a = PETSC_TRUE;
    }
    b->i[0] = 0;
    for (i = 1; i < B->rmap->n + 1; i++) b->i[i] = b->i[i - 1] + b->imax[i - 1];
  } else {
    b->free_a  = PETSC_FALSE;
    b->free_ij = PETSC_FALSE;
  }

  if (b->ipre && nnz != b->ipre && b->imax) {
    /* reserve user-requested sparsity */
    PetscCall(PetscArraycpy(b->ipre, b->imax, B->rmap->n));
  }

  b->nz               = 0;
  b->maxnz            = nz;
  B->info.nz_unneeded = (double)b->maxnz;
  if (realalloc) PetscCall(MatSetOption(B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
  B->was_assembled = PETSC_FALSE;
  B->assembled     = PETSC_FALSE;
  /* We simply deem preallocation has changed nonzero state. Updating the state
     will give clients (like AIJKokkos) a chance to know something has happened.
  */
  B->nonzerostate++;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PETSC_INTERN PetscErrorCode MatSetSeqAIJKokkosWithCSRMatrix_mine(Mat A, Mat_SeqAIJKokkos *akok)
{
  Mat_SeqAIJ *aseq;
  PetscInt    i, m, n;
  auto       &exec = PetscGetKokkosExecutionSpace();

  PetscFunctionBegin;
  PetscCheck(!A->spptr, PETSC_COMM_SELF, PETSC_ERR_PLIB, "A->spptr is supposed to be empty");

  m = akok->nrows();
  n = akok->ncols();
  PetscCall(MatSetSizes(A, m, n, m, n));
  PetscCall(MatSetType(A, MATSEQAIJKOKKOS));

  /* Set up data structures of A as a MATSEQAIJ */
  PetscCall(MatSeqAIJSetPreallocation_SeqAIJ_mine(A, MAT_SKIP_ALLOCATION, NULL));
  aseq = (Mat_SeqAIJ *)A->data;

  PetscCallCXX(akok->i_dual.sync_host(exec)); /* We always need sync'ed i, j on host */
  PetscCallCXX(akok->j_dual.sync_host(exec));
  PetscCallCXX(exec.fence());

  aseq->i       = akok->i_host_data();
  aseq->j       = akok->j_host_data();
  aseq->a       = akok->a_host_data();
  aseq->nonew   = -1; /*this indicates that inserting a new value in the matrix that generates a new nonzero is an error*/
  aseq->free_a  = PETSC_FALSE;
  aseq->free_ij = PETSC_FALSE;
  aseq->nz      = akok->nnz();
  aseq->maxnz   = aseq->nz;

  PetscCall(PetscMalloc1(m, &aseq->imax));
  PetscCall(PetscMalloc1(m, &aseq->ilen));
  for (i = 0; i < m; i++) aseq->ilen[i] = aseq->imax[i] = aseq->i[i + 1] - aseq->i[i];

  /* It is critical to set the nonzerostate, as we use it to check if sparsity pattern (hence data) has changed on host in MatAssemblyEnd */
  akok->nonzerostate = A->nonzerostate;
  A->spptr           = akok; /* Set A->spptr before MatAssembly so that A->spptr won't be allocated again there */
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

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
   PetscInt row_ao, col_ao;
   MatType mat_type;
   PetscInt max_nnzs_total, nnzs_local, nnzs_nonlocal;
   PetscInt nnzs_match_local, nnzs_match_nonlocal, colmap_match_size;

   MatGetType(*input_mat, &mat_type);
   // Are we in parallel?
   bool mpi = strcmp(mat_type, MATMPIAIJKOKKOS) == 0;

   Mat_MPIAIJ *mat_mpi;
   Mat mat_local, mat_nonlocal;
   Mat_SeqAIJKokkos *mat_localkok, *mat_nonlocalkok;
   if (mpi)
   {
      mat_mpi = (Mat_MPIAIJ *)(*input_mat)->data;
      mat_local = mat_mpi->A;
      mat_nonlocal = mat_mpi->B;
      MatGetSize(mat_nonlocal, &row_ao, &col_ao);
      mat_localkok = static_cast<Mat_SeqAIJKokkos *>(mat_local->spptr);
      mat_nonlocalkok = static_cast<Mat_SeqAIJKokkos *>(mat_nonlocal->spptr);      
   }
   else
   {
      mat_local = *input_mat;
      mat_localkok = static_cast<Mat_SeqAIJKokkos *>(mat_local->spptr);
   }

   // Get the comm
   PetscObjectGetComm((PetscObject)*input_mat, &MPI_COMM_MATRIX);
   MatGetLocalSize(*input_mat, &local_rows, &local_cols);
   MatGetSize(*input_mat, &global_rows, &global_cols);
   // This returns the global index of the local portion of the matrix
   MatGetOwnershipRange(*input_mat, &global_row_start, &global_row_end_plus_one);
   MatGetOwnershipRangeColumn(*input_mat, &global_col_start, &global_col_end_plus_one);   

   // ~~~~~~~~~~~~
   // Get pointers to the i,j,vals on the device
   // ~~~~~~~~~~~~
   const PetscInt *device_local_i, *device_local_j, *device_nonlocal_i, *device_nonlocal_j;
   PetscMemType mtype;
   PetscScalar *device_local_vals, *device_nonlocal_vals;  
   MatSeqAIJGetCSRAndMemType(mat_local, &device_local_i, &device_local_j, &device_local_vals, &mtype);  
   if (mpi) MatSeqAIJGetCSRAndMemType(mat_nonlocal, &device_nonlocal_i, &device_nonlocal_j, &device_nonlocal_vals, &mtype);          

   // ~~~~~~~~~~~~
   // Get the number of nnzs
   // ~~~~~~~~~~~~
   nnzs_local = mat_localkok->csrmat.nnz();
   nnzs_nonlocal = 0;
   if (mpi) nnzs_nonlocal = mat_nonlocalkok->csrmat.nnz();
   max_nnzs_total = nnzs_local + nnzs_nonlocal;
   nnzs_match_local = 0;
   nnzs_match_nonlocal = 0;
   colmap_match_size = 0;

   // ~~~~~~~~~~~~
   // Create the output mat
   // ~~~~~~~~~~~~
   MatCreate(MPI_COMM_MATRIX, output_mat);
   MatSetSizes(*output_mat, local_rows, local_cols, global_rows, global_cols);
   // Match the output type
   MatSetType(*output_mat, mat_type);
   MatSetUp(*output_mat);

   // Don't set any off processor entries so no need for a reduction when assembling
   MatSetOption(*output_mat, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
   MatSetOption(*output_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
   MatSetOption(*output_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);     

   // Row and column indices for our assembly
   // We know we never have more to do than the original nnzs
   PetscInt *row_indices, *col_indices;
   PetscMalloc2(max_nnzs_total, &row_indices, max_nnzs_total, &col_indices); 

      MatScalarKokkosView a_local;
      MatRowMapKokkosView i_local;   
      MatColIdxKokkosView j_local;

   // Let's calculate the row and column indices on the device
   // then copy them back to the host for MatSetPreallocationCOO
   {
      // Let's create device memory with a dual view
      auto &exec = PetscGetKokkosExecutionSpace();

      MatRowMapKokkosViewHost row_indices_h(row_indices, max_nnzs_total);
      auto row_indices_d = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, exec, row_indices_h);
      auto row_indices_dual = MatRowMapKokkosDualView(row_indices_d, row_indices_h); 
      // By default drop everything
      Kokkos::deep_copy(row_indices_d, -1);      

      MatColIdxKokkosViewHost col_indices_h(col_indices, max_nnzs_total);
      auto col_indices_d = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, exec, col_indices_h);
      auto col_indices_dual = MatColIdxKokkosDualView(col_indices_d, col_indices_h);             
      // By default drop everything
      Kokkos::deep_copy(col_indices_d, -1);       

      // We need the relative row tolerance, let's create some device memory to store it
      PetscScalarKokkosView rel_row_tol_d("rel_row_tol_d", local_rows);    
      // Copy in the tolerance
      Kokkos::deep_copy(rel_row_tol_d, tol);     
      // We need to know how many entries are in each row after our dropping  
      PetscIntKokkosView nnz_match_local_row_d("nnz_match_local_row_d", local_rows);            
      // Initialize to zero
      Kokkos::deep_copy(nnz_match_local_row_d, 0);   
      PetscIntKokkosView nnz_match_nonlocal_row_d("nnz_match_nonlocal_row_d", local_rows);            
      // Initialize to zero
      Kokkos::deep_copy(nnz_match_nonlocal_row_d, 0);         
      
      // Compute the relative row tolerances if needed
      if (relative_max_row_tolerance_int) 
      {       
         Kokkos::parallel_for(
            Kokkos::RangePolicy<>(0, local_rows), KOKKOS_LAMBDA(int i) {

               PetscInt ncols_local = device_local_i[i + 1] - device_local_i[i];
               PetscScalar max_val = -1.0;

               // Should really make these parallel reductions
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

      // These loops just set the row and col indices to not be -1
      // if we are including it in the matrix
      Kokkos::parallel_for( // for each row
         Kokkos::TeamPolicy<>(PetscGetKokkosExecutionSpace(), local_rows, Kokkos::AUTO()), KOKKOS_LAMBDA(const KokkosTeamMemberType &t) {

         PetscInt i   = t.league_rank(); // row i
         PetscInt ncols_local = device_local_i[i + 1] - device_local_i[i];

         // We also have to check if we have an existing diagonal if we are lumping
         if (lump_int)
         {
            PetscInt diag_index = -1;
            for (int j = 0; j < ncols_local; j++)
            {
               if (device_local_j[device_local_i[i] + j] + global_col_start == i + global_row_start) diag_index = j;
            }
            // This is the case where we have no existing diagonal but lumping
            // We have to add in a diagonal entry
            // Not sure this needs to be atomic as we have only one thread per row here
            if (diag_index == -1) Kokkos::atomic_increment(&nnz_match_local_row_d(i));
         }
         
         // scale entries on the row
         Kokkos::parallel_for(
            Kokkos::TeamThreadRange(t, ncols_local), [&](PetscInt j) { 

            // Set the row/col to be included (ie not -1) 
            // if it is bigger than the tolerance
            if (abs(device_local_vals[device_local_i[i] + j]) >= rel_row_tol_d(i))
            {
               // if (i == 89)
               // {
               //    std::cout << "correct " << i + global_row_start << " " << device_local_j[device_local_i[i] + j] + global_col_start << std::endl;
               // }
               row_indices_d(device_local_i[i] + j) = i + global_row_start;
               // Careful here to use global_col_start in case we are rectangular
               col_indices_d(device_local_i[i] + j) = device_local_j[device_local_i[i] + j] + global_col_start;
               Kokkos::atomic_increment(&nnz_match_local_row_d(i));
            }
            // If the entry is small and we are lumping, then add it to the diagonal
            // or if this is the diagonal and it's small but we are not dropping it 
            else if (lump_int || \
                  (!allow_drop_diagonal_int && \
                     device_local_j[device_local_i[i] + j] + global_col_start == i + global_row_start))
            {
               // if (i == 89)
               // {
               //    std::cout << "correct LUMP/DIAG " << i + global_row_start << " " << device_local_j[device_local_i[i] + j] + global_col_start << "val " << device_local_vals[device_local_i[i] + j] << std::endl;
               // }               
               row_indices_d(device_local_i[i] + j) = i + global_row_start;
               col_indices_d(device_local_i[i] + j) = i + global_row_start;
               // Essential to add in a diagonal in the case where the diagonal is present, but it 
               // is too small and would be dropped if we weren't doing lumping
               // Need this guard here as we're not adding new entries if we're doing lumping
               if (device_local_j[device_local_i[i] + j] + global_col_start == i + global_row_start) {
                  Kokkos::atomic_increment(&nnz_match_local_row_d(i));
               }
            }            
         });
      });

      // Do a reduction to get the local nnzs we end up with
      Kokkos::parallel_reduce ("ReductionLocal", local_rows, KOKKOS_LAMBDA (const PetscInt i, PetscInt& update) {
         update += nnz_match_local_row_d(i);
      }, nnzs_match_local);

      // Need to do a scan on nnz_match_local_row_d to get where each row starts
      Kokkos::parallel_scan (local_rows, KOKKOS_LAMBDA (const PetscInt i, PetscInt& update, const bool final) {
         // Inclusive scan
         update += nnz_match_local_row_d(i);         
         if (final) {
            nnz_match_local_row_d(i) = update; // only update array on final pass
         }

      });            

      if (mpi) 
      {
         // We also copy the colmap over to the device as we need it below
         PetscIntKokkosViewHost colmap_h(mat_mpi->garray, col_ao);
         auto colmap_d = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, exec, colmap_h);
         auto colmap_dual = PetscIntKokkosDualView(colmap_d, colmap_h);

         // We've modified the data on the host
         colmap_dual.modify_host();
         // Let's copy the data back to the device
         colmap_dual.sync_device(exec);  

         // These loops just set the row and col indices to not be -1
         // if we are including it in the matrix
         Kokkos::parallel_for( // for each row
            Kokkos::TeamPolicy<>(PetscGetKokkosExecutionSpace(), local_rows, Kokkos::AUTO()), KOKKOS_LAMBDA(const KokkosTeamMemberType &t) {

            PetscInt i   = t.league_rank(); // row i
            PetscInt ncols_nonlocal = device_nonlocal_i[i + 1] - device_nonlocal_i[i];
            
            // scale entries on the row
            Kokkos::parallel_for(
               Kokkos::TeamThreadRange(t, ncols_nonlocal), [&](PetscInt j) { 

               // Set the row/col to be included (ie not -1) 
               // if it is bigger than the tolerance
               if (abs(device_nonlocal_vals[device_nonlocal_i[i] + j]) >= rel_row_tol_d(i))
               {
                  row_indices_d(nnzs_local + device_nonlocal_i[i] + j) = i + global_row_start;
                  // garray is the colmap
                  col_indices_d(nnzs_local + device_nonlocal_i[i] + j) = colmap_d(device_nonlocal_j[device_nonlocal_i[i] + j]);
                  Kokkos::atomic_increment(&nnz_match_nonlocal_row_d(i));
               }
               // Don't need to test if we're the diagonal here as we're in the off-diagonal block
               else if (lump_int)
               {
                  row_indices_d(nnzs_local + device_nonlocal_i[i] + j) = i + global_row_start;
                  col_indices_d(nnzs_local + device_nonlocal_i[i] + j) = i + global_row_start;
               }            
            });
         }); 

         // Do a reduction to get the nonlocal nnzs we end up with
         Kokkos::parallel_reduce ("ReductionLocal", local_rows, KOKKOS_LAMBDA (const PetscInt i, PetscInt& update) {
            update += nnz_match_nonlocal_row_d(i);
         }, nnzs_match_nonlocal);     

         // Need to do a scan on nnz_match_nonlocal_row_d to get where each row starts
         Kokkos::parallel_scan (local_rows, KOKKOS_LAMBDA (const PetscInt i, PetscInt& update, const bool final) {
            // Inclusive scan
            update += nnz_match_nonlocal_row_d(i);         
            if (final) {
               nnz_match_nonlocal_row_d(i) = update; // only update array on final pass
            }
         });               
      }     

      // // We've modified the data on the device
      // row_indices_dual.modify_device();
      // col_indices_dual.modify_device();
      // // Let's copy the data back to the host
      // row_indices_dual.sync_host(exec);      
      // col_indices_dual.sync_host(exec);

      // Now set the sparsity on the host
      //MatSetPreallocationCOO(*output_mat, max_nnzs_total, row_indices, col_indices);  

      // ~~~~~~~~~~~~~~~~~  
      // ~~~~~~~~~~~~~~~~~

      // We need to assemble our i,j, vals - here are the local ones
      MatScalarKokkosDualView a_local_dual = MatScalarKokkosDualView("a_local", nnzs_match_local);
      MatRowMapKokkosDualView i_local_dual = MatRowMapKokkosDualView("i_local", local_rows+1);
      MatColIdxKokkosDualView j_local_dual = MatColIdxKokkosDualView("j_local", nnzs_match_local);

      a_local = a_local_dual.view_device();
      i_local = i_local_dual.view_device();
      // Initialise to zero
      Kokkos::deep_copy(i_local, 0);       
      j_local = j_local_dual.view_device();  

      // And the nonlocal ones
      MatScalarKokkosDualView a_nonlocal_dual;
      MatRowMapKokkosDualView i_nonlocal_dual;
      MatColIdxKokkosDualView j_nonlocal_dual;

      MatScalarKokkosView a_nonlocal;
      MatRowMapKokkosView i_nonlocal;
      MatColIdxKokkosView j_nonlocal;          

      if (mpi) 
      {
         a_nonlocal_dual = MatScalarKokkosDualView("a_nonlocal", nnzs_match_nonlocal);
         i_nonlocal_dual = MatRowMapKokkosDualView("i_nonlocal", local_rows+1);
         j_nonlocal_dual = MatColIdxKokkosDualView("j_nonlocal", nnzs_match_nonlocal);        

         a_nonlocal = a_nonlocal_dual.view_device();
         i_nonlocal = i_nonlocal_dual.view_device();
         // Initialise to zero
         Kokkos::deep_copy(i_nonlocal, 0);         
         j_nonlocal = j_nonlocal_dual.view_device();   
      }          

      Kokkos::parallel_for( // for each row
         Kokkos::TeamPolicy<>(PetscGetKokkosExecutionSpace(), local_rows, Kokkos::AUTO()), KOKKOS_LAMBDA(const KokkosTeamMemberType &t) {

         PetscInt i   = t.league_rank(); // row i
         PetscInt ncols_local = device_local_i[i + 1] - device_local_i[i];

         // The start of our row index comes from the scan
         i_local(i + 1) = nnz_match_local_row_d(i);

         // Really need some parallelism here! Should rewrite this!
         PetscInt counter = 0;
         PetscScalar lump_val = 0;
         for (int j = 0; j < ncols_local; j++)
         {
            // Set the row/col to be included (ie not -1) 
            // if it is bigger than the tolerance
            if (abs(device_local_vals[device_local_i[i] + j]) >= rel_row_tol_d(i))
            {
               // Have to give it the local column indices
               j_local(i_local(i) + counter) = device_local_j[device_local_i[i] + j];
               a_local(i_local(i) + counter) = device_local_vals[device_local_i[i] + j];
               // if (i_local(i) + counter == 759) {
               //    std::cout << i << " " << device_local_j[device_local_i[i] + j] << "STRAIGHT IN existing a_local " << a_local(i_local(i) + counter) << std::endl;
               // }               
               counter++;
               // Diagonal to lump if we need
               if (device_local_j[device_local_i[i] + j] + global_col_start == i + global_row_start) lump_val+= device_local_vals[device_local_i[i] + j];
            }
            // If the entry is small and we are lumping, then add it to the diagonal
            // or if this is the diagonal and it's small but we are not dropping it 
            else if (lump_int || \
                  (!allow_drop_diagonal_int && \
                     device_local_j[device_local_i[i] + j] + global_col_start == i + global_row_start))
            {
               // If we have an existing diagonal and we're either lumping or not getting rid of it
               // then we need to ensure it stays in, but it's value will be replaced after this loop 
               // with that of lump_val
               if (device_local_j[device_local_i[i] + j] + global_col_start == i + global_row_start)
               {
                  // Have to give it the local column indices
                  j_local(i_local(i) + counter) = device_local_j[device_local_i[i] + j];
                  a_local(i_local(i) + counter) = device_local_vals[device_local_i[i] + j];
                  // if (i_local(i) + counter == 759) {
                  //    std::cout << "STRAIGHT IN EXISTING DIAG NOT DROP existing a_local " << a_local(i_local(i) + counter) << std::endl;
                  // }                   
                  counter++;                  
               }
               lump_val+= device_local_vals[device_local_i[i] + j];
               // if (i == 89) {
               //    std::cout << "here is some lumping " << device_local_vals[device_local_i[i] + j] << " lump val " << lump_val << std::endl;
               // }               
            }   
         }
         // We have two cases here with lumping
         // One where we have an existing diagonal and one where we don't
         if (lump_int && ncols_local != 0)
         {

            // Find where the diagonal is - this has to happen after we've put in the existing entries
            // because the number of columns isn't necessary ncols_local, its i_local(i+1) - i_local(i)
            PetscInt diag_index = -1;
            PetscInt before_diag_index = -1;
            for (int j = 0; j < i_local(i+1) - i_local(i); j++)
            {
               // Get the index one before where the diagonal would be
               if (j_local(i_local(i) + j) + global_col_start < i + global_row_start) 
               {
                  before_diag_index = j;
               }
               // Get the diagonal index (if it exists)
               else if (j_local(i_local(i) + j) + global_col_start == i + global_row_start)
               {
                  diag_index = j;
                  //if (i == 89) std::cout << "found diag" << " ncols_local " << ncols_local << std::endl;
               } 
            }

         //if (i == 89) std::cout << "diag_index " << diag_index << " before_diag_index " << before_diag_index << std::endl;

            //if (i == 89) std::cout << "GOT INTO LUMPING PART" << std::endl;
            // If we have an existing diagonal, replace the value with the lumped value
            if (diag_index != -1)
            {
               //std::cout << i << " i_local(i) " << i_local(i) << " diag_index " << diag_index << " indx " << i_local(i) + diag_index << " nnzs_match_local " << nnzs_match_local << std::endl;
               // if (i == 89) {
               //    std::cout << "idx " << i_local(i) + diag_index << "existing a_local " << a_local(i_local(i) + diag_index) << " lump " << lump_val << std::endl;
               // }
               a_local(i_local(i) + diag_index) = lump_val;
            }
            // If we don't have an existing diagonal, add it in where it would be
            else
            {  
               if (before_diag_index == -1) exit(1);
               for (int j = ncols_local-1; j > before_diag_index; j--)
               {
                  j_local(i_local(i) + j + 1) = j_local(i_local(i) + j);
                  a_local(i_local(i) + j + 1) = a_local(i_local(i) + j);

               // if (i_local(i) + j + 1 == 759) {
               //    std::cout << "MOVED existing a_local " << a_local(i_local(i) + j + 1) << std::endl;
               // }                    
               }

               // if (i_local(i) + before_diag_index + 1 == 759) {
               //    std::cout << "NOT EXISTING DIAG existing a_local " << a_local(i_local(i) + before_diag_index + 1) << " lump " << lump_val << std::endl;
               // }               

               // Has to be the local column index
               j_local(i_local(i) + before_diag_index + 1) = i;
               a_local(i_local(i) + before_diag_index + 1) = lump_val;
               counter++;
            }
         }
      });     
      
    auto ckok = new Mat_SeqAIJKokkos(local_rows, local_cols, nnzs_match_local, i_local_dual, j_local_dual, a_local_dual);  
    MatCreate(PETSC_COMM_SELF, output_mat);
    MatSetSeqAIJKokkosWithCSRMatrix_mine(*output_mat, ckok);



   }
   PetscFree2(row_indices, col_indices);   

   // // MatSetValuesCOO all happens on the device
   // {
   //    // Have to copy the local and then the nonlocal data into one array on the device
   //    if (mpi)
   //    {         
   //       // coo_v is stored on the device
   //       PetscScalarKokkosView coo_v("coo_v", max_nnzs_total);    

   //       // Copy the local values - could also do this with views, subviews
   //       // and calling MatSeqAIJGetKokkosView, but we know our device_local_vals 
   //       // and coo_v share the same execution space
   //       Kokkos::parallel_for(
   //          Kokkos::RangePolicy<>(0, nnzs_local), KOKKOS_LAMBDA(int i) {
   //             coo_v(i) = device_local_vals[i];
   //          });    

   //       // Copy the nonlocal values
   //       Kokkos::parallel_for(
   //          Kokkos::RangePolicy<>(0, nnzs_nonlocal), KOKKOS_LAMBDA(int i) {
   //             coo_v(nnzs_local + i) = device_nonlocal_vals[i];
   //          });           

   //       // This should all happen on the gpu
   //       MatSetValuesCOO(*output_mat, coo_v.data(), INSERT_VALUES);             
   //    }
   //    // Don't have to do a copy if we have no non-local data
   //    else
   //    {
   //       // This should all happen on the gpu
   //       MatSetValuesCOO(*output_mat, device_local_vals, INSERT_VALUES);           
   //    }
   // }  

   // MatSeqAIJGetCSRAndMemType(*output_mat, &device_local_i, &device_local_j, &device_local_vals, &mtype);  

   //    for (int i = 0; i < local_rows+1; i++)
   //    {
   //       if (abs(i_local(i) != device_local_i[i]))
   //       {
   //          std::cout << i << " i_local " << i_local(i) << " device_local_i " << device_local_i[i] << std::endl;
   //          exit(1);
   //       }
   //    }      

   //    for (int i = 0; i < nnzs_match_local; i++)
   //    {
   //       if (abs(j_local(i) != device_local_j[i]))
   //       {
   //          std::cout << i << " j_local " << j_local(i) << " device_local_j " << device_local_j[i] << std::endl;
   //          exit(1);
   //       }
   //       if (abs(a_local(i) - device_local_vals[i]) > 1e-14)
   //       {
   //          std::cout << i << " a_local " << a_local(i) << " device_local_vals " << device_local_vals[i] << std::endl;
   //          exit(1);
   //       }
   //    }      

      //exit(1);    


   return;
}