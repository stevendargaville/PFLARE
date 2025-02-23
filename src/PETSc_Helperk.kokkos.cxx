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

struct ReduceData {
   PetscInt count;
   bool found_diagonal;
   
   // Set count to zero and found_diagonal to false
   KOKKOS_INLINE_FUNCTION
   ReduceData() : count(0), found_diagonal(false) {}
   
   // We use this in our parallel reduction
   KOKKOS_INLINE_FUNCTION
   void operator+=(const ReduceData& src) {
      // Add all the counts
      count += src.count;
      // If we have found a diagonal entry at any point in this row
      // found_diagonal becomes true      
      found_diagonal |= src.found_diagonal;
   }

   // Required for Kokkos reduction
   KOKKOS_INLINE_FUNCTION
   static void join(volatile ReduceData& dest, const volatile ReduceData& src) {
      dest.count += src.count;
      dest.found_diagonal |= src.found_diagonal;
   }   
};

namespace Kokkos {
    template<>
    struct reduction_identity<ReduceData> {
        KOKKOS_INLINE_FUNCTION
        static ReduceData sum() {
            return ReduceData();  // Returns {count=0, found_diagonal=false}
        }
    };
}

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

// Drop according to a tolerance but with kokkos - keeping everything on the device
PETSC_INTERN void remove_small_from_sparse_kokkos(Mat *input_mat, PetscReal tol, Mat *output_mat, \
                  int relative_max_row_tolerance_int, int lump_int, int allow_drop_diagonal_int)
{
   MPI_Comm MPI_COMM_MATRIX;
   PetscInt local_rows, local_cols, global_rows, global_cols;
   PetscInt global_row_start, global_row_end_plus_one;
   PetscInt global_col_start, global_col_end_plus_one;
   PetscInt row_ao, col_ao;
   MatType mat_type;
   PetscInt nnzs_match_local, nnzs_match_nonlocal, colmap_match_size;

   MatGetType(*input_mat, &mat_type);
   // Are we in parallel?
   bool mpi = strcmp(mat_type, MATMPIAIJKOKKOS) == 0;

   Mat_MPIAIJ *mat_mpi;
   Mat mat_local, mat_nonlocal;
   if (mpi)
   {
      mat_mpi = (Mat_MPIAIJ *)(*input_mat)->data;
      mat_local = mat_mpi->A;
      mat_nonlocal = mat_mpi->B;
      MatGetSize(mat_nonlocal, &row_ao, &col_ao);    
   }
   else
   {
      mat_local = *input_mat;
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

   // ~~~~~~~~~~~~~~~~~~~~~~~
   // Let's build our i, j, and a on the device
   // ~~~~~~~~~~~~~~~~~~~~~~~
   auto &exec = PetscGetKokkosExecutionSpace();    

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
      // Reduction over all rows
      Kokkos::parallel_for(
         Kokkos::TeamPolicy<>(PetscGetKokkosExecutionSpace(), local_rows, Kokkos::AUTO()),
         KOKKOS_LAMBDA(const KokkosTeamMemberType &t) {

            PetscInt i = t.league_rank();
            PetscInt ncols_local = device_local_i[i + 1] - device_local_i[i];
            PetscScalar max_val = -1.0;

            // Reduce over local columns
            Kokkos::parallel_reduce(
               Kokkos::TeamThreadRange(t, ncols_local),
               [&](const PetscInt j, PetscScalar& thread_max) {
                  // If our current tolerance is bigger than the max value we've seen so far
                  PetscScalar val = abs(device_local_vals[device_local_i[i] + j]);
                  if (val > thread_max) thread_max = val;
               },
               Kokkos::Max<PetscScalar>(max_val)
            );

            if (mpi) {
               PetscInt ncols_nonlocal = device_nonlocal_i[i + 1] - device_nonlocal_i[i];
               
               // Reduce over nonlocal columns
               Kokkos::parallel_reduce(
                  Kokkos::TeamThreadRange(t, ncols_nonlocal),
                  [&](const PetscInt j, PetscScalar& thread_max) {
                     // If our current tolerance is bigger than the max value we've seen so far
                     PetscScalar val = abs(device_nonlocal_vals[device_nonlocal_i[i] + j]);
                     if (val > thread_max) thread_max = val;
                  },
                  Kokkos::Max<PetscScalar>(max_val)
               );
            }

            // Only want one thread in the team to write the result
            Kokkos::single(Kokkos::PerTeam(t), [&]() {
               rel_row_tol_d(i) *= max_val;
            });
      });
   }

   // ~~~~~~~~~~~~
   // Need to count the number of nnzs we end up with, on each row and in total
   // ~~~~~~~~~~~~
   // Reduce over all the rows
   Kokkos::parallel_reduce(
      Kokkos::TeamPolicy<>(PetscGetKokkosExecutionSpace(), local_rows, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const KokkosTeamMemberType &t, PetscInt& thread_total) {

      PetscInt i   = t.league_rank(); // row i
      PetscInt ncols_local = device_local_i[i + 1] - device_local_i[i];
      // We have a custom reduction type defined - ReduceData
      // Which has both a nnz count for this row, but also tracks whether we 
      // found the diagonal
      ReduceData row_result;

      // Reduce over all the columns
      Kokkos::parallel_reduce(
         Kokkos::TeamThreadRange(t, ncols_local),
         [&](const PetscInt j, ReduceData& thread_data) {

            // Is this column the diagonal
            bool is_diagonal = (device_local_j[device_local_i[i] + j] + global_col_start == i + global_row_start);
            
            // We have found a diagonal in this row
            if (is_diagonal) {
               thread_data.found_diagonal = true;
            }
            
            // If the value is bigger than the tolerance, we keep it
            if (abs(device_local_vals[device_local_i[i] + j]) >= rel_row_tol_d(i)) {
               thread_data.count++;
            }
            // Or if it's small but its the diagonal and we're keeping diagonals
            else if (!allow_drop_diagonal_int && is_diagonal) {
               thread_data.count++;
            }
         }, row_result
      );

      // We're finished our parallel reduction for this row
      // If we're lumping but there was no diagonal in this row
      // we'll have to add in a diagonal
      // This will add one for every thread in this team, but all 
      // the threads in this team share the same result after the reduction
      if (lump_int && !row_result.found_diagonal) row_result.count++;
      // Only want one thread in the team to write the result
      Kokkos::single(Kokkos::PerTeam(t), [&]() {      
         nnz_match_local_row_d(i) = row_result.count;
      });
      // thread_total is automatically only updated once per team
      // because it is managed by the kokkos reduction      
      thread_total += row_result.count;
      },
      nnzs_match_local
   );

   // ~~~~~~~~~~~~

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

      // ~~~~~~~~~~~~
      // Need to count the number of nnzs we end up with, on each row and in total
      // ~~~~~~~~~~~~
      // Reduce over all the rows
      Kokkos::parallel_reduce(
         Kokkos::TeamPolicy<>(PetscGetKokkosExecutionSpace(), local_rows, Kokkos::AUTO()),
         KOKKOS_LAMBDA(const KokkosTeamMemberType &t, PetscInt& thread_total) {
            
            PetscInt i = t.league_rank();
            PetscInt ncols_nonlocal = device_nonlocal_i[i + 1] - device_nonlocal_i[i];
            ReduceData row_result;

            // Reduce over all the columns
            Kokkos::parallel_reduce(
               Kokkos::TeamThreadRange(t, ncols_nonlocal),
               [&](const PetscInt j, ReduceData& thread_data) {
                  // If the value is bigger than the tolerance, we keep it
                  if (abs(device_nonlocal_vals[device_nonlocal_i[i] + j]) >= rel_row_tol_d(i)) {
                     thread_data.count++;
                  }
               },
               row_result
            );
            
            // Store row count and update total
            // Only want one thread in the team to write the result
            // (all threads in the team share the same result)
            Kokkos::single(Kokkos::PerTeam(t), [&]() {             
               nnz_match_nonlocal_row_d(i) = row_result.count;
            });
            thread_total += row_result.count;
         },
         nnzs_match_nonlocal
      );

      // ~~~~~~~~~~~~

      // Need to do a scan on nnz_match_nonlocal_row_d to get where each row starts
      Kokkos::parallel_scan (local_rows, KOKKOS_LAMBDA (const PetscInt i, PetscInt& update, const bool final) {
         // Inclusive scan
         update += nnz_match_nonlocal_row_d(i);         
         if (final) {
            nnz_match_nonlocal_row_d(i) = update; // only update array on final pass
         }
      });               
   }       

   // ~~~~~~~~~~~~~~~~~  
   // ~~~~~~~~~~~~~~~~~

   // We need to assemble our i,j, vals - here are the local ones
   MatScalarKokkosDualView a_local_dual = MatScalarKokkosDualView("a_local", nnzs_match_local);
   MatRowMapKokkosDualView i_local_dual = MatRowMapKokkosDualView("i_local", local_rows+1);
   MatColIdxKokkosDualView j_local_dual = MatColIdxKokkosDualView("j_local", nnzs_match_local);

   MatScalarKokkosView a_local = a_local_dual.view_device();
   MatRowMapKokkosView i_local = i_local_dual.view_device();
   // Initialise to zero
   Kokkos::deep_copy(i_local, 0);       
   MatColIdxKokkosView j_local = j_local_dual.view_device();  

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

   // Have to use rangepolicy here rather than teampolicy, as we don't want
   // threads in each team to get their own copies of things like counter
   // Desperately need to rewrite this whole loop to exploit more parallelism!
   Kokkos::parallel_for(
      Kokkos::RangePolicy<>(0, local_rows), KOKKOS_LAMBDA(int i) {
         
      PetscInt ncols_local = device_local_i[i + 1] - device_local_i[i];
      // The start of our row index comes from the scan
      i_local(i + 1) = nnz_match_local_row_d(i);
      PetscInt pos = (i == 0) ? 0 : nnz_match_local_row_d(i-1);

      // Keep track of the lumping
      PetscScalar lump_val = 0;
      
      // Single sequential pass through columns
      for (int j = 0; j < ncols_local; j++) {
         bool keep_col = false;
         bool is_diagonal = (device_local_j[device_local_i[i] + j] + global_col_start == i + global_row_start);

         // If we hit a diagonal put it in the lump'd value
         if (is_diagonal && lump_int) lump_val+= device_local_vals[device_local_i[i] + j];            
         
         // Check if we keep this column because of size
         if (abs(device_local_vals[device_local_i[i] + j]) >= rel_row_tol_d(i)) {
            keep_col = true;
         }
         // Or if we keep it because we're not dropping diagonals
         else if (!allow_drop_diagonal_int && is_diagonal) {
            keep_col = true;
         }
         
         // Add in this column if we have to keep it
         if (keep_col) {
            j_local(pos) = device_local_j[device_local_i[i] + j];
            a_local(pos) = device_local_vals[device_local_i[i] + j];
            pos++;
         }
         // If we're not on the diagonal and we're small enough to lump
         // add in
         else if (lump_int && !is_diagonal) {
            lump_val+= device_local_vals[device_local_i[i] + j];
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
            } 
         }
         // If we have an existing diagonal, replace the value with the lumped value
         if (diag_index != -1)
         {
            a_local(i_local(i) + diag_index) = lump_val;
         }
         // If we don't have an existing diagonal, add it in where it would be
         else
         {  
            for (int j = ncols_local-1; j > before_diag_index; j--)
            {
               j_local(i_local(i) + j + 1) = j_local(i_local(i) + j);
               a_local(i_local(i) + j + 1) = a_local(i_local(i) + j);                 
            }            

            // Has to be the local column index
            j_local(i_local(i) + before_diag_index + 1) = i;
            a_local(i_local(i) + before_diag_index + 1) = lump_val;
         }
      }      
   });

   a_local_dual.modify_device();
   i_local_dual.modify_device();
   j_local_dual.modify_device();

   // We can create our matrix directly on the device
   // See MatSeqAIJKokkosMergeMats for example
   auto akok = new Mat_SeqAIJKokkos(local_rows, local_cols, nnzs_match_local, i_local_dual, j_local_dual, a_local_dual);   

   MatCreate(PETSC_COMM_SELF, output_mat);
   // Why isn't this publically available??
   MatSetSeqAIJKokkosWithCSRMatrix_mine(*output_mat, akok);

   return;
}