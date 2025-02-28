#include <petscvec_kokkos.hpp>
#include <petsc.h>
#include <iostream>
#include <../src/mat/impls/aij/seq/aij.h>
#include <../src/mat/impls/aij/mpi/mpiaij.h>
#include <../src/mat/impls/aij/seq/kokkos/aijkok.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <../src/vec/vec/impls/seq/kokkos/veckokkosimpl.hpp>

using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;
using DefaultMemorySpace    = Kokkos::DefaultExecutionSpace::memory_space;
using intKokkosViewHost = Kokkos::View<int *, Kokkos::HostSpace>;
using intKokkosView = Kokkos::View<int *, Kokkos::DefaultExecutionSpace>;
using boolKokkosView = Kokkos::View<bool *, Kokkos::DefaultExecutionSpace>;

// Define C and F points in the CF marker 
int C_POINT = 1;
int F_POINT = -1;

//------------------------------------------------------------------------------------------------------------------------

// PMISR cf splitting but on the device
PETSC_INTERN void pmisr_kokkos(Mat *strength_mat, int max_luby_steps, int pmis_int, PetscReal *measure_local, int *cf_markers_local, int zero_measure_c_point_int)
{

   MPI_Comm MPI_COMM_MATRIX;
   PetscInt local_rows, local_cols, global_rows, global_cols;
   PetscInt global_row_start, global_row_end_plus_one;
   PetscInt rows_ao, cols_ao;
   MatType mat_type;

   MatGetType(*strength_mat, &mat_type);
   // Are we in parallel?
   bool mpi = strcmp(mat_type, MATMPIAIJKOKKOS) == 0;

   Mat_MPIAIJ *mat_mpi = nullptr;
   Mat mat_local, mat_nonlocal;

   if (mpi)
   {
      mat_mpi = (Mat_MPIAIJ *)(*strength_mat)->data;
      mat_local = mat_mpi->A;
      mat_nonlocal = mat_mpi->B;
      MatGetSize(mat_nonlocal, &rows_ao, &cols_ao); 
   }
   else
   {
      mat_local = *strength_mat;
   }

   // Get the comm
   PetscObjectGetComm((PetscObject)*strength_mat, &MPI_COMM_MATRIX);
   MatGetLocalSize(*strength_mat, &local_rows, &local_cols);
   MatGetSize(*strength_mat, &global_rows, &global_cols);
   // This returns the global index of the local portion of the matrix
   MatGetOwnershipRange(*strength_mat, &global_row_start, &global_row_end_plus_one);

   // ~~~~~~~~~~~~
   // Get pointers to the i,j,vals on the device
   // ~~~~~~~~~~~~
   const PetscInt *device_local_i = nullptr, *device_local_j = nullptr, *device_nonlocal_i = nullptr, *device_nonlocal_j = nullptr;
   PetscMemType mtype;
   PetscScalar *device_local_vals = nullptr, *device_nonlocal_vals = nullptr;  
   MatSeqAIJGetCSRAndMemType(mat_local, &device_local_i, &device_local_j, &device_local_vals, &mtype);  
   if (mpi) MatSeqAIJGetCSRAndMemType(mat_nonlocal, &device_nonlocal_i, &device_nonlocal_j, &device_nonlocal_vals, &mtype);          

   // Host and device memory for the cf_markers - be careful these aren't petsc ints
   intKokkosViewHost cf_markers_local_h(cf_markers_local, local_rows);
   intKokkosView cf_markers_local_d("cf_markers_local_d", local_rows);

   // Host and device memory for the measure
   PetscScalarKokkosViewHost measure_local_h(measure_local, local_rows);
   PetscScalarKokkosView measure_local_d("measure_local_d", local_rows);   
   PetscScalarKokkosView measure_nonlocal_d;
   if (mpi) measure_nonlocal_d = PetscScalarKokkosView("measure_nonlocal_d", cols_ao);   

   // Device memory for the mark
   boolKokkosView mark_d("mark_d", local_rows);   

   // Copy the input measure from host to device
   Kokkos::deep_copy(measure_local_d, measure_local_h);  

   // Compute the measure
   Kokkos::parallel_for(
      Kokkos::RangePolicy<>(0, local_rows), KOKKOS_LAMBDA(int i) {
         
      PetscInt ncols_local = device_local_i[i + 1] - device_local_i[i];
      measure_local_d(i) += ncols_local;

      if (mpi)
      {
         PetscInt ncols_nonlocal = device_nonlocal_i[i + 1] - device_nonlocal_i[i];
         measure_local_d(i) += ncols_nonlocal;
      }
      // Flip the sign if pmis
      if (pmis_int == 1) measure_local_d(i) *= -1;
   });      

   // We use a vec to do our scatter
   Vec cf_markers_vec;
   PetscScalarKokkosView cf_markers_local_real_d;
   ConstPetscScalarKokkosView cf_markers_local_real_const_d;
   // Create a kokkos vector
   MatCreateVecs(*strength_mat, &cf_markers_vec, NULL);

   // Get a device pointer for cf_markers_vec
   VecGetKokkosViewWrite(cf_markers_vec, &cf_markers_local_real_d);
   // Copy the measure over
   Kokkos::deep_copy(cf_markers_local_real_d, measure_local_d);  
   VecRestoreKokkosViewWrite(cf_markers_vec, &cf_markers_local_real_d);

   if (mpi)
   {
      // Scatter the measure
      VecScatterBegin(mat_mpi->Mvctx, cf_markers_vec, mat_mpi->lvec, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(mat_mpi->Mvctx, cf_markers_vec, mat_mpi->lvec, INSERT_VALUES, SCATTER_FORWARD);

      // Get a device pointer for cf_markers_vec non local components
      ConstPetscScalarKokkosView cf_markers_nonlocal_real_d;
      VecGetKokkosView(mat_mpi->lvec, &cf_markers_nonlocal_real_d);
      // Copy the non local measure
      Kokkos::deep_copy(measure_nonlocal_d, cf_markers_nonlocal_real_d);    
      VecRestoreKokkosView(mat_mpi->lvec, &cf_markers_nonlocal_real_d);   
   }

   // Initialise the set
   int counter_in_set_start = 0;
   // Count how many in the set to begin with
   Kokkos::parallel_reduce ("Reduction", local_rows, KOKKOS_LAMBDA (const int i, int& update) {
      if (abs(measure_local_d(i) < 1)) update++;
   }, counter_in_set_start);

   VecGetKokkosViewWrite(cf_markers_vec, &cf_markers_local_real_d);

   Kokkos::parallel_for(
      Kokkos::RangePolicy<>(0, local_rows), KOKKOS_LAMBDA(int i) {
         
      if (abs(measure_local_d(i) < 1))
      {
         if (zero_measure_c_point_int == 1) {
            if (pmis_int == 1) {
               // Set as F here but reversed below to become C
               cf_markers_local_real_d(i) = F_POINT;
            }
            else {
               // Becomes C
               cf_markers_local_real_d(i) = C_POINT;
            }  
         }
         else {
            if (pmis_int == 1) {
               // Set as C here but reversed below to become F
               // Otherwise dirichlet conditions persist down onto the coarsest grid
               cf_markers_local_real_d(i) = C_POINT;
            }
            else {
               // Becomes F
               cf_markers_local_real_d(i) = F_POINT;
            }
         }
      }
      else
      {
         cf_markers_local_real_d(i) = 0;
      }
   });  
   VecRestoreKokkosViewWrite(cf_markers_vec, &cf_markers_local_real_d);

   // Check the total number of undecided in parallel
   int counter_undecided, counter_parallel;
   if (max_luby_steps < 0) {
      // Assuming here we don't have more than 2B local rows
      counter_undecided = int(local_rows) - counter_in_set_start;
      // Parallel reduction!
      MPI_Allreduce(&counter_undecided, &counter_parallel, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_MATRIX);
      counter_undecided = counter_parallel;
      
   // If we're doing a fixed number of steps, then we don't care
   // how many undecided nodes we have - have to take care here not to use
   // local_rows for counter_undecided, as we may have zero DOFs on some procs
   // but we have to enter the loop below for the collective scatters 
   }
   else {
      counter_undecided = 1;
   }   

   // ~~~~~~~~~~~~
   // Now go through the outer Luby loop
   // ~~~~~~~~~~~~      

   // Let's keep track of how many times we go through the loops
   int loops_through = -1;
   do 
   {
      // If max_luby_steps is positive, then we only take that many times through this top loop
      // We typically find 2-3 iterations decides >99% of the nodes 
      // and a fixed number of outer loops means we don't have to do any parallel reductions
      // We will do redundant nearest neighbour comms in the case we have already 
      // finished deciding all the nodes, but who cares
      // Any undecided nodes just get turned into C points
      // We can do this as we know we won't ruin Aff by doing so, unlike in a normal multigrid
      if (max_luby_steps > 0 && max_luby_steps+1 == -loops_through) break;

      // ~~~~~~~~~
      // Start the async scatter of the nonlocal cf_markers
      // ~~~~~~~~~
      if (mpi) VecScatterBegin(mat_mpi->Mvctx, cf_markers_vec, mat_mpi->lvec, INSERT_VALUES, SCATTER_FORWARD);

      // This keeps track of which of the candidate nodes can become in the set
      // Only need this because we want to do async comms so we need a way to trigger
      // a node not being in the set due to either strong local neighbours *or* strong offproc neighbours
      Kokkos::deep_copy(mark_d, true);   

      VecGetKokkosView(cf_markers_vec, &cf_markers_local_real_const_d);

      // Any that aren't zero cf marker are already assigned so set to to false
      Kokkos::parallel_for(
         Kokkos::RangePolicy<>(0, local_rows), KOKKOS_LAMBDA(int i) {

            if (cf_markers_local_real_const_d(i) != 0) mark_d(i) = false;
      });

      // ~~~~~~~~
      // Go and do the local component
      // ~~~~~~~~      
      Kokkos::parallel_for(
         Kokkos::TeamPolicy<>(PetscGetKokkosExecutionSpace(), local_rows, Kokkos::AUTO()),
         KOKKOS_LAMBDA(const KokkosTeamMemberType &t) {

            // Row
            PetscInt i = t.league_rank();
            int strong_neighbours = 0;

            // Check this row isn't already marked
            if (cf_markers_local_real_const_d(i) == 0)
            {
               PetscInt i = t.league_rank();
               PetscInt ncols_local = device_local_i[i + 1] - device_local_i[i];

               // Reduce over local columns to get the number of strong neighbours
               Kokkos::parallel_reduce(
                  Kokkos::TeamThreadRange(t, ncols_local),
                  [&](const PetscInt j, int& strong_count) {     

                  // Have to only check active strong neighbours
                  if (measure_local_d(i) >= measure_local_d(device_local_j[device_local_i[i] + j]) && \
                        cf_markers_local_real_const_d(device_local_j[device_local_i[i] + j]) == 0)
                  {
                     strong_count++;
                  }
               
               }, strong_neighbours
               );     

               // Only want one thread in the team to write the result
               Kokkos::single(Kokkos::PerTeam(t), [&]() {                  
                  // If we have any strong neighbours
                  if (strong_neighbours > 0) mark_d(i) = false;     
               });
            }
      });

      // ~~~~~~~~
      // Now go through and do the non-local part of the matrix
      // ~~~~~~~~           
      if (mpi) {

         ConstPetscScalarKokkosView cf_markers_nonlocal_real_d;
         VecGetKokkosView(mat_mpi->lvec, &cf_markers_nonlocal_real_d);         

         // Finish the async scatter
         VecScatterEnd(mat_mpi->Mvctx, cf_markers_vec, mat_mpi->lvec, INSERT_VALUES, SCATTER_FORWARD);

         Kokkos::parallel_for(
            Kokkos::TeamPolicy<>(PetscGetKokkosExecutionSpace(), local_rows, Kokkos::AUTO()),
            KOKKOS_LAMBDA(const KokkosTeamMemberType &t) {

               // Row
               PetscInt i = t.league_rank();
               int strong_neighbours = 0;

               // Check this row isn't already marked
               if (cf_markers_local_real_const_d(i) == 0)
               {
                  PetscInt i = t.league_rank();
                  PetscInt ncols_nonlocal = device_nonlocal_i[i + 1] - device_nonlocal_i[i];

                  // Reduce over nonlocal columns to get the number of strong neighbours
                  Kokkos::parallel_reduce(
                     Kokkos::TeamThreadRange(t, ncols_nonlocal),
                     [&](const PetscInt j, int& strong_count) {     

                     if (measure_nonlocal_d(i) >= measure_nonlocal_d(device_nonlocal_j[device_nonlocal_i[i] + j])  && \
                              cf_markers_nonlocal_real_d(device_nonlocal_j[device_nonlocal_i[i] + j]) == 0)
                     {
                        strong_count++;
                     }
                  
                  }, strong_neighbours
                  );     

                  // Only want one thread in the team to write the result
                  Kokkos::single(Kokkos::PerTeam(t), [&]() {                  
                     // If we have any strong neighbours
                     if (strong_neighbours > 0) mark_d(i) = false;     
                  });
               }
         }); 

         VecRestoreKokkosView(mat_mpi->lvec, &cf_markers_nonlocal_real_d);            
      }

      VecRestoreKokkosView(cf_markers_vec, &cf_markers_local_real_const_d);
      VecGetKokkosViewWrite(cf_markers_vec, &cf_markers_local_real_d);

      // Any that aren't zero cf marker are already assigned so set to to false
      Kokkos::parallel_for(
         Kokkos::RangePolicy<>(0, local_rows), KOKKOS_LAMBDA(int i) {

            if (mark_d(i)) cf_markers_local_real_d(i) = double(loops_through);
      });      
      VecRestoreKokkosViewWrite(cf_markers_vec, &cf_markers_local_real_d);
      VecGetKokkosView(cf_markers_vec, &cf_markers_local_real_const_d);

      if (mpi) 
      {

         // We need a non const view now
         PetscScalarKokkosView cf_markers_nonlocal_real_d;
         VecGetKokkosViewWrite(mat_mpi->lvec, &cf_markers_nonlocal_real_d); 

         // We're going to do an add reverse scatter, so set them to zero
         Kokkos::deep_copy(cf_markers_nonlocal_real_d, 0.0);  

         Kokkos::parallel_for(
            Kokkos::TeamPolicy<>(PetscGetKokkosExecutionSpace(), local_rows, Kokkos::AUTO()),
            KOKKOS_LAMBDA(const KokkosTeamMemberType &t) {

               // Row
               PetscInt i = t.league_rank();

               // Check if this node has been assigned during this top loop
               if (cf_markers_local_real_const_d(i) == loops_through)
               {
                  PetscInt i = t.league_rank();
                  PetscInt ncols_nonlocal = device_nonlocal_i[i + 1] - device_nonlocal_i[i];

                  // For over nonlocal columns
                  Kokkos::parallel_for(
                     Kokkos::TeamThreadRange(t, ncols_nonlocal), [&](const PetscInt j) {

                        // Needs to be atomic as may being set by many threads
                        Kokkos::atomic_store(&cf_markers_nonlocal_real_d(device_nonlocal_j[device_nonlocal_i[i] + j]), double(C_POINT));     
                  });     
               }
         }); 

         VecRestoreKokkosViewWrite(mat_mpi->lvec, &cf_markers_nonlocal_real_d); 

         // We've updated the values in cf_markers_nonlocal, which is a pointer to lvec
         // Calling a reverse scatter add will then update the values of cf_markers_vec
         // Begin the scatter asynchronously
         VecScatterBegin(mat_mpi->Mvctx, cf_markers_vec, mat_mpi->lvec, ADD_VALUES, SCATTER_REVERSE);
      }

      VecGetKokkosViewWrite(cf_markers_vec, &cf_markers_local_real_d);
      // Go and do local
      Kokkos::parallel_for(
         Kokkos::TeamPolicy<>(PetscGetKokkosExecutionSpace(), local_rows, Kokkos::AUTO()),
         KOKKOS_LAMBDA(const KokkosTeamMemberType &t) {

            // Row
            PetscInt i = t.league_rank();

            // Check if this node has been assigned during this top loop
            if (cf_markers_local_real_const_d(i) == loops_through)
            {
               PetscInt i = t.league_rank();
               PetscInt ncols_local = device_local_i[i + 1] - device_local_i[i];

               // For over nonlocal columns
               Kokkos::parallel_for(
                  Kokkos::TeamThreadRange(t, ncols_local), [&](const PetscInt j) {

                     // Needs to be atomic as may being set by many threads
                     Kokkos::atomic_store(&cf_markers_local_real_d(device_local_j[device_local_i[i] + j]), double(C_POINT));     
               });     
            }
      });   
      VecRestoreKokkosViewWrite(cf_markers_vec, &cf_markers_local_real_d);    
      VecRestoreKokkosView(cf_markers_vec, &cf_markers_local_real_const_d);

      if (mpi) 
      {
         // Finish the scatter
         VecScatterEnd(mat_mpi->Mvctx, cf_markers_vec, mat_mpi->lvec, ADD_VALUES, SCATTER_REVERSE);
      }

      // We've done another top level loop
      loops_through = loops_through - 1;

      // ~~~~~~~~~~~~
      // Check the total number of undecided in parallel before we loop again
      // ~~~~~~~~~~~~
      if (max_luby_steps < 0) {

         VecGetKokkosView(cf_markers_vec, &cf_markers_local_real_const_d);

         counter_undecided = 0;  
         Kokkos::parallel_reduce ("ReductionCounter_undecided", local_rows, KOKKOS_LAMBDA (const int i, int& update) {
            if (cf_markers_local_real_const_d(i) == 0) update++;
         }, counter_undecided);  

         VecRestoreKokkosView(cf_markers_vec, &cf_markers_local_real_const_d);       

         // Parallel reduction!
         MPI_Allreduce(&counter_undecided, &counter_parallel, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_MATRIX);
         counter_undecided = counter_parallel;            
      }

   }
   while (counter_undecided != 0);

   // ~~~~~~~~~
   // Now assign our final cf markers
   // ~~~~~~~~~
   VecGetKokkosView(cf_markers_vec, &cf_markers_local_real_const_d);

   Kokkos::parallel_for(
      Kokkos::RangePolicy<>(0, local_rows), KOKKOS_LAMBDA(int i) {
         
      if (cf_markers_local_real_const_d(i) == 0)
      {
         cf_markers_local_d(i) = C_POINT;
      }
      else if (cf_markers_local_real_const_d(i) < 0)
      {
         cf_markers_local_d(i) = F_POINT;
      }
      else
      {
         cf_markers_local_d(i) = C_POINT;
      }
      if (pmis_int) cf_markers_local_d(i) *= -1;
   });  

   VecRestoreKokkosView(cf_markers_vec, &cf_markers_local_real_const_d);   

   // Now copy device cf_markers_local_d back to host
   Kokkos::deep_copy(cf_markers_local_h, cf_markers_local_d);
   // Can now destroy our vec
   VecDestroy(&cf_markers_vec);

   return;
}

//------------------------------------------------------------------------------------------------------------------------