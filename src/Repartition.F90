module repartition

   use petsc
   use c_petsc_interfaces

#include "petsc/finclude/petsc.h"
                
   implicit none

#include "petsc_legacy.h"

   public

   ! -------------------------------------------------------------------------------------------------------------------------------
   ! -------------------------------------------------------------------------------------------------------------------------------
   ! Some helper functions that handle repartitioning of petsc matrices
   ! -------------------------------------------------------------------------------------------------------------------------------
   ! -------------------------------------------------------------------------------------------------------------------------------      

   contains 

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine compute_mat_ratio_local_nonlocal_nnzs(input_mat, no_active_cores, ratio)
      
      ! Computes the local to non local ratios of nnzs in input_mat

      ! ~~~~~~
      type(tMat), target, intent(in)      :: input_mat
      PetscInt, intent(in)                :: no_active_cores
      real, intent(out)                   :: ratio

      ! Local
      PetscInt :: local_rows, local_cols, global_rows, global_cols, i_loc
      PetscInt :: ncols, local_nnzs, off_proc_nnzs
      integer :: errorcode
      PetscErrorCode :: ierr
      MPI_Comm :: MPI_COMM_MATRIX      
      real :: ratio_local_nnzs_off_proc, ratio_parallel
      type(tMat) :: Ad, Ao
      PetscOffset :: iicol
      PetscInt :: icol(1)

      ! ~~~~~~  

      call PetscObjectGetComm(input_mat, MPI_COMM_MATRIX, ierr)       

      ! Let's check the ratio of local to non-local nnzs on the coarse matrix
      ! If it's less than some ratio then we'll process agglomerate
      call MatMPIAIJGetSeqAIJ(input_mat, Ad, Ao, icol, iicol, ierr)
      local_nnzs = 0
      off_proc_nnzs = 0
      call MatGetLocalSize(input_mat, local_rows, local_cols, ierr)
      do i_loc = 1, local_rows         
         call MatGetRow(Ad, i_loc-1, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         local_nnzs = local_nnzs + ncols
         call MatRestoreRow(Ad, i_loc-1, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         call MatGetRow(Ao, i_loc-1, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         off_proc_nnzs = off_proc_nnzs + ncols
         call MatRestoreRow(Ao, i_loc-1, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)               
      end do      
      
      ! ~~~~~~~~~~~
      ! Get the ratio of local to non-local nnzs
      ! ~~~~~~~~~~~
      ! If a processor is entirely local, don't include it in the ratio
      ! It's a bit hard to decide here what to do, as we don't want to comm how many 
      ! processors actually have non local work
      ! So the proxy for that is no_active_cores but
      ! there is no guarantee that that many procs all have nonlocal entries
      if (off_proc_nnzs == 0) then
         ratio = 0
      else
         ratio = real(local_nnzs)/real(off_proc_nnzs)
      end if

      call MPI_Allreduce(ratio, ratio_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_MATRIX, errorcode)

      ratio = ratio_parallel  
      ! Only divide by the number of processors that are active
      ratio = ratio / real(no_active_cores)

   end subroutine      
   
! -------------------------------------------------------------------------------------------------------------------------------

   subroutine calculate_repartition(input_mat, proc_stride, no_active_cores, simple, index)
      
      ! Return an IS that represents a repartitioning of the matrix
      ! Simple being true is just processor aggregation to load balance
      ! partitions

      ! ~~~~~~
      type(tMat), target, intent(in)      :: input_mat
      PetscInt, intent(in)                :: proc_stride
      PetscInt, intent(inout)             :: no_active_cores
      logical, intent(in)                 :: simple
      type(tIS), intent(out)              :: index

      ! Local
      PetscInt :: global_rows, global_cols
      integer :: errorcode, comm_size
      PetscErrorCode :: ierr
      MPI_Comm :: MPI_COMM_MATRIX      
      type(tMat) :: adj, input_transpose
      PetscInt :: local_size_is, start_is
      integer(c_long_long) :: A_array, index_array
      PetscInt, parameter :: one=1, zero=0

      ! ~~~~~~  

      call PetscObjectGetComm(input_mat, MPI_COMM_MATRIX, ierr)  
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)      

      ! If doing simple processor aggregation
      if (simple) then

         call MatGetSize(input_mat, global_rows, global_cols, ierr)
         call GenerateIS_ProcAgglomeration_c(proc_stride, global_rows, local_size_is, start_is)

         ! Specifically on comm_world as we aren't going onto a subcomm
         ! If we're on a processor we wish to accumulate onto the local_size_is will have
         ! other enties, it will be zero on processors that will have all entries removed from
         call ISCreateStride(MPI_COMM_WORLD, local_size_is, start_is, one, index, ierr)      

      ! Else call one of the partitioners through petsc
      else

         ! Number of cores we want dofs on
         no_active_cores = floor(real(comm_size)/real(proc_stride))

         ! Have to symmetrize the input matrix or it won't work in parmetis
         ! as it expects a symmetric graph
         call MatTranspose(input_mat, MAT_INITIAL_MATRIX, input_transpose, ierr)
         call MatAXPY(input_transpose, 1.0, input_mat, DIFFERENT_NONZERO_PATTERN, ierr) 

         ! Compute the adjancency graph of the symmetrized input matrix
         call MatConvert(input_transpose, MATMPIADJ, MAT_INITIAL_MATRIX, adj, ierr)

         A_array = adj%v
         call MatPartitioning_c(A_array, no_active_cores, proc_stride, index_array)

         ! Assign the index to the IS pointer we get back from c
         index%v = index_array

         ! Destroy the adjacency matrices
         call MatDestroy(adj, ierr)
         call MatDestroy(input_transpose, ierr)

      end if

   end subroutine

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine MatMPICreateNonemptySubcomm(input_mat, on_subcomm, output_mat)
      
      ! Just calls MatMPICreateNonemptySubcomm_c

      ! ~~~~~~
      type(tMat), target, intent(in)      :: input_mat
      logical, intent(out)                :: on_subcomm
      type(tMat), target, intent(inout)   :: output_mat

      integer(c_long_long) :: A_array, B_array
      integer(c_int)       :: on_subcomm_int

      ! ~~~~~~  

      A_array = input_mat%v
      call MatMPICreateNonemptySubcomm_c(A_array, on_subcomm_int, B_array)
      if (on_subcomm_int == 1) then
         on_subcomm = .TRUE.
      else
         on_subcomm = .FALSE.
      end if
      ! Assign the index to the mat we get back from c
      output_mat%v = B_array

   end subroutine     

   !-------------------------------------------------------------------------------------------------------------------------------

end module repartition

