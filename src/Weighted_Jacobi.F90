module weighted_jacobi

   use petsc

#include "petsc/finclude/petsc.h"

   implicit none

   public
   
   PetscEnum, parameter :: PFLAREINV_WJACOBI=6
   PetscEnum, parameter :: PFLAREINV_JACOBI=7    
   
   contains

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine calculate_and_build_weighted_jacobi_inverse(matrix, weighted, inv_matrix)


      ! Builds an assembled weighted jacobi approximate inverse
      ! We could just apply the weighted diagonal without storing it in a petsc mat
      ! but its just convenient to have it like that so we don't have to change any code

      ! ~~~~~~
      type(tMat), target, intent(in)                    :: matrix
      logical, intent(in)                               :: weighted
      type(tMat), intent(inout)                         :: inv_matrix

      ! Local variables
      integer :: comm_size, errorcode
      MPI_Comm :: MPI_COMM_MATRIX 
      PetscInt :: local_rows, local_cols, ncols, global_row_start, global_row_end_plus_one
      PetscInt :: global_rows, global_cols, i_loc
      PetscErrorCode :: ierr
      type(tMat) :: temp_mat
      type(tVec) :: rhs_copy, diag_vec
      real :: norm_inf, weight
      PetscInt, parameter :: one=1, zero=0

      ! ~~~~~~    

      ! We might want to call this on a sub communicator
      ! so let's get the comm attached to the matrix and make sure to use that 
      call PetscObjectGetComm(matrix, MPI_COMM_MATRIX, ierr)    
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)    
      
      ! Get the local sizes
      call MatGetLocalSize(matrix, local_rows, local_cols, ierr)
      call MatGetSize(matrix, global_rows, global_cols, ierr)
      ! This returns the global index of the local portion of the matrix
      call MatGetOwnershipRange(matrix, global_row_start, global_row_end_plus_one, ierr)        

      call MatDuplicate(matrix, MAT_COPY_VALUES, temp_mat, ierr)          
      call MatCreateVecs(matrix, rhs_copy, diag_vec, ierr)
      call MatGetDiagonal(matrix, diag_vec, ierr)

      ! If weighting the Jacobi
      if (weighted) then

         ! D^(1/2)
         call VecSqrtAbs(diag_vec, ierr)
         ! D^(-1/2)
         call VecReciprocal(diag_vec, ierr)
         ! D^(-1/2) * A * D^(-1/2)
         call MatDiagonalScale(temp_mat, diag_vec, diag_vec, ierr)          
         ! || D^(-1/2) * A * D^(-1/2) ||_inf
         call MatNorm(temp_mat, NORM_INFINITY, norm_inf, ierr)
         call MatDestroy(temp_mat, ierr)

         call MatGetDiagonal(matrix, diag_vec, ierr)

         ! This is the weight that hypre uses, even in assymetric problems
         ! 3 / ( 4 * || D^(-1/2) * A * D^(-1/2) ||_inf )
         weight = 3.0/(4.0 * norm_inf)

      ! Unweighted
      else

         weight = 1.0

      end if

      call VecReciprocal(diag_vec, ierr)
      call VecScale(diag_vec, weight, ierr)      

      ! Let's create a matrix to represent the inverse diagonal
      ! Can't use matdiagonal as we want to do symbolic matmat products
      ! and don't want to have to define how that is done
      
      ! If not re-using
      if (inv_matrix == PETSC_NULL_MAT) then      

         if (comm_size/=1) then
            call MatCreateAIJ(MPI_COMM_MATRIX, local_rows, local_cols, &
                     global_rows, global_cols, &
                     one, PETSC_NULL_INTEGER, &
                     zero, PETSC_NULL_INTEGER, &
                     inv_matrix, ierr)   
         else
            call MatCreateSeqAIJ(MPI_COMM_MATRIX, local_rows, local_cols, &
                     one, PETSC_NULL_INTEGER, &
                     inv_matrix, ierr)            
         end if 
      end if

      ! Don't set any off processor entries so no need for a reduction when assembling
      call MatSetOption(inv_matrix, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)
      call MatSetOption(inv_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE,  ierr)          

      ! Set the diagonal
      do i_loc = global_row_start, global_row_end_plus_one-1
         call MatSetValue(inv_matrix, i_loc, i_loc, &
               1.0, INSERT_VALUES, ierr)
      end do
      call MatAssemblyBegin(inv_matrix, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(inv_matrix, MAT_FINAL_ASSEMBLY, ierr)    

      ! Set the diagonal to our weighted jacobi
      call MatDiagonalSet(inv_matrix, diag_vec, INSERT_VALUES, ierr)
   
      call VecDestroy(diag_vec, ierr)
      call VecDestroy(rhs_copy, ierr)     

   end subroutine calculate_and_build_weighted_jacobi_inverse


! -------------------------------------------------------------------------------------------------------------------------------

end module weighted_jacobi

