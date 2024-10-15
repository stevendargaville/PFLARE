module approx_inverse_setup

   use gmres_poly
   use gmres_poly_newton
   use neumann_poly
   use weighted_jacobi
   use sai_z
   use repartition
   use matshell

#include "petsc/finclude/petsc.h"
      
   implicit none

   public

   contains

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine calculate_and_build_approximate_inverse(matrix, inverse_type, &
                  poly_order, poly_sparsity_order, &
                  matrix_free, subcomm, &
                  inv_matrix)

      ! Builds an approximate inverse
      ! inverse_type:
      ! PFLAREINV_POWER - GMRES polynomial with the power basis 
      ! PFLAREINV_ARNOLDI - GMRES polynomial with the arnoldi basis 
      ! PFLAREINV_NEWTON - GMRES polynomial with the newton basis - can only be used matrix-free atm      
      ! PFLAREINV_NEUMANN - Neumann polynomial
      ! PFLAREINV_SAI - SAI
      ! PFLAREINV_ISAI - Incomplete SAI (ie a restricted additive schwartz)
      ! PFLAREINV_WJACOBI - Weighted Jacobi with weight 3 / ( 4 * || Dff^(-1/2) * Aff * Dff^(-1/2) ||_inf )
      ! PFLAREINV_JACOBI - Unweighted Jacobi                  
      ! This is just a wrapper around the start and finish routines below
      ! No opportunity to put work between the async comms
      ! If you want to do that, you will have to call the individual routines 

      ! ~~~~~~
      type(tMat), intent(in)                            :: matrix
      integer, intent(in)                               :: inverse_type, poly_order
      integer, intent(in)                               :: poly_sparsity_order
      logical, intent(in)                               :: matrix_free, subcomm
      type(tMat), intent(inout)                         :: inv_matrix

      type(tsqr_buffers)                           :: buffers 
      ! This pointer must be declared contiguous, as we require coefficients
      ! to be contiguous in build_gmres_polynomial_newton_inverse as 
      ! we muck about with the rank and if this pointer is not declared
      ! contiguous it creates a temporary copy which then disappears and
      ! we segfault when trying to apply mf
      real, dimension(:, :), contiguous, pointer   :: coefficients
      real, dimension(poly_order + 1, 1), target   :: coefficients_stack
      type(mat_ctxtype), pointer :: mat_ctx => null()
      PetscErrorCode :: ierr
      type(tMat) :: reuse_mat

      ! ~~~~~~  
      
      ! Have to allocate heap memory if matrix-free as the context in inv_matrix
      ! just points at the coefficients
      if (matrix_free) then
         if (inverse_type == PFLAREINV_NEWTON) then
            ! Newton basis needs storage for real and imaginary roots
            allocate(coefficients(poly_order + 1, 2))
         else 
            allocate(coefficients(poly_order + 1, 1))
         end if
      else
         coefficients => coefficients_stack
      end if

      ! The inverse is always returned on the same comm as the input matrix
      ! but some methods benefit from having their intermediate calculations 
      ! done on a subcomm
      buffers%subcomm = subcomm
      ! Don't do any re-use, if you want reuse call the start/finish yourself
      reuse_mat = PETSC_NULL_MAT

      ! Start the calculation
      call start_approximate_inverse(matrix, inverse_type, &
                  poly_order, &
                  buffers, coefficients)
      ! Finish it
      call finish_approximate_inverse(matrix, inverse_type, &
                  poly_order, poly_sparsity_order, &
                  buffers, coefficients, &
                  matrix_free, &
                  reuse_mat, &
                  inv_matrix)
            
      ! Have to tell the matshell to own the coefficient pointers
      if (matrix_free) then
         call MatShellGetContext(inv_matrix, mat_ctx, ierr)
         mat_ctx%own_coefficients = .TRUE.               
      end if
      if (reuse_mat /= PETSC_NULL_MAT) then
         call MatDestroy(reuse_mat, ierr)
      end if

   end subroutine calculate_and_build_approximate_inverse    

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine start_approximate_inverse(matrix, inverse_type, poly_order, &
                  buffers, coefficients)

      ! Starts the assembly of an approximate inverse
      ! We have different types of inverses we can use
      ! Have to call finish_approximate_inverse before it can be used

      ! ~~~~~~
      type(tMat), intent(in)                            :: matrix
      integer, intent(in)                               :: inverse_type, poly_order
      type(tsqr_buffers), intent(inout)                 :: buffers 
      real, dimension(:, :), intent(inout)              :: coefficients

      ! Local variables
      type(tMat) :: matrix_subcomm
      ! ~~~~~~    

      if (buffers%subcomm .AND. &
            (inverse_type == PFLAREINV_POWER .OR. inverse_type == PFLAREINV_ARNOLDI .OR. inverse_type == PFLAREINV_NEWTON)) then

         ! Create a version of matrix on a subcomm and store it in buffers%matrix
         ! If all ranks have entries, then this just returns a pointer to the original 
         ! matrix (ie buffers%matrix = matrix)
         ! If not, then a copy is taken on a subcomm with ranks that have entries
         ! That is the matrix we use to start our polynomial coefficient calculation
         ! But *importantly* the inverse matrix we create is always on MPI_COMM_WORLD
         ! I did try building a mg hierarchy that progressively moved onto smaller 
         ! subcomms with processor agglomeration, but it ended up expensive
         ! This way, only the reductions done in the coefficient calculation happen on 
         ! the subcomm
         ! On any rank that isn't part of this subcomm, nothing happens and the coefficients
         ! are not comm'd onto those ranks, which is fine as those ranks have no rows anyway
         ! and hence they don't need them         
         call MatMPICreateNonemptySubcomm(matrix, buffers%matrix)        
      else
         buffers%matrix = matrix
      end if

      ! ~~~~~~~~~~~~~~~
      ! ~~~~~~~~~~~~~~~

      ! If we're on the subcomm, we don't need to do anything
      if (buffers%matrix /= PETSC_NULL_MAT) then

         ! Gmres poylnomial with power basis
         if (inverse_type == PFLAREINV_POWER) then

            ! Only does one reduction in parallel
            ! Want to start the coefficient calculation asap, as we have a non-blocking
            ! all reduce in it
            call start_gmres_polynomial_coefficients_power(buffers%matrix, poly_order, &
                  buffers, coefficients(:, 1))         

         ! Gmres polynomial with arnoldi basis
         else if (inverse_type == PFLAREINV_ARNOLDI) then

            ! Does lots of reductions in parallel
            call calculate_gmres_polynomial_coefficients_arnoldi(buffers%matrix, poly_order, coefficients(:, 1))

         ! Gmres polynomial with Newton basis - Can only use matrix-free
         else if (inverse_type == PFLAREINV_NEWTON) then

            ! Does lots of reductions in parallel
            call calculate_gmres_polynomial_roots_newton(buffers%matrix, poly_order, coefficients)
         end if
      end if

   end subroutine start_approximate_inverse

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine finish_approximate_inverse(matrix, inverse_type, &
                  poly_order, poly_sparsity_order, &
                  buffers, coefficients, &
                  matrix_free, reuse_mat, inv_matrix)

      ! Finish the assembly of an approximate inverse

      ! ~~~~~~
      type(tMat), intent(in)                            :: matrix
      integer, intent(in)                               :: inverse_type, poly_order
      integer, intent(in)                               :: poly_sparsity_order
      type(tsqr_buffers), intent(inout)                 :: buffers      
      real, dimension(:, :), contiguous, intent(inout)  :: coefficients
      logical, intent(in)                               :: matrix_free
      type(tMat), intent(inout)                         :: reuse_mat, inv_matrix

      logical :: incomplete
      PetscErrorCode :: ierr
      integer :: errorcode

      ! ~~~~~~    

      ! ~~~~~~~~~~~~
      ! For any calculations started in start_approximate_inverse that happen on 
      ! a subcomm, buffers%request should be null and hence only processors on subcomms finish
      ! ~~~~~~~~~~~~

      ! Gmres poylnomial with power or arnoldi basis
      if (inverse_type == PFLAREINV_POWER .OR. inverse_type == PFLAREINV_ARNOLDI) then

         call build_gmres_polynomial_inverse(matrix, poly_order, &
                  buffers, coefficients(:, 1), &
                  poly_sparsity_order, matrix_free, reuse_mat, inv_matrix)  

      ! Gmres polynomial with newton basis
      else if (inverse_type == PFLAREINV_NEWTON) then

         if (.NOT. matrix_free) then
            print *, "GMRES polynomial with Newton basis must be applied matrix-free"
            call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, errorcode)
         end if

         call build_gmres_polynomial_newton_inverse(matrix, poly_order, &
                           coefficients, &
                           inv_matrix)         

      ! Neumann polynomial
      else if (inverse_type == PFLAREINV_NEUMANN) then

         call calculate_and_build_neumann_polynomial_inverse(matrix, poly_order, &
                     buffers, poly_sparsity_order, matrix_free, reuse_mat, inv_matrix)        
                 
      ! Sparse approximate inverse
      else if (inverse_type == PFLAREINV_SAI .OR. inverse_type == PFLAREINV_ISAI) then

         ! ~~~~~~~~~~~
         ! ~~~~~~~~~~~               

         if (inverse_type == PFLAREINV_SAI) then
            incomplete = .FALSE.
         
         ! Incomplete sparse approximate inverse
         ! This is equivalent to a one-level restricted additive schwartz
         ! where each subdomain is a single unknown with the overlap 
         ! given by neighbouring unknowns            
         else
            incomplete = .TRUE.
         end if

         call calculate_and_build_sai(matrix, poly_sparsity_order, incomplete, reuse_mat, inv_matrix)
         
      ! Weighted jacobi
      else if (inverse_type == PFLAREINV_WJACOBI) then

         call calculate_and_build_weighted_jacobi_inverse(matrix, .TRUE., inv_matrix)          

      ! Unweighted jacobi
      else if (inverse_type == PFLAREINV_JACOBI) then

         call calculate_and_build_weighted_jacobi_inverse(matrix, .FALSE., inv_matrix)           

      end if

      ! If we created a gmres polynomial inverse, we may have created 
      ! a matrix on a subcomm
      if (inverse_type == PFLAREINV_POWER .OR. inverse_type == PFLAREINV_ARNOLDI .OR. inverse_type == PFLAREINV_NEWTON) then

         ! buffers%matrix could just be a pointer to the original matrix, so we 
         ! don't want to destroy it in that case
         if (matrix%v /= buffers%matrix%v) then
            ! If we're not on a subcomm there is no matrix
            if (buffers%matrix /= PETSC_NULL_MAT) then
               call MatDestroy(buffers%matrix, ierr)
               buffers%matrix = PETSC_NULL_MAT
            end if
         end if
      end if

   end subroutine finish_approximate_inverse
   
   ! -------------------------------------------------------------------------------------------------------------------------------

   subroutine reset_inverse_mat(matrix)

      ! Resets the matrix

      ! ~~~~~~
      type(tMat), intent(inout) :: matrix

      PetscErrorCode :: ierr
      MatType:: mat_type
      type(mat_ctxtype), pointer :: mat_ctx
      ! ~~~~~~

      if (matrix /= PETSC_NULL_MAT) then
         call MatGetType(matrix, mat_type, ierr)
         ! If its a matshell, make sure to delete its ctx
         if (mat_type==MATSHELL) then
            call MatShellGetContext(matrix, mat_ctx, ierr)
            if (mat_ctx%own_coefficients) deallocate(mat_ctx%coefficients)
            deallocate(mat_ctx)
         end if               
         call MatDestroy(matrix, ierr)
         matrix = PETSC_NULL_MAT
      end if  
      
   end subroutine reset_inverse_mat

   ! -------------------------------------------------------------------------------------------------------------------------------

end module approx_inverse_setup

