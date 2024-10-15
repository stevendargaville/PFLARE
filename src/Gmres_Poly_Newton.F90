module gmres_poly_newton

   use petsc
   use gmres_poly

#include "petsc/finclude/petsc.h"   

   implicit none

   public 
   
   contains

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine modified_leja(real_roots, imag_roots)

      ! Computes a modified leja ordering of the eigenvalues
      ! and re-orders in place

      ! ~~~~~~
      real, dimension(:), intent(inout)  :: real_roots, imag_roots

      ! Local variables
      integer :: i_loc, k_loc, counter
      integer :: max_loc(1)
      integer, dimension(:), allocatable :: indices
      real, dimension(:), allocatable :: magnitude
      real :: a, b, squares

      ! ~~~~~~    

      allocate(magnitude(size(real_roots)))
      allocate(indices(size(real_roots)))

      ! Compute the magnitudes of the evals
      do i_loc = 1, size(real_roots)
         magnitude(i_loc) = sqrt(real_roots(i_loc)**2 + imag_roots(i_loc)**2)
      end do
      ! Find the biggest
      max_loc = maxloc(magnitude)
      counter = 1

      ! That is our first entry
      indices(counter) = max_loc(1)
      counter = counter + 1
      ! If it was imaginary its complex conjugate is next
      if (imag_roots(indices(counter-1)) /= 0.0) then
         indices(counter) = indices(counter - 1) + 1
         counter = counter + 1
      end if

      ! Do while we still have some sorting to do
      do while (counter-1 < size(real_roots))

         ! For each value compute a product of differences
         do i_loc = 1, size(real_roots)

            magnitude(i_loc) = 1.0

            ! Loop over all those we've sorted so far
            k_loop: do k_loc = 1, counter-1
               
               ! Distance
               a = real_roots(i_loc) - real_roots(indices(k_loc))
               b = imag_roots(i_loc) - imag_roots(indices(k_loc))

               squares = a**2 + b**2
               if (squares > 1e-14) then
                  magnitude(i_loc) = magnitude(i_loc) + &
                     log10(sqrt(squares))
               else
                  magnitude(i_loc) = -huge(0)
                  exit k_loop
               end if
            end do k_loop
         end do

         ! The new index is the biggest in distance 
         max_loc = maxloc(magnitude)
         indices(counter) = max_loc(1)
         counter = counter + 1
         ! If it was imaginary its complex conjugate is next
         if (imag_roots(indices(counter - 1)) /= 0.0) then
            indices(counter) = indices(counter - 1) + 1
            counter = counter + 1
         end if
      end do

      ! Reorder the roots
      real_roots = real_roots(indices)
      imag_roots = imag_roots(indices)

      deallocate(magnitude, indices)

   end subroutine modified_leja   

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine calculate_gmres_polynomial_roots_newton(matrix, poly_order, coefficients)

      ! Computes a fixed order gmres polynomial for the matrix passed in
      ! and outputs the Harmonic Ritz values (ie the roots) which we can use to apply 
      ! a polynomial in the Newton basis
      ! This should be stable at high order, although we don't add extra roots like Loe 
      ! for stability, so we are likely not as stable. The only reason we don't is 
      ! that would change the number of roots.
      ! The cost of using the newton basis is many reductions in parallel.
      ! We don't provide a way to compute the monomial polynomial coefficients from the roots 
      ! (although we could by rearranging the newton polynomial) as I don't think this would be 
      ! stable at high order anyway. The only reason we use the monomial coefficients 
      ! is to easily build an assembled (approximate) matrix inverse, and 
      ! f you want that with these roots you could build one directly (although it would not 
      ! be very sparse at high order, and a fixed sparsity version at higher order than either the 
      ! power basis or Arnoldi basis could compute is likely not any better)

      ! ~~~~~~
      type(tMat), intent(in)                            :: matrix
      integer, intent(in)                               :: poly_order
      real, dimension(:, :), intent(out)                :: coefficients

      ! Local variables
      PetscInt :: global_rows, global_cols, local_rows, local_cols
      integer :: lwork, subspace_size, rank, i_loc, comm_size, comm_rank, errorcode
      PetscErrorCode :: ierr      
      MPI_Comm :: MPI_COMM_MATRIX
      real, dimension(poly_order+2,poly_order+1) :: H_n
      real, dimension(poly_order+1,poly_order+1) :: H_n_square
      real, dimension(poly_order+1) :: e_d, r, c, solution, gamma
      integer, dimension(poly_order+1) :: iwork, pivots
      real, dimension(1) :: ferr, berr
      real, dimension(:), allocatable :: work
      real, dimension(:,:), allocatable :: VL, VR
      real :: beta, rcond, avg_real, avg_imag
      real, dimension(:), pointer :: v_one, v_two
      real, dimension(:, :), allocatable, target :: K_m_plus_1_data, K_m_plus_1_data_not_zero
      real, dimension(:, :), pointer :: K_m_plus_1_pointer
      type(tVec) :: w_j
      type(tVec), dimension(poly_order+2) :: V_n
      character(1) :: equed
      logical :: use_harmonic_ritz = .TRUE.

      ! ~~~~~~    

      ! This is how many columns we have in K_m
      subspace_size = poly_order + 1

      ! We might want to call the gmres poly creation on a sub communicator
      ! so let's get the comm attached to the matrix and make sure to use that 
      call PetscObjectGetComm(matrix, MPI_COMM_MATRIX, ierr)  
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)     
      ! Get the comm rank
      call MPI_Comm_rank(MPI_COMM_MATRIX, comm_rank, errorcode)     
      ! Get the matrix sizes
      call MatGetSize(matrix, global_rows, global_cols, ierr)
      call MatGetLocalSize(matrix, local_rows, local_cols, ierr)

      if (subspace_size > global_rows) then
         print *, "The input subspace size is greater than the matrix size"
         call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, errorcode)
      end if

      ! ~~~~~~~~~~
      ! Allocate space and create random numbers 
      ! There are individual petsc vecs built in V_n that point at 
      ! the data in K_m_plus_1_data
      ! The first vec has random numbers in it
      ! ~~~~~~~~~~ 
      call create_temp_space_box_mueller(MPI_COMM_MATRIX, comm_size, comm_rank, &
               local_rows, global_rows, subspace_size, &
               K_m_plus_1_data, K_m_plus_1_data_not_zero, K_m_plus_1_pointer, &
               V_n)
      
      ! Create an extra vector for storage
      call VecDuplicate(V_n(1), w_j, ierr)      
      
      ! Do the Arnoldi and compute H_n
      ! Use the same lucky tolerance as petsc
      call arnoldi(matrix, poly_order, 1e-30, V_n, w_j, beta, H_n)

      ! ~~~~~~~~~~~
      ! Now the Ritz values are just the eigenvalues of the square part of H_n
      ! We're actually going to use the Harmonic Ritz values
      ! see Embree - Polynomial preconditioned arnoldi with stability control
      ! ~~~~~~~~~~~
      if (use_harmonic_ritz) then

         e_d = 0
         e_d(poly_order + 1) = 1

         allocate(work(4 * (poly_order + 1)))
         ! Compute H_d^-H e_d
         ! Doesn't modify H_n on output
         call dgesvx('N', 'T', poly_order + 1, 1, &
               H_n, size(H_n, 1), H_n_square, poly_order + 1, &
               pivots, equed, r, c, &
               e_d, poly_order + 1, &
               solution, poly_order + 1, &
               rcond, ferr, berr, work, iwork, errorcode)
         deallocate(work)
         if (errorcode /= 0) then
            print *, "Harmonic Ritz solve failed - Try lower order"
            call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, errorcode)
         end if
         ! Rearrange given the row permutations done by the LU
         solution(pivots) = solution

         ! Scale f by H(d+1,d)^2
         solution = solution * H_n(poly_order + 2, poly_order + 1)**2
         ! Add to the last column of H_n
         H_n(1:poly_order + 1, poly_order + 1) = &
            H_n(1:poly_order + 1, poly_order + 1) + solution
      end if

      ! ~~~~~~~~~~~~~~
      ! Now compute the eigenvalues of the square part of H_n
      ! ie either compute the Ritz or Harmonic Ritz values
      ! ~~~~~~~~~~~~~~

      allocate(work(1))
      lwork = -1      

      ! Compute the eigenvalues
      call dgeev('N', 'N', poly_order + 1, H_n, size(H_n, 1), &
                     coefficients(1, 1), coefficients(1, 2), VL, 1, VR, 1, &
                     work, lwork, errorcode)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork)) 
      call dgeev('N', 'N', poly_order + 1, H_n, size(H_n, 1), &
                     coefficients(1, 1), coefficients(1, 2), VL, 1, VR, 1, &
                     work, lwork, errorcode)
      deallocate(work)

      if (errorcode /= 0) then
         print *, "Eig decomposition failed"
         call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, errorcode)
      end if  

      ! ~~~~~~~~~~~~~~
      ! If you want to add roots for stability that would go here
      ! ~~~~~~~~~~~~~~       

      ! ~~~~~~~~~~~~~~
      ! Now compute a modified leja ordering for stability
      ! ~~~~~~~~~~~~~~      
      call modified_leja(coefficients(:,1), coefficients(:,2))

      ! Cleanup
      deallocate(K_m_plus_1_data)
      if (allocated(K_m_plus_1_data_not_zero)) deallocate(K_m_plus_1_data_not_zero)
      do i_loc = 1, subspace_size+1
         call VecDestroy(V_n(i_loc), ierr)
      end do
      call VecDestroy(w_j, ierr)

   end subroutine calculate_gmres_polynomial_roots_newton   

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine petsc_matvec_gmres_newton_mf(mat, x, y)

      ! Applies a gmres polynomial in the newton basis matrix-free as an inverse
      ! The roots are stored in mat_ctx%real_roots, mat_ctx%imag_roots in the input matshell
      ! Based on Loe 2021 Toward efficient polynomial preconditioning for GMRES
      ! y = A x

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! Input
      type(tMat), intent(in)    :: mat
      type(tVec) :: x
      type(tVec) :: y

      ! Local
      integer :: order, errorcode
      PetscErrorCode :: ierr      
      type(mat_ctxtype), pointer :: mat_ctx => null()
      type(tVec) :: prod, temp_result, temp

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~

      call MatShellGetContext(mat, mat_ctx, ierr)
      if (.NOT. associated(mat_ctx%real_roots)) then
         print *, "Roots in context aren't found"
         call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, errorcode)
      end if

      ! prod = x
      call VecDuplicate(x, prod, ierr)     
      call VecDuplicate(x, temp_result, ierr)         
      call VecDuplicate(x, temp, ierr)         
      call VecCopy(x, prod, ierr)
      ! y = 0
      call VecSet(y, 0.0, ierr)

      ! ~~~~~~~~~~~~
      ! Iterate over the order
      ! ~~~~~~~~~~~~
      order = 1
      do while (order .le. size(mat_ctx%real_roots) - 1)

         ! If real this is easy
         if (mat_ctx%imag_roots(order) == 0.0) then

            ! y = y + theta_i * prod
            call VecAXPBY(y, &
                     1.0/mat_ctx%real_roots(order), &
                     1.0, &
                     prod, ierr)    
                     
            ! temp_result = A * prod
            call MatMult(mat_ctx%mat, prod, temp_result, ierr)
            ! prod = prod - theta_i * temp_result
            call VecAXPBY(prod, &
                     -1.0/mat_ctx%real_roots(order), &
                     1.0, &
                     temp_result, ierr) 

            order = order + 1

         ! If imaginary, then have to combine the e'val and its
         ! complex conjugate to keep the arithmetic real
         ! Relies on the complex conjugate being next to each other
         else

            ! temp_result = A * prod
            call MatMult(mat_ctx%mat, prod, temp_result, ierr)    
            ! temp_result = 2 * Re(theta_i) * prod - temp_result
            call VecAXPBY(temp_result, &
                  2 * mat_ctx%real_roots(order), &
                  -1.0, &
                  prod, ierr)

            ! y = y + 1/(Re(theta_i)^2 + Imag(theta_i)^2) * temp_result
            call VecAXPBY(y, &
                     1.0/(mat_ctx%real_roots(order)**2 + mat_ctx%imag_roots(order)**2), &
                     1.0, &
                     temp_result, ierr)  
                     
            if (order .le. size(mat_ctx%real_roots) - 2) then
               ! temp = A * temp_result
               call MatMult(mat_ctx%mat, temp_result, temp, ierr)    

               ! prod = prod - 1/(Re(theta_i)^2 + Imag(theta_i)^2) * temp
               call VecAXPBY(prod, &
                        -1.0/(mat_ctx%real_roots(order)**2 + mat_ctx%imag_roots(order)**2), &
                        1.0, &
                        temp, ierr)               
            end if

            ! Skip two evals
            order = order + 2

         end if
      end do

      ! Final step if last root is real
      if (mat_ctx%imag_roots(size(mat_ctx%real_roots)) == 0.0) then
         ! y = y + theta_i * prod
         call VecAXPBY(y, &
                  1.0/mat_ctx%real_roots(order), &
                  1.0, &
                  prod, ierr)          
      end if

      call VecDestroy(temp_result, ierr)
      call VecDestroy(prod, ierr)
      call VecDestroy(temp, ierr)

   end subroutine petsc_matvec_gmres_newton_mf      

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine build_gmres_polynomial_newton_inverse(matrix, poly_order, &
                  coefficients, &
                  inv_matrix)

      ! Assembles a matrix which is an approximation to the inverse of a matrix using the 
      ! gmres polynomial in the newton basis
      ! Can only be applied matrix-free

      ! ~~~~~~
      type(tMat), intent(in)                                      :: matrix
      integer, intent(in)                                         :: poly_order
      real, dimension(:, :), target, contiguous, intent(inout)    :: coefficients
      type(tMat), intent(inout)                                   :: inv_matrix

      ! Local variables
      PetscInt :: global_rows, global_cols, local_rows, local_cols
      integer :: comm_size, errorcode, order
      PetscErrorCode :: ierr      
      MPI_Comm :: MPI_COMM_MATRIX
      type(mat_ctxtype), pointer :: mat_ctx

      ! ~~~~~~       

      ! We might want to call the gmres poly creation on a sub communicator
      ! so let's get the comm attached to the matrix and make sure to use that 
      call PetscObjectGetComm(matrix, MPI_COMM_MATRIX, ierr)    
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)        

      ! Get the local sizes
      call MatGetLocalSize(matrix, local_rows, local_cols, ierr)
      call MatGetSize(matrix, global_rows, global_cols, ierr)   
      
      ! ~~~~~~~
      ! Just build a matshell that applies our polynomial matrix-free
      ! ~~~~~~~

      ! If not re-using
      if (inv_matrix == PETSC_NULL_MAT) then

         ! Have to dynamically allocate this
         allocate(mat_ctx)

         ! We pass in the polynomial coefficients as the context
         call MatCreateShell(MPI_COMM_MATRIX, local_rows, local_cols, global_rows, global_cols, &
                     mat_ctx, inv_matrix, ierr)
         ! The subroutine petsc_matvec_gmres_newton_mf applies the polynomial inverse
         call MatShellSetOperation(inv_matrix, &
                     MATOP_MULT, petsc_matvec_gmres_newton_mf, ierr)

         call MatAssemblyBegin(inv_matrix, MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(inv_matrix, MAT_FINAL_ASSEMBLY, ierr)         

      ! Reusing 
      else
         call MatShellGetContext(inv_matrix, mat_ctx, ierr)

      end if

      mat_ctx%real_roots => coefficients(:, 1)
      mat_ctx%imag_roots => coefficients(:, 2)
      ! Now because the context reset deallocates the coefficient pointer 
      ! we want to make sure we don't leak memory, so we use pointer remapping here 
      ! to turn the 2D coefficient pointer into a 1D that we can store in mat_ctx%coefficients
      ! and then the deallocate on mat_ctx%coefficients should still delete all the memory
      mat_ctx%coefficients(1:2*size(coefficients,1)) => coefficients(:, :)
      ! This is the matrix whose inverse we are applying (just copying the pointer here)
      mat_ctx%mat = matrix       

   end subroutine build_gmres_polynomial_newton_inverse       

! -------------------------------------------------------------------------------------------------------------------------------


end module gmres_poly_newton

