module gmres_poly_newton

   use petsc
   use gmres_poly

#include "petsc/finclude/petsc.h"   

   implicit none

#include "petsc_legacy.h"   

   public 
   
   contains

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine modified_leja(real_roots, imag_roots, indices)

      ! Computes a modified leja ordering of the eigenvalues
      ! and re-orders in place
      ! The roots must be passed in with complex conjugate pairs next to each other
      ! with positive imag e'vals first

      ! ~~~~~~
      real, dimension(:), intent(inout)  :: real_roots, imag_roots
      integer, dimension(:), allocatable, intent(inout) :: indices

      ! Local variables
      integer :: i_loc, k_loc, counter
      integer :: max_loc(1)
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
      ! If it was imaginary its complex conjugate is next to this one
      if (imag_roots(indices(counter-1)) /= 0.0) then
         ! dgeev returns roots with the positive imaginary part first
         if (imag_roots(indices(counter-1)) > 0) then
            ! So if positive we know the conjugate is ahead one
            indices(counter) = indices(counter - 1) + 1
         else
            ! If negative we know the conjugate is one behind
            indices(counter) = indices(counter - 1) - 1
         end if
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
               ! This is just a test for non zero, as squares can't be negative
               ! This ensures if we have already added in a complex root, we skip its
               ! complex conjugate 
               if (squares > 0) then
                  magnitude(i_loc) = magnitude(i_loc) + &
                     log10(sqrt(squares))
               else
                  ! Be careful here to use huge(0.0) rather than huge(0)!
                  magnitude(i_loc) = -huge(0.0)
                  exit k_loop
               end if
            end do k_loop
         end do

         ! The new index is the biggest in distance 
         max_loc = maxloc(magnitude)
         indices(counter) = max_loc(1)
         counter = counter + 1
         ! If it was imaginary its complex conjugate is next to this one
         if (imag_roots(indices(counter - 1)) /= 0.0) then
            ! dgeev returns roots with the positive imaginary part first
            ! We did actually find a situation where due to rounding differnces, the 
            ! conjugate with negative imaginary root had a bigger product of factors
            ! (by 1e-14) than the positive imaginary root, so we can't guarantee we hit 
            ! the positive one first
            if (imag_roots(indices(counter-1)) > 0) then
               ! So if positive we know the conjugate is ahead one
               indices(counter) = indices(counter - 1) + 1
            else
               ! If negative we know the conjugate is one behind
               indices(counter) = indices(counter - 1) - 1
            end if            
            counter = counter + 1
         end if
      end do

      deallocate(magnitude)

   end subroutine modified_leja   

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine calculate_gmres_polynomial_roots_newton(matrix, poly_order, add_roots, coefficients)

      ! Computes a fixed order gmres polynomial for the matrix passed in
      ! and outputs the Harmonic Ritz values (ie the roots) which we can use to apply 
      ! a polynomial in the Newton basis. This should be stable at high order
      ! The cost of computing the newton basis is many reductions in parallel.
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
      logical, intent(in)                               :: add_roots
      real, dimension(:, :), pointer, intent(inout)     :: coefficients

      ! Local variables
      PetscInt :: global_rows, global_cols, local_rows, local_cols
      integer :: lwork, subspace_size, rank, i_loc, comm_size, comm_rank, errorcode, iwork_size, j_loc
      integer :: total_extra, counter, k_loc
      PetscErrorCode :: ierr      
      MPI_Comm :: MPI_COMM_MATRIX
      real, dimension(poly_order+2,poly_order+1) :: H_n
      real, dimension(poly_order+1,poly_order+2) :: H_n_T
      real, dimension(poly_order+1,poly_order+1) :: H_n_square
      real, dimension(poly_order+1) :: e_d, r, col_scale, solution, gamma, s, pof
      integer, dimension(poly_order+1) :: iwork, pivots, extra_pair_roots, overflow
      integer, dimension(:), allocatable :: iwork_allocated, indices
      real, dimension(1) :: ferr, berr
      real, dimension(:), allocatable :: work
      real, dimension(:,:), allocatable :: VL, VR
      real :: beta, div_real, div_imag, a, b, c, d, div_mag
      real, dimension(:), pointer :: v_one, v_two
      real, dimension(:, :), allocatable, target :: K_m_plus_1_data, K_m_plus_1_data_not_zero
      real, dimension(:, :), allocatable :: coefficients_temp
      real, dimension(:, :), pointer :: K_m_plus_1_pointer
      type(tVec) :: w_j
      type(tVec), dimension(poly_order+2) :: V_n
      character(1) :: equed
      logical :: use_harmonic_ritz = .TRUE.
      real :: rcond = 1e-12

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
      ! We're actually going to use the Harmonic Ritz values which are the reciprocals
      ! of the (ordinary) Ritz values of A^-1
      ! see Embree - Polynomial preconditioned arnoldi with stability control 
      ! and Goossens - Ritz and Harmonic Ritz Values and the Convergence of FOM and GMRES
      ! ~~~~~~~~~~~
      if (use_harmonic_ritz) then

         e_d = 0
         e_d(poly_order + 1) = 1

         ! Now we have to be careful here 
         ! PETSc and Trilinos both just use an LU factorisation here to solve
         ! H_d^-H e_d
         ! But we can have the situation where the user has requested a high 
         ! polynomial order and H_d ends up not full rank
         ! In that case the LU factorisation fails 
         ! Previously we just used to abort and tell the user to use lower order
         ! But that is not super helpful
         ! Instead we use a rank revealing factorisation that gives the 
         ! minimum norm solution
         ! What we find is that when use this to compute eigenvalues we find e-vals 
         ! as we might expect up to the rank
         ! but then we have some eigenvalues that are numerically zero
         ! We keep those and our application of the newton polynomial in 
         ! petsc_matvec_gmres_newton_mf just skips them and hence we don't do any 
         ! extra work in the application phase than we would have done with lower order

         ! ~~~~~~~~~~~
         ! Direct LU
         ! ~~~~~~~~~~~

         ! allocate(work(4 * (poly_order + 1)))
         ! ! Compute H_d^-H e_d
         ! ! Doesn't modify H_n on output
         ! call dgesvx('N', 'T', poly_order + 1, 1, &
         !       H_n, size(H_n, 1), H_n_square, poly_order + 1, &
         !       pivots, equed, r, col_scale, &
         !       e_d, poly_order + 1, &
         !       solution, poly_order + 1, &
         !       rcond, ferr, berr, work, iwork, errorcode)
         ! deallocate(work)

         ! Rearrange given the row permutations done by the LU
         !solution(pivots) = solution        
         
         ! ~~~~~~~~~~~
         ! ~~~~~~~~~~~

         ! Rank revealing min norm solution
         H_n_T = transpose(H_n)

         allocate(work(1))
         allocate(iwork_allocated(1))
         lwork = -1         

         ! We have rcond = 1e-12 which is used to decide what singular values to drop
         ! Matlab uses max(size(A))*eps(norm(A)) in their pinv
         call dgelsd(poly_order + 1, poly_order + 1, 1, H_n_T, size(H_n_T, 1), &
                        e_d, size(e_d), s, rcond, rank, &
                        work, lwork, iwork_allocated, errorcode)
         lwork = int(work(1))
         iwork_size = iwork_allocated(1)
         deallocate(work, iwork_allocated)
         allocate(work(lwork)) 
         allocate(iwork_allocated(iwork_size))
         call dgelsd(poly_order + 1, poly_order + 1, 1, H_n_T, size(H_n_T, 1), &
                        e_d, size(e_d), s, rcond, rank, &
                        work, lwork, iwork_allocated, errorcode)
         deallocate(work, iwork_allocated)         

         ! Copy in the solution
         solution = e_d

         if (errorcode /= 0) then
            print *, "Harmonic Ritz solve failed"
            call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, errorcode)
         end if

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
      ! Add roots for stability
      ! ~~~~~~~~~~~~~~         
      if (add_roots) then 

         ! Compute the product of factors
         pof = 1   
         extra_pair_roots = 0
         overflow = 0
         total_extra = 0
         do k_loc = 1, poly_order + 1

            a = coefficients(k_loc, 1)
            b = coefficients(k_loc, 2)      
            
            ! We have already computed pof for the positive imaginary complex conjugate
            if (b < 0) cycle

            ! Skips eigenvalues that are numerically zero
            if (abs(a) < 1e-12) cycle
            if (a**2 + b**2 < 1e-12) cycle

            ! Compute product(k)_{i, j/=i} * | 1 - theta_j/theta_i|
            do i_loc = 1, poly_order + 1

               ! Skip
               if (k_loc == i_loc) cycle

               c = coefficients(i_loc, 1)
               d = coefficients(i_loc, 2)

               ! Skips eigenvalues that are numerically zero
               if (abs(c) < 1e-12) cycle
               if (c**2 + d**2 < 1e-12) cycle

               ! theta_k/theta_i
               div_real = (a * c + b * d)/(c**2 + d**2)
               div_imag = (b * c - a * d)/(c**2 + d**2)

               ! |1 - theta_k/theta_i|
               div_mag = sqrt((1 - div_real)**2 + div_imag**2)

               ! Pof is about to overflow, store the exponent and 
               ! reset pof back to one
               ! We can hit this for very high order polynomials, where we have to 
               ! add more roots than 22 (ie pof > 1e308)
               if (log10(pof(k_loc)) + log10(div_mag) > 307) then
                  overflow(k_loc) = overflow(k_loc) + log10(pof(k_loc))
                  pof(k_loc) = 1
               end if            

               ! Product
               pof(k_loc) = pof(k_loc) * div_mag

            end do

            ! If pof > 1e4, we add an extra root, plus one extra for every 1e14
            if (log10(pof(k_loc)) > 4 .OR. overflow(k_loc) /= 0) then

               ! if real extra_pair_roots counts each distinct real root we're adding
               ! if imaginary it only counts a pair as one
               extra_pair_roots(k_loc) = ceiling((log10(pof(k_loc)) + overflow(k_loc) - 4.0)/14.0)
               total_extra = total_extra + extra_pair_roots(k_loc)

               ! If imaginary, the pof is the same for the conjugate, let's just set it to -1
               if (b > 0) then
                  ! We know the positive imaginary value is first, so the conjugate follows it
                  pof(k_loc+1) = -1
                  ! We need the conjugates as well
                  total_extra = total_extra + extra_pair_roots(k_loc)

               end if            
            end if
         end do

         ! If we have extra roots we need to resize the coefficients storage
         if (total_extra > 0) then
            allocate(coefficients_temp(size(coefficients, 1), size(coefficients, 2)))
            coefficients_temp(1:size(coefficients, 1), 1:size(coefficients, 2)) = coefficients
            deallocate(coefficients)
            allocate(coefficients(size(coefficients_temp, 1) + total_extra, 2))
            coefficients = 0
            coefficients(1:size(coefficients_temp, 1), :) = coefficients_temp
            deallocate(coefficients_temp)
         end if
      end if

      ! Take a copy of the existing roots
      coefficients_temp = coefficients

      if (add_roots) then
         
         ! Add the extra copies of roots, ensuring conjugate pairs we add 
         ! are next to each other
         counter = size(extra_pair_roots)+1
         do i_loc = 1, size(extra_pair_roots)

            ! For each extra root pair to add
            do j_loc = 1, extra_pair_roots(i_loc)

               coefficients(counter, :) = coefficients(i_loc, :)
               ! Add in the conjugate
               if (coefficients(i_loc, 2) > 0) then
                  coefficients(counter+1, 1) = coefficients(i_loc, 1)
                  coefficients(counter+1, 2) = -coefficients(i_loc, 2)
               end if

               ! Store a perturbed root so we have unique values for the leja sort below
               ! Just peturbing the real value
               coefficients_temp(counter, 1) = coefficients(i_loc, 1) + j_loc * 5e-8
               coefficients_temp(counter, 2) = coefficients(i_loc, 2)
               ! Add in the conjugate
               if (coefficients(i_loc, 2) > 0) then
                  coefficients_temp(counter+1, 1) = coefficients(i_loc, 1) + j_loc * 5e-8
                  coefficients_temp(counter+1, 2) = -coefficients(i_loc, 2)
               end if            

               counter = counter + 1
               if (coefficients(i_loc, 2) > 0) counter = counter + 1
            end do
         end do
      end if

      ! ~~~~~~~~~~~~~~
      ! Now compute a modified leja ordering for stability
      ! ~~~~~~~~~~~~~~      
      ! Called with the peturbed extra roots
      call modified_leja(coefficients_temp(:,1), coefficients_temp(:,2), indices)   

      ! Reorder the (non-peturbed) roots 
      coefficients(:,1) = coefficients(indices,1)
      coefficients(:,2) = coefficients(indices,2)

      ! Cleanup
      deallocate(coefficients_temp)
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

            ! Skips eigenvalues that are numerically zero - see 
            ! the comment in calculate_gmres_polynomial_roots_newton 
            if (abs(mat_ctx%real_roots(order)) < 1e-12) then
               order = order + 1
               cycle
            end if

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

            ! Skips eigenvalues that are numerically zero
            if (mat_ctx%real_roots(order)**2 + mat_ctx%imag_roots(order)**2 < 1e-12) then
               order = order + 2
               cycle
            end if            

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

         ! Skips eigenvalues that are numerically zero
         if (abs(mat_ctx%real_roots(order)) > 1e-12) then

            ! y = y + theta_i * prod
            call VecAXPBY(y, &
                     1.0/mat_ctx%real_roots(order), &
                     1.0, &
                     prod, ierr) 
         end if
      end if

      call VecDestroy(temp_result, ierr)
      call VecDestroy(prod, ierr)
      call VecDestroy(temp, ierr)

   end subroutine petsc_matvec_gmres_newton_mf      

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine build_gmres_polynomial_newton_inverse(matrix, poly_order, &
                  coefficients, &
                  inv_matrix)

      ! Builds a matrix which is an approximation to the inverse of a matrix using the 
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
      if (PetscMatIsNull(inv_matrix)) then

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

