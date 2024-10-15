module cf_splitting

   use petsc
   use pmisr_ddc
   use aggregation

#include "petsc/finclude/petsc.h"

   implicit none

   public   

   PetscEnum, parameter :: CF_PMISR_DDC=0
   PetscEnum, parameter :: CF_PMIS=1
   PetscEnum, parameter :: CF_PMIS_DIST2=2
   PetscEnum, parameter :: CF_AGG=3
   
   contains

! -------------------------------------------------------------------------------------------------------------------------------
   
   subroutine create_cf_is(input_mat, cf_markers_local, is_fine, is_coarse)

      ! Creates IS's for C and F points given a cf_markers_local array
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(in) :: input_mat
      integer, allocatable, dimension(:), intent(in) :: cf_markers_local
      type(tIS), intent(inout) :: is_fine, is_coarse

      PetscInt, dimension(:), allocatable :: fine_indices, coarse_indices
      PetscInt :: global_row_start, global_row_end_plus_one
      PetscInt :: fine_counter, coarse_counter, i_loc
      integer :: n, n_fine, n_coarse
      PetscErrorCode :: ierr
      MPI_Comm :: MPI_COMM_MATRIX      
      
      ! ~~~~~~~~~~

      call PetscObjectGetComm(input_mat, MPI_COMM_MATRIX, ierr)    
      call MatGetOwnershipRange(input_mat, global_row_start, global_row_end_plus_one, ierr)   

      ! This is the local_rows
      n = size(cf_markers_local)
      n_fine = count(cf_markers_local == F_POINT)
      n_coarse = count(cf_markers_local == C_POINT) 

      allocate(fine_indices(n_fine))
      allocate(coarse_indices(n_coarse))

      fine_counter = 0
      coarse_counter = 0
      do i_loc = 1, n
         if (cf_markers_local(i_loc) == F_POINT) then
            fine_counter = fine_counter + 1
            ! Remember petsc indices start at zero
            fine_indices(fine_counter) = global_row_start + i_loc-1
         else
            coarse_counter = coarse_counter + 1
            ! Remember petsc indices start at zero
            coarse_indices(coarse_counter) = global_row_start + i_loc-1
         end if
      end do   
      
      call ISCreateGeneral(MPI_COMM_MATRIX, fine_counter, &
            fine_indices(1:fine_counter), &
            PETSC_COPY_VALUES, is_fine, ierr)  
      call ISCreateGeneral(MPI_COMM_MATRIX, coarse_counter, &
            coarse_indices(1:coarse_counter), &
            PETSC_COPY_VALUES, is_coarse, ierr)   

      deallocate(fine_indices, coarse_indices)
         
   end subroutine create_cf_is
   
!------------------------------------------------------------------------------------------------------------------------
   
   subroutine generate_sabs(input_mat, strong_threshold, symmetrize, square, output_mat, allow_drop_diagonal)
      
      ! Generate strength of connection matrix with absolute value 
      ! Output has no diagonal entries
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(in)     :: input_mat
      type(tMat), intent(inout)  :: output_mat
      real, intent(in)           :: strong_threshold
      logical, intent(in)        :: symmetrize, square
      logical, intent(in), optional :: allow_drop_diagonal
      
      PetscInt :: col, ncols, ifree, max_nnzs
      PetscInt :: local_rows, local_cols, global_rows, global_cols, global_row_start, global_row_end_plus_one
      PetscInt :: global_col_start, global_col_end_plus_one, jfree, diagonal_index
      integer :: counter, errorcode, comm_size
      PetscErrorCode :: ierr
      PetscInt, dimension(:), allocatable :: nnzs_row, onzs_row, cols, cols_mod
      real, dimension(:), allocatable :: vals, vals_copy
      PetscInt, parameter :: nz_ignore = -1, one=1, zero=0
      real :: rel_row_tol, abs_biggest_entry
      MPI_Comm :: MPI_COMM_MATRIX
      type(tMat) :: transpose_mat, temp_mat
      type(tIS) :: zero_diags
      PetscInt, dimension(:), pointer :: zero_diags_pointer
      logical :: drop_diag
      
      ! ~~~~~~~~~~

      drop_diag = .TRUE.
      if (present(allow_drop_diagonal)) drop_diag = allow_drop_diagonal

      call PetscObjectGetComm(input_mat, MPI_COMM_MATRIX, ierr)    
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)

      ! Get the local sizes
      call MatGetLocalSize(input_mat, local_rows, local_cols, ierr)
      call MatGetSize(input_mat, global_rows, global_cols, ierr)
      ! This returns the global index of the local portion of the matrix
      call MatGetOwnershipRange(input_mat, global_row_start, global_row_end_plus_one, ierr)  
      call MatGetOwnershipRangeColumn(input_mat, global_col_start, global_col_end_plus_one, ierr)  
      
      allocate(nnzs_row(local_rows))
      nnzs_row = 0

      max_nnzs = 0
      if (comm_size/=1) then
         allocate(onzs_row(local_rows))
         onzs_row = 0
      end if

      do ifree = global_row_start, global_row_end_plus_one-1                  

         call MatGetRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr)
         if (ncols > max_nnzs) max_nnzs = ncols
         call MatRestoreRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr)
      end do

      allocate(cols(max_nnzs))
      allocate(cols_mod(max_nnzs))
      allocate(vals(max_nnzs)) 
      allocate(vals_copy(max_nnzs))

      ! ~~~~~~~~~~
      ! Get the nnzs
      ! ~~~~~~~~~~
      
      ! Loop over global row indices
      do ifree = global_row_start, global_row_end_plus_one-1                  
         ! Get the row
         call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)

         abs_biggest_entry = -huge(0)
         ! Find the biggest entry in the row thats not the diagonal and the diagonal index
         do jfree = 1, ncols
            if (cols(jfree) == ifree) then
               diagonal_index = jfree
            else if (abs(vals(jfree)) > abs_biggest_entry) then
               abs_biggest_entry = abs(vals(jfree))
            end if 
         end do          
                
         ! Abs biggest entry in the row
         rel_row_tol = strong_threshold * abs_biggest_entry

         do col = 1, ncols

            ! Have to count the diagonal and off-diagonal nnzs in parallel           
            if (cols(col) .ge. global_col_start .AND. cols(col) .le. global_col_end_plus_one - 1) then

               !if (abs(vals(col)) .ge. rel_row_tol) then
               if (abs(vals(col)) .ge. rel_row_tol .OR. (.NOT. drop_diag .AND. cols(col) == ifree)) then
                  ! Convert to local row indices
                  nnzs_row(ifree - global_row_start + 1) = nnzs_row(ifree - global_row_start + 1) + 1
               end if
            else
               if (abs(vals(col)) .ge. rel_row_tol .OR. (.NOT. drop_diag .AND. cols(col) == ifree)) then
                  ! Convert to local row indices
                  onzs_row(ifree - global_row_start + 1) = onzs_row(ifree - global_row_start + 1) + 1
               end if               
            end if
         end do

         ! Must call otherwise petsc leaks memory
         call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)
      end do

      ! Create the output matrix
      if (comm_size/=1) then
         call MatCreateAIJ(MPI_COMM_MATRIX, local_rows, local_cols, &
                  global_rows, global_cols, &
                  nz_ignore, nnzs_row, &
                  nz_ignore, onzs_row, &
                  output_mat, ierr)   
      else
         call MatCreateSeqAIJ(MPI_COMM_MATRIX, global_rows, global_cols, nz_ignore, nnzs_row, &
                  output_mat, ierr)            
      end if   
       
      ! Just in case there are some zeros in the input mat, ignore them
      call MatSetOption(output_mat, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)     
      ! Don't set any off processor entries so no need for a reduction when assembling
      call MatSetOption(output_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)      
      call MatSetOption(output_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE,  ierr)     

      ! ~~~~~~~~~~
      ! Get the nnzs
      ! ~~~~~~~~~~      
      
      ! Now go and fill the new matrix
      ! Loop over global row indices
      do ifree = global_row_start, global_row_end_plus_one-1                  
      
         ! Get the row
         call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)  
         
         abs_biggest_entry = -huge(0)
         ! Find the biggest entry in the row thats not the diagonal and the diagonal index
         do jfree = 1, ncols
            if (cols(jfree) == ifree) then
               diagonal_index = jfree
            else if (abs(vals(jfree)) > abs_biggest_entry) then
               abs_biggest_entry = abs(vals(jfree))
            end if 
         end do            
         
         ! Abs biggest entry in the row
         rel_row_tol = strong_threshold * abs_biggest_entry

         cols_mod(1:ncols) = cols(1:ncols)
         ! Just set one to indicate strength
         vals_copy(1:ncols) = 1.0
                  
         do col = 1, ncols
            ! Get petsc to not insert the value here - drop diagonals
            if (abs(vals(col)) < rel_row_tol .OR. (drop_diag .AND. cols(col) == ifree)) then
               cols_mod(col) = -1
            end if
         end do

         ! Much quicker to call setvalues
         call MatSetValues(output_mat, one, ifree, ncols, cols_mod, &
                  vals_copy, INSERT_VALUES, ierr)

         ! Must call otherwise petsc leaks memory
         call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)   
      end do           
      
      call MatAssemblyBegin(output_mat, MAT_FINAL_ASSEMBLY, ierr)

      deallocate(cols, vals, nnzs_row, cols_mod, vals_copy)
      if (comm_size/=1) deallocate(onzs_row)
      call MatAssemblyEnd(output_mat, MAT_FINAL_ASSEMBLY, ierr)  

      ! Now symmetrize if desired
      if (symmetrize) then

         ! We could just do a symbolic transpose and add the two sets of indices together, 
         ! but its so much simpler to just add the two together - and the symbolic will be the expensive part
         ! anyway
         call MatTranspose(output_mat, MAT_INITIAL_MATRIX, transpose_mat, ierr)
         call MatAXPY(output_mat, 1.0, transpose_mat, DIFFERENT_NONZERO_PATTERN, ierr)     

         ! Don't forget to destroy the explicit transpose
         call MatDestroy(transpose_mat, ierr)

      end if

      ! Square the strength matrix to aggressively coarsen (gives a distance 2 MIS)
      if (square) then

         if (symmetrize) then
            call MatMatMult(output_mat, output_mat, &
                        MAT_INITIAL_MATRIX, 1.0, transpose_mat, ierr)     
         else
            call MatTransposeMatMult(output_mat, output_mat, &
                        MAT_INITIAL_MATRIX, 1.0, transpose_mat, ierr)          
         endif     

         ! Also have to add in the original distance 1 connections to the square
         ! as the dist 1 strength matrix has had the diagonals removed, so the square won't 
         ! have the dist 1 connetions in it
         call MatAXPY(transpose_mat, 1.0, output_mat, DIFFERENT_NONZERO_PATTERN, ierr)     
         call MatDestroy(output_mat, ierr)

         ! Can end up with diagonal entries we have to remove
         ! Let's get the diagonals that are zero or unassigned
         call MatFindZeroDiagonals(transpose_mat, zero_diags, ierr)
         call ISGetIndicesF90(zero_diags, zero_diags_pointer, ierr)
         ! Then let's just set every other row to have a zero diagonal
         ! as we know they're already preallocated
         counter = 1
         do ifree = 1, local_rows

            if (counter .le. size(zero_diags_pointer)) then
               ! Skip over any rows that don't have diagonals or are already zero
               if (zero_diags_pointer(counter) - global_row_start + 1 == ifree) then
                  counter = counter + 1 
                  cycle
               end if
            end if
         
            ! Set the diagonal to 0
            call MatSetValue(transpose_mat, ifree - 1 + global_row_start, ifree - 1 + global_row_start, 0.0, INSERT_VALUES, ierr)
         end do
         
         call MatAssemblyBegin(transpose_mat, MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(transpose_mat, MAT_FINAL_ASSEMBLY, ierr)

         ! Could call MatEliminateZeros in later versions of petsc, but for here
         ! given we know the entries are ==1, we will just create a copy with "small" stuff removed
         ! ie the zero diagonal
         call remove_small_from_sparse(transpose_mat, 1e-100, output_mat, allow_drop_diagonal = .TRUE.) 
         call MatDestroy(transpose_mat, ierr)

      end if   
      
      ! Reset the entries in the strength matrix back to 1
      if (symmetrize .OR. square) then
         max_nnzs = 0
         do ifree = global_row_start, global_row_end_plus_one-1     
            call MatGetRow(output_mat, ifree, ncols, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr)  
            if (ncols > max_nnzs) max_nnzs = ncols
            call MatRestoreRow(output_mat, ifree, ncols, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr)               
         end do
         allocate(vals(max_nnzs))
         vals = 1.0
         ! Set all the values in the matrix to one
         do ifree = global_row_start, global_row_end_plus_one-1     
            ! Set every entry in this row
            call MatSetValuesRow(output_mat, ifree, vals, ierr)
         end do
         deallocate(vals)
         call MatAssemblyBegin(output_mat, MAT_FINAL_ASSEMBLY, ierr)
         call MatAssemblyEnd(output_mat, MAT_FINAL_ASSEMBLY, ierr)          
      end if

   end subroutine generate_sabs     

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine first_pass_splitting(input_mat, symmetric, strong_threshold, max_luby_steps, cf_splitting_type, cf_markers_local)

      ! Compute a symmetrized strength matrix and then call the first pass of CF splitting

      ! ~~~~~~
      type(tMat), target, intent(in)      :: input_mat
      logical, intent(in)                 :: symmetric
      real, intent(in)                    :: strong_threshold
      integer, intent(in)                 :: max_luby_steps, cf_splitting_type
      integer, dimension(:), allocatable, intent(inout) :: cf_markers_local

      ! Local
      PetscInt :: local_c_size, global_row_start, global_row_end_plus_one, i_loc, counter
      PetscInt :: global_row_start_dist2, global_row_end_plus_one_dist2
      integer :: errorcode
      PetscErrorCode :: ierr
      type(tMat) :: strength_mat, strength_mat_fc, prolongators, temp_mat, strength_mat_c
      type(tIS) :: is_fine, is_coarse, zero_diags
      integer, dimension(:), allocatable :: cf_markers_local_c
      PetscInt, dimension(:), allocatable :: aggregates, nnzs_row
      PetscInt, dimension(:), pointer :: is_pointer_coarse, is_pointer_fine, zero_diags_pointer
      type(tVec) :: diag_vec
      PetscInt, parameter :: nz_ignore = -1

      ! ~~~~~~  

      ! Generate the strength matrix - symmetrize
      ! Use classical strength of connection
      if (cf_splitting_type == CF_PMIS_DIST2) then
         ! Square the strength of connection graph - S'S + S
         ! Note we are symmetrizing the strength matrix here
         call generate_sabs(input_mat, strong_threshold, .TRUE., .TRUE., strength_mat)
      else if (cf_splitting_type == CF_PMIS) then
         ! Note we are symmetrizing the strength matrix here
         call generate_sabs(input_mat, strong_threshold, .TRUE., .FALSE., strength_mat)
      ! pmisr ddc and aggregation
      else 
         call generate_sabs(input_mat, strong_threshold, .NOT. symmetric, .FALSE., strength_mat)         
      end if

      ! Generate an independent set
      if (cf_splitting_type == CF_PMISR_DDC) then
         ! Do PMISR
         call pmisr(strength_mat, max_luby_steps, .FALSE., cf_markers_local)

      else if (cf_splitting_type == CF_PMIS) then
         ! Do distance 1 PMIS
         call pmisr(strength_mat, max_luby_steps, .TRUE., cf_markers_local)

      else if (cf_splitting_type == CF_PMIS_DIST2) then

         ! As we have generated S'S + S for the strength matrix, this will do distance 2 PMIS
         ! Need to do as many Luby steps as required or terrible convergence, as missing even a single
         ! C point with dist 2 can leave huge patches of F points together
         call pmisr(strength_mat, -1, .TRUE., cf_markers_local)

         ! ~~~~~~~~~~
         ! The code below instead relies on the distance 1 strength matrix and does
         ! MIS(MIS(S_dist1)), which coarsens quite a bit faster than MIS(S'S)
         ! ~~~~~~~~~~

         ! ! We need IS's for the C and F points
         ! call create_cf_is(strength_mat, cf_markers_local, is_fine, is_coarse) 

         ! ! Get the number of local C points
         ! call ISGetLocalSize(is_coarse, local_c_size, ierr)

         ! ! Extract S_fc
         ! call MatCreateSubMatrix(strength_mat, &
         !          is_fine, is_coarse, MAT_INITIAL_MATRIX, &
         !          strength_mat_fc, ierr)   
                  
         ! call MatGetOwnershipRange(strength_mat, global_row_start, global_row_end_plus_one, ierr)                        
         ! ! Create [S_fc; 0]'
         ! call compute_P_from_W(strength_mat_fc, global_row_start, &
         !          is_fine, is_coarse, .FALSE., &
         !          prolongators)                     

         ! ! Now the strength matrix for C points can be computed from
         ! ! [Scf; 0] * strength_mat * [Sfc; 0]'
         ! ! which is equivalent to [Scf Sff Sfc], but given Aff will be so big 
         ! ! we don't bother to extract Aff, we just explicitly create [Sfc; 0]' 
         ! ! and then do the PtaP on the full S
         ! ! Also S is symmetric so Scf = Sfc'
         ! call MatPtap(strength_mat, prolongators, &
         !          MAT_INITIAL_MATRIX, 1.58, temp_mat, ierr)   
         
         ! call ISGetIndicesF90(is_fine, is_pointer_fine, ierr)
         ! call ISGetIndicesF90(is_coarse, is_pointer_coarse, ierr)           

         ! ! Can end up with diagonal entries we have to remove
         ! ! Let's get the diagonals that are zero or unassigned
         ! call MatFindZeroDiagonals(temp_mat, zero_diags, ierr)
         ! call ISGetIndicesF90(zero_diags, zero_diags_pointer, ierr)
         ! ! Then let's just set every other row to have a zero diagonal
         ! ! as we know they're already preallocated
         ! counter = 1
         ! call MatGetOwnershipRange(temp_mat, global_row_start_dist2, global_row_end_plus_one_dist2, ierr)                    
         ! do i_loc = 1, size(is_pointer_coarse)

         !    if (counter .le. size(zero_diags_pointer)) then
         !       ! Skip over any rows that don't have diagonals or are already zero
         !       if (zero_diags_pointer(counter) - global_row_start_dist2 + 1 == i_loc) then
         !          counter = counter + 1 
         !          cycle
         !       end if
         !    end if
         
         !    ! Set the diagonal to 0
         !    call MatSetValue(temp_mat, i_loc - 1 + global_row_start_dist2, i_loc - 1 + global_row_start_dist2, 0.0, INSERT_VALUES, ierr)
         ! end do
         
         ! call MatAssemblyBegin(temp_mat, MAT_FINAL_ASSEMBLY, ierr)
         ! call MatAssemblyEnd(temp_mat, MAT_FINAL_ASSEMBLY, ierr)

         ! ! Could call MatEliminateZeros in later versions of petsc, but for here
         ! ! given we know the entries are ==1, we will just create a copy with "small" stuff removed
         ! ! ie the zero diagonal
         ! call remove_small_from_sparse(temp_mat, 1e-100, strength_mat_c, allow_drop_diagonal = .TRUE.)                   
         ! call MatDestroy(temp_mat, ierr)
         ! call ISRestoreIndicesF90(zero_diags, zero_diags_pointer, ierr)
         ! call ISDestroy(zero_diags, ierr)

         ! ! Do the second distance 1 PMIS on just the C points
         ! ! Now any C-point that is not strongly connected to any other C-point we want to retain as C-points
         ! ! It is only C-points that are not distance 2 away from other C-points that we want to turn into F points
         ! call pmisr(strength_mat_c, -1, .TRUE., cf_markers_local_c, zero_measure_c_point=.TRUE.)

         ! ! Then make sure we put the new definitions back into the original cf_markers
         ! cf_markers_local(is_pointer_coarse - global_row_start + 1) = cf_markers_local_c

         ! call ISRestoreIndicesF90(is_coarse, is_pointer_coarse, ierr)
         ! call ISRestoreIndicesF90(is_fine, is_pointer_fine, ierr)         
         ! call ISDestroy(is_fine, ierr)
         ! call ISDestroy(is_coarse, ierr)
         ! deallocate(cf_markers_local_c)
         ! call MatDestroy(strength_mat_fc, ierr)
         ! call MatDestroy(strength_mat_c, ierr)
         ! call MatDestroy(prolongators, ierr)

      else if (cf_splitting_type == CF_AGG) then         

         ! Call an aggregation algorithm
         call generate_aggregation(input_mat, strength_mat, cf_markers_local, aggregates)
         deallocate(aggregates)

      else
         ! Just do the same thing as distance 2 multiple times if you want higher distance
         print *, "Unknown CF splitting algorithm"
         call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, errorcode)
      end if

      ! Destroy the strength matrix
      call MatDestroy(strength_mat, ierr)

   end subroutine first_pass_splitting


! -------------------------------------------------------------------------------------------------------------------------------

   subroutine compute_cf_splitting(input_mat, symmetric, &
                     strong_threshold, max_luby_steps, &
                     cf_splitting_type, fraction_swap, &
                     is_fine, is_coarse)

      ! Computes a CF splitting and returns the F and C point ISs

      ! ~~~~~~
      type(tMat), target, intent(in)      :: input_mat
      logical, intent(in)                 :: symmetric
      real, intent(in)                    :: strong_threshold
      integer, intent(in)                 :: max_luby_steps, cf_splitting_type
      real, intent(in)                    :: fraction_swap
      type(tIS), intent(inout)            :: is_fine, is_coarse

      PetscErrorCode :: ierr
      integer, dimension(:), allocatable :: cf_markers_local

      ! ~~~~~~  

      ! Call the first pass CF splitting with a symmetrized strength matrix
      call first_pass_splitting(input_mat, symmetric, strong_threshold, max_luby_steps, cf_splitting_type, cf_markers_local)

      ! Create the IS for the CF splittings
      call create_cf_is(input_mat, cf_markers_local, is_fine, is_coarse)   
      
      ! Only bother doing the second pass if we haven't requested an exact independent set
      ! ie only do if we don't have diagonal Aff
      if (strong_threshold /= 0.0) then

         ! Do the second pass cleanup - this will directly modify the values in cf_markers_local
         call ddc(input_mat, is_fine, fraction_swap, cf_markers_local)

         ! If we did anything in our ddc second pass
         if (fraction_swap /= 0.0) then
         
            ! These are now outdated
            call ISDestroy(is_fine, ierr)
            call ISDestroy(is_coarse, ierr)

            ! Create the new CF ISs
            call create_cf_is(input_mat, cf_markers_local, is_fine, is_coarse) 
         end if
      end if

      deallocate(cf_markers_local)

   end subroutine compute_cf_splitting 
   
! -------------------------------------------------------------------------------------------------------------------------------

end module cf_splitting

