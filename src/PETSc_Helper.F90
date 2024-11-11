module petsc_helper

   use petsc
   use c_petsc_interfaces

#include "petsc/finclude/petsc.h"
                
   implicit none

#include "petsc_legacy.h"

   public

   ! -------------------------------------------------------------------------------------------------------------------------------
   ! -------------------------------------------------------------------------------------------------------------------------------
   ! Some helper functions that involve PETSc matrices, mainly concerning sparsification
   ! -------------------------------------------------------------------------------------------------------------------------------
   ! -------------------------------------------------------------------------------------------------------------------------------      

   contains 
 
   !------------------------------------------------------------------------------------------------------------------------
   
   subroutine remove_small_from_sparse(input_mat, tol, output_mat, relative_max_row_tolerance, lump, allow_drop_diagonal)

      ! Returns a copy of a sparse matrix with entries below abs(val) < tol removed
      ! If relative_max_row_tolerance is true, then the tol is taken to be a relative scaling 
      ! of the max row val on each row
      ! If lumped is true the removed entries are added to the diagonal
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(in) :: input_mat
      type(tMat), intent(inout) :: output_mat
      real, intent(in) :: tol
      logical, intent(in), optional :: relative_max_row_tolerance, lump, allow_drop_diagonal

      PetscInt :: col, ncols, ifree, max_nnzs
      PetscInt :: local_rows, local_cols, global_rows, global_cols, global_row_start, global_row_end_plus_one
      PetscInt :: global_col_start, global_col_end_plus_one, diagonal_index
      PetscErrorCode :: ierr
      integer :: errorcode, comm_size
      PetscInt, dimension(:), allocatable :: nnzs_row, onzs_row, cols, cols_mod
      real, dimension(:), allocatable :: vals, vals_copy
      PetscInt, parameter :: nz_ignore = -1, one=1, zero=0
      logical :: rel_row_tol_logical, lump_entries, drop_diag
      real :: rel_row_tol, lump_sum
      MPI_Comm :: MPI_COMM_MATRIX
      
      ! ~~~~~~~~~~
      ! If the tolerance is 0 we still want to go through this routine and drop the zeros

      call PetscObjectGetComm(input_mat, MPI_COMM_MATRIX, ierr)    
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)

      lump_entries = .FALSE.
      drop_diag = .FALSE.
      if (present(lump)) lump_entries = lump
      if (present(allow_drop_diagonal)) drop_diag = allow_drop_diagonal
      rel_row_tol_logical = .FALSE.
      rel_row_tol = tol
      if (present(relative_max_row_tolerance)) then
         if (relative_max_row_tolerance) then
            rel_row_tol_logical = .TRUE.
            rel_row_tol = 1.0
         end if
      end if

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

         call MatGetRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         if (ncols > max_nnzs) max_nnzs = ncols
         call MatRestoreRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
      end do

      allocate(cols(max_nnzs))
      allocate(cols_mod(max_nnzs))
      allocate(vals(max_nnzs)) 
      allocate(vals_copy(max_nnzs))
      
      if (comm_size/=1) then
         ! Loop over global row indices
         do ifree = global_row_start, global_row_end_plus_one-1                  
            ! Get the row
            call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)
            
            if (rel_row_tol_logical) then
               rel_row_tol = tol * maxval(abs(vals(1:ncols)))
            end if  

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
      else
         ! Loop over global row indices
         do ifree = global_row_start, global_row_end_plus_one-1                  
            ! Get the row
            call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)
   
            if (rel_row_tol_logical) then
               rel_row_tol = tol * maxval(abs(vals(1:ncols)))
            end if  
   
            do col = 1, ncols
               if (abs(vals(col)) .ge. rel_row_tol .OR. (.NOT. drop_diag .AND. cols(col) == ifree)) then
               !if (abs(vals(col)) .ge. rel_row_tol) then
                  ! Convert to local row indices
                  nnzs_row(ifree - global_row_start + 1) = nnzs_row(ifree - global_row_start + 1) + 1
               end if
            end do
   
            ! Must call otherwise petsc leaks memory
            call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)
         end do         
      end if

      ! Create the output matrix
      if (comm_size/=1) then
         call MatCreateAIJ(MPI_COMM_MATRIX, local_rows, local_cols, &
                  global_rows, global_cols, &
                  nz_ignore, nnzs_row, &
                  nz_ignore, onzs_row, &
                  output_mat, ierr)   
      else
         call MatCreateSeqAIJ(PETSC_COMM_SELF, global_rows, global_cols, nz_ignore, nnzs_row, &
                     output_mat, ierr)            
      end if   
       
      ! Just in case there are some zeros in the input mat, ignore them
      call MatSetOption(output_mat, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)     
      ! Don't set any off processor entries so no need for a reduction when assembling
      call MatSetOption(output_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)      
      call MatSetOption(output_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE,  ierr)     
      
      ! Now go and fill the new matrix
      ! Loop over global row indices
      do ifree = global_row_start, global_row_end_plus_one-1                  
      
         ! Get the row
         call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)  
         
         ! Find where the diagonal is
         diagonal_index = -1
         lump_sum = 0
         do col = 1, ncols
            if (cols(col) == ifree) then
               diagonal_index = col
               exit
            end if
         end do            

         if (rel_row_tol_logical) then
            rel_row_tol = tol * maxval(abs(vals(1:ncols)))
         end if 

         cols_mod(1:ncols) = cols(1:ncols)
         vals_copy(1:ncols) = vals(1:ncols)
                  
         ! Allow the dropping of the diagonal
         if (drop_diag) then
            do col = 1, ncols
               if (abs(vals(col)) < rel_row_tol) then
                  cols_mod(col) = -1    
                  if (lump_entries) lump_sum = lump_sum + vals(col)             
               end if
            end do               
         else
            do col = 1, ncols
               if (abs(vals(col)) < rel_row_tol .AND. cols(col) /= ifree) then
                  cols_mod(col) = -1    
                  if (lump_entries) lump_sum = lump_sum + vals(col)             
               end if
            end do
         end if

         ! Add lumped terms to the diagonal
         if (lump_entries) then
            vals_copy(diagonal_index) = vals_copy(diagonal_index) + lump_sum
         end if

         ! Much quicker to call setvalues
         call MatSetValues(output_mat, one, [ifree], ncols, cols_mod, &
                  vals_copy, INSERT_VALUES, ierr)

         ! Must call otherwise petsc leaks memory
         call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)   
      end do           
      
      call MatAssemblyBegin(output_mat, MAT_FINAL_ASSEMBLY, ierr)

      deallocate(cols, vals, nnzs_row, cols_mod, vals_copy)
      if (comm_size/=1) deallocate(onzs_row)

      call MatAssemblyEnd(output_mat, MAT_FINAL_ASSEMBLY, ierr) 
         
   end subroutine remove_small_from_sparse

   !------------------------------------------------------------------------------------------------------------------------
   
   subroutine remove_from_sparse_match_no_lump(input_mat, output_mat)

      ! Returns a copy of a sparse matrix with entries that don't match the sparsity
      ! of the other input matrix dropped
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(in) :: input_mat
      type(tMat), intent(inout) :: output_mat

      PetscInt :: col, ncols, ifree, max_nnzs
      PetscInt :: global_row_start, global_row_end_plus_one
      PetscErrorCode :: ierr
      PetscInt, dimension(:), allocatable :: cols
      real, dimension(:), allocatable :: vals
      PetscInt, parameter :: nz_ignore = -1, one=1, zero=0
      
      ! ~~~~~~~~~~
      ! This returns the global index of the local portion of the matrix
      call MatGetOwnershipRange(input_mat, global_row_start, global_row_end_plus_one, ierr)  

      max_nnzs = 0
      do ifree = global_row_start, global_row_end_plus_one-1                  

         call MatGetRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         if (ncols > max_nnzs) max_nnzs = ncols
         call MatRestoreRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         call MatGetRow(output_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         if (ncols > max_nnzs) max_nnzs = ncols
         call MatRestoreRow(output_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)         
      end do

      allocate(cols(max_nnzs))
      allocate(vals(max_nnzs)) 
       
      ! Just in case there are some zeros in the input mat, ignore them
      call MatSetOption(output_mat, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)     
      ! Don't set any off processor entries so no need for a reduction when assembling
      call MatSetOption(output_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)      
      call MatSetOption(output_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE,  ierr)    
      ! This ensures any entries outside the existing sparsity of output_mat are dropped 
      call MatSetOption(output_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)
      
      ! Now go and fill the new matrix
      ! Loop over global row indices
      do ifree = global_row_start, global_row_end_plus_one-1                  
      
         ! Get the row
         call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)
         call MatSetValues(output_mat, one, [ifree], ncols, cols, &
                  vals, INSERT_VALUES, ierr)
         call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)                    
 
      end do  
            
      call MatAssemblyBegin(output_mat, MAT_FINAL_ASSEMBLY, ierr)
      deallocate(cols, vals)
      call MatAssemblyEnd(output_mat, MAT_FINAL_ASSEMBLY, ierr) 
         
   end subroutine remove_from_sparse_match_no_lump   


   !------------------------------------------------------------------------------------------------------------------------
   
   subroutine remove_from_sparse_match(input_mat, sparsity_mat, output_mat, lump)

      ! Returns a copy of a sparse matrix with entries that don't match the sparsity
      ! of the other input matrix dropped
      ! If lumped is true the removed entries are added to the diagonal
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(in) :: input_mat, sparsity_mat
      type(tMat), intent(inout) :: output_mat
      logical, intent(in), optional :: lump

      PetscInt :: col, ncols, ifree, max_nnzs, ncols_mod, index1, index2
      PetscInt :: local_rows, local_cols, global_rows, global_cols, global_row_start, global_row_end_plus_one
      PetscInt :: global_col_start, global_col_end_plus_one, diagonal_index
      PetscErrorCode :: ierr
      integer :: errorcode, comm_size
      PetscInt, dimension(:), allocatable :: nnzs_row, onzs_row, cols, cols_mod
      real, dimension(:), allocatable :: vals, vals_copy
      PetscInt, parameter :: nz_ignore = -1, one=1, zero=0
      logical :: lump_entries
      real :: lump_sum
      MPI_Comm :: MPI_COMM_MATRIX
      
      ! ~~~~~~~~~~
      ! If the tolerance is 0 we still want to go through this routine and drop the zeros

      call PetscObjectGetComm(input_mat, MPI_COMM_MATRIX, ierr)    
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)

      lump_entries = .FALSE.
      if (present(lump)) lump_entries = lump

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

         call MatGetRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         if (ncols > max_nnzs) max_nnzs = ncols
         call MatRestoreRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         call MatGetRow(sparsity_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         if (ncols > max_nnzs) max_nnzs = ncols
         call MatRestoreRow(sparsity_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)         
      end do

      allocate(cols(max_nnzs))
      allocate(cols_mod(max_nnzs))
      allocate(vals(max_nnzs)) 
      allocate(vals_copy(max_nnzs))

      ! Duplicate the sparsity
      call MatDuplicate(sparsity_mat, MAT_DO_NOT_COPY_VALUES, output_mat, ierr)
       
      ! Just in case there are some zeros in the input mat, ignore them
      call MatSetOption(output_mat, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)     
      ! Don't set any off processor entries so no need for a reduction when assembling
      call MatSetOption(output_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)      
      call MatSetOption(output_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE,  ierr)     
      
      ! Now go and fill the new matrix
      ! Loop over global row indices
      do ifree = global_row_start, global_row_end_plus_one-1                  
      
         ! Get the row
         call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)  
         
         ! Find where the diagonal is
         diagonal_index = -1
         lump_sum = 0
         do col = 1, ncols
            if (cols(col) == ifree) then
               diagonal_index = col
               exit
            end if
         end do            

         ncols_mod = ncols
         cols_mod(1:ncols) = cols(1:ncols)
         vals_copy(1:ncols) = vals(1:ncols)

         ! Must call otherwise petsc leaks memory
         call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)  
         ! Get the sparsity row
         call MatGetRow(sparsity_mat, ifree, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)    
         
         ! Now loop through and do the intersection
         ! Anything that is not in both is not inserted
         ! cols should be a subset of cols_mod
         index1 = 1
         index2 = 1
         do while (index1 .le. ncols_mod .and. index2 .le. ncols) 
            if (cols_mod(index1) == cols(index2)) then
               index1 = index1 + 1
               index2 = index2 + 1
            elseif (cols_mod(index1) < cols(index2)) then
               ! This value won't be inserted into the matrix
               cols_mod(index1) = -1
               if (lump_entries) lump_sum = lump_sum + vals_copy(index1)
               index1 = index1 + 1
            else
               index2 = index2 + 1
            end if  
         end do
         ! Do the rest 
         do col = index1, ncols_mod
            cols_mod(col) = -1
            if (lump_entries) lump_sum = lump_sum + vals_copy(col)
         end do

         ! Must call otherwise petsc leaks memory
         call MatRestoreRow(sparsity_mat, ifree, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)          

         ! Add lumped terms to the diagonal
         if (lump_entries) then
            vals_copy(diagonal_index) = vals_copy(diagonal_index) + lump_sum
         end if

         ! Much quicker to call setvalues
         call MatSetValues(output_mat, one, [ifree], ncols_mod, cols_mod, &
                  vals_copy, INSERT_VALUES, ierr)                  
 
      end do  
            
      call MatAssemblyBegin(output_mat, MAT_FINAL_ASSEMBLY, ierr)

      deallocate(cols, vals, nnzs_row, cols_mod, vals_copy)
      if (comm_size/=1) deallocate(onzs_row)

      call MatAssemblyEnd(output_mat, MAT_FINAL_ASSEMBLY, ierr) 
         
   end subroutine remove_from_sparse_match

  !------------------------------------------------------------------------------------------------------------------------
   
   subroutine mat_duplicate_copy_plus_diag(input_mat, output_mat)

      ! Duplicates and copies the values from input matrix into the output mat, but ensures
      ! there are always diagonal entries present that are set to zero if absent
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(in) :: input_mat
      type(tMat), intent(inout) :: output_mat

      PetscInt :: col, ncols, ifree, max_nnzs
      PetscInt :: local_rows, local_cols, global_rows, global_cols, global_row_start, global_row_end_plus_one
      PetscInt :: global_col_start, global_col_end_plus_one, diagonal_index
      PetscErrorCode :: ierr
      integer :: errorcode, comm_size
      PetscInt, dimension(:), allocatable :: nnzs_row, onzs_row, cols
      real, dimension(:), allocatable :: vals
      PetscInt, parameter :: nz_ignore = -1, one=1, zero=0
      MPI_Comm :: MPI_COMM_MATRIX
      
      ! ~~~~~~~~~~
      ! If the tolerance is 0 we still want to go through this routine and drop the zeros

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

         call MatGetRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         if (ncols > max_nnzs) max_nnzs = ncols
         call MatRestoreRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
      end do

      allocate(cols(max_nnzs))
      allocate(vals(max_nnzs)) 
      
      if (comm_size/=1) then
         ! Loop over global row indices
         do ifree = global_row_start, global_row_end_plus_one-1                  
            ! Get the row
            call MatGetRow(input_mat, ifree, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)

            diagonal_index = -1
            do col = 1, ncols

               if (cols(col) == ifree) then
                  diagonal_index = col
               end if

               ! Have to count the diagonal and off-diagonal nnzs in parallel           
               if (cols(col) .ge. global_col_start .AND. cols(col) .le. global_col_end_plus_one - 1) then
                  ! Convert to local row indices
                  nnzs_row(ifree - global_row_start + 1) = nnzs_row(ifree - global_row_start + 1) + 1
               else
                  ! Convert to local row indices
                  onzs_row(ifree - global_row_start + 1) = onzs_row(ifree - global_row_start + 1) + 1
               end if
            end do

            ! If there was no diagonal, we're going to add one in
            if (diagonal_index == -1) nnzs_row(ifree - global_row_start + 1) = nnzs_row(ifree - global_row_start + 1) + 1

            ! Must call otherwise petsc leaks memory
            call MatRestoreRow(input_mat, ifree, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)
         end do
      else
         ! Loop over global row indices
         do ifree = global_row_start, global_row_end_plus_one-1                  
            ! Get the row
            call MatGetRow(input_mat, ifree, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)
   
            diagonal_index = -1
            do col = 1, ncols

               if (cols(col) == ifree) then
                  diagonal_index = col
               end if

               ! Convert to local row indices
               nnzs_row(ifree - global_row_start + 1) = nnzs_row(ifree - global_row_start + 1) + 1
            end do

            ! If there was no diagonal, we're going to add one in
            if (diagonal_index == -1) nnzs_row(ifree - global_row_start + 1) = nnzs_row(ifree - global_row_start + 1) + 1

            ! Must call otherwise petsc leaks memory
            call MatRestoreRow(input_mat, ifree, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)
         end do         
      end if

      ! Create the output matrix
      if (comm_size/=1) then
         call MatCreateAIJ(MPI_COMM_MATRIX, local_rows, local_cols, &
                  global_rows, global_cols, &
                  nz_ignore, nnzs_row, &
                  nz_ignore, onzs_row, &
                  output_mat, ierr)   
      else
         call MatCreateSeqAIJ(PETSC_COMM_SELF, global_rows, global_cols, nz_ignore, nnzs_row, &
                     output_mat, ierr)            
      end if   
       
      ! Don't set any off processor entries so no need for a reduction when assembling
      call MatSetOption(output_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)      
      call MatSetOption(output_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE,  ierr)     
      
      ! Now go and fill the new matrix
      ! Loop over global row indices
      do ifree = global_row_start, global_row_end_plus_one-1                  
      
         ! Get the row
         call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)  
         
         ! Find where the diagonal is
         diagonal_index = -1
         do col = 1, ncols
            if (cols(col) == ifree) then
               diagonal_index = col
               exit
            end if
         end do            

         ! Much quicker to call setvalues
         call MatSetValues(output_mat, one, [ifree], ncols, cols, &
                  vals, INSERT_VALUES, ierr)

         ! Set the diagonal to zero if there isn't one
         if (diagonal_index == -1) then
            call MatSetValue(output_mat, ifree, ifree, 0.0, INSERT_VALUES, ierr)
         end if

         ! Must call otherwise petsc leaks memory
         call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)   
      end do           
      
      call MatAssemblyBegin(output_mat, MAT_FINAL_ASSEMBLY, ierr)

      deallocate(cols, vals, nnzs_row)
      if (comm_size/=1) deallocate(onzs_row)

      call MatAssemblyEnd(output_mat, MAT_FINAL_ASSEMBLY, ierr) 
         
   end subroutine mat_duplicate_copy_plus_diag   

   !------------------------------------------------------------------------------------------------------------------------
   
   subroutine generate_one_point_from_sparse(input_mat, output_mat)

      ! Returns a copy of a sparse matrix with only the biggest relative entry preserved
      ! We use this to calculate our one point approximate ideal prolongator for example
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(in) :: input_mat
      type(tMat), intent(inout) :: output_mat
      
      PetscInt :: col, ncols, ifree, max_nnzs
      PetscInt :: local_rows, local_cols, global_rows, global_cols, global_row_start, global_row_end_plus_one
      PetscInt :: global_col_start, global_col_end_plus_one
      PetscErrorCode :: ierr
      PetscInt, dimension(:), allocatable :: nnzs_row, onzs_row, cols
      real, dimension(:), allocatable :: vals
      PetscInt, parameter :: nz_ignore = -1, one=1, zero=0
      integer :: max_loc(1)
      integer :: comm_size, errorcode
      MPI_Comm :: MPI_COMM_MATRIX
      
      ! ~~~~~~~~~~

      call PetscObjectGetComm(input_mat, MPI_COMM_MATRIX, ierr)   
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)

      ! Get the local sizes
      call MatGetLocalSize(input_mat, local_rows, local_cols, ierr)
      call MatGetSize(input_mat, global_rows, global_cols, ierr)
      ! This returns the global index of the local portion of the matrix
      call MatGetOwnershipRange(input_mat, global_row_start, global_row_end_plus_one, ierr)  
      call MatGetOwnershipRangeColumn(input_mat, global_col_start, global_col_end_plus_one, ierr)  

      max_nnzs = 0
      do ifree = global_row_start, global_row_end_plus_one-1                  

         call MatGetRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         if (ncols > max_nnzs) max_nnzs = ncols
         call MatRestoreRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
      end do

      allocate(cols(max_nnzs))
      allocate(vals(max_nnzs))       

      ! Can do this quicker in serial 
      if (comm_size/=1) then
      
         allocate(nnzs_row(local_rows))
         allocate(onzs_row(local_rows))
         nnzs_row = 0
         onzs_row = 0      

         ! Loop over global row indices
         do ifree = global_row_start, global_row_end_plus_one-1                  
            ! Get the row
            call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)

            max_loc = maxloc(abs(vals(1:ncols)))

            do col = 1, ncols

               ! Have to count the diagonal and off-diagonal nnzs in parallel           
               if (cols(max_loc(1)) .ge. global_col_start .AND. cols(max_loc(1)) .le. global_col_end_plus_one - 1) then
                  ! Convert to local row indices
                  nnzs_row(ifree - global_row_start + 1) =1
               else
                  ! Convert to local row indices
                  onzs_row(ifree - global_row_start + 1) = 1
               end if
            end do

            ! Must call otherwise petsc leaks memory
            call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)
         end do

         ! Create the output matrix
         call MatCreateAIJ(MPI_COMM_MATRIX, local_rows, local_cols, &
                  global_rows, global_cols, &
                  nz_ignore, nnzs_row, &
                  nz_ignore, onzs_row, &
                  output_mat, ierr)   

      else
         call MatCreateSeqAIJ(PETSC_COMM_SELF, global_rows, global_cols, one, PETSC_NULL_INTEGER_ARRAY, &
                     output_mat, ierr)            
      end if    
      
      ! Just in case there are some zeros in the input mat, ignore them
      call MatSetOption(output_mat, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)   
      ! Don't set any off processor entries so no need for a reduction when assembling
      call MatSetOption(output_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)       
      call MatSetOption(output_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE,  ierr)   
      
      ! Now go and fill the new matrix
      ! Loop over global row indices
      do ifree = global_row_start, global_row_end_plus_one-1                  
      
         ! Get the row
         call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)   
         if (ncols /= 0) then
            max_loc = maxloc(abs(vals(1:ncols)))
            call MatSetValue(output_mat, ifree, cols(max_loc(1)), &
                        vals(max_loc(1)), INSERT_VALUES, ierr)   
         end if         
   
         ! Must call otherwise petsc leaks memory
         call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)   
      end do           

      call MatAssemblyBegin(output_mat, MAT_FINAL_ASSEMBLY, ierr)

      deallocate(cols, vals)
      if (comm_size/=1) deallocate(nnzs_row, onzs_row)

      call MatAssemblyEnd(output_mat, MAT_FINAL_ASSEMBLY, ierr) 
         
   end subroutine generate_one_point_from_sparse    

   !------------------------------------------------------------------------------------------------------------------------
   
   subroutine generate_one_point_with_one_entry_from_sparse(input_mat, output_mat)

      ! Returns a copy of a sparse matrix, but with only one in the spot of the biggest entry
      ! This can be used to generate a classical one point prolongator for example
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(in) :: input_mat
      type(tMat), intent(inout) :: output_mat
      
      PetscInt :: col, row, ncols, ifree, max_nnzs
      PetscInt :: local_rows, local_cols, global_rows, global_cols, global_row_start, global_row_end_plus_one
      PetscInt :: global_col_start, global_col_end_plus_one
      PetscErrorCode :: ierr
      PetscInt, dimension(:), allocatable :: nnzs_row, onzs_row, cols
      real, dimension(:), allocatable :: vals
      PetscInt, parameter :: nz_ignore = -1, one=1, zero=0
      integer :: max_loc(1)
      integer :: comm_size, errorcode
      MPI_Comm :: MPI_COMM_MATRIX
      
      ! ~~~~~~~~~~

      call PetscObjectGetComm(input_mat, MPI_COMM_MATRIX, ierr)  
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)

      ! Get the local sizes
      call MatGetLocalSize(input_mat, local_rows, local_cols, ierr)
      call MatGetSize(input_mat, global_rows, global_cols, ierr)
      ! This returns the global index of the local portion of the matrix
      call MatGetOwnershipRange(input_mat, global_row_start, global_row_end_plus_one, ierr)  
      call MatGetOwnershipRangeColumn(input_mat, global_col_start, global_col_end_plus_one, ierr)  

      max_nnzs = 0
      do ifree = global_row_start, global_row_end_plus_one-1                  

         call MatGetRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         if (ncols > max_nnzs) max_nnzs = ncols
         call MatRestoreRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
      end do

      allocate(cols(max_nnzs))
      allocate(vals(max_nnzs))       

      ! Can do this quicker in serial 
      if (comm_size /= 1) then
      
         allocate(nnzs_row(local_rows))
         allocate(onzs_row(local_rows))
         nnzs_row = 0
         onzs_row = 0      

         ! Loop over global row indices
         do ifree = global_row_start, global_row_end_plus_one-1                  
            ! Get the row
            call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)

            max_loc = maxloc(abs(vals(1:ncols)))

            do col = 1, ncols

               ! Have to count the diagonal and off-diagonal nnzs in parallel           
               if (cols(max_loc(1)) .ge. global_col_start .AND. cols(max_loc(1)) .le. global_col_end_plus_one - 1) then
                  ! Convert to local row indices
                  nnzs_row(ifree - global_row_start + 1) =1
               else
                  ! Convert to local row indices
                  onzs_row(ifree - global_row_start + 1) = 1
               end if
            end do

            ! Must call otherwise petsc leaks memory
            call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)
         end do

         ! Create the output matrix
         call MatCreateAIJ(MPI_COMM_MATRIX, local_rows, local_cols, &
                  global_rows, global_cols, &
                  nz_ignore, nnzs_row, &
                  nz_ignore, onzs_row, &
                  output_mat, ierr)   

      else
         call MatCreateSeqAIJ(PETSC_COMM_SELF, global_rows, global_cols, one, PETSC_NULL_INTEGER_ARRAY, &
                     output_mat, ierr)            
      end if    
      
      ! Just in case there are some zeros in the input mat, ignore them
      call MatSetOption(output_mat, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)   
      ! Don't set any off processor entries so no need for a reduction when assembling
      call MatSetOption(output_mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)        
      call MatSetOption(output_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE,  ierr)    
      
      ! Now go and fill the new matrix
      ! Loop over global row indices
      do ifree = global_row_start, global_row_end_plus_one-1                  
      
         ! Get the row
         call MatGetRow(input_mat, ifree, ncols, cols, vals, ierr)   
         if (ncols /= 0) then
            max_loc = maxloc(abs(vals(1:ncols)))
            call MatSetValue(output_mat, ifree, cols(max_loc(1)), &
                        1.0, INSERT_VALUES, ierr)   
         end if         
   
         ! Must call otherwise petsc leaks memory
         call MatRestoreRow(input_mat, ifree, ncols, cols, vals, ierr)   
      end do           

      call MatAssemblyBegin(output_mat, MAT_FINAL_ASSEMBLY, ierr)

      deallocate(cols, vals)
      if (comm_size /= 1) deallocate(nnzs_row, onzs_row)

      call MatAssemblyEnd(output_mat, MAT_FINAL_ASSEMBLY, ierr) 
         
   end subroutine generate_one_point_with_one_entry_from_sparse       

  !------------------------------------------------------------------------------------------------------------------------
   
   subroutine compute_P_from_W(W, global_row_start, is_fine, is_coarse, identity, reuse, P)

      ! Pass in W and get out P = [W I]' (or [W 0] if identity is false)
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(inout) :: W, P
      type(tIS), intent(in)     :: is_fine, is_coarse
      PetscInt, intent(in)      :: global_row_start
      logical, intent(in) :: identity, reuse

      PetscInt :: global_row_start_W, global_row_end_plus_one_W, global_col_start_W, global_col_end_plus_one_W
      PetscInt :: local_rows_coarse, local_rows, local_cols, local_cols_coarse, max_ncols, i_loc, ncols
      PetscInt :: global_cols, global_rows, global_rows_coarse, global_cols_coarse
      PetscInt :: col, cols_z, rows_z, local_rows_fine
      integer :: errorcode, comm_size, comm_size_world
      PetscErrorCode :: ierr
      MPI_Comm :: MPI_COMM_MATRIX
      PetscInt, dimension(:), allocatable :: nnzs_row, onzs_row, cols
      real, dimension(:), allocatable :: vals
      PetscInt, parameter :: nz_ignore = -1, one=1, zero=0
      PetscInt, dimension(:), pointer :: is_pointer_coarse, is_pointer_fine

      ! ~~~~~~~~~~

      call PetscObjectGetComm(W, MPI_COMM_MATRIX, ierr)    
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)
      call MPI_Comm_size(MPI_COMM_WORLD, comm_size_world, errorcode)      

      call MatGetOwnershipRange(W, global_row_start_W, global_row_end_plus_one_W, ierr)  
      call MatGetOwnershipRangeColumn(W, global_col_start_W, global_col_end_plus_one_W, ierr)                  

      call MatGetSize(W, cols_z, rows_z, ierr) 

      call ISGetIndicesF90(is_fine, is_pointer_fine, ierr)
      call ISGetIndicesF90(is_coarse, is_pointer_coarse, ierr)      

      call IsGetLocalSize(is_coarse, local_rows_coarse, ierr)
      call IsGetLocalSize(is_fine, local_rows_fine, ierr)

      local_cols_coarse = local_rows_coarse
      local_cols = local_rows_coarse + local_rows_fine
      local_rows = local_cols 
      
      global_cols = rows_z + cols_z
      global_rows = global_cols
      global_rows_coarse = rows_z
      global_cols_coarse = rows_z      

      allocate(nnzs_row(local_rows))
      allocate(onzs_row(local_rows))
      nnzs_row = 0
      onzs_row = 0

      ! Get the max number of columns before anything
      max_ncols = -1
      do i_loc = global_row_start_W, global_row_end_plus_one_W-1
         call MatGetRow(W, i_loc, &
                  ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         if (ncols > max_ncols) max_ncols = ncols
         call MatRestoreRow(W, i_loc, &
                  ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)          
      end do
      max_ncols = max_ncols + 1
      allocate(cols(max_ncols))
      allocate(vals(max_ncols))

      ! We may be reusing with the same sparsity, so don't need to preallocate
      if (.NOT. reuse) then
      
         ! Have to get the number of local and off processor nnzs
         do i_loc = global_row_start_W, global_row_end_plus_one_W-1
            call MatGetRow(W, i_loc, &
                     ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)
            do col = 1, ncols
               if (cols(col) .ge. global_col_start_W .AND. cols(col) .le. global_col_end_plus_one_W - 1) then
                  nnzs_row(is_pointer_fine(i_loc - global_row_start_W + 1)-global_row_start +1) = nnzs_row(is_pointer_fine(i_loc - global_row_start_W + 1)-global_row_start +1) + 1
               else
                  onzs_row(is_pointer_fine(i_loc - global_row_start_W + 1)-global_row_start +1) = onzs_row(is_pointer_fine(i_loc - global_row_start_W + 1)-global_row_start +1) + 1
               end if
            end do
            call MatRestoreRow(W, i_loc, &
                     ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)          
         end do

         ! Identity
         if (identity) then
            do i_loc = 1, local_rows_coarse
               nnzs_row(is_pointer_coarse(i_loc)-global_row_start +1) = nnzs_row(is_pointer_coarse(i_loc)-global_row_start +1) + 1
            end do
         end if

         if (comm_size /= 1) then
            call MatCreateAIJ(MPI_COMM_MATRIX, local_rows, local_cols_coarse, global_rows, global_cols_coarse, &
                        nz_ignore, nnzs_row, &
                        nz_ignore, onzs_row, &
                        P, ierr)           
         else
            call MatCreateSeqAIJ(MPI_COMM_MATRIX, local_rows, local_rows_coarse, &
                              nz_ignore, nnzs_row, &
                              P, ierr)          
         end if
      end if

      ! Just in case there are some zeros in the input mat, ignore them
      call MatSetOption(P, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)     
      ! Don't set any off processor entries so no need for a reduction when assembling
      call MatSetOption(P, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)      
      call MatSetOption(P, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE,  ierr)           

      ! Copy in the values of W
      do i_loc = global_row_start_W, global_row_end_plus_one_W-1
         call MatGetRow(W, i_loc, &
                  ncols, cols, vals, ierr)

         call MatSetValues(P,&
                  one, [is_pointer_fine(i_loc - global_row_start_W + 1)], &
                  ncols, cols(1:ncols), &
                  vals(1:ncols), INSERT_VALUES, ierr)                     

         call MatRestoreRow(W, i_loc, &
                  ncols, cols, vals, ierr)          
      end do

      ! If we want the identity block or just leave it zero
      if (identity) then
         do i_loc = 1, local_rows_coarse
            call MatSetValue(P, &
                     is_pointer_coarse(i_loc), &
                     i_loc - 1 + global_col_start_W, &
                     1.0, INSERT_VALUES, ierr)
         end do     
      end if
      deallocate(nnzs_row, onzs_row, cols, vals) 
      
      ! Finish the assembly at the end
      call MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY, ierr)  
      
      call ISRestoreIndicesF90(is_coarse, is_pointer_coarse, ierr)
      call ISRestoreIndicesF90(is_fine, is_pointer_fine, ierr)       
         
   end subroutine compute_P_from_W      


 !------------------------------------------------------------------------------------------------------------------------
   
   subroutine compute_R_from_Z(Z, global_row_start, is_fine, is_coarse, &
                     orig_fine_col_indices, &
                     identity, reuse, &
                     R)

      ! Pass in Z and get out R = [Z I] (or [Z 0] if identity is false)
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(inout) :: Z, R
      PetscInt, intent(in)      :: global_row_start
      type(tIS), intent(in)     :: is_fine, is_coarse
      type(tIS), intent(inout)  :: orig_fine_col_indices
      logical, intent(in) :: identity, reuse

      PetscInt :: global_row_start_Z, global_row_end_plus_one_Z, global_col_start_Z, global_col_end_plus_one_Z
      PetscInt :: local_rows_coarse, local_rows_fine, local_rows, local_cols, local_cols_coarse, max_ncols, i_loc, ncols
      PetscInt :: global_cols, global_rows, global_rows_coarse, global_cols_coarse
      PetscInt :: col, cols_z, rows_z, rows_ao, cols_ao, rows_ad, cols_ad, size_cols
      integer :: comm_size, comm_size_world, errorcode, comm_rank
      PetscErrorCode :: ierr
      MPI_Comm :: MPI_COMM_MATRIX      
      PetscInt, dimension(:), allocatable :: nnzs_row, onzs_row, cols
      real, dimension(:), allocatable :: vals
      PetscOffset :: iicol
      PetscInt :: icol(1)       
      PetscInt, parameter :: nz_ignore = -1, one=1, zero=0
      type(tMat) :: Ad, Ao
      type(c_ptr) :: colmap_c_ptr
      PetscInt, pointer :: colmap_c(:)
      PetscInt, dimension(:), pointer :: col_indices_off_proc_array
      type(tIS) :: col_indices
      PetscInt, dimension(:), pointer :: is_pointer_orig_fine_col, is_pointer_coarse, is_pointer_fine
      integer(c_long_long) :: A_array

      ! ~~~~~~~~~~

      call PetscObjectGetComm(Z, MPI_COMM_MATRIX, ierr)    
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)
      call MPI_Comm_size(MPI_COMM_WORLD, comm_size_world, errorcode)   
      
      call ISGetIndicesF90(is_fine, is_pointer_fine, ierr)
      call ISGetIndicesF90(is_coarse, is_pointer_coarse, ierr)      

      call IsGetLocalSize(is_coarse, local_rows_coarse, ierr)
      call IsGetLocalSize(is_fine, local_rows_fine, ierr)

      local_cols_coarse = local_rows_coarse
      local_cols = local_rows_coarse + local_rows_fine
      local_rows = local_cols

      call MatGetSize(Z, rows_z, cols_z, ierr) 
      
      global_cols = rows_z + cols_z
      global_rows = global_cols
      global_rows_coarse = rows_z
      global_cols_coarse = rows_z      
      call MatGetOwnershipRange(Z, global_row_start_Z, global_row_end_plus_one_Z, ierr)  
      call MatGetOwnershipRangeColumn(Z, global_col_start_Z, global_col_end_plus_one_Z, ierr)   

      ! Get the local non-local components and sizes
      if (comm_size /= 1) then

         call MatMPIAIJGetSeqAIJ(Z, Ad, Ao, icol, iicol, ierr)
         ! We know the col size of Ao is the size of colmap, the number of non-zero offprocessor columns
         call MatGetSize(Ao, rows_ao, cols_ao, ierr)    
         call MatGetSize(Ad, rows_ad, cols_ad, ierr)  
         
      else
         call MatGetSize(Z, rows_ad, cols_ad, ierr)    
         Ad = Z         
      end if

      ! We can reuse the orig_fine_col_indices as they can be expensive to generate in parallel
      if (PetscISIsNull(orig_fine_col_indices)) then
            
         ! Now we need the global off-processor column indices in Z
         if (comm_size /= 1) then 

            A_array = Z%v
            call get_colmap_c(A_array, colmap_c_ptr)
            call c_f_pointer(colmap_c_ptr, colmap_c, shape=[cols_ao])
            
            ! These are the global indices of the columns we want
            allocate(col_indices_off_proc_array(cols_ad + cols_ao))
            size_cols = cols_ad + cols_ao
            ! Local rows (as global indices)
            do i_loc = 1, cols_ad
               col_indices_off_proc_array(i_loc) = global_col_start_Z + i_loc - 1
            end do
            ! Off diagonal rows we want (as global indices)
            do i_loc = 1, cols_ao
               col_indices_off_proc_array(cols_ad + i_loc) = colmap_c(i_loc)
            end do            
            ! These indices are not sorted, we deliberately have the local ones first then the 
            ! off processor ones, as we insert Z into the full matrix in pieces
            ! with Ad, then Ao, so that way we can index into orig_fine_col_indices in that order

         else

            ! These are the global indices of the columns we want
            allocate(col_indices_off_proc_array(cols_ad))
            size_cols = cols_ad
            ! Local rows (as global indices)
            do i_loc = 1, cols_ad
               col_indices_off_proc_array(i_loc) = global_col_start_Z + i_loc - 1
            end do
         end if

         ! Create the IS we want with the cols we want (written as global indices)
         call ISCreateGeneral(MPI_COMM_MATRIX, size_cols, col_indices_off_proc_array, PETSC_USE_POINTER, col_indices, ierr)

         ! Now let's do the comms to get what the original column indices in the full matrix are, given these indices for all 
         ! the columns of Z - ie we need to check in the original fine indices at the positions given by col_indices_off_proc_array
         ! This could be expensive as the number of off-processor columns in Z grows!
         call ISCreateSubIS(is_fine, col_indices, orig_fine_col_indices, ierr)

         ! We've now built the original fine indices
         call ISDestroy(col_indices, ierr)
         deallocate(col_indices_off_proc_array)  

      end if

      ! Get the indices
      call ISGetIndicesF90(orig_fine_col_indices, is_pointer_orig_fine_col, ierr)

      allocate(nnzs_row(local_rows_coarse))
      allocate(onzs_row(local_rows_coarse))
      nnzs_row = 0
      onzs_row = 0
      ! Z
      ! Get the max number of columns before anything
      max_ncols = -1
      do i_loc = global_row_start_Z, global_row_end_plus_one_Z-1
         call MatGetRow(Z, i_loc, &
                  ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         if (ncols > max_ncols) max_ncols = ncols
         call MatRestoreRow(Z, i_loc, &
                  ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)          
      end do
      max_ncols = max_ncols + 1
      allocate(cols(max_ncols))
      allocate(vals(max_ncols))

      ! We may be reusing with the same sparsity, so don't need to preallocate
      if (.NOT. reuse) then

         ! Have to get the number of local and off processor nnzs
         do i_loc = global_row_start_Z, global_row_end_plus_one_Z-1
            call MatGetRow(Z, i_loc, &
                     ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)
            do col = 1, ncols
               if (cols(col) .ge. global_col_start_Z .AND. cols(col) .le. global_col_end_plus_one_Z - 1) then
                  nnzs_row(i_loc - global_row_start_Z + 1) = nnzs_row(i_loc - global_row_start_Z + 1) + 1
               else
                  onzs_row(i_loc - global_row_start_Z + 1) = onzs_row(i_loc - global_row_start_Z + 1) + 1
               end if
            end do
            call MatRestoreRow(Z, i_loc, &
                     ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)          
         end do

         ! Identity
         do i_loc = 1, local_rows_coarse
            nnzs_row(i_loc) = nnzs_row(i_loc) + 1
         end do         

         if (comm_size /= 1) then
            call MatCreateAIJ(MPI_COMM_MATRIX, local_rows_coarse, local_cols, global_rows_coarse, global_cols, &
                        nz_ignore, nnzs_row, &
                        nz_ignore, onzs_row, &
                        R, ierr)           
         else
            call MatCreateSeqAIJ(MPI_COMM_MATRIX, local_rows_coarse, local_cols, &
                                 nz_ignore, nnzs_row, &
                                 R, ierr)        
         end if    
      end if

      ! Just in case there are some zeros in the input mat, ignore them
      call MatSetOption(R, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)     
      ! Don't set any off processor entries so no need for a reduction when assembling
      call MatSetOption(R, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE, ierr)      
      call MatSetOption(R, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE,  ierr) 

      ! Z - do Ad and Ao separately, as that way we have local indices into is_pointer_orig_fine_col
      ! to give us the original column numbers 
      ! Let's start with Ad - remember Ad and Ao are serial
      do i_loc = 1, local_rows_coarse
         call MatGetRow(Ad, i_loc-1, &
                  ncols, cols, vals, ierr)
         call MatSetValues(R,&
                  one, [i_loc -1 + global_row_start_Z], &
                  ncols, is_pointer_orig_fine_col(cols(1:ncols)+1), &
                  vals(1:ncols), INSERT_VALUES, ierr)
         call MatRestoreRow(Ad, i_loc-1, &
                  ncols, cols, vals, ierr)          
      end do
      if (comm_size /= 1) then
         ! Then Ao
         do i_loc = 1, local_rows_coarse
            call MatGetRow(Ao, i_loc-1, &
                     ncols, cols, vals, ierr)
            call MatSetValues(R,&
                     one, [i_loc -1 + global_row_start_Z], &
                     ncols, is_pointer_orig_fine_col(cols(1:ncols)+1 + cols_ad), &
                     vals(1:ncols), INSERT_VALUES, ierr)
            call MatRestoreRow(Ao, i_loc-1, &
                     ncols, cols, vals, ierr)          
         end do
      end if

      ! If we want the identity block or just leave it zero
      if (identity) then
         do i_loc = 1, local_rows_coarse
            call MatSetValue(R, &
                     i_loc - 1 + global_row_start_Z, &
                     is_pointer_coarse(i_loc), &
                     1.0, INSERT_VALUES, ierr)
         end do  
      end if                                   
      deallocate(nnzs_row, onzs_row, cols, vals)

      call MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY, ierr)   

      call ISRestoreIndicesF90(orig_fine_col_indices, is_pointer_orig_fine_col, ierr)
      call ISRestoreIndicesF90(is_coarse, is_pointer_coarse, ierr)
      call ISRestoreIndicesF90(is_fine, is_pointer_fine, ierr)       
         
   end subroutine compute_R_from_Z
   
   !------------------------------------------------------------------------------------------------------------------------
   
   subroutine get_nnzs_petsc_sparse(input_mat, nnzs)
      ! Returns nnzs in a sparse petsc matrix
   
      ! ~~~~~~~~~~
      ! Input 
      type(tMat), intent(in) :: input_mat
      integer(kind=8), intent(out) :: nnzs

      PetscInt :: ncols, i_loc
      PetscInt :: global_row_start, global_row_end_plus_one
      integer :: comm_size, errorcode
      PetscErrorCode :: ierr
      MPI_Comm :: MPI_COMM_MATRIX
      integer(kind=8) :: local_nnzs
      
      ! ~~~~~~~~~~

      call PetscObjectGetComm(input_mat, MPI_COMM_MATRIX, ierr)  
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)

      ! This returns the global index of the local portion of the matrix
      call MatGetOwnershipRange(input_mat, global_row_start, global_row_end_plus_one, ierr)  
      
      local_nnzs = 0

      ! This will be the nnzs associated with the local rows
      do i_loc = global_row_start, global_row_end_plus_one-1                  
         call MatGetRow(input_mat, i_loc, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         local_nnzs = local_nnzs + ncols
         call MatRestoreRow(input_mat, i_loc, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
      end do          

      ! Do an accumulate if in parallel
      if (comm_size/=1) then
         call MPI_Allreduce(local_nnzs, nnzs, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_MATRIX, errorcode)
      else
         nnzs = local_nnzs
      end if    
         
   end subroutine get_nnzs_petsc_sparse

   !-------------------------------------------------------------------------------------------------------------------------------

   subroutine svd(input, U, sigma, VT)

      ! ~~~~~~~~~~~~~~~~

      real, dimension(:, :), intent(in) :: input
      real, dimension(size(input, 1), size(input, 1)), intent(out) :: U
      real, dimension(min(size(input, 1), size(input, 2))), intent(out) :: sigma
      real, dimension(size(input, 2), size(input, 2)), intent(out) :: VT
  
      real, dimension(size(input, 1), size(input, 2)) :: tmp_input
      real, dimension(:), allocatable :: WORK
      integer :: LWORK, M, N, info, errorcode

      ! ~~~~~~~~~~~~~~~~
  
     tmp_input = input
     M = size(input, 1)
     N = size(input, 2)
     LWORK = 2*MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
     allocate(WORK(LWORK))
  
     call DGESVD('A', 'A', M, N, tmp_input, M, &
            sigma, U, M, VT, N, WORK, LWORK, info)
  
      if (info /= 0) then
         print *, "SVD fail"
         call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, errorcode)         
      end if             
     deallocate(WORK)

    end subroutine svd   

   !-------------------------------------------------------------------------------------------------------------------------------

    subroutine pseudo_inv(input, output)

      ! ~~~~~~~~~~~~~~~~

      real, dimension(:, :), intent(in) :: input
      real, dimension(min(size(input, 1), size(input, 2))), intent(out) :: output

      real, dimension(size(input, 1), size(input, 1)) :: U
      real, dimension(min(size(input, 1), size(input, 2))) :: sigma
      real, dimension(size(input, 2), size(input, 2)) :: VT

      integer :: iloc, errorcode

      ! ~~~~~~~~~~~~~~~~

      ! Compute the svd
      call svd(input, U, sigma, VT)

      ! Now the pseudoinverse is V * inv(sigma) * U^T
      ! and sigma is diagonal 
      ! So scale each column of U (given the transpose)
      do iloc = 1, size(input,1)
         if (abs(sigma(iloc)) > 1e-13) then
            U(:, iloc) = U(:, iloc) * 1.0/sigma(iloc)
         else
            U(:, iloc) = 0.0
         end if
      end do

      ! Do the matmatmult, making sure to transpose both
      call dgemm("T", "T", size(input,1), size(input,1), size(input,1), &
               1.0, VT, size(input,1), &
               U, size(input,1), &
               0.0, output, size(input,1))          

      ! nan check
      if (any(output /= output)) then
         print *, "NaN in pseudo inverse"
         call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, errorcode)         
      end if      
  
    end subroutine pseudo_inv   

   !-------------------------------------------------------------------------------------------------------------------------------

end module petsc_helper

