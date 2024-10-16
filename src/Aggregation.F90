module aggregation

   use petsc

#include "petsc/finclude/petsc.h"

   implicit none

#include "petsc_legacy.h"

   public   
   
   contains

! -------------------------------------------------------------------------------------------------------------------------------

   subroutine generate_aggregation(input_mat, strength_mat, cf_markers, aggregates)

      ! Do an aggregation algorithm - this is currently an exact copy of the aggregation 
      ! algorithm in PyAMG
      ! Only in serial

      ! ~~~~~~
      type(tMat), target, intent(in)      :: input_mat, strength_mat
      integer, dimension(:), allocatable, intent(inout) :: cf_markers
      PetscInt, dimension(:), allocatable, intent(inout) :: aggregates

      ! Local
      PetscInt :: local_rows, local_cols, global_rows, global_cols
      PetscInt :: a_global_row_start, a_global_row_end_plus_one, ifree, ncols, kfree
      PetscInt :: max_nnzs, ncols_store, aggregate, max_neighbour_index, max_neighbour_value
      integer :: comm_rank, errorcode, seed_size, jfree, comm_size
      PetscErrorCode :: ierr
      MPI_Comm :: MPI_COMM_MATRIX      
      PetscInt, dimension(:), allocatable :: indices, cols, neighbours
      integer, dimension(:), allocatable :: seed
      real, dimension(:), allocatable :: measure, vals
      logical :: mark, mark_neigh

      ! ~~~~~~   

      call PetscObjectGetComm(input_mat, MPI_COMM_MATRIX, ierr)  
      ! Get the comm rank 
      call MPI_Comm_rank(MPI_COMM_MATRIX, comm_rank, errorcode)       
      ! Get the comm size 
      call MPI_Comm_size(MPI_COMM_MATRIX, comm_size, errorcode)      
      ! Serial only!
      if (comm_size /= 1) then
         print *, "Fix me - aggregation algorithm is serial only"
         call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, errorcode)
      end if

      ! Get the local sizes
      call MatGetLocalSize(strength_mat, local_rows, local_cols, ierr)
      call MatGetSize(strength_mat, global_rows, global_cols, ierr)      
      call MatGetOwnershipRange(strength_mat, a_global_row_start, a_global_row_end_plus_one, ierr)   

      ! Get the number of connections in S
      allocate(measure(local_rows))
      allocate(indices(local_rows))
      allocate(cf_markers(local_rows)) 
      cf_markers = 0 
      allocate(aggregates(local_rows))  
      aggregates = 0

      ! Seed the measure between 0 and 1
      call random_seed(size=seed_size)
      allocate(seed(seed_size))
      ! Ensure we seed the same so subsequent runs get the same random
      ! Seed different on each process so we don't have the same random numbers on each processor      
      do jfree = 1, seed_size
         seed(jfree) = comm_rank + 1 + jfree
      end do   
      call random_seed(put=seed) 
      call random_number(measure)

      ! Add the number of connections in S to the randomly seeded measure
      max_nnzs = 0
      do ifree = a_global_row_start, a_global_row_end_plus_one-1                  
         call MatGetRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
         measure(ifree - a_global_row_start + 1) = measure(ifree - a_global_row_start + 1) + real(ncols)
         if (ncols > max_nnzs) max_nnzs = ncols
         call MatRestoreRow(input_mat, ifree, ncols, PETSC_NULL_INTEGER_ARRAY, PETSC_NULL_SCALAR_ARRAY, ierr)
      end do        

      allocate(cols(max_nnzs))    
      allocate(vals(max_nnzs))  
   
      ! Randomise the ordering - not as good as the sort
      do ifree = 1, local_rows
         indices(ifree) = local_rows - (ifree-1)
      end do

      deallocate(measure, seed)
      aggregate = 1

      ! Loop over all the rows - Step 1 - initial covering
      do ifree = local_rows, 1, -1

         ! Get S_i
         call MatGetRow(strength_mat, a_global_row_start + indices(ifree)-1, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)
         ! Store the distance 1 neighbours
         ncols_store = ncols
         allocate(neighbours(ncols_store))
         neighbours = cols(1:ncols)
         call MatRestoreRow(strength_mat, a_global_row_start + indices(ifree)-1, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)

         if (ncols_store == 0) then

            ! No neighbours it stays fine
            cf_markers(indices(ifree)) = -1

         else

            ! print *, "checking node step 1", indices(ifree)
            ! print *, "strong neighbours", neighbours+1
            ! print *, "aggregates neighbours", aggregates(neighbours+1)   
            ! print *, "cf markers neighbours", cf_markers(neighbours+1)         

            ! Check if this node or any of the strongly connected neighbours are already in an aggregate
            mark_neigh = any(cf_markers(neighbours+1) /= 0 .OR. cf_markers(indices(ifree)) /= 0)      
         
            ! Skip if this node or any strong neighbours are already in an aggregate
            if (.NOT. mark_neigh) then
   
               ! This is the root node
               cf_markers(indices(ifree)) = 1
               aggregates(indices(ifree)) = aggregate
               !print *, "this node in aggregate", aggregate
               ! It's neighbours are fine
               do jfree = 1, size(neighbours)
                  cf_markers(neighbours(jfree)+1) = -1
                  aggregates(neighbours(jfree)+1) = aggregate
               end do         
               ! Advance to a new aggregate
               aggregate = aggregate + 1
            end if            

         end if

         deallocate(neighbours)

      end do

      !do ifree = 1, local_rows
      !   print *, ifree, "aggregates after step 1", aggregates(ifree)
      !end do

      ! Step 2 - enlarging the sets
      do ifree = local_rows, 1, -1

         ! Check if this node has been assigned
         if (cf_markers(indices(ifree)) /= 0) cycle

         ! Get S_i
         call MatGetRow(strength_mat, a_global_row_start + indices(ifree)-1, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)
         ! Store the distance 1 neighbours
         ncols_store = ncols
         allocate(neighbours(ncols_store))
         neighbours = cols(1:ncols)
         call MatRestoreRow(strength_mat, a_global_row_start + indices(ifree)-1, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)

         max_neighbour_index = -1
         max_neighbour_value = 0
         ! Could just keep the size of the strong connections in strength_mat and avoid this
         call MatGetRow(input_mat, a_global_row_start + indices(ifree)-1, ncols, cols, vals, ierr)

         !print *, "checking node step 2", indices(ifree)
         !print *, "strong neighbours", neighbours+1
         !print *, "aggregates neighbours", aggregates(neighbours+1)
         !print *, "cols of input mat not strength", cols(1:ncols)

         ! Loop over strongly connected neighbours
         j_loop: do jfree = 1, size(neighbours)

            ! Have we hit a strongly connected neighbour in an aggregate
            ! that wasn't added as part of stage 2
            if (aggregates(neighbours(jfree)+1) .le. 0) cycle
            ! Get the actual matrix entry
            ! Let's find the actual size of this strong connection
            !print *, "checking neighbour", neighbours(jfree)+1
            k_loop: do kfree = 1, ncols
               if (cols(kfree) == neighbours(jfree)) then
                  ! Store which neighbour has the biggest connection
                  !if (abs(vals(kfree)) > max_neighbour_value) then
                  !   max_neighbour_value = abs(vals(kfree))
                     !print *, "found match", jfree
                     max_neighbour_index = jfree
                     exit j_loop
                  !end if
               end if
            end do k_loop
         end do j_loop
         call MatRestoreRow(input_mat, a_global_row_start + indices(ifree)-1, ncols, cols, vals, ierr)

         ! If we didn't have a strongly connected neighbour, cycle
         if (max_neighbour_index /= -1) then

            ! Otherwise we have identified the strongest connection this unassigned node has
            ! And so this node becomes part of that aggregate (ie an F point)
            cf_markers(indices(ifree)) = -1
            ! Make this negative to record that it was added as part of stage 2
            aggregates(indices(ifree)) = -aggregates(neighbours(max_neighbour_index)+1)
            !print *, "setting as", aggregates(neighbours(max_neighbour_index)+1)
         
            ! Advance to a new aggregate
            aggregate = aggregate + 1
         end if

         deallocate(neighbours)
      end do      

      ! Swap the negative ones back to positive
      do ifree = 1, local_rows
         if (aggregates(ifree) < 0) aggregates(ifree) = aggregates(ifree) * (-1.0)
      end do

      ! Step 3 - any remnants
      do ifree = local_rows, 1, -1

         ! Check if this node has been assigned
         if (cf_markers(indices(ifree)) /= 0) cycle

         ! Get S_i
         call MatGetRow(strength_mat, a_global_row_start + indices(ifree)-1, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)
         ! Store the distance 1 neighbours
         ncols_store = ncols
         allocate(neighbours(ncols_store))
         neighbours = cols(1:ncols)
         call MatRestoreRow(strength_mat, a_global_row_start + indices(ifree)-1, ncols, cols, PETSC_NULL_SCALAR_ARRAY, ierr)

         ! This is the root node
         cf_markers(indices(ifree)) = 1
         aggregates(indices(ifree)) = aggregate

         ! Any unassigned strong neighbours are fine and in this aggregate
         do jfree = 1, size(neighbours)
            if (cf_markers(neighbours(jfree)+1) /= 0) cycle            
            cf_markers(neighbours(jfree)+1) = -1
            aggregates(neighbours(jfree)+1) = aggregate
         end do         
         ! Advance to a new aggregate
         aggregate = aggregate + 1

         deallocate(neighbours)

      end do      

      deallocate(indices, cols, vals)

   end subroutine generate_aggregation   

! -------------------------------------------------------------------------------------------------------------------------------

end module aggregation

