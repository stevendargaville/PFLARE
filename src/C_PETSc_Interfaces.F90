module c_petsc_interfaces

   use iso_c_binding
   use petsc

#include "petsc/finclude/petsc.h"

   implicit none

   public

   ! -------------------------------------------------------------------------------------------------------------------------------
   ! -------------------------------------------------------------------------------------------------------------------------------
   ! Contains interfaces to C functions defined in C_routines.cpp which 
   ! are used to get around a lack of some PETSc fortran interfaces
   ! -------------------------------------------------------------------------------------------------------------------------------
   ! -------------------------------------------------------------------------------------------------------------------------------      

   interface   
       
      subroutine mat_mat_symbolic_c(A_array, B_array, C_array) &
         bind(c, name="mat_mat_symbolic_c")
         use iso_c_binding
         integer(c_long_long) :: A_array
         integer(c_long_long) :: B_array
         integer(c_long_long) :: C_array
      end subroutine mat_mat_symbolic_c         

   end interface        

   interface   
      
      subroutine get_colmap_c(A_array, colmap) &
         bind(c, name="get_colmap_c")
         use iso_c_binding
         integer(c_long_long) :: A_array
         type(c_ptr) :: colmap

      end subroutine get_colmap_c         
 
   end interface  

   interface   
      
      subroutine vecscatter_mat_begin_c(A_array, vec_long, cf_markers_nonlocal) &
         bind(c, name="vecscatter_mat_begin_c")
         use iso_c_binding
         integer(c_long_long) :: A_array, vec_long
         type(c_ptr) :: cf_markers_nonlocal

      end subroutine vecscatter_mat_begin_c         
 
   end interface   
   
   interface   
      
      subroutine vecscatter_mat_end_c(A_array, vec_long, cf_markers_nonlocal) &
         bind(c, name="vecscatter_mat_end_c")
         use iso_c_binding
         integer(c_long_long) :: A_array, vec_long
         type(c_ptr) :: cf_markers_nonlocal

      end subroutine vecscatter_mat_end_c         
 
   end interface     

   interface   
      
      subroutine vecscatter_mat_restore_c(A_array, cf_markers_nonlocal) &
         bind(c, name="vecscatter_mat_restore_c")
         use iso_c_binding
         integer(c_long_long) :: A_array
         type(c_ptr) :: cf_markers_nonlocal

      end subroutine vecscatter_mat_restore_c         
 
   end interface  
   
   interface   
      
      subroutine vecscatter_mat_reverse_begin_c(A_array, vec_long) &
         bind(c, name="vecscatter_mat_reverse_begin_c")
         use iso_c_binding
         integer(c_long_long) :: A_array, vec_long

      end subroutine vecscatter_mat_reverse_begin_c         
 
   end interface   
   
   interface   
      
      subroutine vecscatter_mat_reverse_end_c(A_array, vec_long) &
         bind(c, name="vecscatter_mat_reverse_end_c")
         use iso_c_binding
         integer(c_long_long) :: A_array, vec_long

      end subroutine vecscatter_mat_reverse_end_c         
 
   end interface    
   
   interface   
      
      subroutine MatSeqAIJGetArrayF90_mine(A_array, array) &
            bind(c, name="MatSeqAIJGetArrayF90_mine")
            use iso_c_binding
            integer(c_long_long) :: A_array
            type(c_ptr) :: array

      end subroutine MatSeqAIJGetArrayF90_mine         
 
   end interface     

   interface   
      
      subroutine GenerateIS_ProcAgglomeration_c(proc_stride, global_size, local_size_reduced, start) &
         bind(c, name="GenerateIS_ProcAgglomeration_c")
         use iso_c_binding
         PetscInt, value :: proc_stride
         PetscInt, value :: global_size
         PetscInt :: local_size_reduced, start

      end subroutine GenerateIS_ProcAgglomeration_c         
 
   end interface 
   
   interface   
      
      subroutine MatPartitioning_c(A_array, n_parts, proc_stride, index) &
         bind(c, name="MatPartitioning_c")
         use iso_c_binding
         integer(c_long_long) :: A_array
         PetscInt, value :: n_parts
         PetscInt :: proc_stride
         integer(c_long_long) :: index
      end subroutine MatPartitioning_c         
 
   end interface  
   
   interface   
      
      subroutine MatMPICreateNonemptySubcomm_c(A_array, on_subcomm, B_array) &
         bind(c, name="MatMPICreateNonemptySubcomm_c")
         use iso_c_binding
         integer(c_long_long) :: A_array
         integer(c_int)       :: on_subcomm
         integer(c_long_long) :: B_array
      end subroutine MatMPICreateNonemptySubcomm_c         
 
   end interface
   
   interface   
      
      subroutine c_PCGetStructureFlag(A_array, flag) &
         bind(c, name="c_PCGetStructureFlag")
         use iso_c_binding
         integer(c_long_long) :: A_array
         integer(c_int) :: flag
      end subroutine c_PCGetStructureFlag
 
   end interface
   
   interface   
      
      subroutine PCGetSetupCalled_c(A_array, setupcalled) &
         bind(c, name="PCGetSetupCalled_c")
         use iso_c_binding
         integer(c_long_long) :: A_array
         PetscInt :: setupcalled
      end subroutine PCGetSetupCalled_c         
 
   end interface

   interface   
      
      subroutine generate_identity_is_kokkos(A_array, index, B_array) &
         bind(c, name="generate_identity_is_kokkos")
         use iso_c_binding
         integer(c_long_long) :: A_array
         integer(c_long_long) :: index
         integer(c_long_long) :: B_array
      end subroutine generate_identity_is_kokkos         
 
   end interface     

! -------------------------------------------------------------------------------------------------------------------------------

end module c_petsc_interfaces

