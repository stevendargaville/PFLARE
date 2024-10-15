module matshell_data_type

   use petsc
   use air_data_type

   implicit none

   public

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
   ! This is our context type for matshells
   ! This has to be in a separate file to matshell
   
   type :: mat_ctxtype

      integer :: our_level = -1
      real, dimension(:), pointer :: coefficients => null()
      logical                     :: own_coefficients = .FALSE.
      real, dimension(:), pointer :: real_roots => null()
      real, dimension(:), pointer :: imag_roots => null()
      type(tMat) :: mat
      type(tVec) :: vec
      type(air_multigrid_data), pointer :: air_data => null()

   end type mat_ctxtype
   
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
   
   contains

end module matshell_data_type

