#include <petscvec_kokkos.hpp>
#include <petsc.h>
#include <iostream>
#include <../src/mat/impls/aij/seq/kokkos/aijkok.hpp>
#include <memory>

using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;
using DefaultMemorySpace    = Kokkos::DefaultExecutionSpace::memory_space;
using PetscIntKokkosView    = Kokkos::View<PetscInt *, DefaultMemorySpace>;
using PetscIntConstKokkosViewHost = Kokkos::View<const PetscInt *, Kokkos::HostSpace>;

using ViewPtr = std::shared_ptr<PetscIntKokkosView>;

// Define array of shared pointers representing fine and coarse IS's 
// on each level on the device
ViewPtr* IS_fine_views = nullptr;
ViewPtr* IS_coarse_views = nullptr;
int max_levels = -1;

//------------------------------------------------------------------------------------------------------------------------

// Destroys the data
PETSC_INTERN void destroy_VecISCopyLocal_kokkos()
{
   if (IS_fine_views) {
      // Will automatically call the destructor on each element
      delete[] IS_fine_views;
      IS_fine_views = nullptr;
   }
   if (IS_coarse_views) {
      delete[] IS_coarse_views;
      IS_coarse_views = nullptr;
   }   

    return;
}

//------------------------------------------------------------------------------------------------------------------------

// Creates the data we need to do the equivalent of veciscopy on local data in kokkos
PETSC_INTERN void create_VecISCopyLocal_kokkos(int max_levels_input)
{
   // If not built
   if (!IS_fine_views)
   {
      // Allocate array of pointers
      max_levels = max_levels_input;

      // Initialise fine
      IS_fine_views = new ViewPtr[max_levels];
      // Initialize each element as null until it's set
      // we don't want to accidently call the constructor on any of the views
      for (int i = 0; i < max_levels; i++) {
         IS_fine_views[i] = nullptr;
      }
      // Initialise coarse
      IS_coarse_views = new ViewPtr[max_levels];
      for (int i = 0; i < max_levels; i++) {
         IS_coarse_views[i] = nullptr;
      }      
   }
   // Built but different max size, destroy and rebuild
   else if (max_levels_input != max_levels)
   {
      destroy_VecISCopyLocal_kokkos();
      create_VecISCopyLocal_kokkos(max_levels_input);
   }

   return;
}

//------------------------------------------------------------------------------------------------------------------------

// Copy the input IS's to the device for our_level
PETSC_INTERN void set_VecISCopyLocal_kokkos_our_level(int our_level, IS *index_fine, IS *index_coarse)
{
   // Get the sizes of the local component of the input IS's
   PetscInt fine_local_size, coarse_local_size;
   ISGetLocalSize(*index_fine, &fine_local_size);
   ISGetLocalSize(*index_coarse, &coarse_local_size);

   // Get pointers to the indices on the host
   const PetscInt *fine_indices_ptr, *coarse_indices_ptr;
   ISGetIndices(*index_fine, &fine_indices_ptr);

   // Create a host view of the existing indices
   auto fine_view_h = PetscIntConstKokkosViewHost(fine_indices_ptr, fine_local_size);
   // Create a device view
   IS_fine_views[our_level] = std::make_shared<PetscIntKokkosView>("IS_fine_view_" + std::to_string(our_level), fine_local_size);
   // Copy the indices over to the device
   Kokkos::deep_copy(*IS_fine_views[our_level], fine_view_h);
   ISRestoreIndices(*index_fine, &fine_indices_ptr);

   ISGetIndices(*index_coarse, &coarse_indices_ptr);
   auto coarse_view_h = PetscIntConstKokkosViewHost(coarse_indices_ptr, coarse_local_size);
   // Create a device view
   IS_coarse_views[our_level] = std::make_shared<PetscIntKokkosView>("IS_coarse_view_" + std::to_string(our_level), coarse_local_size);
   // Copy the indices over to the device
   Kokkos::deep_copy(*IS_coarse_views[our_level], coarse_view_h);  
   ISRestoreIndices(*index_coarse, &coarse_indices_ptr); 

   return;
}

//------------------------------------------------------------------------------------------------------------------------

// Do the equivalent of veciscopy on local data using the IS data on the device
PETSC_INTERN void VecISCopyLocal_kokkos(int our_level, int fine_int, Vec *vfull, int mode_int, Vec *vreduced)
{

   //std::cout << "inside our one" << std::endl;
   PetscInt vfull_lock_state, vreduced_lock_state;
   PetscInt vfull_num_locks = 0, vreduced_num_locks = 0;

   // If vfull_lock_state is greater than 0 vfull has been locked for reading!
   // This happens on the top level as we are fed in ksp->vec_rhs
   // We need to turn that off to get access to the pointers
   // we will push a lock on at the end   
   VecLockGet(*vfull, &vfull_lock_state);
   if (vfull_lock_state > 0)
   {
      do 
      {
         vfull_num_locks++;
         VecLockReadPop(*vfull);
         VecLockGet(*vfull, &vfull_lock_state);
      }
      while (vfull_lock_state > 0);
   }
   //std::cout << "vfull_num_locks " << vfull_num_locks << std::endl;

   VecLockGet(*vreduced, &vreduced_lock_state);
   if (vreduced_lock_state > 0)
   {
      do 
      {
         vreduced_num_locks++;
         VecLockReadPop(*vreduced);
         VecLockGet(*vreduced, &vreduced_lock_state);
      }
      while (vreduced_lock_state > 0);
   }
   //std::cout << "vreduced_num_locks " << vreduced_num_locks << std::endl;   

   // Get the device pointers
   PetscScalar *vfull_d, *vreduced_d;
   PetscMemType mtype;
   // @@@ the restore we call after says we modified device memory
   // and that isn't true for both pointers - lets modify to access the pointers directly
   // so we can say when we have modified
   // I think I could get rid of the locks then too because they are only checked in 
   // VecGetArrayAndMemType
   VecGetArrayAndMemType(*vfull, &vfull_d, &mtype);
   VecGetArrayAndMemType(*vreduced, &vreduced_d, &mtype);

   // @@@ think i have to make a shallow copy of the individual view in is_fine_views
   // that is in scope here as kokkos on ese-peak is not compiling due to 
   // error: identifier "IS_fine_views" is undefined in device code
   // or it might actually be the shared pointers!

   PetscIntKokkosView is_d;
   if (fine_int)
   {
      is_d = *IS_fine_views[our_level];
   }
   else
   {
      is_d = *IS_coarse_views[our_level];
   } 

   // SCATTER_REVERSE=1
   // vreduced[i] = vfull[is[i]]
   if (mode_int == 1)
   {
      Kokkos::parallel_for(
         Kokkos::RangePolicy<>(0, is_d.extent(0)), KOKKOS_LAMBDA(int i) {           
            vreduced_d[i] = vfull_d[is_d(i)];
      });
   }        
   // SCATTER_FORWARD=0
   // vfull[is[i]] = vreduced[i]
   else if (mode_int == 0)
   {
      Kokkos::parallel_for(
         Kokkos::RangePolicy<>(0, is_d.extent(0)), KOKKOS_LAMBDA(int i) {           
            vfull_d[is_d(i)] = vreduced_d[i];
      });         
   }

   VecRestoreArrayAndMemType(*vfull, &vfull_d);
   VecRestoreArrayAndMemType(*vreduced, &vreduced_d);   

   // Push as many locks back on as we found when we started
   for (int i = 0; i < vfull_num_locks; i++) VecLockReadPush(*vfull);
   for (int i = 0; i < vreduced_num_locks; i++) VecLockReadPush(*vreduced);

   return;
}

//------------------------------------------------------------------------------------------------------------------------