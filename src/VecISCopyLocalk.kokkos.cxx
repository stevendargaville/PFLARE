#include <petscvec_kokkos.hpp>
#include <petsc.h>
#include <iostream>
#include <../src/mat/impls/aij/seq/kokkos/aijkok.hpp>
#include <Kokkos_DualView.hpp>
#include <../src/vec/vec/impls/seq/kokkos/veckokkosimpl.hpp>

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

 // Copied from veckok_kokkos.cxx
 // Sync a Kokkos::DualView to  in execution space 
 // If  is HostSpace, fence the exec so that the data on host is immediately available.
 template <class MemorySpace, typename Type>
 static PetscErrorCode KokkosDualViewSync_mine(Kokkos::DualView<Type *> &v_dual, const Kokkos::DefaultExecutionSpace &exec)
 {
   size_t bytes = v_dual.extent(0) * sizeof(Type);

   PetscFunctionBegin;
   PetscCall(PetscLogGpuTimeBegin());
   if (std::is_same_v<MemorySpace, HostMirrorMemorySpace>) {
     if (v_dual.need_sync_host()) {
       PetscCallCXX(v_dual.sync_host(exec));
       PetscCall(PetscLogGpuToCpu(bytes));
     }
     // even if v_d and v_h share the same memory (as on AMD MI300A) and thus we don't need to sync_host,
     // we still need to fence the execution space as v_d might being populated by some async kernel,
     // and we need to finish it.
     PetscCallCXX(exec.fence());
   } else {
     if (v_dual.need_sync_device()) {
       PetscCallCXX(v_dual.sync_device(exec));
       PetscCall(PetscLogCpuToGpu(bytes));
     }
   }
   PetscCall(PetscLogGpuTimeEnd());
   PetscFunctionReturn(PETSC_SUCCESS);
 }

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
   // Get the device pointers
   PetscScalar *vfull_d, *vreduced_d;

   // ~~~~~~~~~~
   // This is just copied from VecGetArrayAndMemType_SeqKokkos, the reason we don't call 
   // the normal VecGetArrayAndMemType is that it checks if there is a read lock 
   // first and for some reason petsc passes in ksp->vec_rhs in the mg which has a read lock on it
   // ~~~~~~~~~~
   Vec_Kokkos *veckok_full = static_cast<Vec_Kokkos *>((*vfull)->spptr);
   Vec_Kokkos *veckok_reduced = static_cast<Vec_Kokkos *>((*vreduced)->spptr);

   /* Always return up-to-date in the default memory space */
#if (PETSC_VERSION_LT(3,22,2))
   veckok_full->v_dual.sync_device(PetscGetKokkosExecutionSpace());
   veckok_reduced->v_dual.sync_device(PetscGetKokkosExecutionSpace());
#else
   KokkosDualViewSync_mine<DefaultMemorySpace>(veckok_full->v_dual, PetscGetKokkosExecutionSpace());
   KokkosDualViewSync_mine<DefaultMemorySpace>(veckok_reduced->v_dual, PetscGetKokkosExecutionSpace());
#endif
   // ~~~~~~~~~~

   vfull_d = veckok_full->v_dual.view_device().data();
   vreduced_d = veckok_reduced->v_dual.view_device().data();   

   // Get the start range in parallel
   PetscInt global_row_start, global_row_end_plus_one;
   VecGetOwnershipRange(*vfull, &global_row_start, &global_row_end_plus_one);

   // Can't use the shared pointer directly within the parallel 
   // regions on the device
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
            vreduced_d[i] = vfull_d[is_d(i) - global_row_start];
      });
   }        
   // SCATTER_FORWARD=0
   // vfull[is[i]] = vreduced[i]
   else if (mode_int == 0)
   {
      Kokkos::parallel_for(
         Kokkos::RangePolicy<>(0, is_d.extent(0)), KOKKOS_LAMBDA(int i) {           
            vfull_d[is_d(i) - global_row_start] = vreduced_d[i];
      });         
   }

   // Restore doesn't actually do anything except call modify_host or modify_device
   // Only want to do that for the vectors we know we've modified
   // Don't have to have called the get before doing this
   if (mode_int == 0) VecRestoreArrayAndMemType(*vfull, &vfull_d);
   if (mode_int == 1) VecRestoreArrayAndMemType(*vreduced, &vreduced_d);   

   return;
}

//------------------------------------------------------------------------------------------------------------------------