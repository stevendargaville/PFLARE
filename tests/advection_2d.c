/*   DMDA/KSP solving a system of linear equations.
     Steady advection equation in 2D with upwinded finite difference
     Gives the same system as PyAMG advection_2d but we don't eliminate the dirichlet dofs
     Modified from ex50.c by Michael Boghosian <boghmic@iit.edu>, 2008,
*/

static char help[] = "Solves 2D advection equation.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>

#include "pflare.h"

extern PetscErrorCode ComputeMat(KSP,Mat,Mat,void*);
extern PetscErrorCode ComputeRHS(KSP,Vec,void*);

typedef struct {
  PetscScalar theta;
} UserContext;

int main(int argc,char **argv)
{
  KSP            ksp;
  DM             da;
  UserContext    user;
  PetscErrorCode ierr;
  PetscInt its;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //       Let's use PFLARE
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

  // Register the pflare types
  PCRegister_PFLARE();

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,11,11,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da);CHKERRQ(ierr);
  ierr = DMSetFromOptions(da);CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);
  ierr = KSPSetDM(ksp,(DM)da);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,&user);CHKERRQ(ierr);

  // Default theta to pi/4 - the advection direction is [cos(theta), sin(theta)]
  PetscReal pi = 4*atan(1.0);
  user.theta = pi/4.0;
  PetscOptionsGetReal(NULL, NULL, "-theta", &user.theta, NULL);
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR >= 17)
      PetscCheck(user.theta <= pi/2.0 && user.theta >= 0.0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Theta must be between 0 and pi/2");
#else
      if (user.theta > pi/2.0 || user.theta < 0.0) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Theta must be between 0 and pi/2");
#endif  

  ierr = KSPSetComputeRHS(ksp,ComputeRHS,&user);CHKERRQ(ierr);
  ierr = KSPSetComputeOperators(ksp,ComputeMat,&user);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,NULL,NULL);CHKERRQ(ierr);

  KSPGetIterationNumber(ksp,&its);
  printf("iterations %d \n", its);

  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

PetscErrorCode ComputeRHS(KSP ksp,Vec b,void *ctx)
{
  UserContext    *user = (UserContext*)ctx;
  PetscErrorCode ierr;
  PetscRandom rctx;

  PetscFunctionBeginUser;

  // Random rhs
  PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
  VecSetRandom(b,rctx);
  PetscRandomDestroy(&rctx);  

  PetscFunctionReturn(0);
}

PetscErrorCode ComputeMat(KSP ksp,Mat J, Mat jac,void *ctx)
{
  UserContext    *user = (UserContext*)ctx;
  PetscErrorCode ierr;
  PetscInt       i, j, M, N, xm, ym, xs, ys, num, numi, numj;
  PetscScalar    v[5], Hx, Hy, HydHx, HxdHy, w1, w2;
  MatStencil     row, col[5];
  DM             da;
  MatNullSpace   nullspace;

  PetscFunctionBeginUser;

  // Coefficients for 2d advection
  w1 = cos(user->theta);
  w2 = sin(user->theta);

  ierr  = KSPGetDM(ksp,&da);CHKERRQ(ierr);
  ierr  = DMDAGetInfo(da,0,&M,&N,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  Hx    = 1.0 / (PetscReal)(M);
  Hy    = 1.0 / (PetscReal)(N);
  HxdHy = Hx/Hy;
  HydHx = Hy/Hx;
  ierr  = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      row.i = i; row.j = j;

      // Boundary values
      if (i==0 || j==0 || i==M-1 || j==N-1) {
         
         // Dirichlets left and bottom
         if (i==0 || j==0) {

            v[0] = 1.0; col[0].i = i;   col[0].j = j;        
            ierr = MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);      

         // Top or right - just gets the normal upwinded stencil
         } else {

            // Upwind advection
            // left
            v[0] = -w1;              col[0].i = i;   col[0].j = j-1;
            // bottom
            v[1] = -w2;              col[1].i = i-1; col[1].j = j;
            // centre
            v[2] = w1 + w2; col[2].i = i;   col[2].j = j;   
            ierr = MatSetValuesStencil(jac,1,&row,3,col,v,INSERT_VALUES);CHKERRQ(ierr);            

         }

      // interior stencil
      } else {

        // Upwind advection
        // left
        v[0] = -w1;              col[0].i = i;   col[0].j = j-1;
        // bottom
        v[1] = -w2;              col[1].i = i-1; col[1].j = j;
        // centre
        v[2] = w1 + w2; col[2].i = i;   col[2].j = j;   
        ierr = MatSetValuesStencil(jac,1,&row,3,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // Eliminate the zeros that the 5 point stencil will have in it 
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=18)
  MatFilter(jac, 0.0, PETSC_TRUE, PETSC_TRUE);
#endif

  PetscFunctionReturn(0);
}