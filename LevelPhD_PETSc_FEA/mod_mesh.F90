module mod_mesh
    use mod_globals
    ! Include PETSc headers for Vector and Matrix operations
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
    use petscvec
    use petscmat
    implicit none

    ! --- Mesh Arrays (Standard Fortran Arrays) ---
    ! We still keep geometry and properties on all processors to simplify assembly
    real(dp), allocatable :: X(:,:)          ! Nodal coordinates (x, y)
    integer(i4), allocatable :: connectivity(:,:) ! Element connectivity (node1, node2)
    
    ! --- Element Properties ---
    real(dp), allocatable :: A(:)      ! Cross-sectional Area
    real(dp), allocatable :: E(:)      ! Elastic Modulus
    real(dp), allocatable :: alpha(:)  ! Thermal expansion coefficient
    real(dp), allocatable :: heat(:)   ! Temperature change (Delta T)

    ! --- Boundary Conditions ---
    integer(i4), allocatable :: fixdof(:)    ! 0 = Free, 1 = Fixed
    real(dp), allocatable :: Uf(:)           ! Prescribed displacement values
    real(dp), allocatable :: P_load(:)       ! External Nodal Forces (Standard Array)

    ! --- PETSc Objects (The Distributed Data Structures) ---
    Mat :: K_mat        ! Global Stiffness Matrix (Parallel)
    Vec :: F_vec        ! Global Force Vector / RHS (Parallel)
    Vec :: D_vec        ! Global Displacement Vector / Solution (Parallel)
    
    ! --- Auxiliary Variables ---
    integer(i4), allocatable :: freedof(:)   ! Indices of free degrees of freedom

end module mod_mesh