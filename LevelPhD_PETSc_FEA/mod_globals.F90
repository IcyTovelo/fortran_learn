module mod_globals
    ! Include PETSc header for global definitions and MPI types
#include <petsc/finclude/petscsys.h>
    use petscsys
    implicit none

    ! --- Standard Fortran Kinds ---
    integer, parameter :: i4 = selected_int_kind(9)
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! --- PETSc MPI Variables ---
    PetscErrorCode :: ierr
    PetscMPIInt    :: rank, size  ! Rank of current processor and total number of processors

    ! --- FEA Global Variables ---
    integer(i4) :: num_elems      ! Total number of elements
    integer(i4) :: num_nodes      ! Total number of nodes
    integer(i4) :: num_mat_props  ! Number of material properties
    integer(i4) :: ndof           ! Degrees of freedom per node (usually 2 for 2D truss)

end module mod_globals