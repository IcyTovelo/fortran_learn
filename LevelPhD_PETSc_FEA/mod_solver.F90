module mod_solver
    use mod_globals
    use mod_mesh
    
    ! Include KSP header
#include <petsc/finclude/petscksp.h>
    use petscksp
    
    implicit none

contains

    ! ==========================================
    ! Subroutine: Setup System
    ! ==========================================
    subroutine setup_system()
        PetscInt :: n_total

        n_total = num_nodes * ndof

        if (rank == 0) print *, "Initializing PETSc Objects (Size:", n_total, ")..."

        call VecCreate(PETSC_COMM_WORLD, D_vec, ierr)
        call VecSetSizes(D_vec, PETSC_DECIDE, n_total, ierr)
        call VecSetFromOptions(D_vec, ierr)
        
        call VecDuplicate(D_vec, F_vec, ierr)

        call MatCreate(PETSC_COMM_WORLD, K_mat, ierr)
        call MatSetSizes(K_mat, PETSC_DECIDE, PETSC_DECIDE, n_total, n_total, ierr)
        call MatSetType(K_mat, MATMPIAIJ, ierr)
        call MatSetUp(K_mat, ierr)

    end subroutine setup_system

    ! ==========================================
    ! Subroutine: Apply Loads
    ! ==========================================
    subroutine apply_loads()
        integer(i4) :: i
        ! [FIX] Use Arrays for VecSetValues (even for single value)
        PetscInt    :: idx_arr(1)
        PetscScalar :: val_arr(1)
        PetscInt    :: one_idx = 1 

        if (rank == 0) print *, "Applying Loads..."

        do i = 1, num_nodes * ndof
            ! Convert to 0-based index and store in array
            idx_arr(1) = i - 1 
            val_arr(1) = P_load(i)
            
            if (abs(val_arr(1)) > 1.0e-9_dp) then
                 ! [FIX] Pass arrays (idx_arr, val_arr) instead of scalars
                 call VecSetValues(F_vec, one_idx, idx_arr, val_arr, ADD_VALUES, ierr)
            end if
        end do
        
        call VecAssemblyBegin(F_vec, ierr)
        call VecAssemblyEnd(F_vec, ierr)

    end subroutine apply_loads

    ! ==========================================
    ! Subroutine: Solve PETSc
    ! ==========================================
    subroutine solve_petsc()
        KSP :: ksp
        PC  :: pc
        IS  :: is_fixed 
        integer(i4) :: i, count
        PetscInt, allocatable :: fixed_indices(:)
        
        PetscScalar :: diag_one = 1.0

        if (rank == 0) print *, "Preparing Solver & Applying BCs..."

        ! Count fixed DOFs
        count = 0
        do i = 1, num_nodes * ndof
            if (fixdof(i) == 1) count = count + 1
        end do
        
        allocate(fixed_indices(count))
        
        count = 1
        do i = 1, num_nodes * ndof
            if (fixdof(i) == 1) then
                fixed_indices(count) = i - 1 
                count = count + 1
            end if
        end do

        call ISCreateGeneral(PETSC_COMM_WORLD, count-1, fixed_indices, PETSC_COPY_VALUES, is_fixed, ierr)

        ! [FIX] Use MatZeroRowsColumnsIS when passing an IS object
        call MatZeroRowsColumnsIS(K_mat, is_fixed, diag_one, D_vec, F_vec, ierr)

        ! Setup Solver
        call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
        call KSPSetOperators(ksp, K_mat, K_mat, ierr)
        
        call KSPGetPC(ksp, pc, ierr)
        call PCSetType(pc, PCJACOBI, ierr)
        
        call KSPSetFromOptions(ksp, ierr)

        if (rank == 0) print *, "Solving Linear System..."
        call KSPSolve(ksp, F_vec, D_vec, ierr)

        if (rank == 0) print *, "Solver Converged!"

        call KSPDestroy(ksp, ierr)
        call ISDestroy(is_fixed, ierr)

    end subroutine solve_petsc

end module mod_solver