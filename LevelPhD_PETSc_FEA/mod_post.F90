module mod_post
    use mod_globals
    use mod_mesh
    
#include <petsc/finclude/petscvec.h>
    use petscvec
    
    implicit none

contains

    ! ==========================================
    ! Subroutine: Compute and Write Stress
    ! Description: Gathers displacements to Rank 0, calculates stress, writes to file.
    ! ==========================================
    subroutine compute_and_write_stress()
        Vec :: D_seq        ! A sequential vector to hold ALL data on Rank 0
        VecScatter :: ctx   ! The "Messenger" object
        PetscScalar, pointer :: sol_array(:) ! Pointer to access raw numbers
        
        integer(i4) :: ielem, n1, n2
        real(dp) :: dx, dy, L, E_val
        real(dp) :: u1, v1, u2, v2
        real(dp) :: extension, strain, stress, force
        integer(i4) :: idx_n1_x, idx_n1_y, idx_n2_x, idx_n2_y

        ! --- 1. Gather Data to Rank 0 ---
        ! We use VecScatterCreateToZero to send everything to Rank 0
        if (rank == 0) print *, "Post-processing: Calculating Stresses..."
        
        call VecScatterCreateToZero(D_vec, ctx, D_seq, ierr)
        call VecScatterBegin(ctx, D_vec, D_seq, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(ctx, D_vec, D_seq, INSERT_VALUES, SCATTER_FORWARD, ierr)

        ! --- 2. Rank 0 Calculations ---
        if (rank == 0) then
            
            ! Get access to the raw Fortran array inside the PETSc Vector
            call VecGetArrayRead(D_seq, sol_array, ierr)

            ! [FIX] Changed 'newunit=10' to 'unit=10'
            open(unit=10, file='results_stress.txt', status='replace')
            
            write(10, *) "ElementID,  Stress(MPa),  Force(kN),  Extension(mm)"
            print *, "---------------------------------------------------------"
            print *, " Elem |   Stress (MPa)  |    Force (kN)   | Extension (mm)"
            print *, "---------------------------------------------------------"

            do ielem = 1, num_elems
                n1 = connectivity(ielem, 1)
                n2 = connectivity(ielem, 2)

                ! Get original length
                dx = X(n2, 1) - X(n1, 1)
                dy = X(n2, 2) - X(n1, 2)
                L = sqrt(dx**2 + dy**2)

                ! Get Displacements from the gathered array
                ! Node n1 indices (1-based for Fortran array)
                idx_n1_x = (n1 - 1) * 2 + 1
                idx_n1_y = (n1 - 1) * 2 + 2
                ! Node n2 indices
                idx_n2_x = (n2 - 1) * 2 + 1
                idx_n2_y = (n2 - 1) * 2 + 2

                u1 = real(sol_array(idx_n1_x))
                v1 = real(sol_array(idx_n1_y))
                u2 = real(sol_array(idx_n2_x))
                v2 = real(sol_array(idx_n2_y))

                ! Calculate Extension
                extension = ((u2 - u1) * dx + (v2 - v1) * dy) / L

                ! Strain & Stress
                strain = extension / L
                E_val = E(ielem)
                stress = E_val * strain
                force = stress * A(ielem)

                ! Output to Screen and File
                write(10, '(I5, 3(1X, F14.4))') ielem, stress, force, extension
                print '(I5, " |", F14.4, " |", F14.4, " |", F14.4)', ielem, stress, force, extension

            end do

            close(10)
            print *, "---------------------------------------------------------"
            print *, "Results saved to 'results_stress.txt'"

            ! Clean up
            call VecRestoreArrayRead(D_seq, sol_array, ierr)
            call VecDestroy(D_seq, ierr)
        end if

        call VecScatterDestroy(ctx, ierr)

    end subroutine compute_and_write_stress

end module mod_post