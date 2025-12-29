module mod_stiffness
    use mod_globals
    use mod_mesh
    
    ! Include PETSc headers
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscmat.h>
    use petscsys
    use petscmat
    
    implicit none

contains

    ! ==========================================
    ! Subroutine: Compute Global Stiffness
    ! ==========================================
    subroutine compute_global_stiffness()
        integer(i4) :: ielem, n1, n2
        
        ! [FIX 1] Define a PetscInt variable for the size (4)
        ! Compilers hate it when we pass the literal '4' to a PetscInt argument
        PetscInt    :: n_size = 4
        PetscInt    :: idx(4) 
        
        real(dp) :: dx, dy, L, c, s, c2, s2, cs
        
        PetscScalar :: ke(4, 4)
        ! [FIX 2] Create a 1D array to pass to PETSc (safer than 2D)
        PetscScalar :: ke_flat(16) 
        PetscScalar :: val_to_add

        if (rank == 0) print *, "Assembling Stiffness Matrix..."

        do ielem = 1, num_elems
            
            ! 1. Get Node IDs
            n1 = connectivity(ielem, 1)
            n2 = connectivity(ielem, 2)

            ! 2. Calculate Geometry
            dx = X(n2, 1) - X(n1, 1)
            dy = X(n2, 2) - X(n1, 2)
            L = sqrt(dx**2 + dy**2)

            if (L <= 0.0_dp) then
                print *, "Error: Zero length element", ielem
                call MPI_Abort(PETSC_COMM_WORLD, 1, ierr)
            end if

            c = dx / L 
            s = dy / L 
            
            c2 = c*c; s2 = s*s; cs = c*s
            
            ! 3. Calculate Local K Matrix
            val_to_add = (E(ielem) * A(ielem)) / L
            
            ke(1,1) = c2;   ke(1,2) = cs;   ke(1,3) = -c2;  ke(1,4) = -cs
            ke(2,1) = cs;   ke(2,2) = s2;   ke(2,3) = -cs;  ke(2,4) = -s2
            ke(3,1) = -c2;  ke(3,2) = -cs;  ke(3,3) = c2;   ke(3,4) = cs
            ke(4,1) = -cs;  ke(4,2) = -s2;  ke(4,3) = cs;   ke(4,4) = s2
            
            ke = ke * val_to_add

            ! [FIX 3] Flatten the 4x4 matrix to 1D array of size 16
            ke_flat = reshape(ke, (/16/))

            ! 4. Determine Global Indices
            idx(1) = (n1 - 1) * 2
            idx(2) = (n1 - 1) * 2 + 1
            idx(3) = (n2 - 1) * 2
            idx(4) = (n2 - 1) * 2 + 1

            ! 5. Insert into PETSc Matrix
            ! Now we pass 'n_size' (PetscInt) instead of '4' (Integer)
            ! And we pass 'ke_flat' (1D array) which always works
            call MatSetValues(K_mat, n_size, idx, n_size, idx, ke_flat, ADD_VALUES, ierr)
            
        end do

        ! 6. Matrix Assembly
        call MatAssemblyBegin(K_mat, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(K_mat, MAT_FINAL_ASSEMBLY, ierr)

        if (rank == 0) print *, "Stiffness Matrix Assembled."

    end subroutine compute_global_stiffness

end module mod_stiffness