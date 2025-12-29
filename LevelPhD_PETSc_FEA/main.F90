program main
    use mod_globals
    use mod_mesh
    use mod_io
    use mod_stiffness
    use mod_solver
    use mod_post  ! <--- [新增] 引入后处理模块
    implicit none

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

    if (rank == 0) then
        print *, "=========================================="
        print *, "   PETSc FINITE ELEMENT SOLVER (Level PhD)"
        print *, "   Running on ", size, " processors"
        print *, "=========================================="
    end if

    ! 1. Read Input
    call read_input_file('Y3-224-PETSc.txt')

    ! 2. Setup System
    call setup_system()

    ! 3. Apply Loads
    call apply_loads()

    ! 4. Assemble Stiffness
    call compute_global_stiffness()

    ! 5. Solve
    call solve_petsc()

    ! 6. [新增] Compute Stress (Post-processing)
    call compute_and_write_stress()

    ! 7. Cleanup
    if (rank == 0) print *, "DONE."

    call MatDestroy(K_mat, ierr)
    call VecDestroy(F_vec, ierr)
    call VecDestroy(D_vec, ierr)
    
    call PetscFinalize(ierr)

end program main