program main
    use mod_globals
    use mod_mesh
    use mod_io
    use mod_stiffness
    use mod_solver
    implicit none

    ! Initialize PETSc
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

    if (rank == 0) then
        print *, "=========================================="
        print *, "   PETSc FINITE ELEMENT SOLVER (Level PhD)"
        print *, "   Running on ", size, " processors"
        print *, "=========================================="
    end if

    ! 1. Read Input File (Geometry & Properties)
    ! call read_input_file('Y3-224.txt')
    ! call read_input_file('petsc_input.txt')
    call read_input_file('Y3-224-PETSc.txt')
    
    

    ! 2. Setup PETSc System (Allocate K, F, D)
    call setup_system()

    ! 3. Apply Loads to Vector F
    call apply_loads()

    ! 4. Assemble Stiffness Matrix K
    call compute_global_stiffness()

    ! 5. Solve the System (Apply BCs & KSP Solve)
    call solve_petsc()

    ! 6. View Results
    if (rank == 0) print *, "------------------------------------------"
    if (rank == 0) print *, "FINAL RESULTS (Displacement Vector D):"
    
    ! View the vector to standard output (screen)
    call VecView(D_vec, PETSC_VIEWER_STDOUT_WORLD, ierr)

    if (rank == 0) print *, "------------------------------------------"
    if (rank == 0) print *, "DONE."

    ! Cleanup Global PETSc Objects
    call MatDestroy(K_mat, ierr)
    call VecDestroy(F_vec, ierr)
    call VecDestroy(D_vec, ierr)
    
    call PetscFinalize(ierr)

end program main