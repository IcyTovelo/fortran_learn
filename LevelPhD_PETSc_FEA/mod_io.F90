module mod_io
    use mod_globals
    use mod_mesh
    implicit none

contains

    ! ==========================================
    ! Subroutine: Read Input File (Smart Version)
    ! ==========================================
    subroutine read_input_file(filename)
        character(len=*), intent(in) :: filename
        integer(i4) :: fid, i, stat, mid
        integer(i4) :: temp_id, n1, n2, mat_id, type_id
        real(dp) :: x_val, y_val, val1, val2
        
        ! Temporary arrays for mapping
        integer(i4), allocatable :: elem_mat_id_map(:)
        real(dp), allocatable :: mat_lib(:,:) ! Stores (Area, E, Alpha) for each material type

        if (rank == 0) print *, "Reading file: ", filename

        open(newunit=fid, file=filename, status='old', action='read', iostat=stat)
        if (stat /= 0) then
            print *, "Error: Cannot open file on Rank", rank
            call MPI_Abort(PETSC_COMM_WORLD, 1, ierr)
        end if

        ! --- 1. Read Control Parameters ---
        read(fid, *) num_nodes
        read(fid, *) num_elems
        read(fid, *) num_mat_props

        ! --- 2. Allocate Mesh Arrays ---
        call allocate_mesh_arrays()
        
        ! Allocate temp mapping arrays
        allocate(elem_mat_id_map(num_elems))
        allocate(mat_lib(num_mat_props, 3)) ! 3 columns: Area, E, Alpha

        ! --- 3. Read Node Coordinates ---
        do i = 1, num_nodes
            read(fid, *) temp_id, x_val, y_val
            X(i, 1) = x_val
            X(i, 2) = y_val
        end do

        ! --- 4. Read Connectivity AND Material ID ---
        do i = 1, num_elems
            read(fid, *) temp_id, n1, n2, mat_id
            connectivity(i, 1) = n1
            connectivity(i, 2) = n2
            
            ! [FIX] Store which material this element uses
            elem_mat_id_map(i) = mat_id
        end do

        ! --- 5. Read Material Properties into Library ---
        do i = 1, num_mat_props
            read(fid, *) temp_id, val1, val2, x_val
            ! Store in our temporary library
            ! val1=Area, val2=E, x_val=Alpha
            if (temp_id <= num_mat_props) then
                mat_lib(temp_id, 1) = val1
                mat_lib(temp_id, 2) = val2
                mat_lib(temp_id, 3) = x_val
            end if
        end do
        
        ! [FIX] Now Distribute Properties to All Elements
        do i = 1, num_elems
            mid = elem_mat_id_map(i) ! Which material does this element use?
            
            ! Check for safety
            if (mid > 0 .and. mid <= num_mat_props) then
                A(i)     = mat_lib(mid, 1)
                E(i)     = mat_lib(mid, 2)
                alpha(i) = mat_lib(mid, 3)
            else
                print *, "Error: Invalid Material ID", mid, " for Element", i
            end if
        end do

        ! --- 6. Read BCs & Loads ---
        fixdof = 0       
        Uf = 0.0_dp      
        P_load = 0.0_dp 

        do
            read(fid, *, iostat=stat) temp_id, type_id, val1, val2
            if (stat /= 0) exit 

            ! Type 1: Boundary Constraint
            if (type_id == 1) then
                if (int(val1) == 1) fixdof((temp_id-1)*2 + 1) = 1
                if (int(val2) == 1) fixdof((temp_id-1)*2 + 2) = 1
            end if

            ! Type 2: Nodal Force
            if (type_id == 2) then
                P_load((temp_id-1)*2 + 1) = val1
                P_load((temp_id-1)*2 + 2) = val2
            end if
            
            if (type_id == 3) heat(temp_id) = val1
        end do

        close(fid)
        
        ! Cleanup temp arrays
        deallocate(elem_mat_id_map)
        deallocate(mat_lib)

        if (rank == 0) print *, "File read successfully (With Material Mapping)."

    end subroutine read_input_file

    ! ==========================================
    ! Subroutine: Allocate Arrays
    ! ==========================================
    subroutine allocate_mesh_arrays()
        ndof = 2
        allocate(X(num_nodes, 2))
        allocate(connectivity(num_elems, 2))
        allocate(A(num_elems))
        allocate(E(num_elems))
        allocate(alpha(num_elems))
        allocate(heat(num_elems))
        allocate(fixdof(num_nodes * ndof))
        allocate(Uf(num_nodes * ndof))
        allocate(P_load(num_nodes * ndof))
        
        A = 0.0_dp; E = 0.0_dp; alpha = 0.0_dp; heat = 0.0_dp
        fixdof = 0; Uf = 0.0_dp; P_load = 0.0_dp
    end subroutine allocate_mesh_arrays

end module mod_io