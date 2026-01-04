module mod_mesh
  use mod_types, only: dp, i4
  use mod_globals
  implicit none

contains

  subroutine ReadMesh(filename)
    character(len=*), intent(in) :: filename
    integer(i4) :: unit_num, ios, i, temp_id, k
    character(len=256) :: line
    
    print *, "--> [Mesh] Opening file: ", trim(filename)
    open(newunit=unit_num, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
       print *, "Error: Could not open file: ", filename
       stop
    end if

    ! --- 1. ANALYSIS_TYPE ---
    call FindAndReadInt(unit_num, "ANALYSIS_TYPE", analysis_type)
    print *, "    Analysis Type:", analysis_type

    ! --- 2. ELEMENTS ---
    call FindAndReadInt(unit_num, "ELEMENTS", nelem)
    print *, "    Number of Elements:", nelem
    
    allocate(lnods(3, nelem))
    do i = 1, nelem
       read(unit_num, *) temp_id, lnods(1,i), lnods(2,i), lnods(3,i)
    end do

    ! --- 3. NODE_COORDINATES ---
    call FindAndReadInt(unit_num, "NODE_COORDINATES", npoin)
    print *, "    Number of Nodes:", npoin
    
    allocate(coord(2, npoin))
    do i = 1, npoin
       read(unit_num, *) temp_id, coord(1,i), coord(2,i)
    end do

    ! --- 4. THICKNESS ---
    call FindAndReadDouble(unit_num, "THICKNESS", thickness)
    print *, "    Thickness:", thickness

    ! --- 5. PRESCRIBED DISPLACEMENTS ---
    call FindAndReadInt(unit_num, "NODES_WITH_PRESCRIBED_DISPLACEMENTS", nfix)
    print *, "    Number of Constraints:", nfix
    
    allocate(bc_node_ids(nfix))
    allocate(bc_types(nfix))
    allocate(bc_values(2, nfix))
    do i = 1, nfix
       read(unit_num, *) bc_node_ids(i), bc_types(i), bc_values(1,i), bc_values(2,i)
    end do

    ! --- 6. MATERIAL PROPERTIES ---
    call FindAndReadDouble(unit_num, "YOUNG_MODULUS", young)
    call FindAndReadDouble(unit_num, "POISSON_RATIO", poisson)
    print *, "    Material: E =", young, " Nu =", poisson

    ! --- 7. POINT LOADS ---
    call FindAndReadInt(unit_num, "POINT_LOADS", nload)
    print *, "    Number of Point Loads:", nload
    
    allocate(load_node_ids(nload))
    allocate(load_values(2, nload))
    do i = 1, nload
       read(unit_num, *) load_node_ids(i), load_values(1,i), load_values(2,i)
    end do

    close(unit_num)
    print *, "--> [Mesh] Data ingestion complete."

  end subroutine ReadMesh

  ! --- Helper: Find line containing keyword and read one Integer from it ---
  subroutine FindAndReadInt(unit, keyword, out_val)
    integer(i4), intent(in) :: unit
    character(len=*), intent(in) :: keyword
    integer(i4), intent(out) :: out_val
    character(len=256) :: line
    integer :: pos_eq
    
    rewind(unit)
    do
       read(unit, '(A)', iostat=out_val) line ! utilize out_val as temporary io status
       if (out_val /= 0) stop "Error: Keyword not found"
       if (index(line, keyword) > 0) exit
    end do
    
    ! Replace '=' with space so list-directed read works
    pos_eq = index(line, '=')
    if (pos_eq > 0) line(pos_eq:pos_eq) = ' '
    
    ! Read the integer from the modified string (skipping the keyword text)
    ! We assume the number is the last thing on the line
    read(line(pos_eq+1:), *) out_val
  end subroutine FindAndReadInt

  ! --- Helper: Find line containing keyword and read one Double from it ---
  subroutine FindAndReadDouble(unit, keyword, out_val)
    integer(i4), intent(in) :: unit
    character(len=*), intent(in) :: keyword
    real(dp), intent(out) :: out_val
    character(len=256) :: line
    integer :: pos_eq, ios
    
    rewind(unit)
    do
       read(unit, '(A)', iostat=ios) line
       if (ios /= 0) stop "Error: Keyword not found"
       if (index(line, keyword) > 0) exit
    end do
    
    pos_eq = index(line, '=')
    if (pos_eq > 0) line(pos_eq:pos_eq) = ' '
    read(line(pos_eq+1:), *) out_val
  end subroutine FindAndReadDouble

end module mod_mesh