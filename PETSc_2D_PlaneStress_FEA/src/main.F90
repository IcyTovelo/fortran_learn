program main
  !---------------------------------------------------------------------------
  ! PROGRAM: PETSc_2D_PlaneStress_FEA
  ! PHASE:   3 - Physics Verification (Single Element Test)
  !---------------------------------------------------------------------------
  use mod_types
  use mod_globals
  use mod_mesh
  use mod_element  ! <--- Added the new Physics Module
  implicit none

  character(len=100) :: input_file
  
  ! Variables for testing Element 1
  integer(i4) :: n1, n2, n3
  real(dp)    :: x_local(3), y_local(3)
  real(dp)    :: Ke_local(6,6)
  integer(i4) :: i, j

  print *, "=========================================================="
  print *, "            2D Linear Elastic FEA Solver        "
  print *, "=========================================================="
  
  ! 1. Data Ingestion
  input_file = "data/input.txt"
  call ReadMesh(input_file)

  ! 2. PHYSICS TEST: Calculate Stiffness for Element #1
  print *, " "
  print *, "--> [Test] Testing Physics Engine on Element 1..."
  
  ! Get Node IDs for Element 1
  n1 = lnods(1, 1)
  n2 = lnods(2, 1)
  n3 = lnods(3, 1)
  
  print *, "    Element 1 Nodes:", n1, n2, n3
  
  ! Get Coordinates
  x_local(1) = coord(1, n1); y_local(1) = coord(2, n1)
  x_local(2) = coord(1, n2); y_local(2) = coord(2, n2)
  x_local(3) = coord(1, n3); y_local(3) = coord(2, n3)
  
  print *, "    Node 1 Coords:", x_local(1), y_local(1)
  print *, "    Node 2 Coords:", x_local(2), y_local(2)
  print *, "    Node 3 Coords:", x_local(3), y_local(3)

  ! CALL THE NEW PHYSICS ROUTINE
  call ComputeKe_CST(x_local, y_local, Ke_local)
  
  ! 3. Print the Result (So we can compare with MATLAB)
  print *, " "
  print *, "--> [Result] Element 1 Stiffness Matrix (Ke):"
  do i = 1, 6
     write(*, '(6(F12.4, 1X))') (Ke_local(i, j), j=1, 6)
  end do
  print *, "----------------------------------------------------------"
  print *, "--> [Status] Phase 3 Single Element Test COMPLETE."

end program main