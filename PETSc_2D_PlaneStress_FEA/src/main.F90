program main
  !---------------------------------------------------------------------------
  ! PROGRAM: PETSc_2D_PlaneStress_FEA
  ! PHASE:   2 - Data Ingestion
  !---------------------------------------------------------------------------
  use mod_types
  use mod_globals
  use mod_mesh   ! <--- Added this line
  implicit none

  character(len=100) :: input_file

  print *, "=========================================================="
  print *, "           2D Linear Elastic FEA Solver       "
  print *, "=========================================================="
  
  ! 1. Initialize
  input_file = "data/input.txt"

  ! 2. Call the Mesh Reader
  print *, "--> [Status] Phase 2: Starting Data Ingestion..."
  call ReadMesh(input_file)

  ! 3. Verification (Double Check)
  print *, "----------------------------------------------------------"
  print *, "    VERIFICATION REPORT:"
  print *, "    Total Nodes Read:    ", npoin
  print *, "    Total Elements Read: ", nelem
  print *, "    Node 1 Coords:       ", coord(1,1), coord(2,1)
  print *, "    Element 1 Connect:   ", lnods(1,1), lnods(2,1), lnods(3,1)
  print *, "----------------------------------------------------------"
  
  print *, "--> [Status] Phase 2 COMPLETE. Ready for Physics Engine."

end program main