program main
  !---------------------------------------------------------------------------
  ! PROGRAM: PETSc_2D_PlaneStress_FEA
  ! AUTHOR:  [Shen Ma]
  ! DATE:    2026-01-04
  ! DESCRIPTION: Main driver for the Finite Element Analysis Solver.
  !              Currently in Phase 1 (Initialization).
  !---------------------------------------------------------------------------
  use mod_types, only: dp, i4
  use mod_globals
  implicit none

  ! --- Main Execution Start ---
  print *, "=========================================================="
  print *, "              2D Linear Elastic FEA Solver       "
  print *, "=========================================================="
  
  ! 1. Precision Verification
  print *, "--> [Status] Checking Precision..."
  print *, "    Double Precision (dp) Epsilon: ", epsilon(1.0_dp)
  
  ! 2. Module Verification
  ! We haven't allocated arrays yet, but if the code compiles, 
  ! it means mod_globals is correctly linked.
  print *, "--> [Status] Modules linked successfully."
  print *, "--> [Status] Ready for Phase 2: Data Ingestion."
  
  print *, "=========================================================="
  print *, "       Phase 1: Infrastructure Setup COMPLETE.            "
  print *, "=========================================================="

end program main