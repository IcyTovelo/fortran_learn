module mod_types
  !---------------------------------------------------------------------------
  ! MODULE: mod_types
  ! DESCRIPTION: Defines the floating point precision for the entire project.
  !              Using standard 64-bit Double Precision for FEA accuracy.
  !---------------------------------------------------------------------------
  implicit none
  
  ! Define Integer Parameter for standard 4-byte integers
  integer, parameter :: i4 = selected_int_kind(9)
  
  ! Define Real Parameter for Double Precision (15 digits of accuracy)
  ! This corresponds to 'double' in C or 'float64' in MATLAB.
  integer, parameter :: dp = selected_real_kind(15, 307)

end module mod_types