module mod_element
  !---------------------------------------------------------------------------
  ! MODULE: mod_element
  ! DESCRIPTION: Contains the physics core for element stiffness calculations.
  !              Implements CST (Constant Strain Triangle) for Plane Stress.
  !---------------------------------------------------------------------------
  use mod_types, only: dp, i4
  use mod_globals, only: young, poisson, thickness, analysis_type
  implicit none

contains

  !---------------------------------------------------------------------------
  ! SUBROUTINE: ComputeKe_CST
  ! PURPOSE: Computes the 6x6 Element Stiffness Matrix (Ke) for a CST element.
  ! INPUT:   x, y (Real arrays of size 3, containing nodal coordinates)
  ! OUTPUT:  ke   (Real array 6x6, the stiffness matrix)
  !---------------------------------------------------------------------------
  subroutine ComputeKe_CST(x, y, ke)
    real(dp), intent(in)  :: x(3), y(3)
    real(dp), intent(out) :: ke(6,6)
    
    ! Local variables for B matrix, D matrix, and Geometry
    real(dp) :: B(3,6), D(3,3), BT(6,3), DB(3,6)
    real(dp) :: area, detJ, const_D
    integer(i4) :: i, j, k

    ! --- 1. Compute Geometry (Area of Triangle) ---
    ! 2*Area = determinant of Jacobian
    ! formula: x1(y2-y3) + x2(y3-y1) + x3(y1-y2)
    detJ = x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2))
    area = 0.5_dp * abs(detJ)

    ! Safety check: Prevent zero area elements
    if (area < 1.0e-14_dp) then
       print *, "Error: Element has zero or negative area!"
       stop
    end if

    ! --- 2. Compute B Matrix (Strain-Displacement) ---
    ! B = (1/2A) * [ y23  0   y31  0   y12  0 ]
    !              [  0  x32   0  x13   0  x21]
    !              [ x32 y23  x13 y31  x21 y12]
    
    B = 0.0_dp
    
    ! Node 1 derivatives
    B(1,1) = y(2) - y(3)        ! y23
    B(2,2) = x(3) - x(2)        ! x32
    B(3,1) = x(3) - x(2)        ! x32
    B(3,2) = y(2) - y(3)        ! y23
    
    ! Node 2 derivatives
    B(1,3) = y(3) - y(1)        ! y31
    B(2,4) = x(1) - x(3)        ! x13
    B(3,3) = x(1) - x(3)        ! x13
    B(3,4) = y(3) - y(1)        ! y31
    
    ! Node 3 derivatives
    B(1,5) = y(1) - y(2)        ! y12
    B(2,6) = x(2) - x(1)        ! x21
    B(3,5) = x(2) - x(1)        ! x21
    B(3,6) = y(1) - y(2)        ! y12
    
    ! Scale B by 1/(2*Area) = 1/detJ
    B = B / detJ

    ! --- 3. Compute D Matrix (Constitutive Matrix) ---
    ! For Plane Stress (analysis_type = 1)
    D = 0.0_dp
    if (analysis_type == 1) then
       ! Factor E / (1 - nu^2)
       const_D = young / (1.0_dp - poisson**2)
       
       D(1,1) = 1.0_dp
       D(1,2) = poisson
       D(2,1) = poisson
       D(2,2) = 1.0_dp
       D(3,3) = (1.0_dp - poisson) / 2.0_dp
       
       D = D * const_D
    else
       ! Plane Strain logic (Placeholder for future)
       print *, "Error: Only Plane Stress (type=1) is implemented so far."
       stop
    end if

    ! --- 4. Compute Ke = B^T * D * B * Thickness * Area ---
    ! Step 4a: Calculate temp product D * B
    ! DB is 3x6
    DB = matmul(D, B)
    
    ! Step 4b: Calculate B^T (Transpose of B)
    BT = transpose(B)
    
    ! Step 4c: Calculate B^T * D * B (Result is 6x6)
    ke = matmul(BT, DB)
    
    ! Step 4d: Scale by Volume factor (Area * Thickness)
    ke = ke * (area * thickness)

  end subroutine ComputeKe_CST

end module mod_element