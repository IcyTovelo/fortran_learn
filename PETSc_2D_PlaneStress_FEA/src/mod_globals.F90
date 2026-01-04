module mod_globals
  !---------------------------------------------------------------------------
  ! MODULE: mod_globals
  ! DESCRIPTION: Global variables container. Stores the mesh, material properties,
  !              boundary conditions, and loads read from the input file.
  !---------------------------------------------------------------------------
  use mod_types, only: i4, dp
  implicit none
  
  ! --- SCALARS: Single values defining model properties ---
  integer(i4) :: analysis_type ! 1 = Plane Stress, 2 = Plane Strain
  integer(i4) :: nelem         ! Total number of elements (triangles)
  integer(i4) :: npoin         ! Total number of nodes (points)
  integer(i4) :: nfix          ! Number of nodes with prescribed displacements (BCs)
  integer(i4) :: nload         ! Number of point loads
  
  real(dp)    :: thickness     ! Plate thickness (for Plane Stress)
  real(dp)    :: young         ! Young's Modulus (E)
  real(dp)    :: poisson       ! Poisson's Ratio (nu)
  real(dp)    :: amp_factor    ! Plotting amplification factor
  
  ! --- ARRAYS: Allocatable tables for mesh data ---
  
  ! 1. Connectivity Matrix (corresponds to 'lnods' in MATLAB)
  !    Dimensions: (3, nelem) -> Stores 3 node IDs for each element
  integer(i4), allocatable :: lnods(:,:) 
  
  ! 2. Nodal Coordinates (corresponds to 'coord' in MATLAB)
  !    Dimensions: (2, npoin) -> Stores X and Y coordinates for each node
  real(dp), allocatable :: coord(:,:)
  
  ! 3. Boundary Conditions (Prescribed Displacements)
  !    Input file format: NodeID  Type(11/10/01)  Val1  Val2
  integer(i4), allocatable :: bc_node_ids(:)    ! Stores the Node ID
  integer(i4), allocatable :: bc_types(:)       ! Stores the constraint type (e.g. 11, 01)
  real(dp),    allocatable :: bc_values(:,:)    ! Dimensions: (2, nfix) -> Stores value 1 and value 2
  
  ! 4. Point Loads
  !    Input file format: NodeID  ForceX  ForceY
  integer(i4), allocatable :: load_node_ids(:)  ! Stores the loaded Node ID
  real(dp),    allocatable :: load_values(:,:)  ! Dimensions: (2, nload) -> Stores Fx and Fy
  
end module mod_globals