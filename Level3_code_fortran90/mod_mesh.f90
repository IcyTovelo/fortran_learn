module mod_mesh
    use mod_globals
    implicit none

    ! --- 标量数据 (Scalars) ---
    integer(i4) :: ndof             ! 对应 MATLAB: ndof (节点自由度)
    integer(i4) :: num_nodes        ! 节点数量 (用于分配数组大小)
    integer(i4) :: num_elems        ! 单元数量
    integer(i4) :: num_mat_props    ! 材料属性数量
    
    ! --- 动态数组 (对应 MATLAB 的矩阵) ---
    ! allocatable 表示这些数组的大小还没定，要读了文件才知道
    
    ! 1. 节点坐标矩阵 X
    ! MATLAB: X (N行, 2列)
    real(dp), allocatable :: X(:, :) 
    
    ! 2. 单元连接矩阵 connectivity
    ! MATLAB: connectivity (N行, 2列) -> 存储节点编号
    integer(i4), allocatable :: connectivity(:, :)
    
    ! 3. 材料属性 E (杨氏模量) 和 A (截面积)
    ! MATLAB: E, A (向量)
    real(dp), allocatable :: E(:)
    real(dp), allocatable :: A(:)
    
    ! 4. 热膨胀相关
    ! MATLAB: alpha, heat
    real(dp), allocatable :: alpha(:)
    real(dp), allocatable :: heat(:)

    ! 5. 边界条件与荷载
    ! MATLAB: fixdof (约束自由度索引), freedof (自由自由度索引)
    integer(i4), allocatable :: fixdof(:)
    integer(i4), allocatable :: freedof(:)
    
    ! MATLAB: P (节点荷载), Uf (已知位移)
    real(dp), allocatable :: P(:)
    real(dp), allocatable :: Uf(:)

    ! ... (接在 Uf(:) 下面)
    
    ! 全局刚度矩阵 K (Global Stiffness Matrix)
    real(dp), allocatable :: K_global(:,:)  ! 二维数组

contains

    ! 这里以后可以放一些清理内存的子程序

end module mod_mesh