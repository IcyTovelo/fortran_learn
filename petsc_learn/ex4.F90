program main
! 注意：我们要引用 ksp 的头文件了，因为它包含了 Vec 和 Mat
#include <petsc/finclude/petscksp.h>
    use petscksp
    implicit none

    PetscErrorCode :: ierr
    PetscMPIInt    :: rank
    PetscInt       :: n, i
    PetscScalar    :: v
    
    ! 定义三个核心对象
    Vec            :: x, b      ! x是解(位移), b是右端项(力)
    Mat            :: A         ! 矩阵
    KSP            :: ksp       ! 求解器对象

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

    n = 4  ! 还是用 4x4 的系统
    if (rank == 0) print *, "Solving Ax = b (Expect x = 2.0)..."

    ! ------------------------------------------
    ! 1. 准备数据 (矩阵 A 和 向量 b)
    ! ------------------------------------------
    
    ! (A) 创建向量 x 和 b
    call VecCreate(PETSC_COMM_WORLD, x, ierr)
    call VecSetSizes(x, PETSC_DECIDE, n, ierr)
    call VecSetFromOptions(x, ierr)
    ! 把 b 复制一份 x 的结构 (这样它们大小就一样了)
    call VecDuplicate(x, b, ierr)

    ! (B) 给 b 赋值：全部设为 20.0
    v = 20.0
    call VecSet(b, v, ierr)
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)

    ! (C) 创建矩阵 A 并赋值：对角线设为 10.0
    call MatCreate(PETSC_COMM_WORLD, A, ierr)
    call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
    call MatSetType(A, MATMPIAIJ, ierr)
    call MatSetUp(A, ierr)

    do i = 0, n-1
        v = 10.0
        call MatSetValue(A, i, i, v, INSERT_VALUES, ierr)
    end do
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

    ! ------------------------------------------
    ! 2. 配置求解器 (KSP) - 关键步骤！
    ! ------------------------------------------
    
    ! (A) 创建求解器
    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)

    ! (B) 告诉求解器我们要解哪个方程 (A * x = b)
    ! 这里的两个 A，意思是“系统矩阵”和“预处理矩阵”都是 A
    call KSPSetOperators(ksp, A, A, ierr)

    ! (C) 设置默认选项 (让 PETSc 自动选择最好的算法)
    call KSPSetFromOptions(ksp, ierr)

    ! ------------------------------------------
    ! 3. 开始求解！ (Solver)
    ! ------------------------------------------
    ! Solve(ksp, b, x) -> 结果存进 x
    call KSPSolve(ksp, b, x, ierr)

    ! ------------------------------------------
    ! 4. 查看结果
    ! ------------------------------------------
    if (rank == 0) print *, "--- Solution Vector x ---"
    call VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr)

    ! 5. 清理 (记得把所有对象都销毁)
    call KSPDestroy(ksp, ierr)
    call VecDestroy(x, ierr)
    call VecDestroy(b, ierr)
    call MatDestroy(A, ierr)
    call PetscFinalize(ierr)

end program main