program main
#include <petsc/finclude/petscmat.h>
    use petscmat
    implicit none

    PetscErrorCode :: ierr
    PetscMPIInt    :: rank
    Mat            :: A       ! 声明一个矩阵对象
    PetscInt       :: n, i
    PetscScalar    :: v

    ! 1. 初始化
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

    ! 我们要造一个 4x4 的矩阵
    n = 4

    if (rank == 0) print *, "Creating a 4x4 Parallel Matrix..."

    ! 2. 创建矩阵 (Buy the notebook)
    call MatCreate(PETSC_COMM_WORLD, A, ierr)
    
    ! 3. 设定大小 (Set Sizes)
    ! PETSC_DECIDE 表示让系统自动决定本地分多少行，我们只指定总共是 n x n
    call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
    
    ! 4. 设定类型 (Set Type)
    ! MATMPIAIJ 是 PETSc 最常用的稀疏矩阵格式 (适合并行)
    call MatSetType(A, MATMPIAIJ, ierr)
    
    ! 5. 真正分配内存 (Setup)
    call MatSetUp(A, ierr)

    ! 6. 填数 (Fill values)
    ! 我们要在对角线上填入 10.0 (即 (0,0), (1,1), (2,2), (3,3) = 10.0)
    ! 注意：这里所有 CPU 都会跑这个循环，重复填没关系，PETSc 会处理
    do i = 0, n-1
        v = 10.0
        ! MatSetValue(矩阵, 行号, 列号, 数值, 模式, 错误码)
        ! INSERT_VALUES 表示“写入新值”
        call MatSetValue(A, i, i, v, INSERT_VALUES, ierr)
    end do

    ! 7. 组装 (Assembly) - 这一步绝对不能少！
    ! 哪怕只是填几个数，也必须告诉所有 CPU：“填完啦，同步一下数据！”
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

    ! 8. 打印出来看看
    if (rank == 0) print *, "--- Matrix Content ---"
    call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)

    ! 9. 清理
    call MatDestroy(A, ierr)
    call PetscFinalize(ierr)

end program main
