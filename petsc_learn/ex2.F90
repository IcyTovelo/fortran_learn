program main
#include <petsc/finclude/petscvec.h>
    use petscvec
    implicit none

    PetscErrorCode :: ierr
    PetscMPIInt    :: rank
    Vec            :: x  ! 声明一个向量对象 (只要一个句柄，不需要定义数组大小)
    PetscInt       :: N_global
    PetscScalar    :: value

    ! 1. 初始化 PETSc
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

    N_global = 10  ! 我们定义一个长度为 10 的向量

    if (rank == 0) print *, "Processor 0 is creating the vector..."

    ! 2. 创建向量 (Create)
    call VecCreate(PETSC_COMM_WORLD, x, ierr)
    call VecSetSizes(x, PETSC_DECIDE, N_global, ierr) ! PETSC_DECIDE: 让系统自动决定怎么分
    call VecSetFromOptions(x, ierr)

    ! 3. 赋值 (Set Values)
    ! 我们把向量的所有元素都设为 10.0
    value = 10.0
    call VecSet(x, value, ierr)

    ! 4. 组装 (Assembly) 
    ! 就像投完信要等邮递员收信一样，赋值后必须组装
    call VecAssemblyBegin(x, ierr)
    call VecAssemblyEnd(x, ierr)

    ! 5. 查看向量 (View)
    ! PETSc 会自动收集所有 CPU 的数据，整整齐齐地打印出来
    if (rank == 0) print *, "--- Vector Content ---"
    call VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr)

    ! 6. 清理内存 (Destroy)
    call VecDestroy(x, ierr)
    call PetscFinalize(ierr)

end program main
