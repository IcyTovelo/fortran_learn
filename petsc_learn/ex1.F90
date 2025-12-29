program main
#include <petsc/finclude/petscsys.h>
    use petscsys
    implicit none

    PetscErrorCode :: ierr
    PetscMPIInt    :: rank, size

    ! 1. 初始化 PETSc (自动启动 MPI 并行环境)
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

    ! 2. 只有初始化成功后，才能调用 MPI 函数
    ! 查询我是第几个进程 (rank)，总共有几个进程 (size)
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

    ! 3. 每个进程都会打印这句话
    print *, "Hello! I am processor", rank, "out of", size

    ! 4. 退出 PETSc，清理内存
    call PetscFinalize(ierr)

end program main
