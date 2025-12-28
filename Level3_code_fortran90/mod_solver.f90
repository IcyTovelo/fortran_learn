module mod_solver
    use mod_globals
    use mod_mesh
    implicit none

contains

    ! ==========================================
    ! 核心求解流程 (Solve K * D = P)
    ! ==========================================
    subroutine solve_system()
        ! 1. 变量声明区 (必须在最前面！)
        integer(i4) :: i, j, k, n_total, n_free
        integer(i4), allocatable :: free_idx(:) 
        real(dp), allocatable :: K_reduced(:,:) 
        real(dp), allocatable :: F_reduced(:)   
        real(dp), allocatable :: D_free(:)      
        
        ! --- 之前报错的变量，现在移到这里 ---
        integer(i4) :: r_idx, c_idx 
        
        ! 2. 执行区
        n_total = num_nodes * ndof
        
        ! 分配全局 D 向量
        if (allocated(D)) deallocate(D)
        allocate(D(n_total))
        D = 0.0_dp 
        
        ! 先把已知位移 (Fixed DOFs) 填进去
        do i = 1, n_total
            if (fixdof(i) == 1) then
                D(i) = Uf(i)
            end if
        end do
        
        ! 统计 Free 自由度
        n_free = count(fixdof == 0)
        allocate(free_idx(n_free))
        
        k = 1
        do i = 1, n_total
            if (fixdof(i) == 0) then
                free_idx(k) = i
                k = k + 1
            end if
        end do
        
        print *, "Solving System... (Free DOFs:", n_free, ")"
        
        ! 组装缩减方程组
        allocate(K_reduced(n_free, n_free))
        allocate(F_reduced(n_free))
        
        do i = 1, n_free
            r_idx = free_idx(i) ! 这里直接赋值，不再声明
            
            ! 初始力 = 外力 P
            F_reduced(i) = P(r_idx) 
            
            ! 减去已知位移产生的影响
            do j = 1, n_total
                if (fixdof(j) == 1) then 
                    F_reduced(i) = F_reduced(i) - K_global(r_idx, j) * D(j)
                end if
            end do
            
            ! 填充缩减刚度矩阵 K_ff
            do j = 1, n_free
                c_idx = free_idx(j) ! 这里直接赋值，不再声明
                K_reduced(i, j) = K_global(r_idx, c_idx)
            end do
        end do
        
        ! 调用高斯消去法求解
        allocate(D_free(n_free))
        call gaussian_elimination(n_free, K_reduced, F_reduced, D_free)
        
        ! 把解填回全局 D 向量
        do i = 1, n_free
            D(free_idx(i)) = D_free(i)
        end do
        
        print *, "Solution Converged. Max Displacement:", maxval(abs(D))
        
    end subroutine solve_system

    ! ==========================================
    ! 高斯消去法求解器
    ! ==========================================
    subroutine gaussian_elimination(n, A, b, x)
        integer(i4), intent(in) :: n
        real(dp), intent(inout) :: A(n, n) 
        real(dp), intent(inout) :: b(n)    
        real(dp), intent(out) :: x(n)
        
        integer(i4) :: i, j, k
        real(dp) :: factor, pivot
        
        ! 前向消元
        do k = 1, n-1
            pivot = A(k,k)
            if (abs(pivot) < 1.0e-12_dp) then
                print *, "Error: Singular Matrix at index", k
                stop
            end if
            
            do i = k+1, n
                factor = A(i,k) / pivot
                do j = k, n
                    A(i,j) = A(i,j) - factor * A(k,j)
                end do
                b(i) = b(i) - factor * b(k)
            end do
        end do
        
        ! 回代求解
        x(n) = b(n) / A(n,n)
        do i = n-1, 1, -1
            do j = i+1, n
                b(i) = b(i) - A(i,j) * x(j)
            end do
            x(i) = b(i) / A(i,i)
        end do
        
    end subroutine gaussian_elimination

end module mod_solver