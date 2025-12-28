module mod_stiffness
    use mod_globals
    use mod_mesh
    implicit none

contains

    subroutine compute_global_stiffness()
        ! 局部变量声明
        integer(i4) :: ielem
        integer(i4) :: n1, n2, idof_start, jdof_start
        real(dp) :: dx, dy, L, c, s
        real(dp) :: coef
        real(dp) :: Ke_sub(2, 2) ! 2x2 的子刚度矩阵
        integer(i4) :: idx(2), jdx(2) ! 自由度索引

        ! 1. 初始化全局刚度矩阵为 0
        ! 总自由度 = 节点数 * ndof
        allocate(K_global(num_nodes * ndof, num_nodes * ndof))
        K_global = 0.0_dp

        print *, "Computing Stiffness Matrix..."

        ! 2. 循环遍历所有单元
        do ielem = 1, num_elems
            
            ! (1) 获取节点编号
            n1 = connectivity(ielem, 1)
            n2 = connectivity(ielem, 2)

            ! (2) 计算几何信息 (长度 L, sin, cos)
            ! MATLAB: Xij = X(nodej,:) - X(nodei,:)
            dx = X(n2, 1) - X(n1, 1)
            dy = X(n2, 2) - X(n1, 2)
            
            ! MATLAB: L = norm(Xij, 2)
            L = sqrt(dx**2 + dy**2)
            
            if (L <= 0.0_dp) then
                print *, "Error: Element length is zero!", ielem
                stop
            end if

            c = dx / L
            s = dy / L

            ! (3) 计算刚度系数 E*A/L
            coef = (E(ielem) * A(ielem)) / L

            ! (4) 计算子矩阵 (T' * T) * coef
            ! [ c*c  c*s ]
            ! [ s*c  s*s ]
            Ke_sub(1, 1) = coef * c * c
            Ke_sub(1, 2) = coef * c * s
            Ke_sub(2, 1) = coef * s * c
            Ke_sub(2, 2) = coef * s * s

            ! (5) 组装到全局矩阵 K_global
            ! 确定全局自由度索引 (Fortran 也是从 1 开始，算法同 MATLAB)
            ! Node 1 的自由度 (x, y)
            idx(1) = (n1 - 1) * ndof + 1
            idx(2) = (n1 - 1) * ndof + 2
            
            ! Node 2 的自由度 (x, y)
            jdx(1) = (n2 - 1) * ndof + 1
            jdx(2) = (n2 - 1) * ndof + 2

            ! MATLAB logic:
            ! K(idof,idof) += Ke
            ! K(idof,jdof) -= Ke
            ! K(jdof,idof) -= Ke
            ! K(jdof,jdof) += Ke
            
            ! Block 1: K(ii) += Ke
            K_global(idx, idx) = K_global(idx, idx) + Ke_sub
            
            ! Block 2: K(ij) -= Ke
            K_global(idx, jdx) = K_global(idx, jdx) - Ke_sub
            
            ! Block 3: K(ji) -= Ke
            K_global(jdx, idx) = K_global(jdx, idx) - Ke_sub
            
            ! Block 4: K(jj) += Ke
            K_global(jdx, jdx) = K_global(jdx, jdx) + Ke_sub

        end do
        
        print *, "Stiffness Matrix Assembled. Size:", size(K_global,1), "x", size(K_global,2)

    end subroutine compute_global_stiffness

end module mod_stiffness