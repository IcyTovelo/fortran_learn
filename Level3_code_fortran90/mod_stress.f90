module mod_stress
    use mod_globals
    use mod_mesh
    implicit none

contains

    subroutine compute_stress()
        ! --- 1. 变量声明区 ---
        integer(i4) :: ielem, n1, n2
        integer(i4) :: idx1, idx2 ! 节点1和节点2的自由度起始索引
        real(dp) :: dx, dy, L, c, s
        real(dp) :: u1, v1, u2, v2
        real(dp) :: elongation, strain
        
        ! --- 2. 执行区 ---
        
        ! 分配内存
        if (allocated(element_stress)) deallocate(element_stress)
        allocate(element_stress(num_elems))
        element_stress = 0.0_dp

        print *, "Computing Element Stresses..."

        do ielem = 1, num_elems
            
            ! 获取节点号
            n1 = connectivity(ielem, 1)
            n2 = connectivity(ielem, 2)

            ! 获取几何信息
            dx = X(n2, 1) - X(n1, 1)
            dy = X(n2, 2) - X(n1, 2)
            L = sqrt(dx**2 + dy**2)
            
            if (L <= 0.0_dp) stop "Error: Zero length element in stress calc"

            c = dx / L
            s = dy / L

            ! 获取节点位移
            ! 节点 n1 的位移 (u1, v1)
            idx1 = (n1 - 1) * ndof + 1
            u1 = D(idx1)
            v1 = D(idx1 + 1)
            
            ! 节点 n2 的位移 (u2, v2)
            idx2 = (n2 - 1) * ndof + 1
            u2 = D(idx2)
            v2 = D(idx2 + 1)

            ! 计算伸长量 (Delta L) = [-c -s c s] * {u}
            ! 物理意义：投影在杆轴线方向的相对位移
            elongation = -c * u1 - s * v1 + c * u2 + s * v2
            
            ! 计算应变 (Strain)
            strain = elongation / L
            
            ! 计算应力 (Stress) = E * Strain
            element_stress(ielem) = E(ielem) * strain
            
        end do

        print *, "Stress calculation completed."

    end subroutine compute_stress

end module mod_stress