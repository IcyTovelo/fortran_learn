
program main
    use mod_globals
    use mod_mesh
    use mod_io
    use mod_stiffness
    use mod_solver
    use mod_stress    ! <--- 1. 新增引用
    implicit none

    integer(i4) :: i, fid_out

    print *, "--- Starting Fortran FEA ---"
    
    ! 1. 读取数据
    call read_input_file('Y3-224.txt')
    
    ! 2. 计算刚度矩阵
    call compute_global_stiffness()
    
    ! 3. 求解位移
    call solve_system()
    
    ! 4. 计算应力
    call compute_stress()  ! <--- 2. 新增调用
    
    ! ==========================================
    ! 5. 导出所有结果
    ! ==========================================
    
    ! (A) 导出位移
    open(newunit=fid_out, file='Results_D.txt', status='replace')
    write(fid_out, *) "Node_ID   Disp_X         Disp_Y"
    do i = 1, num_nodes
        write(fid_out, '(I5, 2F15.8)') i, D((i-1)*2+1), D((i-1)*2+2)
    end do
    close(fid_out)
    
    ! (B) 导出应力 (新增)
    open(newunit=fid_out, file='Results_Stress.txt', status='replace')
    write(fid_out, *) "Elem_ID   Stress (MPa/ForceUnit)"
    do i = 1, num_elems
        ! 格式: 整数占5位，小数占15位
        write(fid_out, '(I5, F20.8)') i, element_stress(i)
    end do
    close(fid_out)
    
    print *, "------------------------------------------------"
    print *, "ALL DONE!"
    print *, "Results saved to: Results_D.txt & Results_Stress.txt"
    print *, "------------------------------------------------"

end program main