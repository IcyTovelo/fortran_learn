program main
    use mod_globals
    use mod_mesh
    use mod_io
    use mod_stiffness
    use mod_solver 
    implicit none

    integer(i4) :: i, fid_out

    print *, "--- Starting Fortran FEA ---"
    
    ! 1. 读取数据
    call read_input_file('Y3-224.txt')
    
    ! 2. 计算刚度矩阵
    call compute_global_stiffness()
    
    ! (已移除) 计算热荷载 (Q0)
    
    ! 3. 求解位移 (D)
    call solve_system()
    
    ! 验证：打印第一个节点的位移
    print *, "--------------------------"
    print *, "Displacement at Node 1 (DX, DY):"
    print *, D(1), D(2)
    print *, "--------------------------"

    ! 导出结果到文件
    open(newunit=fid_out, file='Results_D.txt', status='replace')
    write(fid_out, *) "Node_ID   Disp_X         Disp_Y"
    do i = 1, num_nodes
        ! 格式化输出: 整数占5位，小数占15位
        write(fid_out, '(I5, 2F15.8)') i, D((i-1)*2+1), D((i-1)*2+2)
    end do
    close(fid_out)
    print *, "Displacement results saved to Results_D.txt"

end program main