program main
    use mod_globals
    use mod_mesh
    use mod_io
    use mod_stiffness
    implicit none

    integer(i4) :: i, fid_out ! 新增的临时变量

    print *, "--- Starting Fortran FEA ---"
    
    ! 1. 读取数据
    call read_input_file('Y3-224.txt')
    
    ! 2. 计算刚度矩阵
    call compute_global_stiffness()
    
    ! 3. 验证 K(1,1)
    print *, "--------------------------"
    print *, "Check K Matrix Value (K(1,1)):"
    print *, K_global(1, 1)
    print *, "--------------------------"

    ! ====================================================
    ! 4. 【新增】把完整矩阵写入 K_matrix.txt 文件
    ! ====================================================
    print *, "Exporting K matrix to K_matrix.txt..."
    
    open(newunit=fid_out, file='K_matrix.txt', status='replace')
    
    do i = 1, size(K_global, 1)
        ! 每一行写一遍，用默认格式 (*) 分隔，方便阅读
        write(fid_out, *) K_global(i, :)
    end do
    
    close(fid_out)
    print *, "Done! You can open K_matrix.txt to see the full matrix."
    ! ====================================================

end program main