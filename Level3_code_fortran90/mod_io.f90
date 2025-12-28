module mod_io
    use mod_globals
    use mod_mesh
    implicit none

contains

    subroutine read_input_file(filename)
        character(len=*), intent(in) :: filename
        
        ! --- 1. 变量声明区 ---
        integer(i4) :: fid, iostat, i
        character(len=256) :: buffer
        character(len=50) :: dummy  ! 加长一点，防止读长单词报错
        integer(i4) :: temp_int
        
        ! 临时数组
        real(dp), allocatable :: temp_props(:,:)
        
        ! 计数器
        integer(i4) :: n_fix, n_point_load, n_elem_load
        
        ! 索引变量
        integer(i4) :: node_idx, dof_x, dof_y, type_fix
        
        ! 数值变量
        real(dp) :: val1, val2
        integer(i4) :: id, n1, n2, mat_id, nid, pid, mat_idx, elem_idx
        real(dp) :: area_val, load_val

        ! --- 2. 执行区 ---
        
        ! 打开文件
        open(newunit=fid, file=filename, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            print *, "Error: Cannot open file ", filename
            stop
        end if
        
        print *, "Reading file: ", filename

        ! 跳过标题 (这两行通常比较长，还是用整行读比较安全)
        read(fid, '(A)') buffer 
        read(fid, '(A)') buffer 
        
        ! 读取控制参数 (会自动跳过空行)
        read(fid, *) num_elems, num_nodes, num_mat_props, n_fix, n_point_load, n_elem_load
        
        ! 读取 ndof
        ! NODE DEGREE OF FREEDOM = 2
        read(fid, *) dummy, dummy, dummy, dummy, dummy, ndof
        
        ! ---------------------------------------------------------
        ! 【修改重点】这里改用 * 读取，自动跳过空行，防止读错位
        ! ---------------------------------------------------------
        ! 读取 LENGTHUNITS = mm
        read(fid, *) dummy, dummy, dummy 
        
        ! 读取 AXIALUNITS = KN
        read(fid, *) dummy, dummy, dummy 
        ! ---------------------------------------------------------
        
        ! 读取 Elements
        ! ELEMENTS = 21
        read(fid, *) dummy, dummy, temp_int 
        
        allocate(connectivity(num_elems, 2))
        allocate(A(num_elems))
        allocate(E(num_elems))     
        allocate(alpha(num_elems)) 
        allocate(heat(num_elems))
        allocate(freedof(num_elems)) 
        
        do i = 1, num_elems
            read(fid, *) id, n1, n2, mat_id, area_val
            connectivity(i, 1) = n1
            connectivity(i, 2) = n2
            A(i) = area_val
            alpha(i) = real(mat_id, dp) 
            heat(i) = 0.0_dp
        end do

        ! 读取节点坐标
        read(fid, *) dummy, dummy, temp_int
        allocate(X(num_nodes, 2))
        
        do i = 1, num_nodes
            read(fid, *) nid, X(i, 1), X(i, 2)
        end do
        
        ! 读取材料
        read(fid, *) dummy, dummy, dummy 
        allocate(temp_props(num_mat_props, 3)) 
        
        do i = 1, num_mat_props
            read(fid, *) pid, temp_props(i, 1), temp_props(i, 2), temp_props(i, 3)
        end do
        
        ! 映射材料
        do i = 1, num_elems
            mat_idx = int(alpha(i))
            if (mat_idx <= num_mat_props .and. mat_idx > 0) then
                E(i)     = temp_props(mat_idx, 1)
                alpha(i) = temp_props(mat_idx, 2) 
            else
                print *, "Error: Material ID out of range for element", i
            end if
        end do
        deallocate(temp_props)

        ! 读取约束
        read(fid, *) dummy, dummy, dummy 
        
        allocate(fixdof(num_nodes * ndof)) 
        allocate(Uf(num_nodes * ndof))
        fixdof = 0 
        Uf = 0.0_dp
        
        do i = 1, n_fix
            read(fid, *) node_idx, type_fix, val1, val2
            dof_x = (node_idx - 1) * ndof + 1
            dof_y = (node_idx - 1) * ndof + 2
            
            if (type_fix == 11) then
                fixdof(dof_x) = 1; Uf(dof_x) = val1
                fixdof(dof_y) = 1; Uf(dof_y) = val2
            else if (type_fix == 10) then
                fixdof(dof_x) = 1; Uf(dof_x) = val1
            else if (type_fix == 1) then 
                fixdof(dof_y) = 1; Uf(dof_y) = val2
            end if
        end do
        
        ! 读取 Point Load
        read(fid, *) dummy, dummy, dummy 
        allocate(P(num_nodes * ndof))
        P = 0.0_dp
        
        if (n_point_load > 0) then
            do i = 1, n_point_load
                read(fid, *) node_idx, val1, val2
                dof_x = (node_idx - 1) * ndof + 1
                dof_y = (node_idx - 1) * ndof + 2
                P(dof_x) = val1
                P(dof_y) = val2
            end do
        end if
        
        ! 读取 Element Load
        read(fid, *) dummy, dummy, dummy 
        if (n_elem_load > 0) then
             do i = 1, n_elem_load
                read(fid, *) elem_idx, dummy, load_val 
                heat(elem_idx) = load_val
             end do
        end if
        
        close(fid)
        print *, "File read successfully. Nodes:", num_nodes, " Elements:", num_elems

    end subroutine read_input_file

end module mod_io