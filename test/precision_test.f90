program precision_test
    implicit none
    
    ! 1. 定义一个默认实数 (通常是单精度)
    real :: x_single
    
    ! 2. 定义一个强制双精度实数 (kind=8)
    real(kind=8) :: x_double
    
    ! --- 赋值 ---
    ! 注意：在 Fortran 里写 3.0 默认是单精度，3.0d0 才是双精度
    x_single = 1.0 / 3.0
    x_double = 1.0d0 / 3.0d0 
    
    ! --- 打印结果 (强制显示小数点后 20 位) ---
    print *, "----------------------------------------"
    print *, "Comparison of 1/3 (0.3333...)"
    print *, "----------------------------------------"
    
    ! F30.20 意思是：总宽30，小数点后显示20位
    write(*, '(A, F30.20)') "Default (Single): ", x_single
    write(*, '(A, F30.20)') "Defined (Double): ", x_double
    
end program precision_test