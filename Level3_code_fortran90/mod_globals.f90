module mod_globals
    implicit none
    
    ! 1. 定义双精度 (Double Precision)
    ! selected_real_kind(15, 307) 对应标准的 64位浮点数
    ! 这相当于 MATLAB 的默认数值类型
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    ! 2. 定义整数类型 (以防万一需要长整型)
    integer, parameter :: i4 = selected_int_kind(9)

end module mod_globals