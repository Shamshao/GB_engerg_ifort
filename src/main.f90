program main
    use GB5DOF_utils
    implicit none

    ! 声明变量
    real(8), dimension(3,3) :: P, Q
    real(8), dimension(:,:), allocatable :: result_geom
    character(len=3) :: whichaxes
    real(8) :: dismax

    ! 初始化输入数据
    P = reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0], [3,3])
    Q = reshape([0.0d0, -1.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0], [3,3])
    whichaxes = '110'
    dismax = 1.0d0

    ! 调用函数
    result_geom = distances_to_set(P, Q, whichaxes, dismax)

    ! 打印结果
    print *, "Result Geometry:"
    print *, result_geom

end program main
