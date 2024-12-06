module GB5DOF_utils
    implicit none
  real(8), parameter :: pi = 3.141592653589793d0  ! 手动定义pi
contains

    ! 计算到一组轴的距离的函数
      function distances_to_set(P, Q, whichaxes, dismax) result(geom)
        ! 输入参数声明
        real(8), dimension(3,3), intent(in) :: P, Q         ! 旋转矩阵
        character(len=*), intent(in) :: whichaxes            ! 轴的类型（'100', '110', '111'）
        real(8), intent(in), optional :: dismax              ! 最大距离（可选）

        ! 局部变量声明
        real(8), dimension(:,:), allocatable :: geom         ! 输出几何参数
        integer :: naxes, i, j, thisindex
        real(8) :: dis, dismax_local, theta, phi, ksi, eta, theta1, theta2, period
        real(8), dimension(:,:), allocatable :: axes, dirs  ! 可分配数组
        real(8), dimension(24) :: distances, ksis, etas, phis
        real(8), dimension(3, 3) :: rotX90, rotY90, rotZ90, rotZ90m, RA
        real(8), dimension(3, 3) :: V(24, 3, 3)  ! 固定维度以存储24个3x3的矩阵
        real(8), dimension(3) :: ax, dir, dir2, n1, n2, m1, m2, axi
        real(8) :: psi, dotp
        real(8), dimension(3, 3) :: R  ! 错位矩阵

        ! 初始化变量
        if (.not. present(dismax)) then
            dismax_local = 0.999999
        else
            dismax_local = dismax
        end if

        ! 根据输入的轴类型定义 axes
        select case (whichaxes)
            case ('110')
                axes = reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0], [3, 3]) / sqrt(2.0d0)
            case ('111')
                axes = reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0], [3, 3]) / sqrt(3.0d0)
            case ('100')
                axes = reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0], [3, 3])
            case default
                print *, "Error: 未定义的轴集"
                stop
        end select

        ! 进行进一步的计算
        naxes = size(axes, 2)
        period = pi * naxes / 6.0d0

        ! 初始化数组
        distances = 0.0d0
        ksis = 0.0d0
        etas = 0.0d0
        phis = 0.0d0
        thisindex = 0

        ! 进行主要计算循环
        do j = 1, 12
            ! 在这里加入基于 P, Q 等的计算逻辑
        end do

        ! 输出几何计算结果
        geom = reshape([distances, ksis, etas, phis], [size(distances), 4])


    end function distances_to_set

end module GB5DOF_utils

