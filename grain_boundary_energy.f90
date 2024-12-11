module weightedmeanenergy_module                             
        implicit none
        real(8), parameter :: pi = 3.141592653589793d0  ! Define pi
contains

        function weightedmeanenergy(geom100, geom110, geom111, pars) result(en)
                ! Inputs
                real(8), intent(in), dimension(4, *) :: geom100, geom110, geom111
                real(8), intent(in), dimension(43) :: pars
                ! Output
                real(8) :: en
                ! Local variables
                real(8) :: eRGB, d0100, d0110, d0111
                real(8) :: weight100, weight110, weight111, offset
                real(8), dimension(:), allocatable :: e100, e110, e111
                real(8), dimension(:), allocatable :: d100, d110, d111
                real(8), dimension(:), allocatable :: s100, s110, s111
                real(8), dimension(:), allocatable :: w100, w110, w111
                integer :: i, n

                ! Extract parameters
                eRGB = pars(1)
                d0100 = pars(2)
                d0110 = pars(3)
                d0111 = pars(4)
                weight100 = pars(5)
                weight110 = pars(6)
                weight111 = pars(7)
                offset = 1.0d-5

                ! Calculate energy contributions from each set
                n = size(geom100, 2)
                allocate(e100(n), e110(n), e111(n))
                allocate(d100(n), d110(n), d111(n))
                allocate(s100(n), s110(n), s111(n))
                allocate(w100(n), w110(n), w111(n))

                e100 = set100(geom100, pars)
                e110 = set110(geom110, pars)
                e111 = set111(geom111, pars)

                d100 = geom100(1, :)
                d110 = geom110(1, :)
                d111 = geom111(1, :)

                ! Calculate weights
                do i = 1, n
                s100(i) = sin(pi/2 * d100(i) / d0100)
                if (d100(i) > d0100) then
                        s100(i) = 1.0d0
                else if (d100(i) < offset * d0100) then
                        s100(i) = offset * pi / 2
                end if
                w100(i) = (1.0d0 / (s100(i) * (1.0d0 - 0.5d0 * log(s100(i)))) - 1.0d0) * weight100
                end do

                do i = 1, n
                s110(i) = sin(pi/2 * d110(i) / d0110)
                if (d110(i) > d0110) then
                        s110(i) = 1.0d0
                else if (d110(i) < offset * d0110) then
                        s110(i) = offset * pi / 2
                end if
                w110(i) = (1.0d0 / (s110(i) * (1.0d0 - 0.5d0 * log(s110(i)))) - 1.0d0) * weight110
                end do

                do i = 1, n
                s111(i) = sin(pi/2 * d111(i) / d0111)
                if (d111(i) > d0111) then
                        s111(i) = 1.0d0
                else if (d111(i) < offset * d0111) then
                        s111(i) = offset * pi / 2
                end if
                w111(i) = (1.0d0 / (s111(i) * (1.0d0 - 0.5d0 * log(s111(i)))) - 1.0d0) * weight111
                end do

                ! Final energy calculation
                en = eRGB * (sum(e100 * w100) + sum(e110 * w110) + sum(e111 * w111) + 1.0d0) / &
                        (sum(w100) + sum(w110) + sum(w111) + 1.0d0)

                ! Deallocate local arrays
                deallocate(e100, e110, e111)
                deallocate(d100, d110, d111)

        end function weightedmeanenergy

        function set100(geom100, pars) result(en)
                ! Inputs
                real(8), intent(in), dimension(4, *) :: geom100
                real(8), intent(in), dimension(43) :: pars

                ! Output
                real(8), dimension(:), allocatable :: en    ! Corrected: en is allocatable
                real(8), dimension(:), allocatable :: ksi, eta, phi
                real(8), dimension(:), allocatable :: entwist, entilt, x

                ! Local variables
                real(8) :: pwr1, pwr2
                integer :: i, n

                ! Extract parameters
                pwr1 = pars(8)
                pwr2 = pars(9)

                ! Get the size of the second dimension of geom100
                n = size(geom100, 2)

                ! Allocate all arrays in a single statement
                allocate(ksi(n), eta(n), phi(n), entwist(n), entilt(n), x(n), en(n))

                ! Assign values from geom100 to the arrays
                ksi = geom100(2, :)
                eta = geom100(3, :)
                phi = geom100(4, :)

                ! Calculate twist and tilt contributions
                entwist = twist100(ksi, pars)
                entilt = atgb100(eta, ksi, pars)

                ! Combine twist and tilt contributions into en
                do i = 1, n
                x(i) = phi(i) / (pi / 2.0d0)
                en(i) = entwist(i) * (1.0d0 - x(i))**pwr1 + entilt(i) * x(i)**pwr2
                end do

                ! Deallocate local arrays
                deallocate(ksi, eta, phi)
                deallocate(entwist, entilt, x)

        end function set100

end module weightedmeanenergy_module
!  另一个模块开始


module grain_boundary_energy
    implicit none
    private
    public :: set100, twist100, atgb100, stgb100

contains

    ! 主函数：计算 set100 的能量
    function set100(geom100, pars) result(en)
        real(8), intent(in), dimension(4, *) :: geom100
        real(8), intent(in), dimension(:) :: pars
        real(8), dimension(size(geom100, 2)) :: en

        real(8) :: pwr1, pwr2
        real(8), dimension(size(geom100, 2)) :: ksi, eta, phi
        real(8), dimension(size(geom100, 2)) :: entwist, entilt, x

        ! 提取参数
        pwr1 = pars(8)
        pwr2 = pars(9)
        ksi = geom100(2, :)
        eta = geom100(3, :)
        phi = geom100(4, :)

        ! 计算 Twist 和 Tilt 能量
        entwist = twist100(ksi, pars)
        entilt = atgb100(eta, ksi, pars)

        ! 结合 phi 权重计算总能量
        x = phi / (pi / 2.0d0)
        en = entwist * (1.0d0 - x)**pwr1 + entilt * x**pwr2
    end function set100

    ! 100 Twist 能量函数
    function twist100(ksi, pars) result(en)
        real(8), intent(in), dimension(:) :: ksi
        real(8), intent(in), dimension(:) :: pars
        real(8), dimension(size(ksi)) :: en

        real(8) :: a, b, perio
        real(8), dimension(size(ksi)) :: sins, xlogx

        ! 提取参数
        a = pars(10)
        b = pars(10) * pars(11)
        perio = pi / 2.0d0

        ! 处理对称性
        ksi = mod(abs(ksi), perio)
        ksi = min(ksi, perio - ksi)

        ! 计算能量
        sins = sin(2.0d0 * ksi)
        xlogx = sins * log(max(sins, 1.0d-20))  ! 避免 log(0) 的问题
        en = a * sins - b * xlogx
    end function twist100

    ! 100 Tilt 能量函数
    function atgb100(eta, ksi, pars) result(en)
        real(8), intent(in), dimension(:) :: eta, ksi
        real(8), intent(in), dimension(:) :: pars
        real(8), dimension(size(ksi)) :: en

        real(8) :: pwr, period
        real(8), dimension(size(ksi)) :: en1, en2
        logical, dimension(size(ksi)) :: select

        pwr = pars(12)
        period = pi / 2.0d0

        ! 计算两端的 Tilt 能量
        en1 = stgb100(ksi, pars)
        en2 = stgb100(period - ksi, pars)

        ! 根据 eta 插值 Tilt 能量
        select = en1 >= en2
        en = 0.0d0
        en(select) = en1(select) - (en1(select) - en2(select)) * (eta(select) / period)**pwr
        en(.not. select) = en2(.not. select) - (en2(.not. select) - en1(.not. select)) * &
                           (1.0d0 - eta(.not. select) / period)**pwr
    end function atgb100

    ! 100 Symmetric Tilt 能量函数
    function stgb100(ksi, pars) result(en)
        real(8), intent(in), dimension(:) :: ksi
        real(8), intent(in), dimension(:) :: pars
        real(8), dimension(size(ksi)) :: en

        ! 参数定义
        real(8) :: en1, en2, en3, en4, en5, en6, en7
        real(8) :: th1, th2, th3, th4, th5, th6, th7
        real(8) :: a12, a23, a34, a45, a56, a67
        logical, dimension(size(ksi)) :: select

        en1 = 0.0d0
        en7 = 0.0d0
        en2 = pars(13)
        en3 = pars(14)
        en4 = pars(15)
        en5 = pars(16)
        en6 = pars(17)

        th1 = 0.0d0
        th2 = pars(18)
        th3 = acos(4.0d0 / 5.0d0)
        th4 = pars(19)
        th5 = acos(3.0d0 / 5.0d0)
        th6 = 2.0d0 * acos(5.0d0 / sqrt(34.0d0))
        th7 = pi / 2.0d0

        a12 = 0.5d0
        a23 = a12
        a34 = a12
        a45 = a12
        a56 = a12
        a67 = a12

        ! 分段能量计算
        en = 0.0d0
        select = (ksi <= th2)
        en(select) = en1 + (en2 - en1) * rsw(ksi(select), th1, th2, a12)

        select = (ksi > th2 .and. ksi <= th3)
        en(select) = en3 + (en2 - en3) * rsw(ksi(select), th3, th2, a23)

        select = (ksi > th3 .and. ksi <= th4)
        en(select) = en3 + (en4 - en3) * rsw(ksi(select), th3, th4, a34)

        select = (ksi > th4 .and. ksi <= th5)
        en(select) = en5 + (en4 - en5) * rsw(ksi(select), th5, th4, a45)

        select = (ksi > th5 .and. ksi <= th6)
        en(select) = en6 + (en5 - en6) * rsw(ksi(select), th6, th5, a56)

        select = (ksi > th6 .and. ksi <= th7)
        en(select) = en7 + (en6 - en7) * rsw(ksi(select), th7, th6, a67)
    end function stgb100

    ! 通用 RSW 函数
    function rsw(x, th1, th2, a) result(rs)
        real(8), intent(in) :: x(:), th1, th2, a
        real(8), dimension(size(x)) :: rs
        rs = sin((x - th1) / (th2 - th1) * pi)**a
    end function rsw

end module grain_boundary_energy
module energy_module
  implicit none
contains

  function set110(geom110, pars) result(en)
    ! Calculate dimensionless contribution to energy from <110> rotations.
    implicit none
    real(dp), intent(in) :: geom110(4, :)   ! Input geometry
    real(dp), intent(in) :: pars(34)        ! Input parameters
    real(dp) :: en(size(geom110, 2))        ! Output energy
    real(dp) :: ksi(size(geom110, 2)), eta(size(geom110, 2)), phi(size(geom110, 2))
    real(dp) :: entwist(size(geom110, 2)), entilt(size(geom110, 2))
    real(dp) :: x
    integer :: i, n

    n = size(geom110, 2)
    ksi = geom110(2, :)
    eta = geom110(3, :)
    phi = geom110(4, :)

    ! Compute twist and tilt contributions
    entwist = twists110(ksi, pars)
    entilt = atgbs110(eta, ksi, pars)

    ! Combine contributions
    do i = 1, n
      x = phi(i) / (pi / 2.0_dp)
      en(i) = entwist(i) * (1.0_dp - x)**pars(20) + entilt(i) * x**pars(21)
    end do
  end function set110

  function atgbs110(eta, ksi, pars) result(en)
    ! Compute 110 asymmetric tilt boundary energy.
    implicit none
    real(dp), intent(in) :: eta(:), ksi(:), pars(:)
    real(dp) :: en(size(eta))
    real(dp) :: en1(size(eta)), en2(size(eta)), period, a
    integer :: i, n
    logical :: select(size(eta))

    n = size(eta)
    period = pi
    a = pars(26)

    ! Calculate tilt energies at the endpoints
    en1 = stgbs110(ksi, pars)
    en2 = stgbs110(period - ksi, pars)

    ! Compute interpolation
    do i = 1, n
      select(i) = en1(i) >= en2(i)
      if (select(i)) then
        en(i) = en2(i) + (en1(i) - en2(i)) * rsw(eta(i), pi, 0.0_dp, a)
      else
        en(i) = en1(i) + (en2(i) - en1(i)) * rsw(eta(i), 0.0_dp, pi, a)
      end if
    end do
  end function atgbs110

  function stgbs110(th, pars) result(en)
    ! Compute 110 symmetric tilt boundary energy.
    implicit none
    real(dp), intent(in) :: th(:), pars(:)
    real(dp) :: en(size(th))
    real(dp) :: en1, en2, en3, en4, en5, en6, en7
    real(dp) :: th1, th2, th3, th4, th5, th6, th7
    real(dp) :: a12, a23, a34, a45, a56, a67
    integer :: i, n

    n = size(th)

    ! Unpack parameters
    en1 = 0.0_dp
    en2 = pars(27)
    en3 = pars(28)
    en4 = pars(29)
    en5 = pars(30)
    en6 = pars(31)
    en7 = 0.0_dp

    th1 = 0.0_dp
    th2 = pars(32)
    th3 = acos(1.0_dp / 3.0_dp)
    th4 = pars(33)
    th5 = acos(-7.0_dp / 11.0_dp)
    th6 = pars(34)
    th7 = pi

    a12 = 0.5_dp; a23 = 0.5_dp; a34 = 0.5_dp
    a45 = 0.5_dp; a56 = 0.5_dp; a67 = 0.5_dp

    ! Adjust symmetry for th
    th = pi - th

    ! Piecewise energy computation
    do i = 1, n
      if (th(i) <= th2) then
        en(i) = en1 + (en2 - en1) * rsw(th(i), th1, th2, a12)
      else if (th(i) <= th3) then
        en(i) = en3 + (en2 - en3) * rsw(th(i), th3, th2, a23)
      else if (th(i) <= th4) then
        en(i) = en3 + (en4 - en3) * rsw(th(i), th3, th4, a34)
      else if (th(i) <= th5) then
        en(i) = en5 + (en4 - en5) * rsw(th(i), th5, th4, a45)
      else if (th(i) <= th6) then
        en(i) = en5 + (en6 - en5) * rsw(th(i), th5, th6, a56)
      else
        en(i) = en7 + (en6 - en7) * rsw(th(i), th7, th6, a67)
      end if
    end do
  end function stgbs110

  function twists110(th, pars) result(en)
    ! Compute 110 twist energy.
    implicit none
    real(dp), intent(in) :: th(:), pars(:)
    real(dp) :: en(size(th))
    real(dp) :: th1, th2, th3
    real(dp) :: en1, en2, en3, a01, a12, a23
    integer :: i, n

    n = size(th)

    th1 = pars(22)
    en1 = pars(23)
    en2 = pars(24)
    en3 = pars(25)
    a01 = 0.5_dp
    a12 = 0.5_dp
    a23 = 0.5_dp

    th2 = acos(1.0_dp / 3.0_dp)
    th3 = pi / 2.0_dp

    ! Handle rotation symmetry
    th = mod(abs(th), pi)
    th(th > pi / 2.0_dp) = pi - th(th > pi / 2.0_dp)

    ! Piecewise energy computation
    do i = 1, n
      if (th(i) <= th1) then
        en(i) = en1 * rsw(th(i), 0.0_dp, th1, a01)
      else if (th(i) <= th2) then
        en(i) = en2 + (en1 - en2) * rsw(th(i), th2, th1, a12)
      else
        en(i) = en3 + (en2 - en3) * rsw(th(i), th3, th2, a23)
      end if
    end do
  end function twists110

  function rsw(x, x1, x2, a) result(rs)
    ! Compute RSW function for interpolation.
    implicit none
    real(dp), intent(in) :: x, x1, x2, a
    real(dp) :: rs
    rs = sin((x - x1) / (x2 - x1) * pi)**a
  end function rsw

end module energy_module
