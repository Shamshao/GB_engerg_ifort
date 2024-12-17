! Fortran implementation of the GB5DOF energy computation
module GB5DOF_module
  implicit none
  private
  public :: GB5DOF, distances_to_set, mat2quat, quat2mat, makeparvec, weightedmeanenergy,rsw
  real(8), parameter :: pi = 3.141592653589793d0  ! Declare pi as a constant
contains
   function cross(a, b) result(c)
   implicit none
        real(8), dimension(3), intent(in) :: a(3), b(3)
        real(8), dimension(3) :: c
        ! Cross product calculation
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
       end function cross
    real(8) function rsw(theta, theta1, theta2, a)
    implicit none
    real(8), intent(in) :: theta, theta1, theta2, a
    real(8) :: dtheta, sins, xlogx,theta_norm
    intrinsic :: sin, log
    ! Calculate the interval of angles
    dtheta = theta2 - theta1
    ! Normalize the angle
    theta_norm = (theta - theta1) / dtheta * (3.141592653589793d0 / 2.0d0)

    ! Compute sin(theta)
    sins = sin(theta_norm)

    ! Handle small values to avoid numerical issues
    if (sins > 1.0d-6) then
        xlogx = sins * log(sins)
    else
        xlogx = 0.0d0
    end if

    ! Compute the RSW function value
    rsw = sins - a * xlogx

end function rsw
  function GB5DOF(P, Q, AlCuParameter, eRGB) result(en)
    implicit none
    real(8), intent(in) :: P(3, 3), Q(3, 3)
    real(8), intent(in), optional :: eRGB
    character(len=*), intent(in) :: AlCuParameter
    real(8) :: en

    ! Local variables
    real(8), allocatable :: geom100(:,:), geom110(:,:), geom111(:,:)
    real(8), allocatable :: parvec(:)

    print *, "Calculating geometry parameters for 100 set"
    geom100 = distances_to_set(P, Q, '100')
    print *, "Calculating geometry parameters for 110 set"
    geom110 = distances_to_set(P, Q, '110')
    print *, "Calculating geometry parameters for 111 set"
    geom111 = distances_to_set(P, Q, '111')

    ! Generate parameter vector
    if (present(eRGB)) then
      print *, "eRGB provided, generating parameter vector with custom values"
      parvec = makeparvec(AlCuParameter, eRGB)
    else
      print *, "Using default eRGB, generating parameter vector"
      parvec = makeparvec(AlCuParameter)
    end if

    ! Calculate energy
    print *, "Calculating weighted mean energy"
    en = weightedmeanenergy(geom100, geom110, geom111, parvec)
    print *, "Energy calculation complete: ", en

    return
  end function GB5DOF
  function distances_to_set(P, Q, whichaxes, dismax) result(geom)
    real(8), intent(in) :: P(3, 3), Q(3, 3)
    character(len=*), intent(in) :: whichaxes
    real(8), intent(in), optional :: dismax
    real(8), allocatable :: geom(:,:)

    ! Local variables
    real(8), dimension(:,:), allocatable :: axes, dirs
    integer :: i, naxes, thisindex, j
    real(8) :: ax(3), dir(3), dir2(3), dis, theta, phi, ksi, eta, period
    real(8), dimension(24, 3, 3) :: V
    real(8), dimension(:), allocatable :: distances, ksis, etas, phis
    real(8) :: rotX90(3, 3), rotY90(3, 3), rotZ90(3, 3), rotZ90m(3, 3)
    real(8) :: local_dismax  ! Local variable for dismax

    ! Default dismax value
    if (.not. present(dismax)) then
        local_dismax = 0.999999d0  ! Default value if dismax is not passed
    else
        local_dismax = dismax  ! Use provided value of dismax
    end if
    ! Define axes and directions based on whichaxes
    select case (whichaxes)
    case ('100')
      allocate(axes(3, 3), dirs(3, 3))
      axes = reshape([1.0d0, 0.0d0, 0.0d0, &
                      0.0d0, 1.0d0, 0.0d0, &
                      0.0d0, 0.0d0, 1.0d0], [3, 3])
      dirs = reshape([0.0d0, 0.0d0, 1.0d0, &
                      1.0d0, 0.0d0, 0.0d0, &
                      0.0d0, 1.0d0, 0.0d0], [3, 3])
    case ('110')
      allocate(axes(3, 6), dirs(3, 6))
      axes = reshape([1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, &
                      1.0d0, -1.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, &
                      0.0d0, 0.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0], [3, 6]) / sqrt(2.0d0)
      dirs = reshape([0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, &
                      0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, &
                      1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], [3, 6])
    case ('111')
      allocate(axes(3, 4), dirs(3, 4))
      axes = reshape([1.0d0, 1.0d0, -1.0d0, -1.0d0, &
                      1.0d0, -1.0d0, 1.0d0, -1.0d0, &
                      1.0d0, -1.0d0, -1.0d0, 1.0d0], [3, 4]) / sqrt(3.0d0)
      dirs = reshape([1.0d0, 1.0d0, 1.0d0, 1.0d0, &
                      -1.0d0, 1.0d0, 1.0d0, -1.0d0, &
                      0.0d0, 0.0d0, 0.0d0, 0.0d0], [3, 4]) / sqrt(2.0d0)
    case default
      stop "Undefined axis set"
    end select

    print *, "Defined axes and directions for ", whichaxes

    naxes = size(axes, 2)
    period = pi * naxes / 6.0d0

    ! Define rotation matrices
    rotX90 = reshape([1.0d0, 0.0d0, 0.0d0, &
                      0.0d0, 0.0d0, -1.0d0, &
                      0.0d0, 1.0d0, 0.0d0], [3, 3])
    rotY90 = reshape([0.0d0, 0.0d0, 1.0d0, &
                      0.0d0, 1.0d0, 0.0d0, &
                     -1.0d0, 0.0d0, 0.0d0], [3, 3])
    rotZ90 = reshape([0.0d0, -1.0d0, 0.0d0, &
                      1.0d0, 0.0d0, 0.0d0, &
                      0.0d0, 0.0d0, 1.0d0], [3, 3])
    rotZ90m = reshape([0.0d0, 1.0d0, 0.0d0, &
                      -1.0d0, 0.0d0, 0.0d0, &
                      0.0d0, 0.0d0, 1.0d0], [3, 3])

    ! Generate 24 symmetry-related variants of Q
    V(1, :, :) = Q
    V(2, :, :) = matmul(V(1, :, :), rotX90)
    V(3, :, :) = matmul(V(2, :, :), rotX90)
    V(4, :, :) = matmul(V(3, :, :), rotX90)

    do j = 1, 12
      V(j + 4, :, :) = matmul(V(j, :, :), rotY90)
    end do

    do j = 1, 4
      V(j + 16, :, :) = matmul(V(j, :, :), rotZ90)
      V(j + 20, :, :) = matmul(V(j, :, :), rotZ90m)
    end do

    print *, "Symmetry variants generated"

    ! Preallocate parameter arrays
    allocate(distances(24 * naxes), ksis(24 * naxes), etas(24 * naxes), phis(24 * naxes))
    distances = 0.0d0; ksis = 0.0d0; etas = 0.0d0; phis = 0.0d0
    thisindex = 0

    ! Iterate through axes and coset elements
    do i = 1, naxes
      ax = axes(:, i)
      dir = dirs(:, i)
      dir2 = cross(ax, dir)

      do j = 1, 24
        print *, "Processing axis: ", i, " symmetry variant: ", j
        call compute_distance(V(j, :, :), ax, dir, dir2, dismax, distances, ksis, etas, phis, thisindex, period)
      end do
    end do

    ! Sort and eliminate redundant values
    print *, "Sorting and cleaning results"
    call sort_and_clean(distances, ksis, etas, phis, thisindex)
    allocate(geom(4, thisindex))
    geom(1, :) = distances(:)
    geom(2, :) = ksis(:)
    geom(3, :) = etas(:)
    geom(4, :) = phis(:)
    !geom = transpose([distances(:), ksis(:), etas(:), phis(:)])

    print *, "Distances_to_set computation completed"

    return
  end function distances_to_set

  function mat2quat(m) result(q)
    implicit none
    real(8), intent(in) :: m(3, 3)
    real(8) :: q(4)

    real(8) :: t, e0, e1, e2, e3
    t = m(1, 1) + m(2, 2) + m(3, 3)
    e0 = sqrt(1.0d0 + t) / 2.0d0

    print *, "Converting matrix to quaternion"

    if (t > -0.999999999d0) then
      q(1) = e0
      q(2) = (m(2, 3) - m(3, 2)) / (4.0d0 * e0)
      q(3) = (m(3, 1) - m(1, 3)) / (4.0d0 * e0)
      q(4) = (m(1, 2) - m(2, 1)) / (4.0d0 * e0)
    else
      e0 = 0.0d0
      e3 = sqrt(-(m(1, 1) + m(2, 2)) / 2.0d0)
      if (abs(e3) > 2.0d-8) then
        q = [0.0d0, m(1, 3) / (2.0d0 * e3), m(2, 3) / (2.0d0 * e3), e3]
      else
        e1 = sqrt((m(1, 1) + 1.0d0) / 2.0d0)
        if (e1 /= 0.0d0) then
          q = [0.0d0, e1, m(2, 1) / (2.0d0 * e1), 0.0d0]
        else
          q = [0.0d0, 0.0d0, 1.0d0, 0.0d0]
        end if
      end if
    end if

    print *, "Quaternion: ", q

    return
  end function mat2quat

  function quat2mat(q) result(m)
    implicit none
    real(8), intent(in) :: q(4)
    real(8) :: m(3, 3)

    real(8) :: e0, e1, e2, e3, norm_q
    e0 = q(1)
    e1 = q(2)
    e2 = q(3)
    e3 = q(4)

    print *, "Converting quaternion to matrix"

    norm_q = e0**2 + e1**2 + e2**2 + e3**2

    m = reshape([e0**2 + e1**2 - e2**2 - e3**2, 2.0d0 * (e1 * e2 - e0 * e3), 2.0d0 * (e1 * e3 + e0 * e2), &
                     2.0d0 * (e1 * e2 + e0 * e3), e0**2 - e1**2 + e2**2 - e3**2, 2.0d0 * (e2 * e3 - e0 * e1), &
                     2.0d0 * (e1 * e3 - e0 * e2), 2.0d0 * (e2 * e3 + e0 * e1), e0**2 - e1**2 - e2**2 + e3**2], [3, 3])
        m = m / norm_q
    !end 


    return
  end function quat2mat

  function makeparvec(AlCuParameter, eRGB, par42Al, par42Cu) result(par43)
    implicit none
    character(len=*), intent(in) :: AlCuParameter
    real(8), intent(in), optional :: eRGB
    real(8), intent(in), optional :: par42Al(:), par42Cu(:)
    real(8), allocatable :: par43(:)

    real(8) :: default_par42Al(42), default_par42Cu(42), default_eRGB
    real(8) :: AlCuParamValue

    ! Default values for parameters
    default_par42Al = [0.405204179289160d0, 0.738862004021890d0, 0.351631012630026d0, 2.40065811939667d0, &
                        1.34694439281655d0, 0.352260396651516d0, 0.602137375062785d0, 1.58082498976078d0, &
                        0.596442399566661d0, 1.30981422643602d0, 3.21443408257354d0, 0.893016409093743d0, &
                        0.835332505166333d0, 0.933176738717594d0, 0.896076948651935d0, 0.775053293192055d0, &
                        0.391719619979054d0, 0.782601780600192d0, 0.678572601273508d0, 1.14716256515278d0, &
                        0.529386201144101d0, 0.909044736601838d0, 0.664018011430602d0, 0.597206897283586d0, &
                        0.200371750006251d0, 0.826325891814124d0, 0.111228512469435d0, 0.664039563157148d0, &
                        0.241537262980083d0, 0.736315075146365d0, 0.514591177241156d0, 1.73804335876546d0, &
                        3.04687038671309d0, 1.48989831680317d0, 0.664965104218438d0, 0.495035051289975d0, &
                        0.495402996460658d0, 0.468878130180681d0, 0.836548944799803d0, 0.619285521065571d0, &
                        0.844685390948170d0, 1.02295427618256d0]

    default_par42Cu = [0.405204179289160d0, 0.738862004021890d0, 0.351631012630026d0, 2.40065811939667d0, &
                        1.34694439281655d0, 3.37892632736175d0, 0.602137375062785d0, 1.58082498976078d0, &
                        0.710489498577995d0, 0.737834049784765d0, 3.21443408257354d0, 0.893016409093743d0, &
                        0.835332505166333d0, 0.933176738717594d0, 0.896076948651935d0, 0.775053293192055d0, &
                        0.509781056492307d0, 0.782601780600192d0, 0.762160812499734d0, 1.10473084066580d0, &
                        0.529386201144101d0, 0.909044736601838d0, 0.664018011430602d0, 0.597206897283586d0, &
                        0.200371750006251d0, 0.826325891814124d0, 0.0226010533470218d0, 0.664039563157148d0, &
                        0.297920289861751d0, 0.666383447163744d0, 0.514591177241156d0, 1.73804335876546d0, &
                        2.69805148576400d0, 1.95956771207484d0, 0.948894352912787d0, 0.495035051289975d0, &
                        0.301975031994664d0, 0.574050577702240d0, 0.836548944799803d0, 0.619285521065571d0, &
                        0.844685390948170d0, 0.0491040633104212d0]

    default_eRGB = 1.03669431227427d0

    ! Assign default values if not provided
    if (.not. present(eRGB)) eRGB = default_eRGB
    if (.not. present(par42Al)) allocate(par42Al(42)); par42Al = default_par42Al
    if (.not. present(par42Cu)) allocate(par42Cu(42)); par42Cu = default_par42Cu

    ! Parse AlCuParameter
    select case (AlCuParameter)
    case ('Ni')
      eRGB = 1.44532834613925d0
      AlCuParamValue = 0.767911805073948d0
    case ('Al')
      eRGB = 0.547128733614891d0
      AlCuParamValue = 0.0d0
    case ('Au')
      eRGB = 0.529912885175204d0
      AlCuParamValue = 0.784289766313152d0
    case ('Cu')
      eRGB = default_eRGB
      AlCuParamValue = 1.0d0
    case default
      stop "Undefined element"
    end select

    ! Compute the 43-parameter vector
    allocate(par43(43))
    par43(1) = eRGB
    par43(2:43) = par42Al + AlCuParamValue * (par42Cu - par42Al)

    print *, "Generated 43-parameter vector"

    return
  end function makeparvec
  function weightedmeanenergy(geom100, geom110, geom111, pars) result(en)
    implicit none
    real(8), intent(in) :: geom100(:,:), geom110(:,:), geom111(:,:)
    real(8), intent(in) :: pars(:)
    real(8) :: en

    ! Local variables
    real(8) :: eRGB, d0100, d0110, d0111
    real(8) :: weight100, weight110, weight111
    real(8), allocatable :: e100(:), e110(:), e111(:)
    real(8), allocatable :: d100(:), d110(:), d111(:)
    real(8), allocatable :: w100(:), w110(:), w111(:)
    real(8), allocatable :: s100(:), s110(:), s111(:)
    real(8) :: offset

    print *, "Starting weighted mean energy calculation"

    ! Extract relevant parameters from pars
    eRGB = pars(1)
    d0100 = pars(2)
    d0110 = pars(3)
    d0111 = pars(4)
    weight100 = pars(5)
    weight110 = pars(6)
    weight111 = pars(7)

    ! Offset for numerical stability
    offset = 0.00001d0

    ! Calculate energies for each set
    e100 = set100(geom100, pars)
    e110 = set110(geom110, pars)
    e111 = set111(geom111, pars)

    ! Extract distances
    d100 = geom100(1, :)
    d110 = geom110(1, :)
    d111 = geom111(1, :)

    ! Calculate weights
    allocate(s100(size(d100)), s110(size(d110)), s111(size(d111)))
    allocate(w100(size(d100)), w110(size(d110)), w111(size(d111)))

    s100 = sin(pi / 2.0d0 * d100 / d0100)
    s100 = merge(s100, 1.0d0, d100 > d0100)
    s100 = merge(s100, offset * pi / 2.0d0, d100 < offset * d0100)
    w100 = (1.0d0 / (s100 * (1.0d0 - 0.5d0 * log(s100))) - 1.0d0) * weight100

    s110 = sin(pi / 2.0d0 * d110 / d0110)
    s110 = merge(s110, 1.0d0, d110 > d0110)
    s110 = merge(s110, offset * pi / 2.0d0, d110 < offset * d0110)
    w110 = (1.0d0 / (s110 * (1.0d0 - 0.5d0 * log(s110))) - 1.0d0) * weight110

    s111 = sin(pi / 2.0d0 * d111 / d0111)
    s111 = merge(s111, 1.0d0, d111 > d0111)
    s111 = merge(s111, offset * pi / 2.0d0, d111 < offset * d0111)
    w111 = (1.0d0 / (s111 * (1.0d0 - 0.5d0 * log(s111))) - 1.0d0) * weight111

    ! Calculate weighted mean energy
    en = eRGB * (sum(e100 * w100) + sum(e110 * w110) + sum(e111 * w111) + 1.0d0) / &
         (sum(w100) + sum(w110) + sum(w111) + 1.0d0)

    print *, "Weighted mean energy calculated: ", en

    ! Deallocate arrays
    deallocate(s100, s110, s111, w100, w110, w111)

    return
  end function weightedmeanenergy
  !Fortran implementation of set100 and related functions

function set100(geom100, pars,n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: geom100(4, n), pars(43)
    real(8) :: en(n)
    real(8) :: ksi(n), eta(n), phi(n), entwist(n), entilt(n), x(n)
    integer :: i

    ! Extract geometrical parameters
    ksi = geom100(2, :)
    eta = geom100(3, :)
    phi = geom100(4, :)

    ! Calculate twist and tilt contributions
    entwist=twist100(ksi, pars, n)
    entilt=atgb100(eta, ksi, pars, n)

    ! Calculate combined contribution
    x = phi / (3.141592653589793d0 / 2.0d0)

    do i = 1, n
        en(i) = entwist(i) * (1.0d0 - x(i))**pars(8) + entilt(i) * x(i)**pars(9)
    end do

end function set100

function twist100(ksi, pars,n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: ksi(n), pars(43)
    real(8) :: en(n)
    real(8) :: a, b, perio, sins(n), xlogx(n)
    integer :: i

    ! Extract parameters
    a = pars(10)
    b = pars(10) * pars(11)
    perio = 3.141592653589793d0 / 2.0d0

    ! Apply rotation symmetry
    do i = 1, n
        ksi(i) = mod(abs(ksi(i)), perio)
        if (ksi(i) > perio / 2.0d0) then
            ksi(i) = perio - ksi(i)
        end if
    end do

    ! Implement rsw function of ksi
    do i = 1, n
        sins(i) = sin(2.0d0 * ksi(i))
        if (sins(i) > 1.0d-6) then
            xlogx(i) = sins(i) * log(sins(i))
        else
            xlogx(i) = 0.0d0
        end if
        en(i) = a * sins(i) - b * xlogx(i)
    end do

end function twist100

function atgb100(eta, ksi, pars,n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: eta(n), ksi(n), pars(43)
    real(8) :: en(n)
    real(8) :: en1(n), en2(n), pwr, period
    integer :: i

    ! Extract parameters
    pwr = pars(12)
    period = 3.141592653589793d0 / 2.0d0

    ! Symmetric tilt boundary energies
    en1=stgb100(ksi, pars,n)
    en2=stgb100(period - ksi, pars,n)

    ! Interpolate based on eta
    do i = 1, n
        if (en1(i) >= en2(i)) then
            en(i) = en1(i) - (en1(i) - en2(i)) * (eta(i) / period)**pwr
        else
            en(i) = en2(i) - (en2(i) - en1(i)) * (1.0d0 - eta(i) / period)**pwr
        end if
    end do

end function atgb100
function stgb100(ksi, pars, n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: ksi(n), pars(43)
    real(8) :: en(n)
    real(8) :: th1, th2, th3, th4, th5, th6, th7
    real(8) :: en1, en2, en3, en4, en5, en6, en7
    real(8) :: a12, a23, a34, a45, a56, a67
    integer :: i

    ! Define piecewise energy levels
    en1 = 0.0d0
    en2 = pars(13)
    en3 = pars(14)
    en4 = pars(15)
    en5 = pars(16)
    en6 = pars(17)
    en7 = 0.0d0

    ! Define angle breaks
    th1 = 0.0d0
    th2 = pars(18)
    th3 = acos(4.0d0 / 5.0d0)
    th4 = pars(19)
    th5 = acos(3.0d0 / 5.0d0)
    th6 = 2.0d0 * acos(5.0d0 / sqrt(34.0d0))
    th7 = 3.141592653589793d0 / 2.0d0

    ! Define rsw shape factors
    a12 = 0.5d0
    a23 = 0.5d0
    a34 = 0.5d0
    a45 = 0.5d0
    a56 = 0.5d0
    a67 = 0.5d0

    ! Piecewise rsw function
    do i = 1, n
        if (ksi(i) <= th2) then
            en(i) = en1 + (en2 - en1) * rsw(ksi(i), th1, th2, a12)
        else if (ksi(i) <= th3) then
            en(i) = en3 + (en2 - en3) * rsw(ksi(i), th3, th2, a23)
        else if (ksi(i) <= th4) then
            en(i) = en3 + (en4 - en3) * rsw(ksi(i), th3, th4, a34)
        else if (ksi(i) <= th5) then
            en(i) = en5 + (en4 - en5) * rsw(ksi(i), th5, th4, a45)
        else if (ksi(i) <= th6) then
            en(i) = en5 + (en6 - en5) * rsw(ksi(i), th5, th6, a56)
        else
            en(i) = en7 + (en6 - en7) * rsw(ksi(i), th7, th6, a67)
        end if
    end do

end function stgb100

! Fortran implementation of set110 and related functions

function set110(geom110, pars,n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: geom110(4, n), pars(43)
    real(8) :: en(n)
    real(8) :: ksi(n), eta(n), phi(n), entwist(n), entilt(n), x(n)
    integer :: i

    ksi = geom110(2, :)
    eta = geom110(3, :)
    phi = geom110(4, :)

    entwist=twists110(ksi, pars, n)
    entilt=atgbs110(eta, ksi, pars,n)

    x = phi / (0.5d0 * 3.141592653589793d0)

    do i = 1, n
        en(i) = entwist(i) * (1.0d0 - x(i))**pars(20) + entilt(i) * x(i)**pars(21)
    end do

end function set110

function atgbs110(eta, ksi, pars, n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: eta(n), ksi(n), pars(43)
    real(8) :: en(n)
    real(8) :: period, en1(n), en2(n), a
    integer :: i

    period = 3.141592653589793d0
    a = pars(26)

    en1=stgbs110(ksi, pars, en1, n)
    en2=stgbs110(period - ksi, pars, en2, n)

    do i = 1, n
        if (en1(i) >= en2(i)) then
            en(i) = en2(i) + (en1(i) - en2(i)) * rsw(eta(i), period, 0.0d0, a)
        else
            en(i) = en1(i) + (en2(i) - en1(i)) * rsw(eta(i), 0.0d0, period, a)
        end if
    end do

end function atgbs110

function stgbs110(th, pars, n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: th(n), pars(43)
    real(8) :: en(n)
    real(8) :: th1, th2, th3, th4, th5, th6, th7
    real(8) :: en1, en2, en3, en4, en5, en6, en7
    real(8) :: a12, a23, a34, a45, a56, a67
    real(8) :: th_temp(n)
    integer :: i

    th1 = 0.0d0
    th3 = acos(1.0d0 / 3.0d0)
    th5 = acos(-7.0d0 / 11.0d0)
    th7 = 3.141592653589793d0

    th2 = pars(32)
    th4 = pars(33)
    th6 = pars(34)

    en1 = 0.0d0
    en2 = pars(27)
    en3 = pars(28)
    en4 = pars(29)
    en5 = pars(30)
    en6 = pars(31)
    en7 = 0.0d0

    a12 = 0.5d0
    a23 = 0.5d0
    a34 = 0.5d0
    a45 = 0.5d0
    a56 = 0.5d0
    a67 = 0.5d0

    do i = 1, n
        th_temp(i) = 3.141592653589793d0 - th(i)
    end do

    do i = 1, n
        if (th_temp(i) <= th2) then
            en(i) = en1 + (en2 - en1) * rsw(th_temp(i), th1, th2, a12)
        else if (th_temp(i) <= th3) then
            en(i) = en3 + (en2 - en3) * rsw(th_temp(i), th3, th2, a23)
        else if (th_temp(i) <= th4) then
            en(i) = en3 + (en4 - en3) * rsw(th_temp(i), th3, th4, a34)
        else if (th_temp(i) <= th5) then
            en(i) = en5 + (en4 - en5) * rsw(th_temp(i), th5, th4, a45)
        else if (th_temp(i) <= th6) then
            en(i) = en5 + (en6 - en5) * rsw(th_temp(i), th5, th6, a56)
        else
            en(i) = en7 + (en6 - en7) * rsw(th_temp(i), th7, th6, a67)
        end if
    end do

end function stgbs110

function twists110(th, pars,n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: th(n), pars(43)
    real(8) :: en(n)
    real(8) :: th1, th2, th3, en1, en2, en3, a01, a12, a23
    real(8) :: th_mod(n)
    integer :: i

    th1 = pars(22)
    en1 = pars(23)
    en2 = pars(24)
    en3 = pars(25)

    a01 = 0.5d0
    a12 = 0.5d0
    a23 = 0.5d0

    th2 = acos(1.0d0 / 3.0d0)
    th3 = 3.141592653589793d0 / 2.0d0

    do i = 1, n
        th_mod(i) = mod(abs(th(i)), 3.141592653589793d0)
        if (th_mod(i) > 3.141592653589793d0 / 2.0d0) then
            th_mod(i) = 3.141592653589793d0 - th_mod(i)
        end if
    end do
      do i = 1, n
        if (th_mod(i) <= th1) then
            en(i) = en1 * rsw(th_mod(i), 0.0d0, th1, a01)
        else if (th_mod(i) <= th2) then
            en(i) = en2 + (en1 - en2) * rsw(th_mod(i), th2, th1, a12)
        else
            en(i) = en3 + (en2 - en3) * rsw(th_mod(i), th3, th2, a23)
        end if
    end do

end function twists110
function set111(geom111, pars,n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: geom111(4, n), pars(43)
    real(8) :: en(n)
    real(8) :: ksi(n), eta(n), phi(n), entwist(n), entilt(n), x(n), a, b
    integer :: i

    a = pars(35)
    b = a - 1.0d0

    ksi = geom111(2, :)
    eta = geom111(3, :)
    phi = geom111(4, :)

    entwist=twists111(ksi, pars,n)
    enatgb=atgbs111(eta, ksi, pars,n)

    x = phi / (0.5d0 * 3.141592653589793d0)

    do i = 1, n
        en(i) = entwist(i) + (entilt(i) - entwist(i)) * (a * x(i) - b * x(i)**2)
    end do

end function set111

function twists111(theta, pars,n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: theta(n), pars(43)
    real(8) :: en(n)
    real(8) :: thd, enm, en2, a1, a2, theta_mod(n)
    integer :: i

    thd = pars(37)
    enm = pars(38)
    en2 = pars(28)
    a1 = pars(36)
    a2 = a1

    do i = 1, n
        if (theta(i) > 3.141592653589793d0 / 3.0d0) then
            theta_mod(i) = 2.0d0 * 3.141592653589793d0 / 3.0d0 - theta(i)
        else
            theta_mod(i) = theta(i)
        end if
    end do

    do i = 1, n
        if (theta_mod(i) <= thd) then
            en(i) = enm * rsw(theta_mod(i), 0.0d0, thd, a1)
        else
            en(i) = en2 + (enm - en2) * rsw(theta_mod(i), 3.141592653589793d0 / 3.0d0, thd, a2)
        end if
    end do

end function twists111

function atgbs111(eta, ksi, pars,n) result(en)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: eta(n), ksi(n), pars(43)
    real(8) :: en(n)
    real(8) :: ksim, enmax, enmin, encnt, a1, a2, etascale, chi(n)
    integer :: i

    ksim = pars(39)
    enmax = pars(40)
    enmin = pars(41)
    encnt = pars(42)
    a1 = 0.5d0
    a2 = 0.5d0
    etascale = pars(43)

    do i = 1, n
        if (ksi(i) > 3.141592653589793d0 / 3.0d0) then
            ksi(i) = 2.0d0 * 3.141592653589793d0 / 3.0d0 - ksi(i)
        end if
        if (eta(i) > 3.141592653589793d0 / 3.0d0) then
            eta(i) = 2.0d0 * 3.141592653589793d0 / 3.0d0 - eta(i)
        end if
    end do

    do i = 1, n
        if (ksi(i) <= ksim) then
            en(i) = enmax * rsw(ksi(i), 0.0d0, ksim, a1)
        else
            chi(i) = enmin + (encnt - enmin) * rsw(eta(i), 0.0d0, 3.141592653589793d0 / (2.0d0 * etascale), 0.5d0)
            en(i) = chi(i) + (enmax - chi(i)) * rsw(ksi(i), 3.141592653589793d0 / 3.0d0, ksim, a2)
        end if
    end do

end function atgbs111
end module GB5DOF_module
