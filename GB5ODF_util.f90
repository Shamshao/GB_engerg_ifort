module GB5DOF_utils
    implicit none
contains

    ! Function GB5DOF
    ! Function to calculate distances for the given plane direction
      function distances_to_set(P, Q, plane) result(distances)
        implicit none
        real(8), intent(in) :: P(3,3), Q(3,3)
        character(len=*) :: plane
        real(8) ,dimension(3) :: distances(3)

        ! Implement the actual logic for distances_to_set
        ! Placeholder implementation for demonstration
        distances = [0.0d0, 0.0d0, 0.0d0]  ! Modify according to actual logic
    end function distances_to_set

    ! Function to generate parameter vector based on the material type and optional eRGB
    function makeparvec(AlCuParameter, eRGB) result(parvec)
        implicit none
        real(8), intent(in) :: AlCuParameter, eRGB
        real(8), allocatable :: parvec(:)

        ! Allocate and fill the parameter vector
        allocate(parvec(2))
        parvec(1) = AlCuParameter
        parvec(2) = eRGB

        ! Return the parameter vector

    end function makeparvec
  function weightedmeanenergy(geom100, geom110, geom111, parvec) result(energy)
        implicit none
        real(8), intent(in) :: geom100, geom110, geom111
        real(8), intent(in) :: parvec(:)
        real(8) :: energy

        ! Implement the onighted mean energy calculation
        ! Placeholder implementation for demonstration
        energy = (geom100 + geom110 + geom111) / 3.0d0  ! Simplified example

    end function weightedmeanenergy

end module GB5DOF_utils
     !接下来要实现 to distances_to_set主程序 81-256
     module GB5DOF_utils
    implicit none
contains

    function distances_to_set(P, Q, whichaxes, dismax) result(geom)
        implicit none
        real(8), intent(in) :: P(3,3), Q(3,3)         ! Rotation matrices
        character(len=*) :: whichaxes                 ! Axis type ('100', '110', '111')
        real(8), intent(in), optional :: dismax       ! Maximum distance for inclusion
        real(8), dimension(:,:), allocatable :: geom  ! Output geometry parameters

        integer :: naxes, i, j, thisindex
        real(8) :: dis, dismax_local, theta, phi, ksi, eta, theta1, theta2, period
        real(8), dimension(:, :), allocatable :: axes, dirs
        real(8), dimension(24) :: distances, ksis, etas, phis
        real(8), dimension(3, 3) :: rotX90, rotY90, rotZ90, rotZ90m, RA
        real(8), dimension(3, 3) :: V(24)  ! For storing symmetrically rotated variants
        real(8), dimension(3) :: ax, dir, dir2, n1, n2, m1, m2, axi
        real(8) :: psi, dotp
        real(8), parameter :: pi = 3.141592653589793d0
        ! Default dismax value if not provided
        if (.not. present(dismax)) then
            dismax_local = 0.999999
        else
            dismax_local = dismax
        end if

        ! Select the correct axes set
        select case (whichaxes)
        case ('110')
            axes = reshape([ &
                1.0d0, 1.0d0, 1.0d0, &
                1.0d0, -1.0d0, 0.0d0, &
                0.0d0, 0.0d0, 1.0d0, &
                1.0d0, 0.0d0, 0.0d0], [3, 4]) / sqrt(2.0d0)
            dirs = reshape([ &
                0.0d0, 0.0d0, 0.0d0, &
                1.0d0, 0.0d0, 1.0d0, &
                0.0d0, 1.0d0, 0.0d0], [3, 4])
        case ('111')
            axes = reshape([ &
                1.0d0, 1.0d0, -1.0d0, &
                1.0d0, -1.0d0, 1.0d0, &
                1.0d0, -1.0d0, -1.0d0, &
                1.0d0], [3, 4]) / sqrt(3.0d0)
            dirs = reshape([ &
                1.0d0, 1.0d0, 1.0d0, &
                -1.0d0, 1.0d0, 1.0d0, &
                0.0d0, 0.0d0, 0.0d0], [3, 4]) / sqrt(2.0d0)
        case ('100')
            axes = reshape([ &
                1.0d0, 0.0d0, 0.0d0, &
                0.0d0, 1.0d0, 0.0d0, &
                0.0d0, 0.0d0, 1.0d0], [3, 3])
            dirs = reshape([ &
                0.0d0, 0.0d0, 1.0d0, &
                1.0d0, 0.0d0, 0.0d0, &
                0.0d0, 1.0d0, 0.0d0], [3, 3])
        case default
            print *, "Error: Undefined axis set"
            stop
        end select

        ! Initialize constants
        axes = size(axes, 2)
        period = pi * naxes / 6.0d0

        ! Define the 90 degree rotation matrices
        rotX90 = reshape([ &
            1.0d0, 0.0d0, 0.0d0, &
            0.0d0, 0.0d0, -1.0d0, &
            0.0d0, 1.0d0, 0.0d0], [3, 3])
        rotY90 = reshape([ &
            0.0d0, 0.0d0, 1.0d0, &
            0.0d0, 1.0d0, 0.0d0, &
            -1.0d0, 0.0d0, 0.0d0], [3, 3])
        rotZ90 = reshape([ &
            0.0d0, -1.0d0, 0.0d0, &
            1.0d0, 0.0d0, 0.0d0, &
            0.0d0, 0.0d0, 1.0d0], [3, 3])
        rotZ90m = reshape([ &
            0.0d0, 1.0d0, 0.0d0, &
            -1.0d0, 0.0d0, 0.0d0, &
            0.0d0, 0.0d0, 1.0d0], [3, 3])

        ! Initialize arrays for storing geometry parameters
        distances = 0.0d0
        ksis = 0.0d0
        etas = 0.0d0
        phis = 0.0d0
        thisindex = 0

        ! Generate 24 symmetry-equivalent variants of Q
        V(1) = Q
        V(2) = matmul(V(1), rotX90)
        V(3) = matmul(V(2), rotX90)
        V(4) = matmul(V(3), rotX90)

        do j = 1, 12
            V(j+4) = matmul(V(j), rotY90)
        end do

        do j = 1, 4
            V(j+16) = matmul(V(j), rotZ90)
            V(j+20) = matmul(V(j), rotZ90m)
        end do

        ! Main loop for calculating distances
        do i = 1, naxes
            ax = axes(:, i)  ! High-symmetry axis
            dir = dirs(:, i)
            dir2 = cross(ax, dir)

            do j = 1, 24
                R = transpose(V(j)) * P  ! Misorientation matrix
                call mat2quat(R, axi, psi)  ! Convert matrix to quaternion
                
                dotp = dot_product(axi, ax)
                dis = 2.0d0 * sqrt(abs(1.0d0 - dotp**2)) * sin(psi / 2.0d0)

                if (dis < dismax_local) then
                    thisindex = thisindex + 1
                    theta = 2.0d0 * atan(dotp * tan(psi / 2.0d0))

                    ! Compute the normal of the best-fitting GB in grain 1
                    n1 = P(1, :)
                    n2 = Q(1, :)

                    ! Rotate vector n2 into the grain frame using the idealized rotation RA
                    RA = quat2mat([cos(theta / 2.0d0), sin(theta / 2.0d0) * ax])
                    m1 = n1 + matmul(transpose(RA), n2)

                    ! Avoid numerical issues when m1 is too small
                    if (norm(m1) < 1.0e-6) then
                        print *, "m1 is singular"
                    end if

                    m1 = m1 / norm(m1)
                    m2 = matmul(RA, m1)

                    phi = acos(abs(dot_product(m1, ax)))

                    if (abs(dot_product(ax, m1)) > 0.9999) then
                        theta1 = -theta / 2.0d0
                        theta2 = theta / 2.0d0
                    else
                        theta1 = atan2(dot_product(dir2, m1), dot_product(dir, m1))
                        theta2 = atan2(dot_product(dir2, m2), dot_product(dir, m2))
                    end if

                    ! Ensure theta1 and theta2 are within the range (-period/2, period/2]
                    theta1 = mod(theta1 + period / 2.0d0, period) - period / 2.0d0
                    theta2 = mod(theta2 + period / 2.0d0, period) - period / 2.0d0

                    ksi = abs(theta2 - theta1)
                    eta = abs(theta2 + theta1)

                    distances(thisindex) = dis
                    ksis(thisindex) = ksi
                    etas(thisindex) = eta
                    phis(thisindex) = phi
                end if
            end do
        end do

        ! Trim excess array entries and sort results by distance
        call sort_arrays(distances, ksis, etas, phis, thisindex)

        ! Return the unique geometry results
        geom = transpose(merge([distances(:), ksis(:), etas(:), phis(:)], 0.0d0, thisindex))

    end function distances_to_set

    ! Additional utility functions can be added here like mat2quat, quat2mat, and sort_arrays

end module GB5DOF_utils
function mat2quat(m) result(q)
    ! Converts a 3x3 rotation matrix m to a quaternion q
    implicit none
    real(8), dimension(3, 3), intent(in) :: m
    real(8), dimension(4) :: q
    real(8) :: t, e0, e3
    real(8), dimension(3) :: e

    t = m(1, 1) + m(2, 2) + m(3, 3)

    if (t > -0.999999999) then
        e0 = sqrt(1.0d0 + t) / 2.0d0
        e(1) = (m(2, 3) - m(3, 2)) / (4.0d0 * e0)
        e(2) = (m(3, 1) - m(1, 3)) / (4.0d0 * e0)
        e(3) = (m(1, 2) - m(2, 1)) / (4.0d0 * e0)
    else
        e0 = 0.0d0
        e3 = sqrt(-(m(1, 1) + m(2, 2)) / 2.0d0)
        if (abs(e3) > 2.0d-8) then
            e(1) = m(1, 3) / (2.0d0 * e3)
            e(2) = m(2, 3) / (2.0d0 * e3)
            e(3) = e3
        else
            e0 = sqrt((m(1, 1) + 1.0d0) / 2.0d0)
            if (e0 /= 0.0d0) then
                e(1) = e0
                e(2) = m(2, 1) / (2.0d0 * e0)
                e(3) = 0.0d0
            else
                e(1) = 0.0d0
                e(2) = 1.0d0
                e(3) = 0.0d0
            end if
        end if
    end if

    q(1) = e0
    q(2) = -e(1)
    q(3) = -e(2)
    q(4) = -e(3)

end function mat2quat
function quat2mat(q) result(m)
    ! Converts a quaternion q to a 3x3 rotation matrix m
    implicit none
    real(8), dimension(4), intent(in) :: q
    real(8), dimension(3, 3) :: m
    real(8) :: e0, e1, e2, e3

    e0 = q(1)
    e1 = q(2)
    e2 = q(3)
    e3 = q(4)

    m(1, 1) = e0**2 + e1**2 - e2**2 - e3**2
    m(1, 2) = 2.0d0 * (e1 * e2 - e0 * e3)
    m(1, 3) = 2.0d0 * (e1 * e3 + e0 * e2)

    m(2, 1) = 2.0d0 * (e1 * e2 + e0 * e3)
    m(2, 2) = e0**2 - e1**2 + e2**2 - e3**2
    m(2, 3) = 2.0d0 * (e2 * e3 - e0 * e1)

    m(3, 1) = 2.0d0 * (e1 * e3 - e0 * e2)
    m(3, 2) = 2.0d0 * (e2 * e3 + e0 * e1)
    m(3, 3) = e0**2 - e1**2 - e2**2 + e3**2

    m = m / (e0**2 + e1**2 + e2**2 + e3**2)

end function quat2mat(q)
    end program main
