module const
    implicit none
    integer, parameter :: sp = selected_int_kind(9)
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: one = 1.0_dp, two = 2.0_dp, four = 4.0_dp
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: half = 0.5_dp, tby4 = 0.25_dp
end module const

module AngularMomentsModule
    use const
    implicit none


end module AngularMomentsModule


module MathFunction
    use const
    implicit none
    private
    type :: AngularMoments
        integer(sp) :: lx, ly, lz
    end type AngularMoments

    integer(sp), parameter :: max_shell_type = 15, max_shell_index = 10
    type(AngularMoments), dimension(max_shell_type, max_shell_index) :: moments

    contains
    subroutine initializeAngularMoments()
        integer(sp) :: shell_type, shell_index

        ! Define the Cartesian angular moments for each shell type and index
        moments = reshape([ &
        AngularMoments(0, 0, 0), &
        AngularMoments(1, 0, 0), &
        AngularMoments(0, 1, 0), &
        AngularMoments(0, 0, 1), &
        AngularMoments(2, 0, 0), &
        AngularMoments(0, 2, 0), &
        AngularMoments(0, 0, 2), &
        AngularMoments(1, 1, 0), &
        AngularMoments(1, 0, 1), &
        AngularMoments(0, 1, 1), &
        AngularMoments(3, 0, 0), &
        AngularMoments(0, 3, 0), &
        AngularMoments(0, 0, 3), &
        AngularMoments(2, 1, 0), &
        AngularMoments(2, 0, 1), &
        AngularMoments(1, 2, 0), &
        AngularMoments(0, 2, 1), &
        AngularMoments(1, 0, 2), &
        AngularMoments(0, 1, 2), &
        AngularMoments(1, 1, 1) &
        ], [max_shell_type, max_shell_index])

    end subroutine initializeAngularMoments

    subroutine get_cart_anglr_moment(shell_type, shell_index, lx, ly, lz)
        integer(sp), intent(in)  :: shell_type, shell_index
        integer(sp), intent(out) :: lx, ly, lz

        ! Assign the Cartesian angular moments based on shell type and index
        lx = moments(shell_type, shell_index)%lx
        ly = moments(shell_type, shell_index)%ly
        lz = moments(shell_type, shell_index)%lz

    end subroutine get_cart_anglr_moment

    function factii(num) result(res)
        integer(sp), intent(in) :: num
        real(dp) :: res
        select case (num)
        case (1); res = 1.0_dp
        case (2); res = 2.0_dp
        case (3); res = 3.0_dp
        case (4); res = 8.0_dp
        case (5); res = 15.0_dp
        case (6); res = 48.0_dp
        case (7); res = 105.0_dp
        case (8); res = 384.0_dp
        case (9); res = 945.0_dp
        case (10); res = 3840.0_dp
        case (11); res = 10395.0_dp
        case (12); res = 46080.0_dp
        case (13); res = 135135.0_dp
        case (14); res = 645120.0_dp
        case (15); res = 2027025.0_dp
        case default; res = 1.0_dp
        end select
    end function factii

    function boysn(X, N) result(res)
        real(dp), intent(in) :: X, N
        real(dp) :: res, val2, twox
        integer(sp) :: IX
        twox = two * X
        val2 = dexp(-X) / twox
        if (int(N) <= 0 .and. dabs(X - zero) < 1.0d-16) then
            res = one
        else
            res = sqrt(pi / (four * X)) * derf(sqrt(X))
        end if
        do IX = 1, floor(N)
            if (X == zero) then
                res = one / (two * IX + one)
            else
                res = ((two * (IX - one) + one) * res) / twox - val2
            end if
        end do
    end function boysn

end module MathFunction

module Gaussian_Integrals
    use const
    use utils
    implicit none
contains



    function ovr_nrm(l, m, n, alp) result(res)
        integer(sp), intent(in) :: l, m, n
        real(dp), intent(in) :: alp
        real(dp) :: res, val, fctlmn
        val = (four * alp) ** ((l + m + n) * one / two)
        fctlmn = sqrt(factii(2 * l - 1) * factii(2 * m - 1) * factii(2 * n - 1))
        res = (((two * alp) / pi) ** tby4) * (val / fctlmn)
    end function ovr_nrm

    function quad_1dir(la, lb, alp, bta, ab, a, b, pa, d1dir) result(res)
        integer(sp), intent(in) :: la, lb
        real(dp), intent(in) :: alp, bta, ab, a, b, pa, d1dir
        real(dp) :: res
        real(dp) :: sgma_1dir

        if (la <= 0 .and. lb <= 0) then
            res = sgma_1dir(1, 1, alp, bta, ab, pa) + b * pa + a * d1dir
        else
            res = sgma_1dir(la + 1, lb + 1, alp, bta, ab, pa) + &
                    sgma_1dir(la + 1, lb, alp, bta, ab, pa) + a * d1dir
        end if
    end function quad_1dir

    function dpl_1dir(la, lb, alp, bta, ab, b, pa, s1dir) result(res)
        integer(sp), intent(in) :: la, lb
        real(dp), intent(in) :: alp, bta, ab, b, pa, s1dir
        real(dp) :: res
        real(dp) :: sgma_1dir

        if (la <= 0 .and. lb <= 0) then
            res = sgma_1dir(0, 1, alp, bta, ab, pa) + b
        else
            res = sgma_1dir(la, lb + 1, alp, bta, ab, pa) + b * s1dir
        end if
    end function dpl_1dir

    function vel_1dir(la, lb, expnt, alp, bta, ab, pa) result(res)
        integer(sp), intent(in) :: la, lb
        real(dp), intent(in) :: expnt, alp, bta, ab, pa
        real(dp) :: res
        real(dp) :: sgma_1dir, twexp

        twexp = two * expnt
        if (la <= 0 .and. lb <= 0) then
            res = -twexp * sgma_1dir(0, 1, alp, bta, ab, pa)
        elseif (la > 0 .and. lb <= 0) then
            res = -twexp * sgma_1dir(la, 1, alp, bta, ab, pa)
        else
            res = lb * sgma_1dir(la, lb - 1, alp, bta, ab, pa) - &
                    twexp * sgma_1dir(la, lb + 1, alp, bta, ab, pa)
        end if
    end function vel_1dir

    function kin_1dir(la, lb, alp, bta, ab, pa) result(res)
        integer(sp), intent(in) :: la, lb
        real(dp), intent(in) :: alp, bta, ab, pa
        real(dp) :: res, twalbt
        real(dp) :: sgma_1dir

        twalbt = two * alp * bta
        if (la == 0 .and. lb == 0) then
            res = twalbt * sgma_1dir(1, 1, alp, bta, ab, pa)
        elseif (la > 0 .and. lb == 0) then
            res = -la * bta * sgma_1dir(la - 1, 1, alp, bta, ab, pa) + &
                    twalbt * sgma_1dir(la + 1, 1, alp, bta, ab, pa)
        elseif (la == 0 .and. lb > 0) then
            res = -lb * alp * sgma_1dir(1, lb - 1, alp, bta, ab, pa) + &
                    twalbt * sgma_1dir(1, lb + 1, alp, bta, ab, pa)
        else
            res = half * (la * lb * sgma_1dir(la - 1, lb - 1, alp, bta, ab, pa) - &
                    two * alp * lb * sgma_1dir(la + 1, lb - 1, alp, bta, ab, pa) - &
                    two * bta * la * sgma_1dir(la - 1, lb + 1, alp, bta, ab, pa) + &
                    two * twalbt * sgma_1dir(la + 1, lb + 1, alp, bta, ab, pa))
        end if
    end function kin_1dir

    function gauss_sgma_1dir(la, lb, alp, bta, gama, ab, qa) result(res)
        integer(sp), intent(in) :: la, lb
        real(dp), intent(in) :: alp, bta, gama, ab, qa
        real(dp) :: res
        res = sgma_1dir(la, lb, alp + bta, gama, ab, qa)
    end function gauss_sgma_1dir

    function sgma_1dir(la, lb, alp, bta, ab, pa) result(res)
        use const
        implicit none
        integer(sp), intent(in) :: la, lb
        real(dp), intent(in) :: alp, bta, ab, pa
        real(dp) :: res, term1, term2, texp, exp_factor

        texp = dabs(pa - ab)
        exp_factor = exp(-alp * bta * (texp ** 2) / (alp + bta))

        if (la == 0 .and. lb == 0) then
            res = exp_factor
        else
            term1 = (la * one) / (two * alp)
            term2 = (lb * one) / (two * bta)
            if (la > 0 .and. lb == 0) then
                res = term1 * sgma_1dir(la - 1, lb, alp, bta, ab, pa)
            elseif (la == 0 .and. lb > 0) then
                res = term2 * sgma_1dir(la, lb - 1, alp, bta, ab, pa)
            else
                res = term1 * sgma_1dir(la - 1, lb, alp, bta, ab, pa) + &
                        term2 * sgma_1dir(la, lb - 1, alp, bta, ab, pa)
            end if
        end if
    end function sgma_1dir

    function sgma_1dir2(la, lb, alp, bta, ab, pa) result(res)
        use const
        implicit none
        integer(sp), intent(in) :: la, lb
        real(dp), intent(in) :: alp, bta, ab, pa
        real(dp) :: res
        real(dp) :: texp, exp_factor
        real(dp), allocatable :: dp_table(:,:)
        integer :: i, j

        ! Calculate the exponential factor
        texp = dabs(pa - ab)
        exp_factor = exp(-alp * bta * (texp ** 2) / (alp + bta))

        ! Allocate memory for the table to store intermediate results
        allocate(dp_table(0:la, 0:lb))

        ! Initialize the base case
        dp_table(0, 0) = exp_factor

        ! Fill the table in a bottom-up manner
        do i = 0, la
            do j = 0, lb
                if (i == 0 .and. j == 0) cycle
                dp_table(i, j) = 0.0_dp
                if (i > 0) dp_table(i, j) = dp_table(i, j) + (i * one) / (two * alp) * dp_table(i - 1, j)
                if (j > 0) dp_table(i, j) = dp_table(i, j) + (j * one) / (two * bta) * dp_table(i, j - 1)
            end do
        end do

        ! The result is stored in dp_table(la, lb)
        res = dp_table(la, lb)

        ! Deallocate the table
        deallocate(dp_table)
    end function sgma_1dir2

end module Gaussian_Integrals
