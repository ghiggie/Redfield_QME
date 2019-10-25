module bath

    use constants
    use variables
    use helper

    implicit none

 contains

    subroutine bc_coeff()

        ! This subroutine will allocate a vector with matsu+1 entries, where
        ! matsu is the number of Matsubara terms. The vector will look like
        !    coeff(0:matsu)
        ! The coefficient for the delta term will be contained in cinf. coeff(0)
        ! will correspond to the coefficient of the exponential gamma term.

        ! We can also edit this subroutine to use the Pade approximation.

        real(kind=DP) :: cotan
        cotan = cos(0.5*gamma/temp)/sin(0.5*gamma/temp)

        coeff(0) = lambda*gamma*CMPLX(cotan,-1)

        do i = 1, matsu
            nu = 2 * PI * temp * i
            coeff(i) = 4*lambda*gamma*temp*nu/(nu**2 - gamma**2)
        end do

        tmp_r = 0
        do i = 1, matsu
            nu = 2 * PI * temp * i
            tmp_r = tmp_r + 1 / (nu**2 - gamma**2)
        end do
        tmp_r = 8*lambda*gamma*temp*tmp_r
        cinf = 4*lambda*temp/gamma - 2*lambda*cotan - tmp_r
    end subroutine bc_coeff

    function lambda_bc(A, B, t)
        real(kind=DP), intent(in) :: t
        real(kind=DP) :: tmp_dt ! I allow for the possibility that this should be shifted
        complex(kind=DP), dimension(:,:), intent(in) :: A, B
        complex(kind=DP), dimension(:,:), allocatable :: lambda_bc
        integer :: M, j
        real(kind=DP) :: t2j, t2j1, t2j2
        complex(kind=DP), dimension(:,:), allocatable :: tmp, tmp1, tmp2
        complex(kind=DP), dimension(:,:), allocatable :: f2j, f2j1, f2j2

        allocate(lambda_bc(ss,ss))

        c3 = 4*lambda*temp/gamma - 2*lambda*cotan - (8*lambda*gamma/temp)/denom

        tmp_dt = dt
        tmp = 0
        M = nint(t / tmp_dt)
        if (mod(M,2) .ne. 0) then
            M = M + 1
            tmp_dt = t / M
        end if
        do j = 0, M/2 - 1
            t2j = 2*j * tmp_dt
            t2j1 = t2j + tmp_dt
            t2j2 = t2j1 + tmp_dt

            tmp1 = (t - t2j) * A
            tmp2 = I2S(tmp1, B)
            f2j = bc(temp, gamma, lambda, t, t2j) * tmp2

            tmp1 = (t - t2j1) * A
            tmp2 = I2S(tmp1, B)
            f2j1 = bc(temp, gamma, lambda, t, t2j1) * tmp2

            tmp1 = (t - t2j2) * A
            tmp2 = I2S(tmp1, B)
            f2j2 = bc(temp, gamma, lambda, t, t2j2) * tmp2

            tmp = tmp + (tmp_dt / 3) * (f2j + 4*f2j1 + f2j2)
        end do
        lambda_bc = tmp + 0.5 * cinf * B
    end function lambda_bc

end module bath
