module bath

    use constants
    use variables
    use helper

    implicit none

    real(kind=DP) :: denom
    complex(kind=DP) :: c1, c2, c3

    denom = (2*PI)**2 - (gamma/temp)**2

    c2 = 8*PI*lambda*gamma/denom

 contains

    subroutine bc_coeff()

        ! This subroutine will operate on a vector with K+2 entries, where
        ! K is the number of Matsubara terms. The vector will look like
        !    coeff(0:K+1)

        real(kind=DP) :: cotan
        cotan = cos(0.5*gamma/temp)/sin(0.5*gamma/temp)

        coeff(0) = lambda*gamma*CMPLX(cotan,-1)

        do i = 1, K
            nu = 2 * PI * temp * i
            coeff(i) = 4*lambda*gamma*temp*nu/(nu**2 - gamma**2)
        end do

        tmp_r = 0
        do i = 1, K
            nu = 2 * PI * temp * i
            tmp_r = tmp_r + 1 / (nu**2 - gamma**2)
        end do
        tmp_r = 8*lambda*gamma*temp*tmp_r
        coeff(K+1) = 4*lambda*temp/gamma - 2*lambda*cotan - tmp_r
    end subroutine bc_coeff

    function bc(tau)
        real(kind=DP), intent(in) ::  tau
        complex(kind=DP) :: bc
        complex(kind=DP) :: tmp

        tmp = c1*exp(-gamma*tau) + c2*exp(-2*PI*temp*tau)

        bc = tmp
    end function bc

    function lambda_bc(A, B, t, dt)
        real(kind=DP), intent(in) :: t
        real(kind=DP) :: dt ! I allow for the possibility that this should be shifted
        complex(kind=DP), dimension(:,:), intent(in) :: A, B
        complex(kind=DP), dimension(:,:), allocatable :: lambda_bc
        integer :: M, j
        real(kind=DP) :: t2j, t2j1, t2j2
        complex(kind=DP), dimension(:,:), allocatable :: tmp, tmp1, tmp2
        complex(kind=DP), dimension(:,:), allocatable :: f2j, f2j1, f2j2

        allocate(lambda_bc(S,S))

        c3 = 4*lambda*temp/gamma - 2*lambda*cotan - (8*lambda*gamma/temp)/denom

        tmp = 0
        M = nint(t / dt)
        if (mod(M,2) .ne. 0) then
            M = M + 1
            dt = t / M
        end if
        if (mod(m,2) .ne. 0) STOP "Value of M is still not even; bath.f95, line 57"
        do j = 0, M/2 - 1
            t2j = 2*j * dt
            t2j1 = t2j + dt
            t2j2 = t2j1 + dt

            tmp1 = (t - t2j) * A
            tmp2 = I2S(tmp1, B)
            f2j = bc(temp, gamma, lambda, t, t2j) * tmp2

            tmp1 = (t - t2j1) * A
            tmp2 = I2S(tmp1, B)
            f2j1 = bc(temp, gamma, lambda, t, t2j1) * tmp2

            tmp1 = (t - t2j2) * A
            tmp2 = I2S(tmp1, B)
            f2j2 = bc(temp, gamma, lambda, t, t2j2) * tmp2

            tmp = tmp + (dt / 3) * (f2j + 4*f2j1 + f2j2)
        end do
        lambda_bc = tmp + 0.5 * c3 * B
    end function lambda_bc

end module bath
