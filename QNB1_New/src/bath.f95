module bath

    use constants
    use parameters
    use helper

    implicit none

 contains

    function bc(temp, gamma, lambda, t1, t2)
        real(kind=DP), intent(in) :: temp, gamma, lambda, t1, t2
        complex(kind=DP) :: bc
        real(kind=DP) :: cotan, denom, tau
        complex(kind=DP) :: c1, c2, c3, tmp

        tau = t1 - t2 ! Note that t1 > t2
        cotan = cos(0.5*gamma/temp)/sin(0.5*gamma/temp)
        denom = (2*PI)**2 - (gamma/temp)**2

        c1 = 2*lambda*gamma*(REAL1*cotan + (-1)*IMAG1)
        c2 = 8*PI*lambda*gamma/denom
        c3 = 4*lambda*temp/gamma - 2*lambda*cotan - (8*lambda*gamma/temp)/denom

        if (tau .eq. 0) then
            tmp = c1 + c2 + c3
        else
            tmp = c1*exp(-gamma*tau) + c2*exp(-2*PI*temp*tau)
        end if
        bc = tmp
    end function bc

    function lambda_bc(A, B, temp, gamma, lambda, t, dt)
        real(kind=DP), intent(in) :: temp, gamma, lambda, t, dt
        complex(kind=DP), dimension(:,:), intent(in) :: A, B
        complex(kind=DP), dimension(:,:), allocatable :: lambda_bc
        integer :: n, M, j
        real(kind=DP) :: t2j, t2j1, t2j2
        complex(kind=DP), dimension(:,:), allocatable :: tmp, tmp1, tmp2
        complex(kind=DP), dimension(:,:), allocatable :: f2j, f2j1, f2j2

        n = size(A, dim=1)
        allocate(lambda_bc(n,n))
        allocate(tmp(n,n))
        allocate(tmp1(n,n))
        allocate(tmp2(n,n))
        allocate(f2j(n,n))
        allocate(f2j1(n,n))
        allocate(f2j2(n,n))

        tmp = 0
        M = int(t / dt)
        do j = 0, M/2 - 1
            t2j = 2*j * dt
            t2j1 = t2j + dt
            t2j2 = t2j1 + dt

            tmp1 = (t - t2j) * A
            tmp2 = ExpOp(tmp1, B)
            f2j = bc(temp, gamma, lambda, t, t2j) * tmp2

            tmp1 = (t - t2j1) * A
            tmp2 = ExpOp(tmp1, B)
            f2j1 = bc(temp, gamma, lambda, t, t2j1) * tmp2

            tmp1 = (t - t2j2) * A
            tmp2 = ExpOp(tmp1, B)
            f2j2 = bc(temp, gamma, lambda, t, t2j2) * tmp2

            tmp = tmp + (dt / 3) * (f2j + 4*f2j1 + f2j2)
        end do
        lambda_bc = tmp
    end function lambda_bc

end module bath
