module bath

    use constants
    use parameters
    use helper

    implicit none

 contains

    ! This function computes the bath correlation function. Right now, it only
    ! uses a single Matsubara exponent.
    function bc(temp, gamma, lambda, t1, t2)
        real(kind=DP), intent(in) :: temp, gamma, lambda, t1, t2
        complex(kind=DP) :: bc
        real(kind=DP) :: cotan, denom, tau
        complex(kind=DP) :: c1, c2, c3

        tau = t1 - t2 ! Note that t1 > t2
        cotan = cos(0.5*gamma/temp)/sin(0.5*gamma/temp)
        denom = (2*PI)**2 - (gamma/temp)**2

        c1 = 2*lambda*gamma*CMPLX(cotan, -1)
        c2 = 8*PI*lambda*gamma/denom
        c3 = 4*lambda*temp/gamma - 2*lambda*cotan - (8*lambda*gamma/temp)/denom

        if (tau .eq. 0) then
            bc = c1 + c2 + c3
        else
            bc = c1*exp(-gamma*tau) + c2*exp(-2*PI*temp*tau)
        end if
    end function bc

    function lambda_bc(A, B, temp, gamma, lambda, t, dt)
        real(kind=DP), intent(in) :: temp, gamma, lambda, t, dt
        complex(kind=DP), dimension(:,:), intent(in) :: A, B
        complex(kind=DP), dimension(:,:), allocatable :: lambda_bc
        integer :: n, M, j
        real(kind=DP) :: tj, tc
        complex(kind=DP), dimension(:,:), allocatable :: j1, j2, j3, j4

        n = size(A, dim=1)
        allocate(lambda_bc(n,n))
        allocate(j1(n,n))
        allocate(j2(n,n))
        allocate(j3(n,n))
        allocate(j4(n,n))

        lambda_bc = 0
        M = int(t / dt)
        do j = 0, M - 1
            tj = j * dt

            ! First RK4 step
            tc = tj
            j1 = bc(temp, gamma, lambda, t, tc) * ExpOp((t - tc) * A, B)

            ! Second RK4 Step
            tc = tj + dt / 2
            j2 = bc(temp, gamma, lambda, t, tc) * ExpOp((t - tc) * A, B)

            ! Because the integral doesn't itself depend on the output function,
            ! the third RK4 step is the same as the second RK4 step.
            j3 = j2

            ! Fourth RK4 step
            tc = tj + dt
            j4 = bc(temp, gamma, lamba, t, tc) * ExpOp((t - tc) * A, B)

            ! Finish the RK4 iteration
            lambda_bc = lambda_bc + (dt / 6) * (j1 + 2*j2 + 2*j3 + j4)
        end do
    end function lambda_bc

end module bath
