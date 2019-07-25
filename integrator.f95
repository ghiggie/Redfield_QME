module integrator

    use constants
    use parameters
    use helper

    implicit none

 contains

    function intop(t, A)
        ! Later, rewrite this so that it uses the eigensystem expansion to
        ! compute the interaction picture operators.
        real(kind=DP) :: t
        complex(kind=DP), dimension(:,:), intent(in) :: A
        complex(kind=DP), dimension(:,:) :: intop
        integer :: n

        n = size(A, dim=1)
        allocate(intop(n,n))
        op = BCH(IMAG1*t*HS,A,4)
    end function intop

    function bath_corr(temp, gamma, lambda, t1, t2)
        real(kind=DP), intent(in) :: temp, gamma, lambda, t1, t2
        complex(kind=DP) :: bath_corr
        real(kind=DP) :: cotan, denom, tau
        complex(kind=DP) :: c1, c2, c3, tmp

        tau = t1 - t2 ! Note that t1 > t2
        cotan = cos(0.5*gamma/temp)/sin(0.5*gamma/temp)
        denom = (2*PI)**2 - (gamma/temp)**2

        c1 = 2*lambda*gamma*(REAL1*cotan + (-1)*IMAG1)
        c2 = 8*PI**lambda*gamma/denom
        c3 = 4*lambda*temp/gamma - 2*lambda*cotan - (8*lambda*gamma/temp)/denom

        if (tau .eq. 0) then
            tmp = c1 + c2 + c3
        else
            tmp = c1*exp(-gamma*tau) + c2*exp(-2*PI*temp*tau)
        end if
        bath_corr = tmp
    end function bath_corr

    function rk4_int(t, M, A, temp, lambda, gamma)
        real(kind=DP), intent(in) :: t, temp, lambda, gamma
        integer, intent(in) :: M
        complex(kind=DP), dimension(:,:), intent(in) :: A
        complex(kind=DP), dimension(:,:) :: rk4_int
        complex(kind=DP), dimension(:,:) :: tmp, j1, j2, j3, j4
        integer :: n, j
        real(kind=DP) :: dtj, t_j

        n = size(A, dim=1)
        allocate(rk4_int(n,n))
        allocate(tmp(n,n))

        tmp = 0
        dtj = t / M
        do j = 0, M - 1
            t_j = j * dtj
            j1=(intop(t,A)*dag(intop(t_j,A))*rho0-dag(intop(t_j,A))*rho0*intop(t,A))*bath_corr(temp,gamma,lambda,t,t_j)
            j2=(intop(t,A)*dag(intop(t_j+dtj/2,A))*rho0-dag(intop(t_j+dtj/2,A))*rho0*intop(t,A))*bath_corr(temp,gamma,lambda,t,t_j+dtj/2)
            j3 = j2 ! In this case, j2 and j3 happen to be the same
            j4=(intop(t,A)*dag(intop(t_j+dtj,A))*rho0-dag(intop(t_j+dtj,A))*rho0*intop(t,A))*bath_corr(temp,gamma,lambda,t,t_j+dtj)
            tmp = tmp + (dtj/6)*(j1+2*j2+2*j3+j4)
        end do
        rk4_int = tmp
    end function rk4_int

end module integrator
