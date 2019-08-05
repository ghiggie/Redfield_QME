module integrator

    use constants
    use parameters
    use helper

    implicit none

 contains

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

    function rk4_int(t, dt, A, temp, lambda, gamma)
        real(kind=DP), intent(in) :: t, dt, temp, lambda, gamma
        complex(kind=DP), dimension(:,:), allocatable, intent(in) :: A
        complex(kind=DP), dimension(:,:), allocatable :: rk4_int
        complex(kind=DP), dimension(:,:), allocatable :: tmp, j1, j2, j3, j4, op1, op2, op3, op4
        integer :: n, j, M
        real(kind=DP) :: t_j
        complex(kind=DP) :: bc1, bc2, bc3

        n = size(A, dim=1)
        allocate(rk4_int(n,n))
        allocate(tmp(n,n))
        allocate(j1(n,n))
        allocate(j2(n,n))
        allocate(j3(n,n))
        allocate(j4(n,n))
        allocate(op1(n,n))
        allocate(op2(n,n))
        allocate(op3(n,n))
        allocate(op4(n,n))

        op1 = intop(t*HS,A)

        tmp = 0
        M = int(t / dt)
        do j = 0, M - 1
            t_j = j * dt
            op2 = intop(t_j*HS,A)
            op3 = intop((t_j + dt/2)*HS,A)
            op4 = intop((t_j+dt)*HS,A)
            bc1 = bath_corr(temp,gamma,lambda,t,t_j)
            bc2 = bath_corr(temp,gamma,lambda,t,t_j+dt/2)
            bc3 = bath_corr(temp,gamma,lambda,t,t_j+dt)
            j1=(-matmul(matmul(op1,dag(op2)),rho0)+matmul(matmul(dag(op2),rho0),op1))*bc1
            j2=(-matmul(matmul(op1,dag(op3)),rho0)+matmul(matmul(dag(op3),rho0),op1))*bc2
            j3 = j2 ! In this case, j2 and j3 happen to be the same
            j4=(-matmul(matmul(op1,dag(op4)),rho0)+matmul(matmul(dag(op4),rho0),op1))*bc3
            tmp = tmp + (dt/6)*(j1+2*j2+2*j3+j4)
        end do
        rk4_int = tmp
    end function rk4_int

end module integrator
