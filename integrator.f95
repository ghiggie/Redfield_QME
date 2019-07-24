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

    function rk4_int(t, M, A)
        real(kind=DP), intent(in) :: t
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
            j1=intop(t,A)*dag(intop(t_j,A))*rho0-dag(intop(t_j,A))*rho0*intop(t,A)
            j2=intop(t,A)*dag(intop(t_j+dtj/2,A))*rho0-dag(intop(t_j+dtj/2,A))*rho0*intop(t,A)
            j3 = j2 ! In this case, j2 and j3 happen to be the same
            j4=intop(t,A)*dag(intop(t_j+dtj,A))*rho0-dag(intop(t_j+dtj,A))*rho0*intop(t,A)
            tmp = tmp + (dtj/6)*(j1+2*j2+2*j3+j4)
        end do
        rk4_int = tmp
    end function rk4_int

end module integrator
