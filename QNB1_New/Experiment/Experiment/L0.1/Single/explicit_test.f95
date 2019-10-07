program explicit_test

    use constants
    use parameters
    use helper
    use bath

    implicit none

    complex(kind=DP), dimension(:,:,:), allocatable :: rho
    integer :: i, j, k, N, M
    complex(kind=DP) :: j1, j2, j3, j4
    complex(kind=DP), dimension(2,2) :: gamma_m, k1, k2, k3, k4, tmp1, tmp2
    logical :: tmpl
    real(kind=DP) :: ti, tj, tin, tc

    allocate(rho0(2,2))
    allocate(HS(2,2))
    allocate(VI(2,2))

    rho0(1, :) = (/cmplx(1., 0), cmplx(0, 0)/)
    rho0(2, :) = (/cmplx(0, 0), cmplx(0, 0)/)
    HS(1, :) = (/cmplx(0.5, 0), cmplx(0, 0)/)
    HS(2, :) = (/cmplx(0, 0), cmplx(-0.5, 0)/)
    VI(1, :) = (/cmplx(0, 0), cmplx(1., 0)/)
    VI(2, :) = (/cmplx(1., 0), cmplx(0, 0)/)

    time_limit = 100
    dt1 = 0.1
    dt2 = 0.1
    temp = 1
    gamma = 0.1
    lambda = 0.1

    N = int(time_limit / dt1)
    allocate(rho(0 : N, 2, 2))
    rho(0, :, :) = rho0

    open(10,file='rho.dat')
    write(10,'(f10.3,8(e15.6))') 0.0,((rho(0,j,k),j=1,2),k=1,2)

    do i = 0, N - 1
        ti = i * dt1

        ! Calculate the gamma matrix for ti
        tin = ti
        gamma_m(1, 1) = 0
        M = int(tin / dt2)
        do j = 0, M - 1
            tj = j * dt2

            j1 = bc(temp, gamma, lambda, tin, tj)
            j2 = bc(temp, gamma, lambda, tin, tj + dt2 / 2)
            j3 = j2
            j4 = bc(temp, gamma, lambda, tin, tj + dt2)
            gamma_m(1,1) = gamma_m(1,1) + (dt2/6)*(j1 + 2*j2 + 2*j3 + j4)
        end do
        gamma_m(1, 2) = 0
        M = int(ti / dt2)
        do j = 0, M - 1
            tj = j * dt2

            tc = tj
            j1 = bc(temp, gamma, lambda, tin, tc)*cmplx(cos(tin-tc),sin(tin-tc))
            tc = tj + dt2 / 2
            j2 = bc(temp, gamma, lambda, tin, tc)*cmplx(cos(tin-tc),sin(tin-tc))
            j3 = j2
            tc = tj + dt2
            j4 = bc(temp, gamma, lambda, tin, tc)*cmplx(cos(tin-tc),sin(tin-tc))
            gamma_m(1,2) = gamma_m(1,2) + (dt2/6)*(j1 + 2*j2 + 2*j3 + j4)
        end do
        gamma_m(2,1) = conjg(gamma_m(1,2))
        gamma_m(2, 2) = gamma_m(1, 1)

        ! Now calulate the first step of RK4 for rho
        tmp1 = rho(i,:,:)
        tmp2 = matmul(VI,matmul(gamma_m,tmp1)) - matmul(tmp1,matmul(dag(gamma_m),VI))
        k1 = -IMAG1*comm(HS, tmp1) - comm(VI, tmp2)

        ! Calculate the gamma matrix for ti + dt1 / 2
        tin = ti + dt1 / 2
        gamma_m(1, 1) = 0
        M = int(tin / dt2)
        do j = 0, M - 1
            tj = j * dt2

            j1 = bc(temp, gamma, lambda, tin, tj)
            j2 = bc(temp, gamma, lambda, tin, tj + dt2 / 2)
            j3 = j2
            j4 = bc(temp, gamma, lambda, tin, tj + dt2)
            gamma_m(1,1) = gamma_m(1,1) + (dt2/6)*(j1 + 2*j2 + 2*j3 + j4)
        end do
        gamma_m(1, 2) = 0
        M = int(ti / dt2)
        do j = 0, M - 1
            tj = j * dt2

            tc = tj
            j1 = bc(temp, gamma, lambda, tin, tc)*cmplx(cos(tin-tc),sin(tin-tc))
            tc = tj + dt2 / 2
            j2 = bc(temp, gamma, lambda, tin, tc)*cmplx(cos(tin-tc),sin(tin-tc))
            j3 = j2
            tc = tj + dt2
            j4 = bc(temp, gamma, lambda, tin, tc)*cmplx(cos(tin-tc),sin(tin-tc))
            gamma_m(1,2) = gamma_m(1,2) + (dt2/6)*(j1 + 2*j2 + 2*j3 + j4)
        end do
        gamma_m(2,1) = conjg(gamma_m(1,2))
        gamma_m(2, 2) = gamma_m(1, 1)

        ! Now calulate the second step of RK4 for rho
        tmp1 = rho(i,:,:) + k1 * dt1 / 2
        tmp2 = matmul(VI,matmul(gamma_m,tmp1)) - matmul(tmp1,matmul(dag(gamma_m),VI))
        k2 = -IMAG1*comm(HS, tmp1) - comm(VI, tmp2)

        ! Now calulate the third step of RK4 for rho
        tmp1 = rho(i,:,:) + k2 * dt1 / 2
        tmp2 = matmul(VI,matmul(gamma_m,tmp1)) - matmul(tmp1,matmul(dag(gamma_m),VI))
        k3 = -IMAG1*comm(HS, tmp1) - comm(VI, tmp2)

        ! Calculate the gamma matrix for ti + dt1
        tin = ti + dt1
        gamma_m(1, 1) = 0
        M = int(tin / dt2)
        do j = 0, M - 1
            tj = j * dt2

            j1 = bc(temp, gamma, lambda, tin, tj)
            j2 = bc(temp, gamma, lambda, tin, tj + dt2 / 2)
            j3 = j2
            j4 = bc(temp, gamma, lambda, tin, tj + dt2)
            gamma_m(1,1) = gamma_m(1,1) + (dt2/6)*(j1 + 2*j2 + 2*j3 + j4)
        end do
        gamma_m(1, 2) = 0
        M = int(ti / dt2)
        do j = 0, M - 1
            tj = j * dt2

            tc = tj
            j1 = bc(temp, gamma, lambda, tin, tc)*cmplx(cos(tin-tc),sin(tin-tc))
            tc = tj + dt2 / 2
            j2 = bc(temp, gamma, lambda, tin, tc)*cmplx(cos(tin-tc),sin(tin-tc))
            j3 = j2
            tc = tj + dt2
            j4 = bc(temp, gamma, lambda, tin, tc)*cmplx(cos(tin-tc),sin(tin-tc))
            gamma_m(1,2) = gamma_m(1,2) + (dt2/6)*(j1 + 2*j2 + 2*j3 + j4)
        end do
        gamma_m(2,1) = conjg(gamma_m(1,2))
        gamma_m(2, 2) = gamma_m(1, 1)

        ! Now calulate the fourth step of RK4 for rho
        tmp1 = rho(i,:,:) + k3 * dt1
        tmp2 = matmul(VI,matmul(gamma_m,tmp1)) - matmul(tmp1,matmul(dag(gamma_m),VI))
        k4 = -IMAG1*comm(HS, tmp1) - comm(VI, tmp2)

        rho(i+1,:,:) = rho(i,:,:) + (dt1/6) * (k1+2*k2+2*k3+k4)
        write(10,'(f10.3,8(e15.6))') ti+dt1,((rho(i+1,j,k),j=1,2),k=1,2)
    end do

    close(10)

    open(10, file='hermitian.dat')
    do i = 0, N
        ti = i * dt1
        tmpl = test_hermitian(rho(i,:,:))
        write(10, '(f10.3,L2)') ti, tmpl
    end do
    close(10)

    open(10, file='trace.dat')
    do i = 0, N
        ti = i * dt1
        tmpl = test_trace(rho(i,:,:))
        write(10, '(f10.3,L2,2e15.6)') ti, tmpl, trace(rho(i,:,:))
    end do
    close(10)

    open(10, file='positivity.dat')
    do i = 0, N
        ti = i * dt1
        tmpl = test_positivity(rho(i,:,:))
        write(10, '(f10.3,L2)') ti, tmpl
    end do

end program explicit_test
