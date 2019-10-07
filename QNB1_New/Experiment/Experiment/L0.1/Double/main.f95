program main

    use constants
    use parameters
    use helper
    use bath

    implicit none

    complex(kind=DP), dimension(:,:,:), allocatable :: rho
    complex(kind=DP), dimension(:,:), allocatable :: k1, k2, k3, k4
    complex(kind=DP), dimension(:,:), allocatable :: tmp1, tmp2, tmp3
    character(len=40) :: filename, arg
    integer :: S, N, i, j, k
    real(kind=DP) :: ti, tc
    logical :: tmpl

    namelist/params/dt1,dt2,time_limit,temp,gamma,lambda,rho0,HS,VI

    ! Get the dimension of the system Hilbert space
    call getarg(1,arg)
    read(arg,*) S

    allocate(rho0(S,S))
    allocate(HS(S,S))
    allocate(VI(S,S))
    allocate(k1(S,S))
    allocate(k2(S,S))
    allocate(k3(S,S))
    allocate(k4(S,S))
    allocate(tmp1(S,S))
    allocate(tmp2(S,S))
    allocate(tmp3(S,S))

    ! Get the name of the parameter file
    call getarg(2, arg)
    filename = trim(arg)

    open(10, file=trim(filename))
    read(10,nml=params)
    close(10)

    N = int(time_limit / dt1)
    allocate(rho(0:N,S,S)) !Sets up the storage for the data points
    rho(0,:,:) = rho0 ! Initialize the storage

    open(10,file='rho.dat')
    write(10,'(f10.3,32(e15.6))') 0.0,((rho(0,j,k),j=1,S),k=1,S)

    do i = 0, N - 1
        ti = i * dt1

        tc = ti
        tmp1 = rho(i,:,:)
        tmp2 = lambda_bc(HS, VI, temp, gamma, lambda, tc, dt2)
        tmp3 = matmul(tmp2, tmp1) - matmul(tmp1, dag(tmp2))
        k1 = -IMAG1*comm(HS, tmp1) - comm(VI, tmp3)

        tc = ti + dt1 / 2
        tmp1 = rho(i,:,:) + k1 * dt1 / 2
        tmp2 = lambda_bc(HS, VI, temp, gamma, lambda, tc, dt2)
        tmp3 = matmul(tmp2, tmp1) - matmul(tmp1, dag(tmp2))
        k2 = -IMAG1*comm(HS, tmp1) - comm(VI, tmp3)

        tc = ti + dt1 / 2
        tmp1 = rho(i,:,:) + k2 * dt1 / 2
        tmp2 = lambda_bc(HS, VI, temp, gamma, lambda, tc, dt2)
        tmp3 = matmul(tmp2, tmp1) - matmul(tmp1, dag(tmp2))
        k3 = -IMAG1*comm(HS, tmp1) - comm(VI, tmp3)

        tc = ti + dt1
        tmp1 = rho(i,:,:) + k3 * dt1
        tmp2 = lambda_bc(HS, VI, temp, gamma, lambda, tc, dt2)
        tmp3 = matmul(tmp2, tmp1) - matmul(tmp1, dag(tmp2))
        k4 = -IMAG1*comm(HS, tmp1) - comm(VI, tmp3)

        rho(i+1,:,:) = rho(i,:,:) + (dt1 / 6) * (k1 + 2*k2 + 2*k3 + k4)
        write(10,'(f10.3,32(e15.6))') ti+dt1,((rho(i+1,j,k),j=1,S),k=1,S)
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
    close(10)

    open(10, file='entropy.dat')
    do i = 0, N
        ti = i * dt1
        write(10, '(f10.3,2(e15.6))') ti, Entropy(rho(i,:,:))
    end do
    close(10)

    ! This line isn't perfect.
    open(10, file='systemB.dat')
    do i = 0, N
        ti = i * dt1
        tmp1 = rhoB(rho(i,:,:))
        write(10,'(f10.3,8(e15.6))') ti,((tmp1(j,k),j=1,2),k=1,2)
    end do

end program main
