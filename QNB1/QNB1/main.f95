program main

!---------------------------------------------------------------------!
!   This program will plot the Born master equation (without the      !
!   Markov or secular approximations) in the weak coupling limit. It  !
!   is expected that this program does not preserve positivity of the !
!   density matrix.                                                   !
!                                                                     !
!   Version 0.1 (July 5, 2019)                                        !
!   Coded by Garrett Higginbotham at the University of Alabama at     !
!   Birmingham                                                        !
!   Copyright 2019 by Garrett Higginbotham                            !
!---------------------------------------------------------------------!

    use constants
    use parameters
    use integrator
    use helper
    implicit none

    complex(kind=DP), dimension(:,:,:), allocatable :: rho
    complex(kind=DP), dimension(:,:), allocatable :: k1, k2, k3, k4
    integer :: i, j, k, N
    real(kind=DP) :: dti, t_i
    character(len=40) :: filename, arg
    logical :: tmpl

    namelist/params/dt1,dt2,time_limit,temp,gamma,lambda,rho0,HS,VI

    ! Get the dimension of the system Hilbert space
    call getarg(1,arg)
    read(arg,*) k

    allocate(rho0(k,k))
    allocate(HS(k,k))
    allocate(VI(k,k))
    allocate(k1(k,k))
    allocate(k2(k,k))
    allocate(k3(k,k))
    allocate(k4(k,k))

    ! Get the name of the parameter file
    call getarg(2, arg)
    filename = trim(arg)

    open(10, file=trim(filename))
    read(10,nml=params)
    close(10)

    N = int(time_limit / dt1)
    allocate(rho(0:N,k,k)) !Sets up the storage for the data points
    rho(0,:,:) = rho0 ! Initialize the storage

    open(10,file='rho.dat')
    write(10,'(f10.3,8(e15.6))') 0.0,((rho(0,j,k),j=1,2),k=1,2)

    do i = 0, N - 1
        t_i = i * dt1
        ! Calculate k1
        k1 = rk4_int(t_i,t_i,dt2,VI,rho(i,:,:),temp,lambda,gamma)
        ! Calculate k2
        k2 = rk4_int(t_i,t_i+dt1/2,dt2,VI,rho(i,:,:),temp,lambda,gamma)
        ! Calculate k3
        k3 = k2
        ! Calculate k4
        k4 = rk4_int(t_i,t_i+dt1,dt2,VI,rho(i,:,:),temp,lambda,gamma)
        ! Calculate rho(i+1)
        rho(i+1,:,:) = rho(i,:,:) + (dt1/6)*(k1+2*k2+2*k3+k4) + (dt1/6)*dag((k1+2*k2+2*k3+k4))
        write(10,'(f10.3,8(e15.6))') t_i+dt1,((rho(i+1,j,k),j=1,2),k=1,2)
    end do

    close(10)

    open(10, file='hermitian.dat')
    do i = 0, N
        t_i = i * dt1
        tmpl = test_hermitian(rho(i,:,:))
        write(10, '(f10.3,L2)') t_i, tmpl
    end do
    close(10)

    open(10, file='trace.dat')
    do i = 0, N
        t_i = i * dt1
        tmpl = test_trace(rho(i,:,:))
        write(10, '(f10.3,L2,2e15.6)') t_i, tmpl, trace(rho(i,:,:))
    end do
    close(10)

    open(10, file='positivity.dat')
    do i = 0, N
        t_i = i * dt1
        tmpl = test_positivity(rho(i,:,:))
        write(10, '(f10.3,L2)') t_i, tmpl
    end do

end program main
