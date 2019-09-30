program main

    use constants
    use parameters
    use helper
    use bath

    implicit none

    complex(kind=DP), dimension(:,:,:), allocatable :: rho
    complex(kind=DP), dimension(:,:), allocatable :: k1, k2, k3, k4
    character(len=40) :: filename, arg
    integer :: S, N

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
    write(10,'(f10.3,8(e15.6))') 0.0,((rho(0,j,k),j=1,S),k=1,S)

    
