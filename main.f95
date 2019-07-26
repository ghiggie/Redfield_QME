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

    complex(kind=DP), dimension(:,:) :: rho, k1, k2, k3, k4
    integer :: i, j, k
    real(kind=DP) :: dti, t_i
    character(len=40) :: filename, arg

    namelist/params/N,time_limit,time_interval,temp1,temp2,gamma1,gamma2,lambda1,lambda2,V1H,V2C,HS,VS,rho0

    ! Get the name of the parameter file
    call getarg(1, arg)
    filename = trim(arg)

    open(10, file=trim(filename))
    read(10,nml=params)
    close(10)

    k = size(HS, dim=1)
    allocate(rho(0:N,k,k)) !Sets up the storage for the data points
    rho(0,:,:) = rho0 ! Initialize the storage

    dti = time_limit / N
    do i = 0, N - 1
        t_i = i * dti
        ! Calculate k1
        k1 = rk4_int(t_i,N,V1H,temp1,lambda1,gamma1)+rk4_int(t_i,N,V2C,temp2,lambda2,gamma2)
        ! Calculate k2
        k2 = rk4_int(t_i+dti/2,N,V1H,temp1,lambda1,gamma1)+rk4_int(t_i+dti/2,N,V2C,temp2,lambda2,gamma2)
        ! Calculate k3
        k3 = k2
        ! Calculate k4
        k4 = rk4_int(t+dti,N,V1H,temp1,lambda1,gamma1)+rk4_int(t_i+dti,N,V2C,temp2,lambda2,gamma2)
        ! Calculate rho(i+1)
        rho(i+1,:,:) = rho(i,:,:) + (dti/6)*(k1+2*k2+2*k3+k4) + (dti/6)*dag((k1+2*k2+2*k3+k4),4)
    end do

    ! Write the density matrices to a file
    open(10,file='rho.dat')
    do i = 0, N
        t_i = i * dti
        write(10,'(f10.3,32(e15.6e3))') time,((rho(i,j,k),j=1,4),k=1,4)
    end do
    close(10)

end program main
