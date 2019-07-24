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

    namelist/params/N,time_limit,time_interval,temp1,temp2,gamma1,&
        & gamma2,lambda1,lambda2,V1H,V2C,HS,VS,rho0




    dti = time_limit / N
    do i = 0, N - 1
        ! For now, it only calculates for hot bath. Further, it does not include
        ! the bath correlation functions.
        t_i = i * dti
        ! Calculate k1
        k1 = rk4_int(t_i,N,V1H)
        ! Calculate k2
        k2 = rk4_int(t_i+dti/2,N,V1H)
        ! Calculate k3
        k3 = k2
        ! Calculate k4
        k4 = rk4_int(t+dti,N,V1H)
        ! Calculate rho(i+1)
        rho(i+1) = rho(i) + (dti/6)*(k1+2*k2+2*k3+k4) + (dti/6)*dag((k1+2*k2+2*k3+k4),4)
        ! Need to store or print rho(i+1)
        ! Need to be sure to include the result for the cold bath when we compute
        ! this sum
        ! Need to be sure to include the bath correlation functions
    end do




end program main
