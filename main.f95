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
    
    implicit none
    
    real(kind = DP) :: dt, time_limit, time_interval, temp1
    real(kind = DP) :: temp2, gamma1, gamma2, lambda1, lambda2
    
    complex(kind = DP), dimension(4,4) :: rho0, HS, V1H, V2C
    
    
    namelist/params/dt,time_limit,time_interval,temp1,temp2,gamma1,&
        & gamma2,lambda1,lambda2,V1H,V2C,HS,rho0
    
    
    
    
    
    
    
    
    
    
end program main
