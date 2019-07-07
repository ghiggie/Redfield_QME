module parameters

!---------------------------------------------------------------------!
!    This file will declare the variables that the user must set.     !
!    This will also help in unifying the code.                        !
!---------------------------------------------------------------------!

    use constants

    implicit none
    
    real(kind = DP) :: dt, time_limit, time_interval, temp1
    real(kind = DP) :: temp2, gamma1, gamma2, lambda1, lambda2
    
    complex(kind = DP), dimension(4,4) :: rho0, HS, V1H, V2C

end module parameters
