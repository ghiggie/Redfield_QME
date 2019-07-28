module parameters

!---------------------------------------------------------------------!
!    This file will declare the variables that the user must set.     !
!    This will also help in unifying the code.                        !
!---------------------------------------------------------------------!

    use constants
    implicit none

    real(kind = DP) :: N, time_limit, time_interval, trunc, temp1
    real(kind = DP) :: temp2, gamma1, gamma2, lambda1, lambda2

    complex(kind = DP), dimension(:,:), allocatable :: rho0, HS, VS, V1H, V2C

end module parameters
