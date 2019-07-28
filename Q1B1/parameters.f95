module parameters

!---------------------------------------------------------------------!
!    This file will declare the variables that the user must set.     !
!    This will also help in unifying the code.                        !
!---------------------------------------------------------------------!

    use constants
    implicit none

    real(kind = DP) :: N, time_limit, time_interval, trunc, temp
    real(kind = DP) :: gamma, lambda

    complex(kind = DP), dimension(:,:), allocatable :: rho0, HS, VS, VI

end module parameters
