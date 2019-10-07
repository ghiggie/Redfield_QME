module parameters

!---------------------------------------------------------------------!
!    This file will declare the variables that the user must set.     !
!    This will also help in unifying the code.                        !
!---------------------------------------------------------------------!

    use constants
    implicit none

    real(kind = DP) :: dt1, dt2, time_limit, temp, gamma, lambda

    complex(kind = DP), dimension(:,:), allocatable :: rho0, HS, VI

end module parameters
