module parameters

!---------------------------------------------------------------------!
!    This file will declare the variables that the user must set.     !
!    This will also help in unifying the code.                        !
!---------------------------------------------------------------------!

    use constants
    implicit none

    integer :: N, trunc
    real(kind = DP) :: time_limit, temp
    real(kind = DP) :: gamma, lambda

    complex(kind = DP), dimension(:,:), allocatable :: rho0, HS, VI

end module parameters
