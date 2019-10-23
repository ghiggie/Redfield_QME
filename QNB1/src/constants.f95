module constants

!---------------------------------------------------------------------!
!   This file declares some common constants.                         !
!---------------------------------------------------------------------!

    implicit none

    integer, parameter :: DP=selected_real_kind(15,307)
    real(DP), parameter :: PI = 4*atan(1._DP)
    real(DP), parameter :: tol=1.0E-6
    character(len=6), parameter :: version='v1.0.1'

end module constants
