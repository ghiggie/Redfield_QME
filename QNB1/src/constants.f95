module constants

!---------------------------------------------------------------------!
!   This file declares some common constants.                         !
!---------------------------------------------------------------------!

    implicit none

    integer, parameter :: DP=selected_real_kind(15,307)
    real(DP), parameter :: PI = 4*atan(1._DP)
    real(DP), parameter :: tol=1.0E-6
    character(len=10), parameter :: version='v1.1.1'
    character(len=20), parameter :: last_update='November 15, 2019'
    character(len=40), parameter :: coded_by='Garrett Higginbotham'

end module constants
