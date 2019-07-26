module constants

!---------------------------------------------------------------------!
!   This file declares some common constants.                         !
!---------------------------------------------------------------------!

    implicit none

    integer, parameter :: DP=selected_real_kind(15,307)
    real(DP), parameter :: PI = 4*atan(1._DP)
    real(DP), parameter :: tol=1.0E-10
    complex(kind=DP) :: IMAG1=(0.0_DP, 1.0_DP), ZERO=(0.0_DP, 0.0_DP)
    complex(kind=DP) :: REAL1=(1.0_DP,0.0_DP)

end module constants
