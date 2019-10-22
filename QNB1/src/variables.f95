module variables

!---------------------------------------------------------------------!
!    This file will declare the variables that the code will use.     !
!    This will also help in unifying the code.                        !
!---------------------------------------------------------------------!

    use constants
    implicit none

    !!! User Parameters set in params.cfg

    real(kind = DP) :: dt1, dt2, time_limit, temp, gamma, lambda
    complex(kind = DP), dimension(:,:), allocatable :: rho0, HS, VI

    !!! Variables needed by main.f95

    complex(kind=DP), dimension(:,:,:), allocatable :: rho
    complex(kind=DP), dimension(:,:), allocatable :: k1, k2, k3, k4
    complex(kind=DP), dimension(:,:), allocatable :: tmp1, tmp2, tmp3
    character(len=40) :: filename, arg, hostname
    integer :: S, N, i, j, k
    real(kind=DP) :: ti, tc, tmp_val1, tmp_val2, tmp_val3
    logical :: tmpl, halt






end module variables
