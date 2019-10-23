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
    logical :: tmpl_1, halt


    real(kind=DP), dimension(:), allocatable :: eigval
    complex(kind=DP), dimension(:,:), allocatable :: eigvect, diag, tmp_arr



contains

    subroutine array_init()

        allocate(rho0(S,S))
        allocate(HS(S,S))
        allocate(VI(S,S))
        allocate(k1(S,S))
        allocate(k2(S,S))
        allocate(k3(S,S))
        allocate(k4(S,S))
        allocate(tmp1(S,S))
        allocate(tmp2(S,S))
        allocate(tmp3(S,S))

        allocate(rho0(0:N,S,S))

        allocate(eigval(S))
        allocate(eigvect(S,S))
        allocate(diag(S,S))
        allocate(tmp_arr(S,S))

        allocate(tmp(S,S))
        allocate(tmp1(S,S))
        allocate(tmp2(S,S))
        allocate(f2j(S,S))
        allocate(f2j1(S,S))
        allocate(f2j2(S,S))






    end subroutine array_init

end module variables
