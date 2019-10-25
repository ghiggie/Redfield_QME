module variables

!---------------------------------------------------------------------!
!    This file will declare the variables that the code will use.     !
!    This will also help in unifying the code.                        !
!---------------------------------------------------------------------!

    use constants
    implicit none

    !!! User Parameters set in params.cfg

    real(kind = DP) :: dt, time_limit, temp, gamma, lambda
    integer :: matsu
    complex(kind = DP), dimension(:,:), allocatable :: rho0, HS, VI

    !!! Variables needed by main.f95

    complex(kind=DP), dimension(:,:,:), allocatable :: rho
    complex(kind=DP), dimension(:,:), allocatable :: k1, k2, k3, k4
    complex(kind=DP), dimension(:,:), allocatable :: tmp1, tmp2, tmp3
    character(len=40) :: filename, arg, hostname
    integer :: ss, N, i, j, k
    real(kind=DP) :: ti, tc, tmp_val1, tmp_val2, tmp_val3
    logical :: tmpl_1, halt


    real(kind=DP), dimension(:), allocatable :: eigval, coeff
    complex(kind=DP), dimension(:,:), allocatable :: eigvect, diag, tmp_arr



contains

    subroutine array_init()

        allocate(rho0(ss,ss))
        allocate(HS(ss,ss))
        allocate(VI(ss,ss))
        allocate(k1(ss,ss))
        allocate(k2(ss,ss))
        allocate(k3(ss,ss))
        allocate(k4(ss,ss))
        allocate(tmp1(ss,ss))
        allocate(tmp2(ss,ss))
        allocate(tmp3(ss,ss))

        allocate(rho0(0:N,ss,ss))

        allocate(eigval(ss))
        allocate(eigvect(ss,ss))
        allocate(diag(ss,ss))
        allocate(tmp_arr(ss,ss))

        allocate(tmp(ss,ss))
        allocate(tmp1(ss,ss))
        allocate(tmp2(ss,ss))
        allocate(f2j(ss,ss))
        allocate(f2j1(ss,ss))
        allocate(f2j2(ss,ss))

        allocate(coeff(0:matsu))






    end subroutine array_init

end module variables
