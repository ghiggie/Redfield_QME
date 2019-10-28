module variables

!---------------------------------------------------------------------!
!    This file will declare the variables that the code will use.     !
!    This will also help in unifying the code.                        !
!---------------------------------------------------------------------!

    use constants
    implicit none

    !!! User Parameters set in params.cfg

    real(kind = DP) :: dt, time_limit, temp, gamma, lambda
    logical :: pade
    integer :: matsu
    complex(kind = DP), dimension(:,:), allocatable :: rho0, HS, VI

    !!! Variables needed by main.f95

    complex(kind=DP), dimension(:,:,:), allocatable :: rho, bath_VI
    complex(kind=DP), dimension(:,:), allocatable :: k1, k2, k3, k4
    complex(kind=DP), dimension(:,:), allocatable :: tmp1, tmp2, tmp3
    character(len=40) :: filename, arg, hostname
    integer :: ss, N, i, j, k
    real(kind=DP) :: ti, tc, tmp_val1, tmp_val2, tmp_val3
    logical :: tmpl_1


    real(kind=DP), dimension(:), allocatable :: eigval, exp_vec
    complex(kind=DP), dimension(:,:), allocatable :: eigvect, diag, tmp_arr
    complex(kind=DP), dimension(:), allocatable :: coeff



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

        if (pade) then
            matsu = 1
            allocate(coeff(0:matsu))
            allocate(exp_vec(0:matsu))
        else
            allocate(coeff(0:matsu))
            allocate(exp_vec(0:matsu))
        end if






    end subroutine array_init

end module variables
