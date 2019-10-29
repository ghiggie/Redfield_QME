module variables

!---------------------------------------------------------------------!
!    This file will declare the variables that the code will use.     !
!    This will also help in unifying the code.                        !
!---------------------------------------------------------------------!
! Aside from global variables, this file should only be used for
! arrays. No atomic types unless necessarily global.
    use constants
    implicit none

    !!! User Parameters set in params.cfg and through flags
    real(kind = DP) :: dt, time_limit, temp, gamma, lambda
    logical :: pade
    integer :: matsu, ss
    complex(kind = DP), dimension(:,:), allocatable :: rho0, HS, VI

    integer :: n_steps
    complex(kind=DP) :: cinf
    complex(kind=DP), dimension(:,:), allocatable :: k1, k2, k3, k4
    complex(kind=DP), dimension(:,:), allocatable :: tmp_arr1, tmp_arr2, tmp_arr3, tmp_arr4
    real(kind=DP), dimension(:), allocatable :: eigval, exp_vec
    complex(kind=DP), dimension(:,:), allocatable :: eigvect, diag
    complex(kind=DP), dimension(:), allocatable :: coeff

    ! Arrays for reduced density calculations
    complex(kind=DP), dimension(4,2) :: b1, b2
    complex(kind=DP), dimension(4,2) :: a1, a2

contains

    subroutine array_init()

        allocate(rho0(ss,ss))
        allocate(HS(ss,ss))
        allocate(VI(ss,ss))
        allocate(eigval(ss))
        allocate(eigvect(ss,ss))

        allocate(k1(ss,ss))
        allocate(k2(ss,ss))
        allocate(k3(ss,ss))
        allocate(k4(ss,ss))

        allocate(tmp_arr1(ss,ss))
        allocate(tmp_arr2(ss,ss))
        allocate(tmp_arr3(ss,ss))
        allocate(tmp_arr4(ss,ss))

        if (pade) then
            matsu = 1
            allocate(coeff(0:matsu))
            allocate(exp_vec(0:matsu))
        else
            allocate(coeff(0:matsu))
            allocate(exp_vec(0:matsu))
        end if

        ! Set up the unit vectors for the shitty reduced density
        b1 = 0
        b2 = 0
        b1(1,1) = 1
        b1(3,2) = 1
        b2(2,1) = 1
        b2(4,2) = 1

        a1 = 0
        a2 = 0
        a1(1,1) = 1
        a1(2,2) = 1
        a2(3,1) = 1
        a2(4,2) = 1

    end subroutine array_init

end module variables
