MODULE variables

!---------------------------------------------------------------------!
!    This file will declare the variables that the code will use.     !
!    This will also help in unifying the code.                        !
!---------------------------------------------------------------------!
! Aside from global variables, this file should only be used for
! arrays. No atomic types unless necessarily global.
   USE constants

   IMPLICIT NONE

   !!! User Parameters set in params.cfg and through flags
   REAL(KIND=DP) :: dt, time_limit, time_write, temp, gamma, lambda
   LOGICAL       :: pade
   INTEGER       :: matsu, ss
   COMPLEX(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: rho0, HS, VI

   INTEGER                                       :: n_steps
   COMPLEX(KIND=DP)                              :: cinf
   COMPLEX(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: k1, k2, k3, k4
   COMPLEX(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: tmp_arr1, tmp_arr2, tmp_arr3, tmp_arr4
   COMPLEX(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: tmp_i2s1, tmp_i2s2
   REAL   (KIND=DP), DIMENSION(:),   ALLOCATABLE :: eigval, exp_vec
   COMPLEX(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: eigvect, diag
   COMPLEX(KIND=DP), DIMENSION(:),   ALLOCATABLE :: coeff
   COMPLEX(KIND=DP), DIMENSION(2,2)              :: diag2, HS_2, eigvect_2, tmp_arr22
   REAL   (KIND=DP), DIMENSION(2)                :: eigval_2
   COMPLEX(KIND=DP), DIMENSION(2,2)              :: rho_A, rho_B, rho_AE, rho_BE, eigvect_inv2, rho_g2
   COMPLEX(KIND=DP), DIMENSION(4,4)              :: eigvect_inv4

   ! Arrays for reduced density calculations
   COMPLEX(KIND=DP), DIMENSION(4,2) :: b1, b2
   COMPLEX(KIND=DP), DIMENSION(4,2) :: a1, a2

CONTAINS

   SUBROUTINE array_init()

      ALLOCATE(   rho0(ss,ss))
      ALLOCATE(     HS(ss,ss))
      ALLOCATE(     VI(ss,ss))
      ALLOCATE( eigval(ss))
      ALLOCATE(eigvect(ss,ss))
      ALLOCATE(   diag(ss,ss))

      ALLOCATE(k1(ss,ss))
      ALLOCATE(k2(ss,ss))
      ALLOCATE(k3(ss,ss))
      ALLOCATE(k4(ss,ss))

      ALLOCATE(tmp_arr1(ss,ss))
      ALLOCATE(tmp_arr2(ss,ss))
      ALLOCATE(tmp_arr3(ss,ss))
      ALLOCATE(tmp_arr4(ss,ss))

      ALLOCATE(tmp_i2s1(ss,ss))
      ALLOCATE(tmp_i2s2(ss,ss))

      ! Set up the unit vectors for the shitty reduced density
      b1      = 0
      b2      = 0
      b1(1,1) = 1
      b1(3,2) = 1
      b2(2,1) = 1
      b2(4,2) = 1

      a1      = 0
      a2      = 0
      a1(1,1) = 1
      a1(2,2) = 1
      a2(3,1) = 1
      a2(4,2) = 1

   END SUBROUTINE array_init

END MODULE variables
