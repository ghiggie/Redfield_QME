MODULE constants

!---------------------------------------------------------------------!
!   This file declares some common constants.                         !
!---------------------------------------------------------------------!

   IMPLICIT NONE

   INTEGER,           PARAMETER :: DP          = SELECTED_REAL_KIND(15,307)
   REAL(KIND=DP),     PARAMETER :: PI          = 4*ATAN(1._DP)
   REAL(KIND=DP),     PARAMETER :: tol         = 1.0E-6
   CHARACTER(LEN=10), PARAMETER :: version     = 'v1.2.6'
   CHARACTER(LEN=20), PARAMETER :: last_update = 'February 17, 2020'
   CHARACTER(LEN=40), PARAMETER :: coded_by    = 'Garrett Higginbotham'

END MODULE constants
