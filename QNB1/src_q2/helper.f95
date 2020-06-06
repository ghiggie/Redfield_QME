MODULE helper

   USE constants
   USE variables

   IMPLICIT NONE

   ! Work area for zheev
   INTEGER,          SAVE            :: info, lwork
   REAL   (KIND=DP), DIMENSION(6)    :: rwork
   COMPLEX(KIND=DP), DIMENSION(1000) :: work

CONTAINS

      SUBROUTINE eigensystem(X,E,V)
         ! X is the matrix we want the eigensystem of
         ! E will store the eigenvalues
         ! V will store the eigenvectors

         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(IN)  :: X
         REAL   (KIND=DP), DIMENSION(:),   INTENT(OUT) :: E
         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(OUT) :: V
         INTEGER                                       :: n

         n = SIZE(X, DIM=1)
         V = X
         lwork = -1
         CALL ZHEEV('V','U',n,V,n,E,work,lwork,rwork,info)
         lwork = MIN(1000, INT(REAL(work(1))))
         CALL ZHEEV('V','U',n,V,n,E,work,lwork,rwork,info)
      END SUBROUTINE eigensystem

      FUNCTION Entropy(A)
         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(IN) :: A
         REAL   (KIND=DP)                             :: Entropy

         CALL eigensystem(A, eigval, eigvect)

         Entropy = -SUM(eigval * LOG(eigval), MASK=eigval>0.0_DP)
      END FUNCTION Entropy

      FUNCTION fidelity(A,B)
         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(IN) :: A, B
         REAL   (KIND=DP)                             :: fidelity

         INTEGER :: i

         diag = 0.0_DP

         ! Compute SQRT(A)
         CALL eigensystem(A, eigval, eigvect)
         DO i = 1, ss
            diag(i,i) = DSQRT(REAL(eigval(i), KIND=DP))
         END DO
         tmp_arr1 = CONJG(TRANSPOSE(eigvect))
         tmp_arr1 = MATMUL(eigvect, MATMUL(diag, tmp_arr1)) ! SQRT(A)

         tmp_arr1 = MATMUL(tmp_arr1, MATMUL(B, tmp_arr1)) !SQRT(A)*B*SQRT(A)
         CALL eigensystem(tmp_arr1, eigval, eigvect)
         DO i = 1, ss
            diag(i,i) = DSQRT(REAL(eigval(i), KIND=DP))
         END DO
         tmp_arr1 = CONJG(TRANSPOSE(eigvect))
         tmp_arr1 = MATMUL(eigvect, MATMUL(diag, tmp_arr1)) !SQRT(SQRT(A)*B*SQRT(A))

         fidelity = (REAL(trace(tmp_arr1), KIND=DP))**2
      END FUNCTION fidelity

      FUNCTION fidelity2(A,B)
         COMPLEX(KIND=DP), DIMENSION(2,2), INTENT(IN) :: A, B
         REAL   (KIND=DP)                             :: fidelity2

         INTEGER :: i

         diag2 = 0.0_DP
         ! Compute SQRT(A)
         CALL eigensystem(A, eigval_2, eigvect_2)
         DO i = 1, 2
            diag2(i,i) = DSQRT(REAL(eigval_2(i), KIND=DP))
         END DO
         tmp_arr22 = CONJG(TRANSPOSE(eigvect_2))
         tmp_arr22 = MATMUL(eigvect_2, MATMUL(diag2, tmp_arr22)) ! SQRT(A)
         tmp_arr22 = MATMUL(tmp_arr22, MATMUL(B, tmp_arr22)) !SQRT(A)*B*SQRT(A)
         CALL eigensystem(tmp_arr22, eigval_2, eigvect_2)
         DO i = 1, 2
            diag2(i,i) = DSQRT(REAL(eigval_2(i), KIND=DP))
         END DO
         tmp_arr22 = CONJG(TRANSPOSE(eigvect_2))
         tmp_arr22 = MATMUL(eigvect_2, MATMUL(diag2, tmp_arr22)) !SQRT(SQRT(A)*B*SQRT(A))

         fidelity2 = (REAL(trace(tmp_arr22), KIND=DP))**2
      END FUNCTION fidelity2

      SUBROUTINE I2S(A, B, C)
         ! Computes
         !        e^{-iA}Be^{iA}
         ! using singular value decomposition.

         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(IN)  :: A, B
         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(OUT) :: C
         INTEGER                                       :: i

         diag = 0.0_DP

         CALL eigensystem(A,eigval,eigvect)

         ! Construct the diagonal matrix of e^{-iA}
         DO i = 1, ss
            diag(i,i) = CMPLX(DCOS(eigval(i)), DSIN(-eigval(i)))
         END DO

         ! Construct e^{-iA}
         tmp_i2s1 = CONJG(TRANSPOSE(eigvect))
         tmp_i2s1 = MATMUL(eigvect, MATMUL(diag, tmp_i2s1))
         tmp_i2s2 = CONJG(TRANSPOSE(tmp_i2s1))

         ! Construct e^{-iA}Be^{iA}
         C = MATMUL(tmp_i2s1, MATMUL(B, tmp_i2s2))
      END SUBROUTINE I2S

      ! Returns the inverse of a matrix calculated by finding the LU
      ! decomposition.  Depends on LAPACK.
      SUBROUTINE matinv(A, Ainv)
         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(IN)      :: A
         COMPLEX(KIND=DP), DIMENSION(SIZE(A,1), SIZE(A,2)) :: Ainv

         REAL(KIND=DP), DIMENSION(SIZE(A,1)) :: work   ! work array for LAPACK
         INTEGER,       DIMENSION(SIZE(A,1)) :: ipiv   ! pivot indices
         INTEGER                             :: n, info

         ! External procedures defined in LAPACK
         EXTERNAL ZGETRF
         EXTERNAL ZGETRI

         ! Store A in Ainv to prevent it from being overwritten by LAPACK
         Ainv = A
         n    = SIZE(A,1)

         ! ZGETRF computes an LU factorization of a general M-by-N matrix A
         ! using partial pivoting with row interchanges.
         CALL ZGETRF(n, n, Ainv, n, ipiv, info)

         IF (info /= 0) THEN
            stop 'Matrix is numerically singular!'
         END IF

         ! ZGETRI computes the inverse of a matrix using the LU factorization
         ! computed by ZGETRF.
         CALL ZGETRI(n, Ainv, n, ipiv, work, n, info)

         IF (info /= 0) THEN
            stop 'Matrix inversion failed!'
         END IF
      END SUBROUTINE matinv

      SUBROUTINE matinv2(A,B)
         !! Performs a direct calculation of the inverse of a 2x2 matrix.
         COMPLEX(KIND=DP), INTENT(IN) :: A(2,2)   !! Matrix
         COMPLEX(KIND=DP)             :: B(2,2)   !! Inverse matrix
         COMPLEX(KIND=DP)             :: detinv

         ! Calculate the inverse determinant of the matrix
         detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

         ! Calculate the inverse of the matrix
         B(1,1) = +detinv * A(2,2)
         B(2,1) = -detinv * A(2,1)
         B(1,2) = -detinv * A(1,2)
         B(2,2) = +detinv * A(1,1)
      END SUBROUTINE matinv2

      SUBROUTINE matinv4(A,B)
         !! Performs a direct calculation of the inverse of a 4×4 matrix.
         COMPLEX(KIND=DP), INTENT(IN) :: A(4,4)   !! Matrix
         COMPLEX(KIND=DP)             :: B(4,4)   !! Inverse matrix
         COMPLEX(KIND=DP)             :: detinv

         ! Calculate the inverse determinant of the matrix
         detinv = &
            1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)&
            -A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
            -A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)&
            -A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
            +A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)&
            -A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)&
            *(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))&
            +A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

         ! Calculate the inverse of the matrix
         B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*&
            A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
         B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*&
            A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
         B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*&
            A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
         B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*&
            A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
         B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*&
            A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
         B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*&
            A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
         B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*&
            A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
         B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*&
            A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
         B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*&
            A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
         B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*&
            A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
         B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*&
            A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
         B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*&
            A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
         B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*&
            A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
         B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*&
            A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
         B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*&
            A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
         B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*&
            A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
      END SUBROUTINE matinv4

      ! Right now, this only works for a four dimensional Hilbert space
      ! that we want to decompose into two two-dimensional Hilbert spaces.
      SUBROUTINE rhoA(A, B)
         COMPLEX(KIND=DP), DIMENSION(4,4), INTENT(IN)  :: A
         COMPLEX(KIND=DP), DIMENSION(2,2), INTENT(OUT) :: B

         B = MATMUL(TRANSPOSE(CONJG(b1)),MATMUL(A,b1)) + MATMUL(TRANSPOSE(CONJG(b2)),MATMUL(A,b2))
      END SUBROUTINE rhoA

      SUBROUTINE rhoB(A, B)
         COMPLEX(KIND=DP), DIMENSION(4,4), INTENT(IN)  :: A
         COMPLEX(KIND=DP), DIMENSION(2,2), INTENT(OUT) :: B

         B = MATMUL(TRANSPOSE(CONJG(a1)),MATMUL(A,a1)) + MATMUL(TRANSPOSE(CONJG(a2)),MATMUL(A,a2))
      END SUBROUTINE rhoB

      FUNCTION test_hermitian(A)
         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(IN) :: A
         LOGICAL                                      :: test_hermitian
         LOGICAL                                      :: test
         INTEGER                                      :: i, j, n

         n    = SIZE(A, DIM=1)
         test = .TRUE.
         DO i = 1, n
            IF (ABS(AIMAG(A(i,i))) > tol) THEN
               test = .FALSE.
            END IF
            DO j = i+1, n
               IF (ABS(A(i,j)-CONJG(A(j,i)))> tol) THEN
                  test = .FALSE.
               END IF
            END DO
         END DO
         test_hermitian = test
      END FUNCTION test_hermitian

      FUNCTION test_positivity(A)
         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(IN) :: A
         LOGICAL                                      :: test_positivity
         INTEGER                                      :: i, n

         n = SIZE(A, DIM=1)
         CALL eigensystem(A, eigval, eigvect)
         test_positivity = .TRUE.

         DO i = 1, n
            IF (eigval(i) < 0.0_DP) test_positivity = .FALSE.
         END DO
      END FUNCTION test_positivity

      FUNCTION test_trace(A)
         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(IN) :: A
         LOGICAL                                      :: test_trace
         REAL(KIND=DP)                                :: tr

         test_trace = .TRUE.
         tr         = DSQRT(REAL(trace(A)*CONJG(trace(A)), KIND=DP))
         IF (ABS(tr - 1.0_DP) > tol) test_trace = .FALSE.
      END FUNCTION test_trace

      FUNCTION time_stamp()
         !---------------------------------------
         ! this returns current time in a string
         !---------------------------------------
         character(18) :: time_stamp
         character(8)  :: date
         character(10) :: time

         CALL DATE_AND_TIME(date,time)

         time_stamp = date(5:6) // '/' // date(7:8) // '/' //  date(1:4) &
                              // '  ' // time(1:2) // ':' // time(3:4)
      END FUNCTION time_stamp

      FUNCTION trace(A)
         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(IN) :: A
         COMPLEX(KIND=DP)                             :: trace
         INTEGER                                      :: n, i

         n     = SIZE(A, DIM=1)
         trace = 0.0_DP
         DO i = 1, n
            trace = trace + A(i,i)
         END DO
      END FUNCTION trace

      FUNCTION trace_distance(A, B)
         COMPLEX(KIND=DP), DIMENSION(:,:), INTENT(IN) :: A, B
         COMPLEX(KIND=DP)                             :: trace_distance
         INTEGER                                      :: i

         tmp_arr1 = MATMUL(CONJG(TRANSPOSE(A - B)), A - B)
         diag = 0.0_DP

         CALL eigensystem(tmp_arr1, eigval, eigvect)
         DO i = 1, ss
            diag(i,i) = DSQRT(REAL(eigval(i), KIND=DP))
         END DO
         tmp_arr1 = CONJG(TRANSPOSE(eigvect))
         tmp_arr1 = MATMUL(eigvect, MATMUL(diag, tmp_arr1))

         trace_distance = 0.5_DP * trace(tmp_arr1)
      END FUNCTION trace_distance
 
      FUNCTION trace_distance2(A, B)
         COMPLEX(KIND=DP), DIMENSION(2,2), INTENT(IN) :: A, B
         COMPLEX(KIND=DP)                             :: trace_distance2
         INTEGER                                      :: i

         tmp_arr22 = MATMUL(CONJG(TRANSPOSE(A - B)), A - B)
         diag2 = 0.0_DP

         CALL eigensystem(tmp_arr22, eigval_2, eigvect_2)
         DO i = 1, 2
            diag2(i,i) = DSQRT(REAL(eigval_2(i), KIND=DP))
         END DO
         tmp_arr22 = CONJG(TRANSPOSE(eigvect_2))
         tmp_arr22 = MATMUL(eigvect_2, MATMUL(diag2, tmp_arr22))

         trace_distance2 = 0.5_DP * trace(tmp_arr22)
      END FUNCTION trace_distance2

END MODULE helper
