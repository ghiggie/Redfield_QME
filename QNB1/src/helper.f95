module helper

    use constants
    use variables
    implicit none

    ! Work area for zheev
    integer, save :: info, lwork
    real(kind=DP), dimension(6) :: rwork
    complex(kind=DP), dimension(1000) :: work

 contains

     subroutine eigensystem(X,E,V)
         ! X is the matrix we want the eigensystem of
         ! E will store the eigenvalues
         ! V will store the eigenvectors

         complex(kind=DP), dimension(:,:), intent(in) :: X
         real(kind=DP), dimension(:), intent(out) :: E
         complex(kind=DP), dimension(:,:), intent(out) :: V
         integer :: n

         n = size(X, dim=1)
         V = X
         lwork = -1
         call ZHEEV('V','U',n,V,n,E,work,lwork,rwork,info)
         lwork = min(1000, int(real(work(1))))
         call ZHEEV('V','U',n,V,n,E,work,lwork,rwork,info)
     end subroutine eigensystem

     ! This entropy calculation will fail if any eigenvalues are 0
     function Entropy(A)
         complex(kind=DP), dimension(:,:), intent(in) :: A
         real(kind=DP) :: Entropy

         call eigensystem(A, eigval, eigvect)

         Entropy = -SUM(eigval * LOG(eigval), mask=eigval>0.0_DP)
     end function Entropy

     subroutine I2S(A, B, C)
        ! Computes
        !        e^{-iA}Be^{iA}
        ! using singular value decomposition.

        complex(kind=DP), dimension(:,:), intent(in) :: A, B
        complex(kind=DP), dimension(:,:), intent(out) :: C
        integer :: i

        diag = 0

        call eigensystem(A,eigval,eigvect)

        ! Construct the diagonal matrix of e^{-iA}
        do i = 1, ss
            diag(i,i) = cmplx(cos(eigval(i)), sin(-eigval(i)))
        end do

        ! Construct e^{-iA}
        tmp_i2s1 = conjg(transpose(eigvect))
        tmp_i2s1 = matmul(eigvect, matmul(diag, tmp_i2s1))
        tmp_i2s2 = conjg(transpose(tmp_i2s1))

        ! Construct e^{-iA}Be^{iA}
        C = matmul(tmp_i2s1, matmul(B, tmp_i2s2))
    end subroutine I2S

     ! Right now, this only works for a four dimensional Hilbert space
     ! that we want to decompose into two two-dimensional Hilbert spaces.
     subroutine rhoA(A, B)
         complex(kind=DP), dimension(:,:), intent(in) :: A
         complex(kind=DP), dimension(2,2), intent(out) :: B

         B = matmul(transpose(conjg(b1)),matmul(A,b1)) + matmul(transpose(conjg(b2)),matmul(A,b2))
     end subroutine rhoA

     subroutine rhoB(A, B)
         complex(kind=DP), dimension(:,:), intent(in) :: A
         complex(kind=DP), dimension(2,2), intent(out) :: B

         B = matmul(transpose(conjg(a1)),matmul(A,a1)) + matmul(transpose(conjg(a2)),matmul(A,a2))
     end subroutine rhoB

     function test_hermitian(A)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        logical :: test_hermitian
        logical :: test
        integer :: i, j, n

        n = size(A, dim=1)
        test = .true.
        do i=1,n
            if (abs(aimag(A(i,i))) > tol) then
                test = .false.
            end if
            do j=i+1,n
                if(abs(A(i,j)-conjg(A(j,i)))> tol) then
                    test = .false.
                end if
            end do
        end do
        test_hermitian = test
     end function test_hermitian

     function test_positivity(A)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        logical :: test_positivity
        integer :: i, n

        n = size(A, dim=1)
        call eigensystem(A, eigval, eigvect)
        test_positivity = .true.

        do i = 1, n
            if (eigval(i) < 0.0_DP) test_positivity = .false.
        end do
     end function test_positivity

     function test_trace(A)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        logical :: test_trace
        real(kind=DP) :: tr

        test_trace = .true.
        tr = sqrt(trace(A)*conjg(trace(A)))
        if (abs(tr-1.0_DP) > tol) test_trace = .false.
     end function test_trace

     function time_stamp()
         !---------------------------------------
         ! this returns current time in a string
         !---------------------------------------
         implicit none

         character(18) :: time_stamp
         character(8)  :: date
         character(10) :: time

         call date_and_time(date,time)

         time_stamp = date(5:6) // '/' // date(7:8) // '/' //  date(1:4) &
                              // '  ' // time(1:2) // ':' // time(3:4)
      end function time_stamp

      function trace(A)
         complex(kind=DP), dimension(:,:), intent(in) :: A
         complex(kind=DP) :: trace
         integer :: n, i

         n = size(A, dim=1)
         trace = 0
         do i = 1, n
           trace = trace + A(i,i)
         end do
      end function trace

end module helper
