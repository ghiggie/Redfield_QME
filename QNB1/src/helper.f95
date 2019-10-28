module helper

    use constants
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

         Entropy = -SUM(eigval * LOG(eigval))
     end function Entropy

     ! This subroutine is better, but still needs to be improved. Come back later
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
        ! Remember that ss is a global variable for the dimension of the system
        do i = 1, ss
            diag(i,i) = cmplx(cos(-eigval(i)), sin(-eigval(i)))
        end do

        ! Construct e^{-iA}
        tmp_arr = matmul(eigvect, matmul(diag, transpose(conjg(eigvect))))

        ! Construct e^{-iA}Be^{iA}
        C = matmul(tmp_arr, matmul(B, transpose(conjg(tmp_arr))))
    end subroutine I2S

     ! Right now, this only works for a four dimensional Hilbert space
     ! that we want to decompose into two two-dimensional Hilbert spaces.
     function rhoA(A)
         complex(kind=DP), dimension(:,:), intent(in) :: A
         complex(kind=DP), dimension(2,2) :: rhoA
         complex(kind=DP), dimension(4,2) :: b1, b2

         b1 = 0
         b2 = 0
         b1(1,1) = 1
         b1(3,2) = 1
         b2(2,1) = 1
         b2(4,2) = 1

         rhoA = matmul(transpose(conjg(b1),matmul(A,b1)) + matmul(transpose(conjg(b2),matmul(A,b2))
     end function rhoA

     function rhoB(A)
         complex(kind=DP), dimension(:,:), intent(in) :: A
         complex(kind=DP), dimension(2,2) :: rhoB
         complex(kind=DP), dimension(4,2) :: a1, a2

         a1 = 0
         a2 = 0
         a1(1,1) = 1
         a1(2,2) = 1
         a2(3,1) = 1
         a2(4,2) = 1

         rhoB = matmul(transpose(conjg(a1),matmul(A,a1)) + matmul(transpose(conjg(a2),matmul(A,a2))
     end function rhoB

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
        complex(kind=DP), dimension(:,:), allocatable :: eigvect
        real(kind=DP), dimension(:), allocatable :: eigval
        integer :: i, n

        n = size(A, dim=1)
        call eigensystem(A, eigval, eigvect)
        test_positivity = .true.

        do i = 1, n
            if (eigval(i) < 0.) test_positivity = .false.
        end do
     end function test_positivity

     function test_trace(A)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        logical :: test_trace
        logical :: test
        real(kind=DP) :: tr

        test = .true.
        tr = sqrt(trace(A)*conjg(trace(A)))
        if (abs(tr-1.0_DP) > tol) then
            test = .false.
        end if
        test_trace = test
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
        complex(kind=DP) :: tr
        integer :: n, i

        n = size(A, dim=1)
        tr = 0
        do i = 1, n
          tr = tr + A(i,i)
        end do

        trace = tr
     end function trace

end module helper
