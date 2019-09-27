module helper

    use constants
    implicit none

 contains

     function dag(A)
         complex(kind=DP), dimension(:,:), intent(in) :: A
         complex(kind=DP), dimension(:,:), allocatable :: dag
         integer :: n

         n = size(A, dim=1)
         allocate(dag(n,n))
         dag = transpose(conjg(A))
     end function dag

     subroutine eigensystem(X,E,V)
         ! X is the matrix we want the eigensystem of
         ! E will store the eigenvalues
         ! V will store the eigenvectors

         complex(kind=DP), dimension(:,:), intent(in) :: X
         real(kind=DP), dimension(:), allocatable, intent(out) :: E
         complex(kind=DP), dimension(:,:), allocatable, intent(out) :: V

         integer, save :: info, lwork
         real(kind=DP), dimension(6) :: rwork
         complex(kind=DP), dimension(1000) :: work

         integer :: n

         n = size(X, dim=1)
         allocate(E(n))
         allocate(V(n,n))

         V = X
         lwork = -1
         call ZHEEV('V','U',n,V,n,E,work,lwork,rwork,info)
         lwork = min(1000, int(real(work(1))))
         call ZHEEV('V','U',n,V,n,E,work,lwork,rwork,info)
     end subroutine eigensystem

     function ExpOp(A, B)
        ! Computes the interaction picture operator
        !        e^{-iA}Be^{iA}
        ! using singular value decomposition.

        complex(kind=DP), dimension(:,:), intent(in) :: A, B
        complex(kind=DP), dimension(:,:), allocatable :: ExpOp
        real(kind=DP), dimension(:), allocatable :: eigval
        complex(kind=DP), dimension(:,:), allocatable :: eigvect
        complex(kind=DP), dimension(:,:), allocatable :: diag1, diag2, tmp1, tmp2
        integer :: i, n

        n = size(A, dim=1)
        allocate(ExpOp(n,n))
        allocate(diag1(n,n))
        allocate(diag2(n,n))
        allocate(tmp1(n,n))
        allocate(tmp2(n,n))

        diag1 = 0
        diag2 = 0

        call eigensystem(A,eigval,eigvect)

        ! Construct the diagonal matrices of e^{iA} and e^{-iA}
        do i = 1, n
            diag1(i,i) = cmplx(cos(-eigval(i)), sin(-eigval(i)))
            diag2(i,i) = cmplx(cos(eigval(i)), sin(eigval(i)))
        end do

        ! Construct e^{iA} and e^{-iA}
        tmp1 = matmul(eigvect, matmul(diag1, dag(eigvect)))
        tmp2 = matmul(eigvect, matmul(diag2, dag(eigvect)))

        ! Construct e^{-iA}Be^{iA}
        ExpOp = matmul(tmp1, matmul(B, tmp2))
     end function ExpOp

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
