module helper

    use constants
    implicit none

 contains

     function acomm(A, B)
         complex(kind=DP), dimension(:,:), intent(in) :: A, B
         complex(kind=DP), dimension(:,:), allocatable :: acomm
         integer :: n

         n = size(A, dim=1)
         allocate(acomm(n,n))
         acomm = matmul(A,B) + matmul(B,A)
     end function acomm

!---------------------------------------------------------------------!
!   This function will calculate the product exp(A)Bexp(-A), where    !
!   A, B are nxn matrixes, and k is the truncation level.             !
!---------------------------------------------------------------------!

     function BCH(A, B, k)
         integer, intent(in) :: k
         complex(kind=DP), dimension(:,:), intent(in) :: A, B
         complex(kind=DP), dimension(:,:), allocatable :: tmp1, tmp2, BCH
         integer :: l, j, n

         n = size(A, dim=1)
         allocate(tmp1(n,n))
         allocate(tmp2(n,n))
         allocate(BCH(n,n))
         tmp1 = 0
         BCH = 0
         if (k .eq. 0) then
             BCH = B
         else if (k .eq. 1) then
             BCH = B + comm(A, B)
         else
             tmp1 = B
             do l = 1, k
                 tmp2 = B
                 do j = 1, l
                     tmp2 = comm(A,tmp2) / j
                 end do
                 tmp1 = tmp1 + tmp2
             end do
             BCH = tmp1
         end if
     end function BCH

     function comm(A, B)
         complex(kind=DP), dimension(:,:), intent(in) :: A, B
         complex(kind=DP), dimension(:,:), allocatable :: comm
         integer :: n

         n = size(A, dim=1)
         allocate(comm(n,n))
         comm = matmul(A,B) - matmul(B,A)
     end function comm

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

     function euler(x)
         real(kind=DP), intent(in) :: x
         complex(kind=DP) :: euler

         euler = cmplx(cos(x), sin(x))
     end function euler

     function dag(A)
         complex(kind=DP), dimension(:,:), intent(in) :: A
         complex(kind=DP), dimension(:,:), allocatable :: dag
         integer :: n

         n = size(A, dim=1)
         allocate(dag(n,n))
         dag = transpose(conjg(A))
     end function dag

     function gexp(A)
         complex(kind=DP), dimension(:,:), intent(in) :: A
         complex(kind=DP), dimension(:,:), allocatable :: gexp
         integer :: n

         n = size(A, dim=1)
         allocate(gexp(n,n))
         gexp = idm(n) + A
     end function gexp

     function idm(n)
         integer, intent(in) :: n
         complex(kind=DP), dimension(n,n) :: id, idm
         integer :: i

         id = 0
         do i = 1, n
             id(i,i) = REAL1
         end do
         idm = id
     end function idm

     function IntOp(A, B)
        ! Computes the interaction picture operator
        !        e^{iA}Be^{-iA}
        ! using singular value decomposition.

        complex(kind=DP), dimension(:,:), intent(in) :: A, B
        complex(kind=DP), dimension(:,:), allocatable :: IntOp
        real(kind=DP), dimension(:), allocatable :: eigval
        complex(kind=DP), dimension(:,:), allocatable :: eigvect
        complex(kind=DP), dimension(:,:), allocatable :: diag1, diag2, tmp1, tmp2
        integer :: i, n

        n = size(A, dim=1)
        allocate(IntOp(n,n))
        allocate(diag1(n,n))
        allocate(diag2(n,n))
        allocate(tmp1(n,n))
        allocate(tmp2(n,n))

        diag1 = 0
        diag2 = 0

        call eigensystem(A,eigval,eigvect)

        ! Construct the diagonal matrices of e^{iA} and e^{-iA}
        do i = 1, n
            diag1(i,i) = euler(eigval(i))
            diag2(i,i) = euler(-eigval(i))
        end do

        ! Construct e^{iA} and e^{-iA}
        tmp1 = matmul(eigvect, matmul(diag1, dag(eigvect)))
        tmp2 = matmul(eigvect, matmul(diag2, dag(eigvect)))

        ! Construct e^{iA}Be^{-iA}
        IntOp = matmul(tmp1, matmul(B, tmp2))
     end function

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
