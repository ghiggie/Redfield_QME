module helper

    use constants
    implicit none

 contains

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

    function gexp(A)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        complex(kind=DP), dimension(:,:) :: gexp
        integer :: n

        n = size(A, dim=1)
        allocate(gexp(n,n))
        gexp = idm(n) + A
    end function gexp

    function comm(A, B)
        complex(kind=DP), dimension(:,:), intent(in) :: A, B
        complex(kind=DP), dimension(:,:) :: comm
        integer :: n

        n = size(A, dim=1)
        allocate(comm(n,n))
        comm = matmul(A,B) - matmul(B,A)
    end function comm

    function acomm(A, B)
        complex(kind=DP), dimension(:,:), intent(in) :: A, B
        complex(kind=DP), dimension(:,:) :: acomm
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
        complex(kind=DP), dimension(:,:) :: tmp1, tmp2, BCH
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

    function euler(x)
        real(kind=DP), intent(in) :: x
        complex(kind=DP) :: euler

        euler = REAL1*cos(x) + IMAG1*sin(x)
    end function euler

    function dag(A)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        complex(kind=DP), dimension(:,:) :: dag
        integer :: n

        n = size(A, dim=1)
        allocate(dag(n,n))
        dag = transpose(conjg(A))
    end function dag

    function trace(A)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        complex(kind=DP) :: trace
        complex(kind=DP) :: tr
        integer :: n

        n = size(A, dim=1)
        tr = 0
        do i = 1, n
          tr = tr + A(i,i)
        end do

        trace = tr
    end function trace

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
        if (abs(tr-1) > tol) then
            test = .false.
        end if
        test_trace = test
    end function test_trace

end module helper
