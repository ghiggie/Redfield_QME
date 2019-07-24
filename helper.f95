module helper

    use constants
    implicit none

 contains

    complex function bath_corr(temp, gamma)
        real(kind=DP), intent(in) :: temp, gamma
        real(kind=DP), dimension(2) :: gj, cj
        real(kind=DP) :: cotan, delta

        ! Matsubara frequency (n=1)
        gj(1) = gamma
        gj(2) = (2.0_DP*pi)*temp
        cotan = cos(0.5_DP*gj(1)/temp)/sin(0.5*gj(1)/temp)
        delta = 2.0_DP/(gj(1)/temp) - cotan - 4.0_DP*gj(1)*temp/(gj(2)**2-gj(1)**2)
        cj(1) = gj(1) * cotan
        cj(2) = 4.0_DP*gj(1)*temp * gj(2)/(gj(2)**2-gj(1)**2)

        bath_corr = cj(1) * REAL1 + cj(2) * IMAG1
    end function bath_corr

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

    complex function euler(x)
        real(kind=DP), intent(in) :: x

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

    complex function trace(A)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        complex(kind=DP) :: tr
        integer :: n

        n = size(A, dim=1)
        tr = 0
        do i = 1, n
          tr = tr + A(i,i)
        end do

        trace = tr
    end function trace

end module helper