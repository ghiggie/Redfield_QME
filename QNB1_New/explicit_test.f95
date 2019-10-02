program explicit_test

    use constants
    use parameters
    use helper
    use bath

    implicit none

    complex(kind=DP), dimension(:,:,:), allocatable :: rho
    integer :: i, N

    rho0(1, :) = /(cmplx(1., 0), cmplx(0, 0)/)
    rho0(2, :) = /(cmplx(0, 0), cmplx(0, 0)/)
    HS(1, :) = /(cmplx(0.5, 0), cmplx(0, 0)/)
    HS(2, :) = /(cmplx(0, 0), cmplx(-0.5, 0)/)
    VI(1, :) = /(cmplx(0, 0), cmplx(1., 0)/)
    VI(2, :) = /(cmplx(1., 0), cmplx(0, 0)/)

    time_limit = 100
    dt1 = 0.1
    dt2 = 0.1
    temp = 1
    gamma = 0.1
    lambda = 0.1

    N = int(time_limit / dt1)
    allocate(rho(0 : N, 2, 2))
    rho(0, :, :) = rho0

    open(10,file='rho.dat')
    write(10,'(f10.3,8(e15.6))') 0.0,((rho(0,j,k),j=1,2),k=1,2)

    do i = 0, N - 1
        ti = i * dt1

        
