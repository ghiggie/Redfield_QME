program main

    use constants
    use variables
    use helper
    use bath

    implicit none

    namelist/params/dt,time_limit,pade,matsu,temp,gamma,lambda,rho0,HS,VI

    n_arg = iargc()
    if (n_arg .eq. 0) then
        write(*,*) ' System size and configuration file are required.'
        write(*,*) ' Usage: ./BornMarkov1B [-h] [-v] system_size config_file'
    else
        do i = 1, n_arg
            call getarg(i, arg)
            select case(trim(arg))
            case('-v')
                write(*,*) 'Coded by ', trim(coded_by)
                write(*,*) 'Last Updated: ', trim(last_update)
                write(*,*) 'Born Markov with one bath: ', trim(version)
                STOP ''
            case('-h')
                write(*,*) ' Usage: ./BornMarkov1B [-h] [-v] system_size config_file'
                write(*,*) ' Required Parameters:'
                write(*,*) '    system_size: size of the system Hilbert space'
                write(*,*) '    config_file: parameter values in namelist format'
                write(*,*) ' Optional Parameters:'
                write(*,*) '    -h: show this usage information'
                write(*,*) '    -v: show the version number'
                STOP ''
            case default
                if (i .eq. 1) then
                    read(arg,*) ss
                else if (i .eq. 2) then
                    filename = trim(arg)
                else
                    write(*,*) ' Too many parameters'
                    write(*,*) ' Usage: ./BornMarkov1B [-h] [-v] system_size config_file'
                end if
            end select
        end do
    end if

    ! Initialize the single time-step arrays
    call array_init()

    open(10, file=trim(filename))
    read(10,nml=params)
    close(10)

    n_steps = nint(time_limit / dt)

    ! Initialize the multi time-step arrays
    allocate(rho0(0:n_steps,ss,ss))
    rho(0,:,:) = rho0 ! Initialize the storage
    allocate(bath_VI(0:n_steps,ss,ss))
    ! Create the arrays needed for the bath correlation calculation
    call bc_coeff()
    ! Create the array of values for lambda_bc, to be stored in bath_VI
    call lambda_bc(VI, bath_VI)

    ! Set up the summary file
    open(20, file = 'BornMarkov1B.out')
    write(20, '(3a)') '*** Born Markov with One Bath (BornMarkov1B)', version, '***'

    write(20, '(/a,a)') 'Job begun at ', time_stamp()
    call get_environment_variable('HOSTNAME', hostname)
    write(20, '(/a,a)') 'Host: ', trim(hostname)

    write(20,'(/a)') 'System Hamiltonian:'

    tmpl_1 = test_hermitian(HS)

    if (tmpl_1) then
        write(20,'(/a)') 'Self-Adjoint.'
    else
        write(20,'(/a)') '**** NOT Self-Adjoint. ****'
    end if

    do i=1,ss
        write(20,'(4(a,f10.7,a,f10.7,a))') ('(',REAL(HS(i,j)),',',AIMAG(HS(i,j)),'),',j=1,ss)
    end do

    open(10,file='rho.dat')
    write(10,'(f10.3,8(e15.6))') 0.0,((rho(0,j,k),j=1,ss),k=1,ss)

    do i = 0, n_steps - 1
        ti = i * dt

        tc = ti
        tmp1 = rho(i,:,:)
        tmp2 = bath_VI(i,:,:)
        tmp3 = matmul(tmp2, tmp1) - matmul(tmp1, transpose(conjg(tmp2)))
        k1 = -CMPLX(0,1)*(matmul(HS,tmp1)-matmul(tmp1,HS)) - (matmul(VI,tmp3)-matmul(tmp3,VI))

        tc = ti + dt / 2
        tmp1 = rho(i,:,:) + k1 * dt / 2
        tmp2 = 0.5 * (bath_VI(i,:,:) + bath_VI(i+1,:,:)) ! Not happy about this
        tmp3 = matmul(tmp2, tmp1) - matmul(tmp1, transpose(conjg(tmp2)))
        k2 = -CMPLX(0,1)*(matmul(HS,tmp1)-matmul(tmp1,HS)) - (matmul(VI,tmp3)-matmul(tmp3,VI))

        tc = ti + dt / 2
        tmp1 = rho(i,:,:) + k2 * dt / 2
        tmp2 = 0.5 * (bath_VI(i,:,:) + bath_VI(i+1,:,:)) ! Not happy about this
        tmp3 = matmul(tmp2, tmp1) - matmul(tmp1, transpose(conjg(tmp2)))
        k3 = -CMPLX(0,1)*(matmul(HS,tmp1)-matmul(tmp1,HS)) - (matmul(VI,tmp3)-matmul(tmp3,VI))

        tc = ti + dt
        tmp1 = rho(i,:,:) + k3 * dt
        tmp2 = bath_VI(i+1,:,:)
        tmp3 = matmul(tmp2, tmp1) - matmul(tmp1, transpose(conjg(tmp2)))
        k4 = -CMPLX(0,1)*(matmul(HS,tmp1)-matmul(tmp1,HS)) - (matmul(VI,tmp3)-matmul(tmp3,VI))

        rho(i+1,:,:) = rho(i,:,:) + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
        write(10,'(f10.3,8(e15.6))') ti+dt,((rho(i+1,j,k),j=1,ss),k=1,ss)
    end do

    close(10)

    open(10, file='hermitian.dat')
    do i = 0, n_steps
        ti = i * dt
        tmpl_1 = test_hermitian(rho(i,:,:))
        write(10, '(f10.3,L2)') ti, tmpl_1
    end do
    close(10)

    open(10, file='trace.dat')
    do i = 0, n_steps
        ti = i * dt
        tmpl_1 = test_trace(rho(i,:,:))
        tmp_val1 = REAL(trace(rho(i,:,:)), kind=DP)
        write(10, '(f10.3,L2,e15.6)') ti, tmpl_1, tmp_val1
    end do
    close(10)

    open(10, file='positivity.dat')
    do i = 0, n_steps
        ti = i * dt
        tmpl_1 = test_positivity(rho(i,:,:))
        write(10, '(f10.3,L2)') ti, tmpl_1
    end do
    close(10)

    open(10, file='entropy.dat')
    do i = 0, n_steps
        ti = i * dt
        tmp_val1 = REAL(Entropy(rho(i,:,:)), kind=DP)
        write(10, '(f10.3,e15.6)') ti, tmp_val1
    end do
    close(10)

    open(10, file='stats_E.dat')
    do i = 0, n_steps
        ti = i * dt
        tmp_val1 = REAL(trace(matmul(HS, rho(i,:,:))), kind=DP)
        tmp_val2 = REAL(trace(matmul(HS, matmul(HS, rho(i,:,:)))), kind=DP)
        tmp_val3 = REAL(SQRT(tmp_val2 - tmp_val1**2), kind=DP)
        write(10, '(f10.3,2(e15.6))') ti, tmp_val1, tmp_val3
    end do
    close(10)

    open(10, file='stats_XS.dat')
    do i = 0, n_steps
        ti = i * dt
        tmp_val1 = REAL(trace(matmul(VI, rho(i,:,:))), kind=DP)
        tmp_val2 = REAL(trace(matmul(VI, matmul(VI, rho(i,:,:)))), kind=DP)
        tmp_val3 = REAL(SQRT(tmp_val2 - tmp_val1**2), kind=DP)
        write(10, '(f10.3,2(e15.6))') ti, tmp_val1, tmp_val3
    end do
    close(10)

    close(20)
end program main
