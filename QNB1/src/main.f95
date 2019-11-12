program main

    use constants
    use variables
    use helper
    use bath

    implicit none

    character(len=40) :: filename, arg, hostname
    complex(kind=DP), dimension(:,:,:), allocatable :: rho, bath_VI, bath_VI_half
    integer :: n_arg, i, j, k, iargc
    logical :: tmp_l1, halt, found_herm, found_trace, found_pos
    real(kind=DP) :: ti, tc, tmp_r1, tmp_r2, tmp_r3
    real :: cputime0, cputime1, cputime2
    character(len=4) :: tmp_str, form_str1, form_str2

    namelist/params/dt,time_limit,pade,matsu,temp,gamma,lambda,rho0,HS,VI

    halt = .false.

    n_arg = iargc()
    if (n_arg .eq. 0) then
        write(*,*) ' System size and configuration file are required.'
        write(*,*) ' Usage: ./Redfield1B [-h] [-v] system_size config_file'
    else
        do i = 1, n_arg
            call getarg(i, arg)
            select case(trim(arg))
            case('-v')
                write(*,*) 'Coded by ', trim(coded_by)
                write(*,*) 'Last Updated: ', trim(last_update)
                write(*,*) 'Redfield with one bath: ', trim(version)
                STOP ''
            case('-h')
                write(*,*) ' Usage: ./Redfield1B [-h] [-v] system_size config_file'
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
                    write(*,*) ' Usage: ./Redfield1B [-h] [-v] system_size config_file'
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
    allocate(rho(0:n_steps,ss,ss))
    rho(0,:,:) = rho0 ! Initialize the storage
    allocate(bath_VI(0:n_steps,ss,ss))
    allocate(bath_VI_half(n_steps,ss,ss))
    ! Create the arrays needed for the bath correlation calculation
    call bc_coeff()

    ! Set up the summary file
    write(form_str1, '(I4)') ss
    write(form_str2, '(I4)') 2*ss**2

    open(20, file = 'Redfield1B.out')
    write(20, '(3a)') '*** Redfield with One Bath (Redfield1B) ', trim(version), ' ***'

    write(20, '(/a,a)') 'Job begun at ', time_stamp()
    call get_environment_variable('HOSTNAME', hostname)
    write(20, '(/a,a)') 'Host: ', trim(hostname)

    write(20,'(/a)') 'System Hamiltonian'
    write(20,'(a/)') '=================='
    write(20,'(a)') 'HS = '
    do i=1,ss
        write(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(HS(i,j)),',',AIMAG(HS(i,j)),'),',j=1,ss)
    end do
    tmp_l1 = test_hermitian(HS)
    if (tmp_l1) then
        write(20,'(/a)') 'The System Hamiltonian is Self-Adjoint.'
    else
        write(20,'(/a)') '**** The System Hamiltonian is NOT Self-Adjoint. ****'
        halt = .true.
    end if
    if (tmp_l1) then
        call eigensystem(HS,eigval,eigvect)
        write(20,'(/a)') 'Energy Eigenvalues and Eigenvectors (in column)'
        write(20,'('//form_str1//'(7x,f10.7,7x))') eigval
        do i=1,ss
            write(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',real(eigvect(i,j)),',',aimag(eigvect(i,j)),'),',j=1,ss)
        end do
    end if

    write(20, '(/a)') 'Coupling Operator to Environment'
    write(20, '(a/)') '================================'
    write(20,'(a)') 'VI = '
    do i=1,ss
        write(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(VI(i,j)),',',AIMAG(VI(i,j)),'),',j=1,ss)
    end do
    tmp_l1 = test_hermitian(VI)
    if (tmp_l1) then
        write(20,'(/a)') 'The Coupling Operator is Self-Adjoint.'
    else
        write(20,'(/a)') '**** The Couplng Operator is NOT Self-Adjoint. ****'
        halt = .true.
    end if
    if (tmp_l1) then
        call eigensystem(VI,eigval,eigvect)
        write(20,'(/a)') 'Coupling Operator Eigenvalues and Eigenvectors (in column)'
        write(20,'('//form_str1//'(7x,f10.7,7x))') eigval
        do i=1,ss
            write(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',real(eigvect(i,j)),',',aimag(eigvect(i,j)),'),',j=1,ss)
        end do
    end if

    write(20, '(/a)') 'Initial Density Matrix'
    write(20, '(a/)') '======================'
    write(20,'(a)') 'rho0 = '
    do i=1,ss
        write(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(rho0(i,j)),',',AIMAG(rho0(i,j)),'),',j=1,ss)
    end do
    tmp_l1 = test_hermitian(rho(0,:,:))
    if (tmp_l1) then
        write(20,'(/a)') 'The Initial Density is Self-Adjoint.'
    else
        write(20,'(/a)') '**** The Initial Density is NOT Self-Adjoint. ****'
        halt = .true.
    end if
    tmp_l1 = test_trace(rho0)
    if (tmp_l1) then
        write(20,'(a,2f10.6)') 'The trace is OK. Trace = ', trace(rho0)
    else
        write(20,'(a)') '**** The Initial Density is NOT Normalized. ****'
        halt = .true.
    end if
    tmp_l1 = test_positivity(rho0)
    if (tmp_l1) then
        write(20,'(a)') 'The Initial Density is Non-Negative.'
    else
        write(20,'(a)') '**** The Initial Density is NOT Non-Negative. ****'
        halt = .true.
    end if
    if (tmp_l1) then
        call eigensystem(rho0,eigval,eigvect)
        write(20,'(/a)') 'Initial Density Eigenvalues and Eigenvectors (in column)'
        write(20,'('//form_str1//'(7x,f10.7,7x))') eigval
        do i=1,ss
            write(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',real(eigvect(i,j)),',',aimag(eigvect(i,j)),'),',j=1,ss)
        end do
    end if

    if (halt) then
        write(20, '(/a)') '**** Bad Input Data ****'
        write(20, '(a)') '**** Execution is Stopped ****'
        STOP ''
    end if

    write(20, '(/a)') 'Environment Parameters'
    write(20, '(a/)') '======================'
    write(20, '(a,f5.3)') '     T = ', temp
    write(20, '(a,f5.3)') ' gamma = ', gamma
    write(20, '(a,f5.3)') 'lambda = ', lambda

    write(20,'(/a)') 'Correlation functions'
    write(20,'(a/)') '====================='
    if(pade) then
       write(20,'(a/)') '[1/1] Pade approximation is used.'
       write(20,'(a,f10.7,a,f10.7,a)') 'cinf = (', REAL(cinf), ',',AIMAG(cinf),')'
       write(20,'(a,f10.7,a,f10.7,a)') 'c0 = (',REAL(coeff(0)),',',AIMAG(coeff(0)),')'
       write(20,'(a,f10.7)') 'gamma0 = ', abs(exp_vec(0))
       write(20,'(a,f10.7,a,f10.7,a)') 'c1 = (',REAL(coeff(1)),',',AIMAG(coeff(1)),')'
       write(20,'(a,f10.7)') 'gamma1 = ', abs(exp_vec(1))
    else
       write(tmp_str, '(I2)') matsu
       write(20,'(2a/)') trim(adjustl(tmp_str)), ' matsubara freaquencies are used.'
       write(20,'(a,f10.7,a,f10.7,a)') 'cinf = (', REAL(cinf), ',',AIMAG(cinf),')'
       do i = 0, matsu
           write(20,'(a,I1,a,f10.7,a,f10.7,a)') 'c',i,' = (',REAL(coeff(i)),',',AIMAG(coeff(i)),')'
           write(20, '(a,I1,a,f10.7)') 'gamma',i,' = ', abs(exp_vec(i))
       end do
    end if

    write(20, '(/a)') 'Execution Parameters'
    write(20, '(a/)') '===================='
    write(20, '(a,f8.3)') 'Time limit = ', time_limit
    write(20, '(a,f5.3)') '        dt = ', dt

    call CPU_TIME(cputime0)

    ! Create the array of values for lambda_bc, to be stored in bath_VI
    call lambda_bc(VI, bath_VI, bath_VI_half)

    do i = 0, n_steps - 1
        ti = i * dt

        tc = ti
        tmp_arr1 = rho(i,:,:)
        tmp_arr2 = bath_VI(i,:,:)
        tmp_arr3 = matmul(tmp_arr2, tmp_arr1) - matmul(tmp_arr1, transpose(conjg(tmp_arr2)))
        k1 = -CMPLX(0,1)*(matmul(HS,tmp_arr1)-matmul(tmp_arr1,HS)) - (matmul(VI,tmp_arr3)-matmul(tmp_arr3,VI))

        tc = ti + dt / 2
        tmp_arr1 = rho(i,:,:) + k1 * dt / 2
        tmp_arr2 = bath_VI_half(i+1,:,:)
        tmp_arr3 = matmul(tmp_arr2, tmp_arr1) - matmul(tmp_arr1, transpose(conjg(tmp_arr2)))
        k2 = -CMPLX(0,1)*(matmul(HS,tmp_arr1)-matmul(tmp_arr1,HS)) - (matmul(VI,tmp_arr3)-matmul(tmp_arr3,VI))

        tc = ti + dt / 2
        tmp_arr1 = rho(i,:,:) + k2 * dt / 2
        tmp_arr2 = bath_VI_half(i+1,:,:)
        tmp_arr3 = matmul(tmp_arr2, tmp_arr1) - matmul(tmp_arr1, transpose(conjg(tmp_arr2)))
        k3 = -CMPLX(0,1)*(matmul(HS,tmp_arr1)-matmul(tmp_arr1,HS)) - (matmul(VI,tmp_arr3)-matmul(tmp_arr3,VI))

        tc = ti + dt
        tmp_arr1 = rho(i,:,:) + k3 * dt
        tmp_arr2 = bath_VI(i+1,:,:)
        tmp_arr3 = matmul(tmp_arr2, tmp_arr1) - matmul(tmp_arr1, transpose(conjg(tmp_arr2)))
        k4 = -CMPLX(0,1)*(matmul(HS,tmp_arr1)-matmul(tmp_arr1,HS)) - (matmul(VI,tmp_arr3)-matmul(tmp_arr3,VI))

        rho(i+1,:,:) = rho(i,:,:) + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
    end do

    call CPU_TIME(cputime1)

    open(10,file='rho.dat')
    do i = 0, n_steps
        ti = i * dt
        if (mod(ti, 1.0) .eq. 0) write(10,'(f10.3,('//form_str2//'e15.6))') ti,((rho(i,j,k),j=1,ss),k=1,ss)
    end do
    close(10)

    open(10, file='hermitian.dat')
    found_herm = .false.
    do i = 0, n_steps
        tmp_l1 = test_hermitian(rho(i,:,:))
        if (.not. tmp_l1) then
            ti = i * dt
            write(10, '(f10.3,L2)') ti, tmp_l1
            found_herm = .true.
        end if
    end do
    if (.not. found_herm) write(10, '(a)') 'The density remains hermitian at all times.'
    close(10)

    open(10, file='trace.dat')
    found_trace = .false.
    do i = 0, n_steps
        tmp_l1 = test_trace(rho(i,:,:))
        if (.not. tmp_l1) then
            ti = i * dt
            write(10, '(f10.3,L2)') ti, tmp_l1
            found_trace = .true.
        end if
    end do
    if (.not. found_trace) write(10, '(a)') 'The density remains normalized at all times.'
    close(10)

    open(10, file='positivity.dat')
    found_pos = .false.
    do i = 0, n_steps
        tmp_l1 = test_positivity(rho(i,:,:))
        if (.not. tmp_l1) then
            ti = i * dt
            write(10, '(f10.3,L2)') ti, tmp_l1
            found_pos = .true.
        end if
    end do
    if (.not. found_pos) write(10, '(a)') 'The density remains positive at all times.'
    close(10)

    open(10, file='entropy.dat')
    do i = 0, n_steps
        ti = i * dt
        tmp_r1 = REAL(Entropy(rho(i,:,:)), kind=DP)
        if (mod(ti, 1.0) .eq. 0) write(10, '(f10.3,e15.6)') ti, tmp_r1
    end do
    close(10)

    open(10, file='stats_E.dat')
    do i = 0, n_steps
        ti = i * dt
        tmp_r1 = REAL(trace(matmul(HS, rho(i,:,:))), kind=DP)
        tmp_r2 = REAL(trace(matmul(HS, matmul(HS, rho(i,:,:)))), kind=DP)
        tmp_r3 = REAL(SQRT(tmp_r2 - tmp_r1**2), kind=DP)
        if (mod(ti, 1.0) .eq. 0) write(10, '(f10.3,2(e15.6))') ti, tmp_r1, tmp_r3
    end do
    close(10)

    open(10, file='stats_XS.dat')
    do i = 0, n_steps
        ti = i * dt
        tmp_r1 = REAL(trace(matmul(VI, rho(i,:,:))), kind=DP)
        tmp_r2 = REAL(trace(matmul(VI, matmul(VI, rho(i,:,:)))), kind=DP)
        tmp_r3 = REAL(SQRT(tmp_r2 - tmp_r1**2), kind=DP)
        if (mod(ti, 1.0) .eq. 0) write(10, '(f10.3,2(e15.6))') ti, tmp_r1, tmp_r3
    end do
    close(10)

    call CPU_TIME(cputime2)

    write(20, '(/a)') 'Final Density Matrix'
    write(20, '(a/)') '======================'
    write(20,'(a)') 'rho(final) = '
    tmp_arr1 = rho(n_steps,:,:)
    do i=1,ss
        write(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(tmp_arr1(i,j)),',',AIMAG(tmp_arr1(i,j)),'),',j=1,ss)
    end do
    tmp_l1 = test_hermitian(tmp_arr1)
    if (tmp_l1) then
        write(20,'(/a)') 'The Final Density is Self-Adjoint.'
    else
        write(20,'(/a)') '**** The Final Density is NOT Self-Adjoint. ****'
    end if
    tmp_l1 = test_trace(tmp_arr1)
    if (tmp_l1) then
        write(20,'(a,2f10.6)') 'The trace is OK. Trace = ', trace(tmp_arr1)
    else
        write(20,'(a)') '**** The Final Density is NOT Normalized. ****'
    end if
    tmp_l1 = test_positivity(tmp_arr1)
    if (tmp_l1) then
        write(20,'(a)') 'The Final Density is Non-Negative.'
    else
        write(20,'(a)') '**** The Final Density is NOT Non-Negative. ****'
    end if
    if (tmp_l1) then
        call eigensystem(tmp_arr1,eigval,eigvect)
        write(20,'(/a)') 'Final Density Eigenvalues and Eigenvectors (in column)'
        write(20,'('//form_str1//'(7x,f10.7,7x))') eigval
        do i=1,ss
            write(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',real(eigvect(i,j)),',',aimag(eigvect(i,j)),'),',j=1,ss)
        end do
    end if

    write(20, '(/a)') 'Other Information'
    write(20, '(a/)') '================='

    if (found_herm) then
        write(20, '(a)') 'Hermiticity is not preserved at all times. See hermitian.dat'
    else
        write(20, '(a)') 'Hermiticity is preserved at all times.'
    end if

    if (found_trace) then
        write(20, '(a)') 'Trace is not preserved at all times. See trace.dat'
    else
        write(20, '(a)') 'Trace is preserved at all times.'
    end if

    if (found_pos) then
        write(20, '(a)') 'Positivity is not preserved at all times. See positivity.dat'
    else
        write(20, '(a)') 'Positivity is preserved at all times.'
    end if

    write(20,'(/a,f10.5)') 'Time needed for time evolution: ', cputime1 - cputime0
    write(20,'(a,f10.5)') 'Time needed for data output: ', cputime2 - cputime1
    write(20,'(/a,a)') 'Job finished at ', time_stamp()

    close(20)

end program main
