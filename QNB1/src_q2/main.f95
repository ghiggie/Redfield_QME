program main

   use constants
   use variables
   use helper
   use bath

   implicit none

   character(len=40) :: filename, arg, hostname
   complex(kind=DP), dimension(:,:,:), allocatable :: rho, bath_VI, bath_VI_half, eta, rho_E
   integer :: n_arg, i, j, k, iargc
   logical :: tmp_l1, halt, found_herm, found_trace, found_pos
   logical :: gibbs, random, normalize
   real(kind=DP) :: ti, tc, tmp_r1, tmp_r2, tmp_r3, mod
   complex(kind=DP) :: Z
   real :: cputime0, cputime1, cputime2
   character(len=10) :: tmp_str1, tmp_str2, form_str1, form_str2
   complex(kind=DP) :: tmp_c1, tmp_c2
   complex(kind=DP), dimension(:), allocatable :: sys_entropy, energy, heat, sprod
   complex(kind=DP), dimension(:), allocatable :: traced, fid
   complex(kind=DP), dimension(:,:), allocatable :: gibbsrho, eigvect_inv
   complex(kind=DP), dimension(2,2) :: rho_A, rho_B, rho_AE, rho_BE, eigvect_inv2, rho_gibbs2
   complex(kind=DP), dimension(2,2) :: HS_2, eigvect_2
   real(kind=DP), dimension(2) :: eigval_2
   complex(kind=DP), dimension(4,4) :: eigvect_inv4

   namelist/params/dt,time_limit,time_write,pade,matsu,temp,gamma,lambda,rho0,HS,VI

   halt = .false.
   gibbs = .false.
   random = .false.
   normalize = .false.

   n_arg = iargc()
   if (n_arg .gt. 0) then
      do i = 1, n_arg
         call getarg(i, arg)
         select case(trim(arg))
         case('-v')
            write(*,*) 'Coded by ', trim(coded_by)
            write(*,*) 'Last Updated: ', trim(last_update)
            write(*,*) 'Redfield with one bath: ', trim(version)
            STOP ''
         case('-h')
            write(*,*) ' Usage: ./Redfield1B [-h] [-v] [-G|-R] [-N] config_file'
            write(*,*) ' Required Parameters:'
            write(*,*) '    config_file: parameter values in namelist format'
            write(*,*) ' Optional Parameters:'
            write(*,*) '    -h: show this usage information'
            write(*,*) '    -v: show the version number'
            write(*,*) '    -G: use an initial Gibbs density'
            write(*,*) '    -R: use an initial random density'
            write(*,*) '    -N: automatically normalize the input density'
            STOP ''
         case('-G')
            gibbs = .true.
         case('-R')
            random = .true.
         case('-N')
            normalize = .true.
         case default
            filename = trim(arg)
         end select
      end do
   else
      write(*,*) ' Configuration file is required.'
      write(*,*) ' Usage: ./Redfield1B [-h] [-v] [-G|-R] [-N] config_file'
      write(*,*) ' Type ./Redfield1B -h for full usage information.'
      STOP ''
   end if

   if (gibbs .and. random) then
      write(*,*) 'You cannot use -G and -R simultaneously.'
      write(*,*) ' Usage: ./Redfield1B [-h] [-v] [-G|-R] [-N] config_file'
      write(*,*) ' Type ./Redfield1B -h for full usage information.'
      STOP ''
   end if

   ! Set the system size to be ss=4
   ss = 4

   ! Initialize the single time-step arrays
   call array_init()

   open(10, file=trim(filename))
   read(10,nml=params)
   close(10)

   n_steps = nint(time_limit / dt)

   if (gibbs) then ! Generate Gibbs state
      rho0 = 0
      call eigensystem(HS, eigval, eigvect)
      do k=1, ss
         do i=1, ss
            do j=1, ss
               rho0(i,j) = rho0(i,j) + eigvect(i,k)*conjg(eigvect(j,k))*exp(-eigval(k)/temp)
            end do
         end do
      end do
      Z = trace(rho0)
      rho0 = rho0 / Z
   else if (random) then ! Generate random state
      do i=1,ss
         do j=1,ss
            call random_number(tmp_r1)
            call random_number(tmp_r2)
            rho0(i,j) = cmplx(tmp_r1, tmp_r2)
         end do
      end do
      rho0 = matmul(rho0,transpose(conjg(rho0)))
      do i=1,ss
         rho0(i,i) = cmplx(real(rho0(i,i),kind=DP),0.0_DP)
      end do
      Z = trace(rho0)
      rho0 = rho0/Z
   end if

   if (normalize) then
      Z = trace(rho0)
      rho0 = rho0 / Z
   end if

   ! Initialize the multi time-step arrays
   allocate(rho(0:n_steps,ss,ss))
   rho(0,:,:) = rho0 ! Initialize the storage
   allocate(rho_E(0:n_steps,ss,ss))
   allocate(bath_VI(0:n_steps,ss,ss))
   allocate(bath_VI_half(n_steps,ss,ss))
   allocate(sys_entropy(0:n_steps))
   allocate(eta(0:n_steps,ss,ss))
   allocate(energy(0:n_steps))
   allocate(heat(0:n_steps))
   allocate(sprod(0:n_steps))
   allocate(traced(0:n_steps))
   allocate(fid(0:n_steps))
   allocate(gibbsrho(ss,ss))
   allocate(eigvect_inv(ss,ss))

   ! Calculate gibbs state
   gibbsrho = 0
   call eigensystem(HS, eigval, eigvect)
   do k=1, ss
      do i=1, ss
         do j=1, ss
            gibbsrho(i,j) = gibbsrho(i,j) + eigvect(i,k)*conjg(eigvect(j,k))*exp(-eigval(k)/temp)
         end do
      end do
   end do
   Z = trace(gibbsrho)
   gibbsrho = gibbsrho / Z

   ! Create the arrays needed for the bath correlation calculation
   call bc_coeff()

   ! Set up the summary file
   write(form_str1, '(I4)') ss
   write(form_str2, '(I4)') 2*ss**2

   open(20, file = 'Redfield1B.out')
   write(20, '(3a)') '*** Redfield with One Bath (Redfield1B) ', trim(version), ' ***'

   write(20, '(/a)') 'Unless otherwise indicated, all matrices are in the atomic basis.'

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
   write(20, '(/a)') 'The Gibbs state for this Hamiltonian is'
   do i=1,ss
      write(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(gibbsrho(i,j)),',',AIMAG(gibbsrho(i,j)),'),',j=1,ss)
   end do

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
   if (gibbs) then
      write(20,'(/a)') 'A Gibbs state is used for the initial density.'
   else if (random) then
      write(20,'(/a)') 'A random state is used for the initial density.'
   end if
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
      write(tmp_str1, '(I3)') matsu
      write(20,'(2a/)') trim(adjustl(tmp_str1)), ' matsubara frequencies are used.'
      write(tmp_str1, '(E10.3)') REAL(cinf)
      write(20,'(2a)') 'cinf    = ', trim(adjustl(tmp_str1))
      write(20,'(a,E10.3,a,E10.3,a)') 'c0      = (',REAL(coeff(0)),',',AIMAG(coeff(0)),')'
      write(tmp_str1, '(f10.7)') abs(exp_vec(0))
      write(20, '(2a)') 'gamma0  = ', trim(adjustl(tmp_str1))
      do i = 1, matsu
         write(tmp_str1, '(I3)') i
         write(tmp_str2, '(E10.3)') REAL(coeff(i))
         write(20,'(4a)') 'c',trim(adjustl(tmp_str1)),'      = ', trim(adjustl(tmp_str2))
         write(tmp_str2, '(f10.7)') abs(exp_vec(i))
         write(20, '(4a)') 'gamma',trim(adjustl(tmp_str1)),'  = ', trim(adjustl(tmp_str2))
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
      tmp_arr2 = bath_VI(n_steps,:,:)
      tmp_arr3 = matmul(tmp_arr2, tmp_arr1) - matmul(tmp_arr1, transpose(conjg(tmp_arr2)))
      k1 = -CMPLX(0,1)*(matmul(HS,tmp_arr1)-matmul(tmp_arr1,HS)) - (matmul(VI,tmp_arr3)-matmul(tmp_arr3,VI))

      tc = ti + dt / 2
      tmp_arr1 = rho(i,:,:) + k1 * dt / 2
      tmp_arr2 = bath_VI(n_steps,:,:)
      tmp_arr3 = matmul(tmp_arr2, tmp_arr1) - matmul(tmp_arr1, transpose(conjg(tmp_arr2)))
      k2 = -CMPLX(0,1)*(matmul(HS,tmp_arr1)-matmul(tmp_arr1,HS)) - (matmul(VI,tmp_arr3)-matmul(tmp_arr3,VI))

      tc = ti + dt / 2
      tmp_arr1 = rho(i,:,:) + k2 * dt / 2
      tmp_arr2 = bath_VI(n_steps,:,:)
      tmp_arr3 = matmul(tmp_arr2, tmp_arr1) - matmul(tmp_arr1, transpose(conjg(tmp_arr2)))
      k3 = -CMPLX(0,1)*(matmul(HS,tmp_arr1)-matmul(tmp_arr1,HS)) - (matmul(VI,tmp_arr3)-matmul(tmp_arr3,VI))

      tc = ti + dt
      tmp_arr1 = rho(i,:,:) + k3 * dt
      tmp_arr2 = bath_VI(n_steps,:,:)
      tmp_arr3 = matmul(tmp_arr2, tmp_arr1) - matmul(tmp_arr1, transpose(conjg(tmp_arr2)))
      k4 = -CMPLX(0,1)*(matmul(HS,tmp_arr1)-matmul(tmp_arr1,HS)) - (matmul(VI,tmp_arr3)-matmul(tmp_arr3,VI))

      rho(i+1,:,:) = rho(i,:,:) + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
   end do

   ! Calculate entropy
   do i = 0, n_steps
      sys_entropy(i) = Entropy(rho(i,:,:))
   end do

   ! Calculate eta
   do i = 0, n_steps
      eta(i,:,:) = CMPLX(0,1)*(matmul(rho(i,:,:),transpose(conjg(bath_VI(i,:,:))))-matmul(bath_VI(i,:,:),rho(i,:,:)))
   end do

   ! Calculate system energy
   do i = 0, n_steps
      energy(i) = trace(matmul(HS, rho(i,:,:)))
   end do

   ! Calculate heat
   do i = 0, n_steps
      tmp_c1 = energy(i) - energy(0)
      tmp_c2 = trace(matmul(VI,eta(i,:,:))) - trace(matmul(VI,eta(0,:,:)))
      heat(i) = tmp_c1 + tmp_c2
   end do

   ! Calculate entropy production
   do i = 0, n_steps
      sprod(i) = (sys_entropy(i)-sys_entropy(0)) - heat(i) / temp
   end do

   ! Calculate fidelity
   do i = 0, n_steps
      fid(i) = fidelity(rho(i,:,:),gibbsrho)
   end do

   ! Calculate trace distance
   do i = 0, n_steps
      traced(i) = trace_distance(rho(i,:,:),gibbsrho)
   end do

   ! Calculate system density in the energy eigenbasis

   call eigensystem(HS, eigval, eigvect)
   call matinv4(eigvect, eigvect_inv4)
   do i = 0, n_steps
      rho_E(i,:,:) = matmul(eigvect_inv4, matmul(rho(i,:,:), eigvect))
   end do

   call CPU_TIME(cputime1)

   open(10,file='rho.dat')
   do i = 0, n_steps
      ti = i * dt
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10,'(f10.3,('//form_str2//'e15.6))') ti,((rho(i,j,k),j=1,ss),k=1,ss)
   end do
   close(10)

   open(10,file='rho_E.dat')
   do i = 0, n_steps
      ti = i * dt
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10,'(f10.3,('//form_str2//'e15.6))') ti,((rho_E(i,j,k),j=1,ss),k=1,ss)
   end do
   close(10)

   open(10,file='lambda.dat')
   do i = 0, n_steps
      ti = i * dt
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10,'(f10.3,('//form_str2//'e15.6))') ti,((bath_VI(i,j,k),j=1,ss),k=1,ss)
   end do
   close(10)

   open(10, file='eta.dat')
   do i = 0, n_steps
      ti = i * dt
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10,'(f10.3,('//form_str2//'e15.6))') ti,((eta(i,j,k),j=1,ss),k=1,ss)
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
      tmp_r1 = REAL(sys_entropy(i), kind=DP)
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10, '(f10.3,e15.6)') ti, tmp_r1
   end do
   close(10)

   open(10, file='stats_E.dat')
   do i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(energy(i), kind=DP)
      tmp_r2 = REAL(trace(matmul(HS, matmul(HS, rho(i,:,:)))), kind=DP)
      tmp_r3 = REAL(SQRT(tmp_r2 - tmp_r1**2), kind=DP)
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10, '(f10.3,2(e15.6))') ti, tmp_r1, tmp_r3
   end do
   close(10)

   open(10, file='stats_XS.dat')
   do i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(trace(matmul(VI, rho(i,:,:))), kind=DP)
      tmp_r2 = REAL(trace(matmul(VI, matmul(VI, rho(i,:,:)))), kind=DP)
      tmp_r3 = REAL(SQRT(tmp_r2 - tmp_r1**2), kind=DP)
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10, '(f10.3,2(e15.6))') ti, tmp_r1, tmp_r3
   end do
   close(10)

   open(10, file='heat.dat')
   do i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(heat(i), kind=DP)
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10, '(f10.3,e15.6)') ti, tmp_r1
   end do
   close(10)

   open(10, file='sprod.dat')
   do i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(sprod(i), kind=DP)
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10, '(f10.3,e15.6)') ti, tmp_r1
   end do
   close(10)

   open(10, file='fidelity.dat')
   do i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(fid(i), kind=DP)
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10, '(f10.3,e15.6)') ti, tmp_r1
   end do
   close(10)

   open(10, file='trace_distance.dat')
   do i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(traced(i), kind=DP)
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10, '(f10.3,e15.6)') ti, tmp_r1
   end do
   close(10)

   open(10, file='rhoA.dat')
   do i = 0, n_steps
      ti = i * dt
      call rhoA(rho(i,:,:), rho_A)
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10,'(f10.3,8e15.6)') ti, ((rho_A(k,j),k=1,2),j=1,2)
   end do
   close(10)

   open(10, file='rhoB.dat')
   do i = 0, n_steps
      ti = i * dt
      call rhoB(rho(i,:,:), rho_B)
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10,'(f10.3,8e15.6)') ti, ((rho_B(k,j),k=1,2),j=1,2)
   end do
   close(10)

   ! First, calculate Gibbs state for the single qubit Hamiltonian

   HS_2(1,:) = (/CMPLX(0.5,0.), CMPLX(0.,0.)/)
   HS_2(2,:) = (/CMPLX(0.,0.), CMPLX(-0.5,0.)/)

   call eigensystem(HS_2, eigval_2, eigvect_2)
   call matinv2(eigvect_2, eigvect_inv2)
   do k=1, 2
      do i=1, 2
         do j=1, 2
            rho_gibbs2(i,j) = rho_gibbs2(i,j) + eigvect_2(i,k)*conjg(eigvect_2(j,k))*exp(-eigval_2(k)/temp)
         end do
      end do
   end do
   Z = trace(rho_gibbs2)
   rho_gibbs2 = rho_gibbs2 / Z

   open(10, file='rhoA_E.dat')
   do i = 0, n_steps
      ti = i * dt
      call rhoA(rho(i,:,:), rho_A)
      ! Now convert rho_A to energy eigenbasis and store in rho_AE
      rho_AE = matmul(eigvect_inv2, matmul(rho_A, eigvect_2))
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10,'(f10.3,8e15.6)') ti, ((rho_AE(k,j),k=1,2),j=1,2)
   end do
   close(10)

   open(10, file='rhoB_E.dat')
   do i = 0, n_steps
      ti = i * dt
      call rhoB(rho(i,:,:), rho_B)
      rho_BE = matmul(eigvect_inv2, matmul(rho_B, eigvect_2))
      mod = ti - nint(ti/time_write)*time_write
      if (mod .eq. 0) write(10,'(f10.3,8e15.6)') ti, ((rho_BE(k,j),k=1,2),j=1,2)
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

   write(20,'(/a)') 'The final density in the energy eigenbasis is'
   tmp_arr1 = rho_E(n_steps,:,:)
   do i=1,ss
      write(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(tmp_arr1(i,j)),',',AIMAG(tmp_arr1(i,j)),'),',j=1,ss)
   end do

   write(20,'(/a)') 'The reduced density for system A is'
   write(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_A(1,j)),',',AIMAG(rho_A(1,j)),'),',j=1,2)
   write(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_A(2,j)),',',AIMAG(rho_A(2,j)),'),',j=1,2)
   write(20,'(/a)') 'The reduced density for system B is'
   write(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_B(1,j)),',',AIMAG(rho_B(1,j)),'),',j=1,2)
   write(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_B(2,j)),',',AIMAG(rho_B(2,j)),'),',j=1,2)
   write(20,'(/a)') 'The reduced density for system A in the energy eigenbasis is'
   write(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_AE(1,j)),',',AIMAG(rho_AE(1,j)),'),',j=1,2)
   write(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_AE(2,j)),',',AIMAG(rho_AE(2,j)),'),',j=1,2)
   write(20,'(/a)') 'The reduced density for system B in the energy eigenbasis is'
   write(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_BE(1,j)),',',AIMAG(rho_BE(1,j)),'),',j=1,2)
   write(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_BE(2,j)),',',AIMAG(rho_BE(2,j)),'),',j=1,2)
   write(20,'(/a)') 'The Gibbs state for the individual qubits is'
   write(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_gibbs2(1,j)),',',AIMAG(rho_gibbs2(1,j)),'),',j=1,2)
   write(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_gibbs2(2,j)),',',AIMAG(rho_gibbs2(2,j)),'),',j=1,2)

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
