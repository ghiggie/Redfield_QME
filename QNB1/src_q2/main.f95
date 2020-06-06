PROGRAM main

   USE constants
   USE variables
   USE helper
   USE bath

   IMPLICIT NONE

   CHARACTER(LEN=40)                               :: filename, arg, hostname
   COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: rho, bath_VI, bath_VI_half, eta, rho_E
   INTEGER                                         :: n_arg, i, j, k, iargc
   LOGICAL                                         :: tmp_l1, halt, found_herm, found_trace, found_pos
   LOGICAL                                         :: gibbs, random, normalize
   REAL   (KIND=DP)                                :: ti, tc, tmp_r1, tmp_r2, tmp_r3, mod
   REAL                                            :: cputime0, cputime1, cputime2
   CHARACTER(LEN=10)                               :: tmp_str1, tmp_str2, form_str1, form_str2
   COMPLEX(KIND=DP)                                :: Z, tmp_c1, tmp_c2
   COMPLEX(KIND=DP), DIMENSION(:),   ALLOCATABLE   :: sys_entropy, energy, heat, sprod
   COMPLEX(KIND=DP), DIMENSION(:),   ALLOCATABLE   :: traced, fid
   COMPLEX(KIND=DP), DIMENSION(:,:), ALLOCATABLE   :: rho_g, eigvect_inv

   NAMELIST/params/dt,time_limit,time_write,pade,matsu,temp,gamma,lambda,rho0,HS,VI

   halt      = .FALSE.
   gibbs     = .FALSE.
   random    = .FALSE.
   normalize = .FALSE.

   n_arg = IARGC()
   IF (n_arg .GT. 0) THEN
      DO i = 1, n_arg
         CALL GETARG(i, arg)
         SELECT CASE(TRIM(arg))
         CASE('-v')
            WRITE(*,*) 'Coded by ',                TRIM(coded_by)
            WRITE(*,*) 'Last Updated: ',           TRIM(last_update)
            WRITE(*,*) 'Redfield with one bath: ', TRIM(version)
            STOP ''
         CASE('-h')
            WRITE(*,*) ' Usage: ./Redfield1B [-h] [-v] [-G|-R] [-N] config_file'
            WRITE(*,*) ' Required Parameters:'
            WRITE(*,*) '    config_file: parameter values in namelist format'
            WRITE(*,*) ' Optional Parameters:'
            WRITE(*,*) '    -h: show this usage information'
            WRITE(*,*) '    -v: show the version number'
            WRITE(*,*) '    -G: use an initial Gibbs density'
            WRITE(*,*) '    -R: use an initial random density'
            WRITE(*,*) '    -N: automatically normalize the input density'
            STOP ''
         CASE('-G')
            gibbs     = .TRUE.
         CASE('-R')
            random    = .TRUE.
         CASE('-N')
            normalize = .TRUE.
         CASE DEFAULT
            filename  = TRIM(arg)
         END SELECT
      END DO
   ELSE
      WRITE(*,*) ' Configuration file is required.'
      WRITE(*,*) ' Usage: ./Redfield1B [-h] [-v] [-G|-R] [-N] config_file'
      WRITE(*,*) ' Type ./Redfield1B -h for full usage information.'
      STOP ''
   END IF

   IF (gibbs .AND. random) THEN
      WRITE(*,*) 'You cannot use -G and -R simultaneously.'
      WRITE(*,*) ' Usage: ./Redfield1B [-h] [-v] [-G|-R] [-N] config_file'
      WRITE(*,*) ' Type ./Redfield1B -h for full usage information.'
      STOP ''
   END IF

   ! Set the system size to be 4
   ss = 4

   ! Initialize the single time-step arrays
   CALL array_init()

   OPEN(10, FILE=TRIM(filename))
   READ(10, NML=params)
   CLOSE(10)

   n_steps = NINT(time_limit / dt)

   ! Calculate gibbs state
   ALLOCATE(rho_g(ss,ss))
   rho_g = 0
   CALL eigensystem(HS, eigval, eigvect)
   DO k = 1, ss
      DO i = 1, ss
         DO j = 1, ss
            rho_g(i,j) = rho_g(i,j) + eigvect(i,k)*CONJG(eigvect(j,k))*EXP(-eigval(k)/temp)
         END DO
      END DO
   END DO
   Z = trace(rho_g)
   rho_g = rho_g / Z

   IF (gibbs) THEN ! Generate Gibbs state
      rho0 = rho_g
   ELSE IF (random) THEN ! Generate random state
      DO i = 1, ss
         DO j = 1, ss
            CALL RANDOM_NUMBER(tmp_r1)
            CALL RANDOM_NUMBER(tmp_r2)
            rho0(i,j) = CMPLX(tmp_r1, tmp_r2)
         END DO
      END DO
      rho0 = MATMUL(rho0,TRANSPOSE(CONJG(rho0)))
      DO i = 1, ss
         rho0(i, i) = CMPLX(REAL(rho0(i,i),KIND=DP),0.0_DP)
      END DO
      Z    = trace(rho0)
      rho0 = rho0/Z
   END IF

   IF (normalize) THEN
      Z    = trace(rho0)
      rho0 = rho0 / Z
   END IF

   ! Initialize the multi time-step arrays
   ALLOCATE(rho(0:n_steps,ss,ss))
   rho(0,:,:) = rho0 ! Initialize the storage
   ALLOCATE(rho_E(0:n_steps,ss,ss))
   ALLOCATE(bath_VI(0:n_steps,ss,ss))
   ALLOCATE(bath_VI_half(n_steps,ss,ss))
   ALLOCATE(sys_entropy(0:n_steps))
   ALLOCATE(eta(0:n_steps,ss,ss))
   ALLOCATE(energy(0:n_steps))
   ALLOCATE(heat(0:n_steps))
   ALLOCATE(sprod(0:n_steps))
   ALLOCATE(traced(0:n_steps))
   ALLOCATE(fid(0:n_steps))
   ALLOCATE(eigvect_inv(ss,ss))

   ! Create the arrays needed for the bath correlation calculation
   CALL bc_coeff()

   ! Set up the summary file
   WRITE(form_str1, '(I4)') ss
   WRITE(form_str2, '(I4)') 2*ss**2

   OPEN(20, FILE = 'Redfield1B.out')
   WRITE(20, '(3a)') '*** Redfield with One Bath (Redfield1B) ', TRIM(version), ' ***'

   WRITE(20, '(/a)') 'Unless otherwise indicated, all matrices are in the atomic basis.'

   WRITE(20, '(/a,a)') 'Job begun at ', time_stamp()
   CALL GET_ENVIRONMENT_VARIABLE('HOSTNAME', hostname)
   WRITE(20, '(/a,a)') 'Host: ', TRIM(hostname)

   WRITE(20,'(/a)') 'System Hamiltonian'
   WRITE(20,'(a/)') '=================='
   WRITE(20,'(a)') 'HS = '
   DO i=1,ss
      WRITE(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(HS(i,j)),',',AIMAG(HS(i,j)),'),',j=1,ss)
   END DO
   tmp_l1 = test_hermitian(HS)
   IF (tmp_l1) THEN
      WRITE(20,'(/a)') 'The System Hamiltonian is Self-Adjoint.'
   ELSE
      WRITE(20,'(/a)') '**** The System Hamiltonian is NOT Self-Adjoint. ****'
      halt = .TRUE.
   END IF
   IF (tmp_l1) THEN
      CALL eigensystem(HS,eigval,eigvect)
      WRITE(20,'(/a)') 'Energy Eigenvalues and Eigenvectors (in column)'
      WRITE(20,'('//form_str1//'(7x,f10.7,7x))') eigval
      DO i=1,ss
         WRITE(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(eigvect(i,j)),',',AIMAG(eigvect(i,j)),'),',j=1,ss)
      END DO
   END IF
   WRITE(20, '(/a)') 'The Gibbs state for this Hamiltonian is'
   DO i=1,ss
      WRITE(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(rho_g(i,j)),',',AIMAG(rho_g(i,j)),'),',j=1,ss)
   END DO

   WRITE(20, '(/a)') 'Coupling Operator to Environment'
   WRITE(20, '(a/)') '================================'
   WRITE(20,'(a)') 'VI = '
   DO i=1,ss
      WRITE(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(VI(i,j)),',',AIMAG(VI(i,j)),'),',j=1,ss)
   END DO
   tmp_l1 = test_hermitian(VI)
   IF (tmp_l1) THEN
      WRITE(20,'(/a)') 'The Coupling Operator is Self-Adjoint.'
   ELSE
      WRITE(20,'(/a)') '**** The Couplng Operator is NOT Self-Adjoint. ****'
      halt = .TRUE.
   END IF
   IF (tmp_l1) THEN
      CALL eigensystem(VI,eigval,eigvect)
      WRITE(20,'(/a)') 'Coupling Operator Eigenvalues and Eigenvectors (in column)'
      WRITE(20,'('//form_str1//'(7x,f10.7,7x))') eigval
      DO i = 1, ss
         WRITE(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(eigvect(i,j)),',',AIMAG(eigvect(i,j)),'),',j=1,ss)
      END DO
   END IF

   WRITE(20, '(/a)') 'Initial Density Matrix'
   WRITE(20, '(a/)') '======================'
   IF (gibbs) THEN
      WRITE(20,'(/a)') 'A Gibbs state is used for the initial density.'
   ELSE IF (random) THEN
      WRITE(20,'(/a)') 'A random state is used for the initial density.'
   END IF
   WRITE(20,'(a)') 'rho0 = '
   DO i=1,ss
      WRITE(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(rho0(i,j)),',',AIMAG(rho0(i,j)),'),',j=1,ss)
   END DO
   tmp_l1 = test_hermitian(rho(0,:,:))
   IF (tmp_l1) THEN
      WRITE(20,'(/a)') 'The Initial Density is Self-Adjoint.'
   ELSE
      WRITE(20,'(/a)') '**** The Initial Density is NOT Self-Adjoint. ****'
      halt = .TRUE.
   END IF
   tmp_l1 = test_trace(rho0)
   IF (tmp_l1) THEN
      WRITE(20,'(a,2f10.6)') 'The trace is OK. Trace = ', trace(rho0)
   ELSE
      WRITE(20,'(a)') '**** The Initial Density is NOT Normalized. ****'
      halt = .TRUE.
   END IF
   tmp_l1 = test_positivity(rho0)
   IF (tmp_l1) THEN
      WRITE(20,'(a)') 'The Initial Density is Non-Negative.'
   ELSE
      WRITE(20,'(a)') '**** The Initial Density is NOT Non-Negative. ****'
      halt = .TRUE.
   END IF
   IF (tmp_l1) THEN
      CALL eigensystem(rho0,eigval,eigvect)
      WRITE(20,'(/a)') 'Initial Density Eigenvalues and Eigenvectors (in column)'
      WRITE(20,'('//form_str1//'(7x,f10.7,7x))') eigval
      DO i=1,ss
         WRITE(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(eigvect(i,j)),',',AIMAG(eigvect(i,j)),'),',j=1,ss)
      END DO
   END IF

   IF (halt) THEN
      WRITE(20, '(/a)') '**** Bad Input Data ****'
      WRITE(20, '(a)') '**** Execution is Stopped ****'
      STOP ''
   END IF

   WRITE(20, '(/a)') 'Environment Parameters'
   WRITE(20, '(a/)') '======================'
   WRITE(20, '(a,f7.5)') '     T = ', temp
   WRITE(20, '(a,f7.5)') ' gamma = ', gamma
   WRITE(20, '(a,f7.5)') 'lambda = ', lambda

   WRITE(20,'(/a)') 'Correlation functions'
   WRITE(20,'(a/)') '====================='
   IF (pade) THEN
      WRITE(20,'(a/)') '[1/1] Pade approximation is used.'
      WRITE(20,'(a,f10.7,a,f10.7,a)') 'cinf = (', REAL(cinf), ',',AIMAG(cinf),')'
      WRITE(20,'(a,f10.7,a,f10.7,a)') 'c0 = (',REAL(coeff(0)),',',AIMAG(coeff(0)),')'
      WRITE(20,'(a,f10.7)') 'gamma0 = ', ABS(exp_vec(0))
      WRITE(20,'(a,f10.7,a,f10.7,a)') 'c1 = (',REAL(coeff(1)),',',AIMAG(coeff(1)),')'
      WRITE(20,'(a,f10.7)') 'gamma1 = ', ABS(exp_vec(1))
   ELSE
      WRITE(tmp_str1, '(I3)') matsu
      WRITE(20,'(2a/)') TRIM(ADJUSTL(tmp_str1)), ' matsubara frequencies are used.'
      WRITE(tmp_str1, '(E10.3)') REAL(cinf)
      WRITE(20,'(2a)') 'cinf    = ', TRIM(ADJUSTL(tmp_str1))
      WRITE(20,'(a,E10.3,a,E10.3,a)') 'c0      = (',REAL(coeff(0)),',',AIMAG(coeff(0)),')'
      WRITE(tmp_str1, '(f10.7)') ABS(exp_vec(0))
      WRITE(20, '(2a)') 'gamma0  = ', TRIM(ADJUSTL(tmp_str1))
      DO i = 1, matsu
         WRITE(tmp_str1, '(I3)') i
         WRITE(tmp_str2, '(E10.3)') REAL(coeff(i))
         WRITE(20,'(4a)') 'c',TRIM(ADJUSTL(tmp_str1)),'      = ', TRIM(ADJUSTL(tmp_str2))
         WRITE(tmp_str2, '(f10.7)') ABS(exp_vec(i))
         WRITE(20, '(4a)') 'gamma',TRIM(ADJUSTL(tmp_str1)),'  = ', TRIM(ADJUSTL(tmp_str2))
      END DO
   END IF

   WRITE(20, '(/a)') 'Execution Parameters'
   WRITE(20, '(a/)') '===================='
   WRITE(20, '(a,f8.3)') 'Time limit = ', time_limit
   WRITE(20, '(a,f5.3)') '        dt = ', dt

   CALL CPU_TIME(cputime0)

   ! Create the array of values for lambda_bc, to be stored in bath_VI
   CALL lambda_bc(VI, bath_VI, bath_VI_half)

   DO i = 0, n_steps - 1
      ti = i * dt

      tc = ti
      tmp_arr1 = rho(i,:,:)
      tmp_arr2 = bath_VI(n_steps,:,:)
      tmp_arr3 = MATMUL(tmp_arr2, tmp_arr1) - MATMUL(tmp_arr1, TRANSPOSE(CONJG(tmp_arr2)))
      k1 = -CMPLX(0._DP,1._DP)*(MATMUL(HS,tmp_arr1)-MATMUL(tmp_arr1,HS)) - (MATMUL(VI,tmp_arr3)-MATMUL(tmp_arr3,VI))

      tc = ti + dt / 2
      tmp_arr1 = rho(i,:,:) + k1 * dt / 2
      tmp_arr2 = bath_VI(n_steps,:,:)
      tmp_arr3 = MATMUL(tmp_arr2, tmp_arr1) - MATMUL(tmp_arr1, TRANSPOSE(CONJG(tmp_arr2)))
      k2 = -CMPLX(0._DP,1._DP)*(MATMUL(HS,tmp_arr1)-MATMUL(tmp_arr1,HS)) - (MATMUL(VI,tmp_arr3)-MATMUL(tmp_arr3,VI))

      tc = ti + dt / 2
      tmp_arr1 = rho(i,:,:) + k2 * dt / 2
      tmp_arr2 = bath_VI(n_steps,:,:)
      tmp_arr3 = MATMUL(tmp_arr2, tmp_arr1) - MATMUL(tmp_arr1, TRANSPOSE(CONJG(tmp_arr2)))
      k3 = -CMPLX(0._DP,1._DP)*(MATMUL(HS,tmp_arr1)-MATMUL(tmp_arr1,HS)) - (MATMUL(VI,tmp_arr3)-MATMUL(tmp_arr3,VI))

      tc = ti + dt
      tmp_arr1 = rho(i,:,:) + k3 * dt
      tmp_arr2 = bath_VI(n_steps,:,:)
      tmp_arr3 = MATMUL(tmp_arr2, tmp_arr1) - MATMUL(tmp_arr1, TRANSPOSE(CONJG(tmp_arr2)))
      k4 = -CMPLX(0._DP,1._DP)*(MATMUL(HS,tmp_arr1)-MATMUL(tmp_arr1,HS)) - (MATMUL(VI,tmp_arr3)-MATMUL(tmp_arr3,VI))

      rho(i+1,:,:) = rho(i,:,:) + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
   END DO

   ! Calculate entropy
   DO i = 0, n_steps
      sys_entropy(i) = Entropy(rho(i,:,:))
   END DO

   ! Calculate eta
   DO i = 0, n_steps
      eta(i,:,:) = CMPLX(0._DP,1._DP)*(MATMUL(rho(i,:,:),TRANSPOSE(CONJG(bath_VI(i,:,:))))-MATMUL(bath_VI(i,:,:),rho(i,:,:)))
   END DO

   ! Calculate system energy
   DO i = 0, n_steps
      energy(i) = trace(MATMUL(HS, rho(i,:,:)))
   END DO

   ! Calculate heat
   DO i = 0, n_steps
      tmp_c1 = energy(i) - energy(0)
      tmp_c2 = trace(MATMUL(VI,eta(i,:,:))) - trace(MATMUL(VI,eta(0,:,:)))
      heat(i) = tmp_c1 + tmp_c2
   END DO

   ! Calculate entropy production
   DO i = 0, n_steps
      sprod(i) = (sys_entropy(i)-sys_entropy(0)) - heat(i) / temp
   END DO

   ! Calculate fidelity
   DO i = 0, n_steps
      fid(i) = fidelity(rho(i,:,:),rho_g)
   END DO

   ! Calculate trace distance
   DO i = 0, n_steps
      traced(i) = trace_distance(rho(i,:,:),rho_g)
   END DO

   ! Calculate system density in the energy eigenbasis

   CALL eigensystem(HS, eigval, eigvect)
   CALL matinv4(eigvect, eigvect_inv4)
   DO i = 0, n_steps
      rho_E(i,:,:) = MATMUL(eigvect_inv4, MATMUL(rho(i,:,:), eigvect))
   END DO

   CALL CPU_TIME(cputime1)

   OPEN(10, FILE='rho.dat')
   DO i = 0, n_steps
      ti = i * dt
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10,'(f10.3,('//form_str2//'e15.6))') ti,((rho(i,j,k),j=1,ss),k=1,ss)
   END DO
   CLOSE(10)

   OPEN(10, FILE='rho_E.dat')
   DO i = 0, n_steps
      ti = i * dt
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10,'(f10.3,('//form_str2//'e15.6))') ti,((rho_E(i,j,k),j=1,ss),k=1,ss)
   END DO
   CLOSE(10)

   OPEN(10, FILE='lambda.dat')
   DO i = 0, n_steps
      ti = i * dt
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10,'(f10.3,('//form_str2//'e15.6))') ti,((bath_VI(i,j,k),j=1,ss),k=1,ss)
   END DO
   CLOSE(10)

   OPEN(10, FILE='eta.dat')
   DO i = 0, n_steps
      ti = i * dt
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10,'(f10.3,('//form_str2//'e15.6))') ti,((eta(i,j,k),j=1,ss),k=1,ss)
   END DO
   CLOSE(10)

   OPEN(10, FILE='hermitian.dat')
   found_herm = .FALSE.
   DO i = 0, n_steps
      tmp_l1 = test_hermitian(rho(i,:,:))
      IF (.NOT. tmp_l1) THEN
         ti = i * dt
         WRITE(10, '(f10.3,L2)') ti, tmp_l1
         found_herm = .TRUE.
      END IF
   END DO
   IF (.NOT. found_herm) WRITE(10, '(a)') 'The density remains hermitian at all times.'
   CLOSE(10)

   OPEN(10, FILE='trace.dat')
   found_trace = .FALSE.
   DO i = 0, n_steps
      tmp_l1 = test_trace(rho(i,:,:))
      IF (.NOT. tmp_l1) THEN
         ti = i * dt
         WRITE(10, '(f10.3,L2)') ti, tmp_l1
         found_trace = .TRUE.
      END IF
   END DO
   IF (.NOT. found_trace) WRITE(10, '(a)') 'The density remains normalized at all times.'
   CLOSE(10)

   OPEN(10, FILE='positivity.dat')
   found_pos = .FALSE.
   DO i = 0, n_steps
      tmp_l1 = test_positivity(rho(i,:,:))
      IF (.NOT. tmp_l1) THEN
         ti = i * dt
         WRITE(10, '(f10.3,L2)') ti, tmp_l1
         found_pos = .TRUE.
      END IF
   END DO
   IF (.NOT. found_pos) WRITE(10, '(a)') 'The density remains positive at all times.'
   CLOSE(10)

   OPEN(10, FILE='entropy.dat')
   DO i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(sys_entropy(i), KIND=DP)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10, '(f10.3,e15.6)') ti, tmp_r1
   END DO
   CLOSE(10)

   OPEN(10, FILE='stats_E.dat')
   DO i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(energy(i), KIND=DP)
      tmp_r2 = REAL(trace(MATMUL(HS, MATMUL(HS, rho(i,:,:)))), KIND=DP)
      tmp_r3 = REAL(SQRT(tmp_r2 - tmp_r1**2), KIND=DP)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10, '(f10.3,2(e15.6))') ti, tmp_r1, tmp_r3
   END DO
   CLOSE(10)

   OPEN(10, FILE='stats_XS.dat')
   DO i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(trace(MATMUL(VI, rho(i,:,:))), KIND=DP)
      tmp_r2 = REAL(trace(MATMUL(VI, MATMUL(VI, rho(i,:,:)))), KIND=DP)
      tmp_r3 = REAL(SQRT(tmp_r2 - tmp_r1**2), KIND=DP)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10, '(f10.3,2(e15.6))') ti, tmp_r1, tmp_r3
   END DO
   CLOSE(10)

   OPEN(10, FILE='heat.dat')
   DO i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(heat(i), KIND=DP)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10, '(f10.3,e15.6)') ti, tmp_r1
   END DO
   CLOSE(10)

   OPEN(10, FILE='sprod.dat')
   DO i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(sprod(i), KIND=DP)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10, '(f10.3,e15.6)') ti, tmp_r1
   END DO
   CLOSE(10)

   OPEN(10, FILE='fidelity.dat')
   DO i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(fid(i), KIND=DP)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10, '(f10.3,e15.6)') ti, tmp_r1
   END DO
   CLOSE(10)

   OPEN(10, FILE='trace_distance.dat')
   DO i = 0, n_steps
      ti = i * dt
      tmp_r1 = REAL(traced(i), KIND=DP)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10, '(f10.3,e15.6)') ti, tmp_r1
   END DO
   CLOSE(10)

   OPEN(10, FILE='rhoA.dat')
   DO i = 0, n_steps
      ti = i * dt
      CALL rhoA(rho(i,:,:), rho_A)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10,'(f10.3,8e15.6)') ti, ((rho_A(k,j),k=1,2),j=1,2)
   END DO
   CLOSE(10)

   OPEN(10, FILE='rhoB.dat')
   DO i = 0, n_steps
      ti = i * dt
      CALL rhoB(rho(i,:,:), rho_B)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) WRITE(10,'(f10.3,8e15.6)') ti, ((rho_B(k,j),k=1,2),j=1,2)
   END DO
   CLOSE(10)

   ! First, calculate Gibbs state for the single qubit Hamiltonian

   rho_g2(1,:) = (/CMPLX(DEXP(0.5*temp),0.), CMPLX(0.,0.)/)
   rho_g2(2,:) = (/CMPLX(0.,0.), CMPLX(DEXP(-0.5*temp),0.)/)
   Z = trace(rho_g2)
   rho_g2 = rho_g2 / Z

   OPEN(10, FILE='rhoA_E.dat')
   OPEN(11, FILE='rhoA_fid.dat')
   OPEN(12, FILE='rhoA_traced.dat')
   DO i = 0, n_steps
      ti = i * dt
      CALL rhoA(rho_E(i,:,:), rho_AE)
      tmp_r1 = fidelity2(rho_AE, rho_g2)
      tmp_r2 = REAL(trace_distance2(rho_AE, rho_g2), KIND=DP)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) THEN
         WRITE(10,'(f10.3,8e15.6)') ti, ((rho_AE(k,j),k=1,2),j=1,2)
         WRITE(11,'(2f10.3)') ti, tmp_r1
         WRITE(12,'(2f10.3)') ti, tmp_r2
      END IF
   END DO
   CLOSE(12)
   CLOSE(11)
   CLOSE(10)

   OPEN(10, FILE='rhoB_E.dat')
   OPEN(11, FILE='rhoB_fid.dat')
   OPEN(12, FILE='rhoB_traced.dat')
   DO i = 0, n_steps
      ti = i * dt
      CALL rhoB(rho_E(i,:,:), rho_BE)
      tmp_r1 = fidelity2(rho_BE, rho_g2)
      tmp_r2 = REAL(trace_distance2(rho_BE, rho_g2), KIND=DP)
      mod = ti - NINT(ti/time_write)*time_write
      IF (mod .EQ. 0) THEN
         WRITE(10,'(f10.3,8e15.6)') ti, ((rho_BE(k,j),k=1,2),j=1,2)
         WRITE(11,'(2f10.3)') ti, tmp_r1
         WRITE(12,'(2f10.3)') ti, tmp_r2
      END IF
   END DO
   CLOSE(12)
   CLOSE(11)
   CLOSE(10)

   CALL CPU_TIME(cputime2)

   WRITE(20, '(/a)') 'Final Density Matrix'
   WRITE(20, '(a/)') '======================'
   WRITE(20,'(a)') 'rho(final) = '
   tmp_arr1 = rho(n_steps,:,:)
   DO i=1,ss
      WRITE(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(tmp_arr1(i,j)),',',AIMAG(tmp_arr1(i,j)),'),',j=1,ss)
   END DO
   tmp_l1 = test_hermitian(tmp_arr1)
   IF (tmp_l1) THEN
      WRITE(20,'(/a)') 'The Final Density is Self-Adjoint.'
   ELSE
      WRITE(20,'(/a)') '**** The Final Density is NOT Self-Adjoint. ****'
   END IF
   tmp_l1 = test_trace(tmp_arr1)
   IF (tmp_l1) THEN
      WRITE(20,'(a,2f10.6)') 'The trace is OK. Trace = ', trace(tmp_arr1)
   ELSE
      WRITE(20,'(a)') '**** The Final Density is NOT Normalized. ****'
   END IF
   tmp_l1 = test_positivity(tmp_arr1)
   IF (tmp_l1) THEN
      WRITE(20,'(a)') 'The Final Density is Non-Negative.'
   ELSE
      WRITE(20,'(a)') '**** The Final Density is NOT Non-Negative. ****'
   END IF
   WRITE(20,'(a,f10.7)') 'The final fidelity between the system density and the Gibbs state is ', REAL(fid(n_steps))
   WRITE(20,'(a,f10.7)') 'The final trace distance between the system density and the Gibbs state is ', REAL(traced(n_steps))

   tmp_r1 = fidelity(rho(n_steps-1,:,:), rho(n_steps,:,:))
   tmp_r2 = trace_distance(rho(n_steps-1,:,:), rho(n_steps,:,:))
   WRITE(20,'(a,f10.7)') 'The fidelity between the system density at the last two time steps is', tmp_r1
   WRITE(20,'(a,f10.7)') 'The trace distance between the system density at the last two time steps is', tmp_r2
   IF ((ABS(1-tmp_r1) .LT. tol) .AND. (ABS(tmp_r2) .LT. tol)) THEN
      WRITE(20,'(a)') 'The values for the fidelity and the trace distance indicate that the system'
      WRITE(20,'(a)') 'has likely reached equilibrium.'
   END IF

   IF (tmp_l1) THEN
      CALL eigensystem(tmp_arr1,eigval,eigvect)
      WRITE(20,'(/a)') 'Final Density Eigenvalues and Eigenvectors (in column)'
      WRITE(20,'('//form_str1//'(7x,f10.7,7x))') eigval
      DO i=1,ss
         WRITE(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(eigvect(i,j)),',',AIMAG(eigvect(i,j)),'),',j=1,ss)
      END DO
   END IF

   WRITE(20,'(/a)') 'The final density in the energy eigenbasis is'
   tmp_arr1 = rho_E(n_steps,:,:)
   DO i=1,ss
      WRITE(20,'('//form_str1//'(a,f10.7,a,f10.7,a))') ('(',REAL(tmp_arr1(i,j)),',',AIMAG(tmp_arr1(i,j)),'),',j=1,ss)
   END DO

   WRITE(20,'(/a)') 'The reduced density for system A is'
   WRITE(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_A(1,j)),',',AIMAG(rho_A(1,j)),'),',j=1,2)
   WRITE(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_A(2,j)),',',AIMAG(rho_A(2,j)),'),',j=1,2)
   WRITE(20,'(/a)') 'The reduced density for system B is'
   WRITE(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_B(1,j)),',',AIMAG(rho_B(1,j)),'),',j=1,2)
   WRITE(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_B(2,j)),',',AIMAG(rho_B(2,j)),'),',j=1,2)
   WRITE(20,'(/a)') 'The reduced density for system A in the energy eigenbasis is'
   WRITE(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_AE(1,j)),',',AIMAG(rho_AE(1,j)),'),',j=1,2)
   WRITE(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_AE(2,j)),',',AIMAG(rho_AE(2,j)),'),',j=1,2)
   WRITE(20,'(/a)') 'The reduced density for system B in the energy eigenbasis is'
   WRITE(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_BE(1,j)),',',AIMAG(rho_BE(1,j)),'),',j=1,2)
   WRITE(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_BE(2,j)),',',AIMAG(rho_BE(2,j)),'),',j=1,2)
   WRITE(20,'(/a)') 'The Gibbs state for the individual qubits is'
   WRITE(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_g2(1,j)),',',AIMAG(rho_g2(1,j)),'),',j=1,2)
   WRITE(20,'(2(a,f10.7,a,f10.7,a))') ('(',REAL(rho_g2(2,j)),',',AIMAG(rho_g2(2,j)),'),',j=1,2)

   WRITE(20, '(/a)') 'Other Information'
   WRITE(20, '(a/)') '================='

   IF (found_herm) THEN
      WRITE(20, '(a)') 'Hermiticity is NOT preserved at all times. See hermitian.dat'
   ELSE
      WRITE(20, '(a)') 'Hermiticity is preserved at all times.'
   END IF

   IF (found_trace) THEN
      WRITE(20, '(a)') 'Trace is NOT preserved at all times. See trace.dat'
   ELSE
      WRITE(20, '(a)') 'Trace is preserved at all times.'
   END IF

   IF (found_pos) THEN
      WRITE(20, '(a)') 'Positivity is NOT preserved at all times. See positivity.dat'
   ELSE
      WRITE(20, '(a)') 'Positivity is preserved at all times.'
   END IF

   WRITE(20,'(/a,f10.5)') 'Time needed for time evolution: ', cputime1 - cputime0
   WRITE(20,'(a,f10.5)') 'Time needed for data output: ', cputime2 - cputime1
   WRITE(20,'(/a,a)') 'Job finished at ', time_stamp()

   CLOSE(20)

END PROGRAM main
