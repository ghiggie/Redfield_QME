MODULE bath

   USE constants
   USE variables
   USE helper

   IMPLICIT NONE

CONTAINS

   SUBROUTINE bc_coeff()

      ! This subroutine will allocate a vector with matsu+1 entries, where
      ! matsu is the number of Matsubara terms. The vector will look like
      !    coeff(0:matsu)
      ! The coefficient for the delta term will be contained in cinf. coeff(0)
      ! will correspond to the coefficient of the exponential gamma term.
      REAL(KIND=DP) :: cotan, nu, tmp_r
      INTEGER       :: i

      IF (pade) THEN
         matsu = 1
         ALLOCATE(  coeff(0:matsu))
         ALLOCATE(exp_vec(0:matsu))
         cotan          = 2*temp/gamma - gamma/(20._DP*temp)-4.9_DP*gamma*temp/(42._DP*temp**2-gamma**2)
         coeff(0)       = lambda*gamma*CMPLX(cotan, -1._DP)
         exp_vec(0)     = -gamma
         coeff(matsu)   = 4.9_DP*lambda*gamma*temp/(42._DP*temp**2-gamma**2)
         exp_vec(matsu) = -temp*DSQRT(42.0_DP)
         cinf           = 0.1_DP * lambda * gamma / temp
      ELSE
         ALLOCATE(  coeff(0:matsu))
         ALLOCATE(exp_vec(0:matsu))
         cotan = DCOS(0.5_DP*gamma/temp)/DSIN(0.5_DP*gamma/temp)

         ! Set up gamma coefficient
         coeff(0)   = lambda*gamma*CMPLX(cotan,-1._DP)
         exp_vec(0) = -gamma
         ! Set up Matsubara coefficients and compute cinf
         tmp_r = 0
         DO i = 1, matsu
            nu         = 2 * PI * temp * i
            coeff(i)   = 4*lambda*gamma*temp*nu/(nu**2 - gamma**2)
            exp_vec(i) = -nu

            tmp_r = tmp_r + 1 / (nu**2 - gamma**2)
         END DO
            tmp_r = 8*lambda*gamma*temp*tmp_r
            cinf  = 4*lambda*temp/gamma - 2*lambda*cotan - tmp_r
      END IF
   END SUBROUTINE bc_coeff

   SUBROUTINE lambda_bc(A, IA, IA_half)
      COMPLEX(KIND=DP), DIMENSION(:,:),    INTENT(IN)  :: A
      COMPLEX(KIND=DP), DIMENSION(0:,:,:), INTENT(OUT) :: IA
      COMPLEX(KIND=DP), DIMENSION(:,:,:),  INTENT(OUT) :: IA_half

      COMPLEX(KIND=DP), DIMENSION(ss,ss) :: f2j, f2j1, f2j2, tmp_arr

      INTEGER          :: M, i, j
      REAL(KIND=DP)    :: ti, tmp_dt, t2j, t2j1, t2j2, tau
      COMPLEX(KIND=DP) :: corr, tc

      IA        = 0
      IA_half   = 0
      IA(0,:,:) = 0.5 * cinf * A
      ! Calculate the integral at ti
      DO i = 0, n_steps - 1
         ti = i * dt
         M = 20
         ! Begin Simpson Loop
         tmp_arr = 0
         tmp_dt  = dt / M
         DO j = 0, M/2 -1
            ! First calculate f2j
            tc = ti + 2*j * tmp_dt
            CALL I2S(tc*HS,VI,f2j)
            corr = SUM(coeff * DEXP(REAL(tc*exp_vec, KIND=DP)))
            f2j  = corr * f2j
            ! Calculate f2j1
            tc = ti + (2*j + 1) * tmp_dt
            CALL I2S(tc*HS,VI,f2j1)
            corr = SUM(coeff * DEXP(REAL(tc * exp_vec, KIND=DP)))
            f2j1 = corr * f2j1
            ! Calculate f2j2
            tc = ti + (2*j + 2) * tmp_dt
            CALL I2S(tc*HS,VI,f2j2)
            corr = SUM(coeff * DEXP(REAL(tc * exp_vec, KIND=DP)))
            f2j2 = corr * f2j2

            tmp_arr = tmp_arr + (tmp_dt / 3.0_DP) * (f2j+4.0_DP*f2j1+f2j2)
            IF (j .EQ. M/4 -1) THEN
               IA_half(i+1,:,:) = IA(i,:,:) + tmp_arr
            END IF
         END DO
         IA(i+1,:,:) = IA(i,:,:) + tmp_arr
      END DO
   END SUBROUTINE lambda_bc

END MODULE bath
