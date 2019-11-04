module bath

    use constants
    use variables
    use helper

    implicit none

 contains

    subroutine bc_coeff()

        ! This subroutine will allocate a vector with matsu+1 entries, where
        ! matsu is the number of Matsubara terms. The vector will look like
        !    coeff(0:matsu)
        ! The coefficient for the delta term will be contained in cinf. coeff(0)
        ! will correspond to the coefficient of the exponential gamma term.

        ! We can also edit this subroutine to use the Pade approximation.
        real(kind=DP) :: cotan, nu, tmp_r
        integer :: i

        if (pade) then
            matsu = 1
            allocate(coeff(0:matsu))
            allocate(exp_vec(0:matsu))
            cotan = 2*temp/gamma - gamma/(20*temp)-4.9*gamma*temp/(42*temp**2-gamma**2)
            coeff(0) = lambda*gamma*CMPLX(cotan, -1)
            exp_vec(0) = -gamma
            coeff(matsu) = 4.9*lambda*gamma*temp/(42*temp**2-gamma**2)
            exp_vec(matsu) = -temp*SQRT(42.0_DP)
            cinf = 0.1 * lambda * gamma / temp
        else
            allocate(coeff(0:matsu))
            allocate(exp_vec(0:matsu))
            cotan = cos(0.5*gamma/temp)/sin(0.5*gamma/temp)

            ! Set up gamma coefficient
            coeff(0) = lambda*gamma*CMPLX(cotan,-1)
            exp_vec(0) = -gamma
            ! Set up Matsubara coefficients and compute cinf
            tmp_r = 0
            do i = 1, matsu
                nu = 2 * PI * temp * i
                coeff(i) = 4*lambda*gamma*temp*nu/(nu**2 - gamma**2)
                exp_vec(i) = -nu

                tmp_r = tmp_r + 1 / (nu**2 - gamma**2)
            end do
            tmp_r = 8*lambda*gamma*temp*tmp_r
            cinf = 4*lambda*temp/gamma - 2*lambda*cotan - tmp_r
        end if
    end subroutine bc_coeff

    subroutine lambda_bc(A, IA, IA_half)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        complex(kind=DP), dimension(0:,:,:), intent(out) :: IA
        complex(kind=DP), dimension(:,:,:), intent(out) :: IA_half

        complex(kind=DP), dimension(ss,ss) :: f2j, f2j1, f2j2, tmp_arr

        integer :: M, i, j
        real(kind=DP) :: ti, tmp_dt, t2j, t2j1, t2j2, tau
        complex(kind=DP) :: corr, tc

        IA = 0
        IA_half = 0
        IA(0,:,:) = 0.5 * cinf * A
        ! Calculate the integral at ti
        do i = 0, n_steps - 1
            ti = i * dt
            M = 20
            ! Begin Simpson Loop
            tmp_arr = 0
            tmp_dt = dt / M
            do j = 0, M/2 -1
                ! First calculate f2j
                tc = ti + 2*j * tmp_dt
                call I2S(tc*HS,VI,f2j)
                corr =sum(coeff*EXP(tc*exp_vec))
                f2j = corr * f2j
                ! Calculate f2j1
                tc = ti + (2*j + 1) * tmp_dt
                call I2S(tc*HS,VI,f2j1)
                corr = sum(coeff * EXP(tc * exp_vec))
                f2j1 = corr * f2j1
                ! Calculate f2j2
                tc = ti + (2*j + 2) * tmp_dt
                call I2S(tc*HS,VI,f2j2)
                corr = sum(coeff * EXP(tc * exp_vec))
                f2j2 = corr * f2j2

                tmp_arr = tmp_arr + (tmp_dt / 3.0_DP) * (f2j+4.0_DP*f2j1+f2j2)
                if (j .eq. M/4 -1) then
                    IA_half(i+1,:,:) = IA(i,:,:) + tmp_arr
                end if
            end do
            IA(i+1,:,:) = IA(i,:,:) + tmp_arr
        end do
    end subroutine lambda_bc

end module bath
