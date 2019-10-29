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
            cotan = 2*temp/gamma - gamma/(20*temp)-4.9*gamma*temp/(42*temp**2-gamma**2)
            coeff(0) = lambda*gamma*CMPLX(cotan, -1)
            exp_vec(0) = -gamma
            coeff(1) = 4.9*lambda*gamma*temp/(42*temp**2-gamma**2)
            exp_vec(1) = -temp*SQRT(42)
            cinf = 0.1 * lambda * gamma / temp
        else
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

    subroutine lambda_bc(A, IA)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        complex(kind=DP), dimension(:,:,:), intent(out) :: IA

        integer :: M, i, j
        real(kind=DP) :: ti, tmp_dt, t2j, t2j1, t2j2, tau
        complex(kind=DP) :: corr

        IA = 0
        IA(0,:,:) = 0.5 * cinf * A
        M = 10 ! Number of Simpson steps. Be sure it's even

        ! This loop will calculate the integral from 0 to time_limit
        do i = 0, n_steps - 1
            ! Compute the ith. integral
            ti = i * dt
            tmp_dt = dt / M

            ! This loop will compute the individual integral from t_i to t_(i+1)
            tmp_arr1 = 0
            do j = 0, M/2 - 1
                t2j = ti + 2 * j * tmp_dt
                t2j1 = t2j + tmp_dt
                t2j2 = t2j1 + tmp_dt

                tau = (ti+dt) - t2j
                ! Build the non-delta coefficients for G(\tau)
                corr = sum(coeff * EXP(exp_vec*tau))
                ! Construct the matrix exponential
                call I2S(tau * HS, A, tmp_arr2)
                tmp_arr2 = corr * tmp_arr2

                tau = (ti+dt) - t2j1
                ! Build the non-delta coefficients for G(\tau)
                corr = sum(coeff * EXP(exp_vec*tau))
                ! Construct the matrix exponential
                call I2S(tau * HS, A, tmp_arr3)
                tmp_arr3 = corr * tmp_arr3

                tau = (ti+dt) - t2j2
                ! Build the non-delta coefficients for G(\tau)
                corr = sum(coeff * EXP(exp_vec*tau))
                ! Construct the matrix exponential
                call I2S(tau * HS, A, tmp_arr4)
                tmp_arr4 = corr * tmp_arr4

                tmp_arr1 = tmp_arr1 + (tmp_dt/3) * (tmp_arr2+4*tmp_arr3+tmp_arr4)
            end do
            IA(i+1,:,:) = IA(i,:,:) + tmp_arr1
        end do
    end subroutine lambda_bc

end module bath
