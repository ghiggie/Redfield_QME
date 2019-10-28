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

        if (pade) then
            cotan = 2*temp/gamma - gamma/(20*temp)-4.9*gamma*temp/(42*T**2-gamma**2)
            coeff(0) = lambda*gamma*CMPLX(cotan, -1)
            exp_vec(0) = -gamma
            coeff(1) = 4.9*lambda*gamma*temp/(42*temp**2-gamma**2)
            exp_vec(1) = -temp*SQRT(42)
            cinf = 0.1 * lambda * gamma / temp
        else
            real(kind=DP) :: cotan
            cotan = cos(0.5*gamma/temp)/sin(0.5*gamma/temp)

            ! Set up gamma coefficient
            coeff(0) = lambda*gamma*CMPLX(cotan,-1)
            exp_vec(0) = -gamma
            ! Set up Matsubara coefficients
            do i = 1, matsu
                nu = 2 * PI * temp * i
                coeff(i) = 4*lambda*gamma*temp*nu/(nu**2 - gamma**2)
                exp_vec(i) = -nu
            end do
            ! Calculate cinf
            tmp_r = 0
            do i = 1, matsu
                nu = 2 * PI * temp * i
                tmp_r = tmp_r + 1 / (nu**2 - gamma**2)
            end do
            tmp_r = 8*lambda*gamma*temp*tmp_r
            cinf = 4*lambda*temp/gamma - 2*lambda*cotan - tmp_r
        end if
    end subroutine bc_coeff

    subroutine lambda_bc(A, IA)
        complex(kind=DP), dimension(:,:), intent(in) :: A
        complex(kind=DP), dimension(:,:,:) :: IA

        IA = 0
        IA(0,:,:) = 0.5 * cinf * A
        M = 10 ! Number of Simpson steps. Be sure it's even

        ! This loop will calculate the integral from 0 to time_limit
        do i = 0, n_steps - 1
            ! Compute the ith. integral
            ti = i * dt
            tmp_dt = dt / M

            ! This loop will compute the individual integral from t_i to t_(i+1)
            dI = 0
            do j = 0, M/2 - 1
                t2j = ti + 2 * j * tmp_dt
                t2j1 = t2j + tmp_dt
                t2j2 = t2j1 + tmp_dt

                tau = (ti+dt) - t2j
                ! Build the non-delta coefficients for G(\tau)
                corr = sum(coeff * EXP(exp_vec*tau))
                ! Construct the matrix exponential
                tmp_arr = 0
                call I2S(tau * HS, A, tmp_arr)
                tmp_arr1 = corr * tmp_arr

                tau = (ti+dt) - t2j1
                ! Build the non-delta coefficients for G(\tau)
                corr = sum(coeff * EXP(exp_vec*tau))
                ! Construct the matrix exponential
                tmp_arr = 0
                call I2S(tau * HS, A, tmp_arr)
                tmp_arr2 = corr * tmp_arr

                tau = (ti+dt) - t2j2
                ! Build the non-delta coefficients for G(\tau)
                corr = sum(coeff * EXP(exp_vec*tau))
                ! Construct the matrix exponential
                tmp_arr = 0
                call I2S(tau * HS, A, tmp_arr)
                tmp_arr2 = corr * tmp_arr

                dI = dI + (tmp_dt/3) * (tmp_arr1+4*tmp_arr2+tmp_arr3)
            end do
            IA(i+1,:,:) = IA(i,:,:) + dI
        end do

    end subroutine lambda_bc

    !!!! This subroutine is complete garbage
    ! subroutine lambda_bc(A, B, C)
    !     real(kind=DP) :: tmp_dt ! I allow for the possibility that this should be shifted
    !     complex(kind=DP), dimension(:,:), intent(in) :: A, B
    !     complex(kind=DP), dimension(:,:), intent(out) :: C
    !     integer :: M, j
    !     real(kind=DP) :: t2j, t2j1, t2j2
    !     complex(kind=DP), dimension(:,:), allocatable :: tmp, tmp1, tmp2
    !     complex(kind=DP), dimension(:,:), allocatable :: f2j, f2j1, f2j2
    !
    !     allocate(lambda_bc(ss,ss))
    !
    !     tmp_dt = dt
    !     tmp = 0
    !     M = nint(t / tmp_dt)
    !     if (mod(M,2) .ne. 0) then
    !         M = M + 1
    !         tmp_dt = t / M
    !     end if
    !     do j = 0, M/2 - 1
    !         t2j = 2*j * tmp_dt
    !         t2j1 = t2j + tmp_dt
    !         t2j2 = t2j1 + tmp_dt
    !
    !         tmp1 = (t - t2j) * A
    !         tmp2 = I2S(tmp1, B)
    !         f2j = bc(temp, gamma, lambda, t, t2j) * tmp2
    !
    !         tmp1 = (t - t2j1) * A
    !         tmp2 = I2S(tmp1, B)
    !         f2j1 = bc(temp, gamma, lambda, t, t2j1) * tmp2
    !
    !         tmp1 = (t - t2j2) * A
    !         tmp2 = I2S(tmp1, B)
    !         f2j2 = bc(temp, gamma, lambda, t, t2j2) * tmp2
    !
    !         tmp = tmp + (tmp_dt / 3) * (f2j + 4*f2j1 + f2j2)
    !     end do
    !     lambda_bc = tmp + 0.5 * cinf * B
    ! end subroutine lambda_bc

end module bath
