module integrator

    use constants
    use parameters
    
    implicit none
    
    
    
    
    complex function bath_corr(temp, gamma)
    
        real(kind=DP), intent(in) :: temp, gamma
    
        ! Matsubara frequency (n=1)
        gj(1) = gamma
        gj(2) = (2.0_DP*pi)*temp
        cotan = cos(0.5_DP*gj(1)/temp)/sin(0.5*gj(1)/temp)
        delta = 2.0_DP/(gj(1)/temp) - cotan - 4.0_DP*gj(1)*temp/(gj(2)**2-gj(1)**2)
        cj(1) = gj(1) * cotan
        cj(2) = 4.0_DP*gj(1)*temp * gj(2)/(gj(2)**2-gj(1)**2)
        
        bath_corr = cj(1) * REAL1 + cj(2) * IMAG1
    end function bath_corr
    
    subroutine rk4(n_steps)
    
        integer, intent(in) :: n_steps
        
        do n = 1, n_steps











        




        end do
        
    end subroutine rk4

end module integrator
