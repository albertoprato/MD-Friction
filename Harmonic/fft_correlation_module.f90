MODULE fft_correlation_module
  USE kinds, ONLY: wp => dp
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT  NONE
  
  INCLUDE 'fftw3.f03'

CONTAINS

  SUBROUTINE compute_correlation_fft(N, x, y, correlation)

    !==============================================
    ! Calculate the cross-correlation throught FFT
    ! If x == y, autocorrelation
    !==============================================
 
    INTEGER(C_INT), INTENT(IN) :: N
    REAL(KIND=wp), INTENT(IN) :: x(N), y(N)
    REAL(KIND=wp), INTENT(OUT) :: correlation(N)

    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: signal_x(:), signal_y(:)
    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: fft_x(:), fft_y(:)
    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: power_spectrum(:), result_complex(:) 
    
    TYPE(C_PTR) :: plan_fwd_x, plan_fwd_y, plan_bwd

    REAL(KIND=wp) :: mean_x, mean_y, norm
    INTEGER :: i

    ! Allocation
    ALLOCATE(signal_x(N), fft_x(N))
    ALLOCATE(signal_y(N), fft_y(N))
    ALLOCATE(power_spectrum(N), result_complex(N))
    
    ! Subtraction of the mean to obtain fluctuations
    mean_x = SUM(x) / DBLE(N)
    mean_y = SUM(y) / DBLE(N)

    ! Copy real values into complex vector
    DO i = 1, N
      signal_x(i) = CMPLX(x(i)- mean_x, 0.0_wp, KIND=C_DOUBLE_COMPLEX) 
      signal_y(i) = CMPLX(y(i)- mean_y, 0.0_wp, KIND=C_DOUBLE_COMPLEX)
    END DO

    ! Creation Plans FFTW
    plan_fwd_x = fftw_plan_dft_1d(N, signal_x, fft_x, FFTW_FORWARD, FFTW_ESTIMATE)
    plan_fwd_y = fftw_plan_dft_1d(N, signal_y, fft_y, FFTW_FORWARD, FFTW_ESTIMATE)
    plan_bwd   = fftw_plan_dft_1d(N, power_spectrum, result_complex, FFTW_BACKWARD, FFTW_ESTIMATE)
    
    ! FFT
    CALL fftw_execute_dft(plan_fwd_x, signal_x, fft_x)
    CALL fftw_execute_dft(plan_fwd_y, signal_y, fft_y)

    ! Calculate conj(F(f)) * F(g)
    DO i = 1, N
        power_spectrum(i) = CONJG(fft_x(i)) * fft_y(i)
    END DO

    ! IFFT
    CALL fftw_execute_dft(plan_bwd, power_spectrum, result_complex)

    norm = 1.0_wp / (DBLE(N) * DBLE(N))

    DO i = 1, N
       correlation(i) = REAL(result_complex(i), KIND=wp) * norm
    END DO

    ! Clean
    CALL fftw_destroy_plan(plan_fwd_x)
    CALL fftw_destroy_plan(plan_fwd_y)
    CALL fftw_destroy_plan(plan_bwd)
    DEALLOCATE(signal_x, fft_x, signal_y, fft_y, power_spectrum, result_complex)

  END SUBROUTINE compute_correlation_fft

END MODULE fft_correlation_module
