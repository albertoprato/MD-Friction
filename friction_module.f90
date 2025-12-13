MODULE friction_module
  USE kinds, ONLY: wp => dp
  USE fft_correlation_module
  IMPLICIT NONE

  CONTAINS

  !================
  ! Friction Tensor
  !================

  subroutine friction_tensor(n_steps, force_hist, dt)

    INTEGER, INTENT(IN) :: n_steps
    REAL (KIND=wp), DIMENSION(4, 3, n_steps), INTENT(IN) :: force_hist
    REAL (KIND=wp), INTENT(IN) :: dt

    ! Local Variables
    
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: vec_i, vec_j
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: corr_result
    REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: final_results
    REAL (KIND=wp), DIMENSION(3,3) :: Friction

    INTEGER :: a, b, I_sol, t_lag, step
    INTEGER(c_int) :: N_c
    REAL (KIND=wp) :: time_val

    ! Allocate 1D array for FFT
    ALLOCATE(vec_i(n_steps), vec_j(n_steps))
    ALLOCATE(corr_result(n_steps))

    N_c = INT(n_steps, KIND(c_int))

    ! Output Files
    OPEN(UNIT=30, FILE='friction_tensor.dat', status='replace')
    WRITE(30, '(A)') "#    Time  |      Friction_XX Friction_YY Friction_ZZ"
          
    PRINT *, "Calculating Friction Tensor throught FFTW3..."
    
    I_sol = 1
    
    ! Create an array to store the diagonal components      
    ALLOCATE(final_results(n_steps, 3))

    ! Loop on x,y,z
    DO a = 1, 3
      ! Extract component 'a' of the force
      vec_i = force_hist(I_sol, a, :)

      ! Calculate autocorrelation <Fa(0) Fa(t)>
      CALL compute_correlation_fft(N_c, vec_i, vec_i, corr_result)

      ! Save results
      final_results(:, a) = corr_result(:)
    
    END DO
    
    ! We only report the first half (N/2) because the cyclic FFT is symmetric and the correlation decays
    DO step = 1, n_steps / 2
      t_lag = step - 1
      time_val = t_lag * dt

      WRITE(30, '(F10.2, 2X, 3(ES14.6, 1X))') &
                  time_val, final_results(step, 1), final_results(step, 2), final_results(step, 3)
    END DO

    CLOSE(30)

    DEALLOCATE(vec_i, vec_j, corr_result, final_results)

    PRINT *, "Friction calculation completed. Results in friction_tensor.dat"

  END SUBROUTINE friction_tensor

END MODULE friction_module
