!==============================================================================
!                           FRICTION TENSOR MODULE
!
!  Calculates the friction tensor for the solute particles. 
!  It processes the history of forces experienced by the solute during the 
!  simulation and computes their time-autocorrelation functions via FFT method.
!
!  INPUTS:
!   - n_steps: Total number of simulation steps.
!   - force_hist: 3D array containing force history for the 4 solute particles.
!   - dt: Time step of the simulation.
!
!  OUTPUT:
!   - friction_tensor.dat: File containing the time evolution of the friction
!                          tensor components
!==============================================================================

MODULE friction_module
  USE kinds, ONLY: wp => dp
  USE fft_correlation_module
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE friction_tensor(n_steps, force_hist, dt)

    INTEGER, INTENT(IN) :: n_steps
    REAL (KIND=wp), DIMENSION(4, 3, n_steps), INTENT(IN) :: force_hist
    REAL (KIND=wp), INTENT(IN) :: dt

    ! Local Variables
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: vec_i, vec_j
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: corr_result
    REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: final_results

    INTEGER :: a, I_sol, J_sol, t_lag, step, idx
    INTEGER(c_int) :: N_c
    REAL (KIND=wp) :: time_val

    ! Allocate 1D array for FFT
    ALLOCATE(vec_i(n_steps), vec_j(n_steps))
    ALLOCATE(corr_result(n_steps))

    N_c = INT(n_steps, KIND(c_int))

    ! Output Files
    OPEN(UNIT=30, FILE='friction_tensor.dat', status='replace')

    WRITE(30, '(A)', ADVANCE='NO') "#    Time      "
    DO I_sol = 1, 4
      WRITE(30, '(3(A,I1), A)', ADVANCE='NO') &
             "  Sol ", I_sol, "_X    Sol", I_sol, "_Y    Sol", I_sol, "_Z "
    END DO
    WRITE(30, *) ""
          
    ! Create an array to store the diagonal components      
    ALLOCATE(final_results(n_steps, 12))
  
    J_sol = 1    
    idx = 0

    ! Loop on solute particle
    DO I_sol = 1, 4
      ! Loop on x,y,z
      DO a = 1, 3
        idx = idx + 1       

        ! Extract component 'a' of the force on particle I
        vec_i = force_hist(I_sol, a, :)

        ! Calculate autocorrelation
        CALL compute_correlation_fft(N_c, vec_i, vec_i, corr_result)

        ! Save results
        final_results(:, idx) = corr_result(:)
    
      END DO
    END DO

    ! We report the first half only due to FFT symmetry
    DO step = 1, n_steps / 2
      t_lag = step - 1
      time_val = t_lag * dt
      
      WRITE(30, '(F10.4, 2X)', ADVANCE='NO') time_val
      DO idx = 1, 12
        WRITE(30, '(ES14.6, 1X)', ADVANCE='NO') final_results(step, idx)
      END DO
      WRITE(30, *) ""
    END DO

    CLOSE(30)

    DEALLOCATE(vec_i, vec_j, corr_result, final_results)

  END SUBROUTINE friction_tensor

END MODULE friction_module
