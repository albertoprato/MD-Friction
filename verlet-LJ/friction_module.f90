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

  SUBROUTINE friction_tensor(n_steps, force_hist, dt, temp, kb)

    INTEGER, INTENT(IN) :: n_steps
    REAL (KIND=wp), DIMENSION(4, 3, n_steps), INTENT(IN) :: force_hist
    REAL (KIND=wp), INTENT(IN) :: dt

    ! Local Variables
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: vec_i, vec_j
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: corr_result
    REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: final_results
    
    INTEGER :: i_idx, j_idx
    INTEGER :: a, b, I_sol, J_sol, t_lag, step, idx
    INTEGER(c_int) :: N_c
    REAL (KIND=wp) :: time_val
    REAL (KIND=wp) :: temp, kb    

    ! Allocate 1D array for FFT
    ALLOCATE(vec_i(n_steps), vec_j(n_steps))
    ALLOCATE(corr_result(n_steps))

    N_c = INT(n_steps, KIND(c_int))

    ! Output Files
    OPEN(UNIT=30, FILE='friction_tensor.dat', status='replace')
    
    WRITE(30, '(A)', ADVANCE='NO') "#    Time      "
    DO I_sol = 1, 4
      DO a = 1, 3
        ! Calculate the 1-12 index for the first dimension
        i_idx = (I_sol - 1) * 3 + a
        
        DO J_sol = 1, 4
          DO b = 1, 3
            ! Calculate the 1-12 index for the second dimension
            j_idx = (J_sol - 1) * 3 + b
            
            ! Write the pair (i,j) using I0 format to avoid extra spaces
            WRITE(30, '(" (", I0, ",", I0, ") ")', ADVANCE='NO') i_idx, j_idx
            
          END DO
        END DO
      END DO
    END DO
    WRITE(30, *) ""

    ! Create an array to store the diagonal components      
    ALLOCATE(final_results(n_steps, 144))
  
    idx = 0

    ! Loop on solute particle
    DO I_sol = 1, 4
      ! Loop on x,y,z
      DO a = 1, 3
        ! Extract component 'a' of the force on particle I
        vec_i = force_hist(I_sol, a, :)

        ! Loop on solute particle
        DO J_sol = 1, 4
          ! Loop on x,y,z
          DO b = 1, 3
          
            idx = idx + 1
          
            ! Extract component 'a' of the force on particle I
            vec_j = force_hist(J_sol, b, :) 
        
            ! Calculate autocorrelation
            CALL compute_correlation_fft(N_c, vec_i, vec_j, corr_result)

            ! Save results
            final_results(:, idx) = corr_result(:) / (temp * kb)
            
          END DO
        END DO
      END DO
    END DO

    ! We report the first half only due to FFT symmetry
    DO step = 1, n_steps / 2
      t_lag = step - 1
      time_val = t_lag * dt
      
      WRITE(30, '(F10.4, 2X)', ADVANCE='NO') time_val
      DO idx = 1, 144
        WRITE(30, '(ES14.6, 1X)', ADVANCE='NO') final_results(step, idx)
      END DO
      WRITE(30, *) ""
    END DO

    CLOSE(30)

    DEALLOCATE(vec_i, vec_j, corr_result, final_results)

  END SUBROUTINE friction_tensor

END MODULE friction_module
