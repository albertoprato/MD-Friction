!===========================================================================
!                           MINIMIZATION MODULE
!
!  Performs energy minimization of the system using the Conjugate Gradient
!  method (Polak-Ribiere variant). This is used to relax the initial solvent
!  configuration and remove high-energy overlaps before dynamics begin.
!
!  INPUTS:
!   - n_solv: Number of solvent particles.
!   - pos_solv: Initial positions of solvent particles.
!   - pos_solute: Positions of solute particles.
!   - epsilon_ss, sigma_ss: LJ parameters for solvent-solvent interaction.
!   - epsilon_int, sigma_int: LJ parameters for solute-solvent interaction.
!   - box_L: Length of the cubic box.
!   - tol: Convergence tolerance for the gradient.
!   - max_iter: Maximum number of minimization steps.
!
!  OUTPUT:
!   - pos_solv: Updated (relaxed) positions of solvent particles.
!===========================================================================


MODULE minimization_module
  USE kinds, ONLY: wp => dp
  USE force_module
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE conjugate_gradient_minimize(n_solv, pos_solv, pos_solute, &
                                         epsilon_ss, sigma_ss, epsilon_int, sigma_int, &
                                         box_L, tol, max_iter)
  
    ! Inputs
    INTEGER, INTENT(IN) :: n_solv, max_iter
    REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) :: pos_solute
    REAL(KIND=wp), INTENT(IN) :: epsilon_ss, sigma_ss, epsilon_int, sigma_int
    REAL(KIND=wp), INTENT(IN) :: box_L
    REAL(KIND=wp), INTENT(IN) :: tol
    REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) :: pos_solv    ! Input and Output

    ! Local Variables
    REAL(KIND=wp), DIMENSION(n_solv, 3) :: forces_solv, forces_solute_dummy
    REAL(KIND=wp), DIMENSION(n_solv, 3) :: g, g_old, d, pos_new
    REAL(KIND=wp) :: epot, epot_new
    REAL(KIND=wp) :: gamma, gamma_old, beta, alpha
    REAL(KIND=wp) :: g_dot_g, g_diff_dot_g
    INTEGER :: k, iter, i
    LOGICAL :: converged

    ! Clipping
    REAL(KIND=wp) :: disp_vec(3), disp_norm, scale_factor
    REAL(KIND=wp), PARAMETER :: max_disp = 0.1_wp     ! Max Displacement per step

    ! Initial calculation of forces and potential energy
    CALL force_calculation(n_solv, pos_solv, pos_solute, forces_solv, forces_solute_dummy, &
                           epsilon_ss, sigma_ss, epsilon_int, sigma_int, box_L, epot)
    
    ! Initial gradient
    g = -forces_solv
    
    ! Initial descent direction (Steepest Descent)
    d = -g  
    
    g_dot_g = SUM(g * g)
    gamma = g_dot_g
    
    iter = 0
    converged = .FALSE.

    ! Optimization Loop
    DO WHILE (iter < max_iter .AND. .NOT. converged)
      iter = iter + 1

      ! Check convergence on the gradient norm
      IF (SQRT(g_dot_g) < tol) THEN
         converged = .TRUE.
         EXIT
      END IF

      ! Line Search to find step size alpha 
      alpha = 1.0e-4_wp    ! Small initial step

      ! Calculation of new positions with displacement limit
      DO i = 1, n_solv
        ! Displacement vector proposed: alpha * direction
        disp_vec(:) = alpha * d(i, :)
          
        ! Calculate displacement lenght 
        disp_norm = SQRT(SUM(disp_vec**2))
          
        ! If the displacement is too large, scale it down
        IF (disp_norm > max_disp) THEN
          scale_factor = max_disp / disp_norm
          disp_vec(:) = disp_vec(:) * scale_factor
        END IF
          
        pos_new(i, :) = pos_solv(i, :) + disp_vec(:)
      END DO
 
      CALL force_calculation(n_solv, pos_new, pos_solute, forces_solv, forces_solute_dummy, &
                              epsilon_ss, sigma_ss, epsilon_int, sigma_int, box_L, epot_new)
       
      ! Backtracking: decrease alpha if energy increases
      DO k = 1, 10
        IF (epot_new < epot) EXIT
    
        alpha = alpha * 0.5_wp

        DO i = 1, n_solv
          disp_vec(:) = alpha * d(i, :)
          disp_norm = SQRT(SUM(disp_vec**2))
          IF (disp_norm > max_disp) THEN
            scale_factor = max_disp / disp_norm
            disp_vec(:) = disp_vec(:) * scale_factor
          END IF
          pos_new(i, :) = pos_solv(i, :) + disp_vec(:)
        END DO
        
        CALL force_calculation(n_solv, pos_new, pos_solute, forces_solv, forces_solute_dummy, &
                               epsilon_ss, sigma_ss, epsilon_int, sigma_int, box_L, epot_new)
      END DO

      ! Update positions and energy
      pos_solv = pos_new
      epot = epot_new
      g_old = g              ! Save the old gradient 
      gamma_old = gamma      ! Save the old dot product

      ! Calculate new gradient
      g = -forces_solv 
      g_dot_g = SUM(g * g)
      gamma = g_dot_g

      ! Polak-Ribiere formula
      g_diff_dot_g = SUM(g * (g - g_old))
      beta = g_diff_dot_g / gamma_old
       
      ! Reset to Steepest Descent if beta < 0
      IF (beta < 0.0_wp) beta = 0.0_wp

      ! Compute new conjugate direction
      d = -g + beta * d
      
    END DO
    
  END SUBROUTINE conjugate_gradient_minimize

END MODULE minimization_module
