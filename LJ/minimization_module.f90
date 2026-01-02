MODULE minimization_module
  USE kinds, ONLY: wp => dp
  USE force_module
  IMPLICIT NONE

  CONTAINS

  !============================================================
  ! Conjugate Gradient Minimization Routine (Polak-Ribiere)
  ! Performs energy minimization to relax the solvent structure
  !============================================================

  SUBROUTINE conjugate_gradient_minimize(n_solv, pos_solv, pos_solute, &
                                         epsilon_ss, sigma_ss, epsilon_int, sigma_int, &
                                         tol, max_iter)
  
    ! Inputs
    INTEGER, INTENT(IN) :: n_solv, max_iter
    REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) :: pos_solv
    REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)    :: pos_solute
    REAL(KIND=wp), INTENT(IN) :: epsilon_ss, sigma_ss, epsilon_int, sigma_int
    REAL(KIND=wp), INTENT(IN) :: tol

    ! Local Variables
    REAL(KIND=wp), DIMENSION(n_solv, 3) :: forces_solv, forces_solute_dummy
    REAL(KIND=wp), DIMENSION(n_solv, 3) :: g, g_old, d, pos_new
    REAL(KIND=wp) :: epot, epot_new
    REAL(KIND=wp) :: gamma, gamma_old, beta, alpha
    REAL(KIND=wp) :: g_dot_g, g_diff_dot_g
    INTEGER :: k, iter
    LOGICAL :: converged

    ! Initial calculation of forces and potential energy
    CALL force_calculation(n_solv, pos_solv, pos_solute, forces_solv, forces_solute_dummy, &
                           epsilon_ss, sigma_ss, epsilon_int, sigma_int, epot)
    
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
      ! We verify that E(pos + alpha*d) < E(pos)  
      alpha = 1.0e-4_wp ! Small initial step to ensure stability
      pos_new = pos_solv + alpha * d
       
      CALL force_calculation(n_solv, pos_new, pos_solute, forces_solv, forces_solute_dummy, &
                              epsilon_ss, sigma_ss, epsilon_int, sigma_int, epot_new)
       
      ! Backtracking: decrease alpha if energy increases
      DO k = 1, 10
        IF (epot_new < epot) EXIT
        alpha = alpha * 0.5_wp
        pos_new = pos_solv + alpha * d
        CALL force_calculation(n_solv, pos_new, pos_solute, forces_solv, forces_solute_dummy, &
                               epsilon_ss, sigma_ss, epsilon_int, sigma_int, epot_new)
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
