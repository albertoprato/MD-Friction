MODULE force_module                                                       
  USE kinds, ONLY: wp => dp                                           
  IMPLICIT NONE                                                                                                  
                                                                      
  CONTAINS                                                           
                                                                      
    SUBROUTINE force_calculation(n_solv, pos_solv, pos_solute, forces_solv, forces_solute, &
                                 mu_B, omega, d0, mu_int, Omega_big, D0_big, epot)

      ! Inputs

      INTEGER, INTENT(IN) :: n_solv
      REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: pos_solv ! (N,3)
      REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: pos_solute ! (4,3)
      
      REAL (KIND=wp), INTENT(IN) :: mu_B, omega, d0
      REAL (KIND=wp), INTENT(IN) :: mu_int, Omega_big, D0_big
      
      ! Outputs

      REAL (KIND=wp), DIMENSION(:,:), INTENT(OUT) ::  forces_solv    ! (N,3)
      REAL (KIND=wp), DIMENSION(:,:), INTENT(OUT) ::  forces_solute  ! (4,3)
      REAL (KIND=wp), INTENT(OUT) :: epot

      ! Local Variables

      INTEGER :: i, k, dim, I_solute  
      REAL (KIND=wp) :: dist_sq, dist
      REAL (KIND=wp) :: diff(3)
      REAL (KIND=wp) :: force_factor, const_B, const_int

      ! Force Constants
      const_B   = 2.0_wp * mu_B * (omega**2)
      const_int = 2.0_wp * mu_int * (Omega_big**2)

      ! Initialize outputs
      forces_solv = 0.0_wp
      forces_solute = 0.0_wp
      epot = 0.0_wp
      

      ! ===========================     
      ! Solvent-Solvent Interaction
      ! ===========================


      DO i = 1, n_solv - 1
        DO k = i + 1, n_solv

          ! Calculate the Distance Vector r_i - r_k
      
          dist_sq = 0.0_wp
          DO dim = 1, 3
            diff(dim) = pos_solv(i, dim) - pos_solv(k, dim)
            dist_sq = dist_sq + diff(dim)**2
          END DO

          dist = sqrt(dist_sq)

          IF (dist > 1.0e-10_wp) THEN ! Avoid division by zero

            epot = epot + mu_B * (omega**2) * ((dist - d0)**2)
      
            force_factor = -const_B * (1.0_wp - (d0 / dist))
        
            DO dim=1, 3

              ! Force on i-th particle
              forces_solv(i, dim) = forces_solv(i, dim) + (force_factor * diff(dim))

              ! Newton's 3rd law
              forces_solv(k, dim) = forces_solv(k, dim) - (force_factor * diff(dim))
            END DO
          ENDIF
        END DO
      END DO             


      ! ===========================
      ! Solute-Solvent Interaction
      ! ===========================

      DO i = 1, n_solv
        DO I_solute = 1, 4

          ! Calculate the Distance Vector r_i - R_k
          dist_sq = 0.0_wp
          DO dim = 1, 3
            diff(dim) = pos_solv(i, dim) - pos_solute(I_solute, dim)
            dist_sq = dist_sq + diff(dim)**2
          END DO

          dist = sqrt(dist_sq)

          IF (dist > 1.0e-10_wp) THEN ! Avoid division by zero

            epot = epot + mu_int * (Omega_big**2) * ((dist - D0_big)**2)

            force_factor = -const_int * (1.0_wp - (D0_big / dist))

            DO dim=1, 3

              ! Force on i-th bath particle
              forces_solv(i, dim) = forces_solv(i, dim) + (force_factor * diff(dim))
              
              ! Force on I-th solute particle
              forces_solute(I_solute, dim) = forces_solute(I_solute, dim) - (force_factor * diff(dim))
            END DO
          ENDIF
        END DO
      END DO          

    END SUBROUTINE force_calculation     
                                                                      
END MODULE force_module
