!=========================================================
! Force Calculation Module
! Computes interaction forces and potential energy between
! solvent-solvent and solute-solvent particles using 
! Lennard-Jones potentials.
!=========================================================

MODULE force_module                                                       
  USE kinds, ONLY: wp => dp                                           
  IMPLICIT NONE                                                                                                  
                                                                      
  CONTAINS                                                           
                                                                      
    SUBROUTINE force_calculation(n_solv, pos_solv, pos_solute, forces_solv, forces_solute, &
                                 epsilon_ss, sigma_ss, epsilon_int, sigma_int, box_L, epot)
      ! Inputs
      INTEGER, INTENT(IN) :: n_solv
      REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: pos_solv ! (N,3)
      REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: pos_solute ! (4,3)
      
      ! Lennard-Jones Parameters
      REAL (KIND=wp), INTENT(IN) :: epsilon_ss   ! epsilon solvent-solvent
      REAL (KIND=wp), INTENT(IN) :: sigma_ss     ! sigma solvent-solvent
      REAL (KIND=wp), INTENT(IN) :: epsilon_int  ! epsilon solute-solvent
      REAL (KIND=wp), INTENT(IN) :: sigma_int    ! sigma solute-solvent
      REAL (KIND=wp), INTENT(IN) :: box_L        ! Cubic Box Length
      
      ! Outputs
      REAL (KIND=wp), DIMENSION(:,:), INTENT(OUT) ::  forces_solv    ! (N,3)
      REAL (KIND=wp), DIMENSION(:,:), INTENT(OUT) ::  forces_solute  ! (4,3)
      REAL (KIND=wp), INTENT(OUT) :: epot
      
      ! Local Variables
      INTEGER :: i, k, dim, I_solute  
      REAL (KIND=wp) :: dist_sq, dist, r_inv, r2_inv, r6_inv, r12_inv
      REAL (KIND=wp) :: diff(3)
      REAL (KIND=wp) :: sigma6, sigma12
      REAL (KIND=wp) :: e_lj, f_lj_factor
      REAL (KIND=wp) :: inv_box, r_cut_sq, half_box 
      REAL(KIND=wp), PARAMETER :: r_cut = 12.0_wp   ! Cutoff

      forces_solv = 0.0_wp
      forces_solute = 0.0_wp
      epot = 0.0_wp
      
      inv_box = 1.0_wp / box_L
      half_box = box_L * 0.5_wp
      
      ! We use the physical cutoff unless the box is too small.
      ! In that case, we use L/2
      IF (half_box < r_cut) THEN
        r_cut_sq = half_box**2
      ELSE
        r_cut_sq = r_cut**2
      END IF
 
      ! ===========================     
      ! Solvent-Solvent Interaction
      ! ===========================
      
      sigma6 = sigma_ss**6
      sigma12 = sigma6**2
      
      DO i = 1, n_solv - 1
        DO k = i + 1, n_solv
          
          dist_sq = 0.0_wp
          DO dim = 1, 3
            diff(dim) = pos_solv(i, dim) - pos_solv(k, dim)
            
            ! Minimum Image Convention
            diff(dim) = diff(dim) - box_L * ANINT(diff(dim) * inv_box)
            
            dist_sq = dist_sq + diff(dim)**2
          END DO
          
          IF (dist_sq > 1.0e-10_wp .AND. dist_sq < r_cut_sq) THEN
            
            r_inv = 1.0_wp / dist
            r2_inv = 1.0_wp / dist_sq
            r6_inv = r2_inv**3
            r12_inv = r6_inv**2
            
            ! Potential energy
            e_lj = 4.0_wp * epsilon_ss * (sigma12 * r12_inv - sigma6 * r6_inv)
            epot = epot + e_lj
            
            ! Force 
            f_lj_factor = 24.0_wp * epsilon_ss * r2_inv * (2.0_wp * sigma12 * r12_inv - sigma6 * r6_inv)
            
            DO dim = 1, 3
              ! Force on the i-th particle
              forces_solv(i, dim) = forces_solv(i, dim) + f_lj_factor * diff(dim)
              ! Newton's law
              forces_solv(k, dim) = forces_solv(k, dim) - f_lj_factor * diff(dim)
            END DO
          ENDIF
        END DO
      END DO             
      
      ! ===========================
      ! Solute-Solvent Interaction
      ! ===========================
      
      sigma6 = sigma_int**6
      sigma12 = sigma6**2

      DO i = 1, n_solv
        DO I_solute = 1, 4
          
          dist_sq = 0.0_wp
          DO dim = 1, 3
            diff(dim) = pos_solv(i, dim) - pos_solute(I_solute, dim)
            
            diff(dim) = diff(dim) - box_L * ANINT(diff(dim) * inv_box)

            dist_sq = dist_sq + diff(dim)**2
          END DO
          
          IF (dist_sq > 1.0e-10_wp .AND. dist_sq < r_cut_sq) THEN
            
            r_inv = 1.0_wp / dist
            r2_inv = 1.0_wp / dist_sq
            r6_inv = r2_inv**3
            r12_inv = r6_inv**2
            
            e_lj = 4.0_wp * epsilon_int * (sigma12 * r12_inv - sigma6 * r6_inv)
            epot = epot + e_lj
            
            f_lj_factor = 24.0_wp * epsilon_int * r2_inv * (2.0_wp * sigma12 * r12_inv - sigma6 * r6_inv)
            
            DO dim = 1, 3
              ! Force on the i-th bath particle
              forces_solv(i, dim) = forces_solv(i, dim) + f_lj_factor * diff(dim)
              
              ! Force on the I-th solute
              forces_solute(I_solute, dim) = forces_solute(I_solute, dim) - f_lj_factor * diff(dim)
            END DO
          ENDIF
        END DO
      END DO          
    END SUBROUTINE force_calculation     
                                                                      
END MODULE force_module
