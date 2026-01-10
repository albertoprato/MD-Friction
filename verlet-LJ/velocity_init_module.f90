!===============================================================
!                VELOCITY INITIALIZATION MODULE
!
! Generates initial velocities for solvent particles according
! to the Maxwell-Boltzmann distribution at a target Temperature.
!
! INPUTS:
!   - n_solv: Number of solvent particles.
!   - mass: Mass of a solvent particle.
!   - temp: Target temperature.
!   - kb: Boltzmann constant.
!
! OUTPUT:
!   - vel_solv: Array of initialized particle velocities scaled 
!               to the target temperature.
!===============================================================

MODULE velocity_init_module
  USE kinds, ONLY: wp => dp
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE init_velocities(n_solv, vel_solv, mass, temp, kb)
    
    ! Inputs
    INTEGER, INTENT(IN) :: n_solv
    REAL(KIND=wp), INTENT(IN) :: mass  ! Solvent mass
    REAL(KIND=wp), INTENT(IN) :: temp  ! Temperature
    REAL(KIND=wp), INTENT(IN) :: kb    ! Boltzmann Constant

    ! Output
    REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT) :: vel_solv ! (N, 3)

    ! Local Variables
    INTEGER :: i, k
    REAL(KIND=wp) :: sigma_v, r1, r2, v_cm(3), kin_en, computed_temp, scale_f
    REAL(KIND=wp) :: pi

    pi = 4.0_wp * ATAN(1.0_wp)

    ! Calculate standard deviation for Gaussian distribution
    sigma_v = SQRT((kb * temp) / mass)

    ! Seed the random number generator
    CALL RANDOM_SEED()

    ! Assign random velocities (Box-Muller transform)
    DO i = 1, n_solv
      DO k = 1, 3
        CALL RANDOM_NUMBER(r1)    ! Generates a uniform number between 0 and 1
        CALL RANDOM_NUMBER(r2)    ! Generates a uniform number between 0 and 1
          
        ! Gaussian number with mean 0 and std dev sigma_v
        vel_solv(i, k) = sigma_v * SQRT(-2.0_wp * LOG(r1)) * COS(2.0_wp * pi * r2)
      END DO
    END DO

    ! Remove Center of Mass Drift
    v_cm = 0.0_wp
    DO i = 1, n_solv
      v_cm(:) = v_cm(:) + vel_solv(i, :)     ! Vector sum of all velocities
    END DO
    v_cm = v_cm / DBLE(n_solv)    ! Average system velocity

    DO i = 1, n_solv
      vel_solv(i, :) = vel_solv(i, :) - v_cm(:)    ! Subtraction 
    END DO

    ! Scale to exact Temperature
    kin_en = 0.0_wp
    DO i = 1, n_solv
      kin_en = kin_en + SUM(vel_solv(i, :)**2)
    END DO
    kin_en = 0.5_wp * mass * kin_en

    computed_temp = (2.0_wp * kin_en) / (3.0_wp * DBLE(n_solv - 1) * kb)
    
    scale_f = SQRT(temp / computed_temp)
    vel_solv = vel_solv * scale_f

  END SUBROUTINE init_velocities

END MODULE velocity_init_module
