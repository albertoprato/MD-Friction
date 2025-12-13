PROGRAM main
  USE kinds, ONLY: wp => dp                                           
  USE force_module
  USE friction_module
                                                        
  IMPLICIT NONE

  ! Input Variables
  
  INTEGER :: nk, n_solv, i, step
  REAL (KIND=wp) :: dt
  REAL (KIND=wp) :: mass_solv
  REAL (KIND=wp) :: mu_B, omega, d0
  REAL (KIND=wp) :: mu_int, Omega_big, D0_big

  ! Array

  REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: pos_solv, vel_solv, force      ! (N, 3)
  REAL (KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: traj_history                 ! (N, 3, Steps)
  REAL (KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: force_history                ! (4, 3, Steps)
  REAL(KIND=wp), DIMENSION(4, 3) :: pos_solute, force_solute 
  
  REAL(KIND=wp) :: e_pot
  
  REAL(KIND=wp) :: inp_x, inp_y, inp_z, inp_vx, inp_vy, inp_vz

  ! Read Input

  OPEN(UNIT=10, FILE='input.txt', STATUS='old')
  
  READ(10, *) nk, dt
  READ(10, *) mass_solv, mu_B, omega, d0
  READ(10, *) mu_int, Omega_big, D0_big
 
  DO i = 1, 4
    READ(10, *) pos_solute(i, 1), pos_solute(i, 2), pos_solute(i, 3)
  END DO
  
  READ(10, *) n_solv

  ALLOCATE(pos_solv(n_solv, 3), vel_solv(n_solv, 3), force(n_solv, 3))
  ALLOCATE(traj_history(n_solv, 3, nk))
  ALLOCATE(force_history(4, 3, nk))

  DO i = 1, n_solv
    READ(10, *) inp_x, inp_y, inp_z, inp_vx, inp_vy, inp_vz
    pos_solv(i, 1) = inp_x
    pos_solv(i, 2) = inp_y
    pos_solv(i, 3) = inp_z
    vel_solv(i, 1) = inp_vx
    vel_solv(i, 2) = inp_vy
    vel_solv(i, 3) = inp_vz
  END DO
  CLOSE(10)
  
  PRINT *, "System initialized with", n_solv, "solvent particles and", 4, "solute particles"
  PRINT *, "Starting Molecular Dynamics..."

  ! Calculate the forces at t = 0
  CALL force_calculation(n_solv, pos_solv, pos_solute, force, force_solute, mu_B, omega, d0, mu_int, Omega_big, D0_big, e_pot) 

  ! Output Files
  OPEN(UNIT=20, FILE='trajectory.xyz', STATUS='replace')

  ! Verlet Algorithm
  
  DO step = 1, nk

    ! Half-kick
    vel_solv = vel_solv + 0.5_wp * dt * (force / mass_solv)

    ! Drift
    pos_solv = pos_solv + dt * vel_solv

    ! Save Trajectory
    traj_history(:, :, step) = pos_solv

    ! Calculate new forces 
    CALL force_calculation(n_solv, pos_solv, pos_solute, force, force_solute, mu_B, omega, d0, mu_int, Omega_big, D0_big, e_pot)
   
    ! Save Forces
    force_history(:, :, step) = force_solute

    ! Half-kick
    vel_solv = vel_solv + 0.5_wp * dt * (force / mass_solv)
  
    ! Output Results (e.g. every 100 steps)
    IF (MOD(step, 100) == 0) THEN
      
      ! Trajectory
      WRITE(20, *) n_solv + 4
      WRITE(20, '(A, I8)') "Step: ", step
      
      ! 4 solute particles
      DO i = 1, 4
        WRITE(20, '(A, 3F15.8)') "C ", pos_solute(i, :)
      END DO

      ! Solvent particles
      DO i = 1, n_solv
        WRITE(20, '(A, 3F15.8)') "O ", pos_solv(i, :)
      END DO
 
    END IF


  END DO

  CLOSE(20)
  CLOSE(21)
  
  PRINT *, "Simulation Completed."


  ! Calculate Friction Tensor
  CALL friction_tensor(nk, force_history, dt)

  DEALLOCATE(pos_solv, vel_solv, force, traj_history, force_history)

END PROGRAM main
