PROGRAM main
  USE kinds, ONLY: wp => dp                                           
  USE force_module
  USE friction_module
                                                        
  IMPLICIT NONE

  ! Input Variables
  
  INTEGER :: nk, n_solv, i, step, n_total_atoms, atom_idx
  REAL (KIND=wp) :: dt
  REAL (KIND=wp) :: mass_solv
  REAL (KIND=wp) :: mu_B, omega, d0
  REAL (KIND=wp) :: mu_int, Omega_big, D0_big
  REAL (KIND=wp) :: r_cut

  ! Array

  REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: pos_solv, vel_solv, force      ! (N, 3)
  REAL (KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: traj_history                 ! (N, 3, Steps)
  REAL (KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: force_history                ! (4, 3, Steps)
  REAL(KIND=wp), DIMENSION(4, 3) :: pos_solute, force_solute 
  
  REAL(KIND=wp) :: e_pot
  CHARACTER(LEN=2) :: atom_label  
  REAL(KIND=wp) :: inp_x, inp_y, inp_z

  ! Read Input

  OPEN(UNIT=10, FILE='input.txt', STATUS='old')
  
  READ(10, *) nk, dt
  READ(10, *) mass_solv, mu_B, omega, d0, r_cut
  READ(10, *) mu_int, Omega_big, D0_big
 
  CLOSE(10)

  OPEN(UNIT=11, FILE='system.xyz', STATUS='old')

  READ(11, *) n_total_atoms
  READ(11, *)
  
  n_solv = n_total_atoms - 4

  PRINT *, "System initialized with", n_solv, "solvent particles and", 4, "solute particles"
 
  ALLOCATE(pos_solv(n_solv, 3), vel_solv(n_solv, 3), force(n_solv, 3))
  ALLOCATE(traj_history(n_solv, 3, nk))
  ALLOCATE(force_history(4, 3, nk))
  
  DO i = 1, n_total_atoms
    READ(11, *) atom_label, inp_x, inp_y, inp_z
       
    IF (i <= 4) THEN
      ! Solute
      pos_solute(i, 1) = inp_x
      pos_solute(i, 2) = inp_y
      pos_solute(i, 3) = inp_z
    ELSE
      ! Solvent
      atom_idx = i - 4
      pos_solv(atom_idx, 1) = inp_x
      pos_solv(atom_idx, 2) = inp_y
      pos_solv(atom_idx, 3) = inp_z
    END IF
  END DO
  CLOSE(11)
    
  PRINT *, "Starting Molecular Dynamics..."

  vel_solv = 0.0_wp

  ! Calculate the forces at t = 0
  CALL force_calculation(n_solv, pos_solv, pos_solute, force, force_solute, mu_B, omega, d0, &
                         mu_int, Omega_big, D0_big, r_cut, e_pot) 

  ! Output Files
  OPEN(UNIT=20, FILE='trajectory.xyz', STATUS='replace')
  WRITE(20, *) n_total_atoms
  WRITE(20, '(A)') "Initial Configuration"
  DO i = 1, 4
    WRITE(20, '(A, 3F15.8)') "C ", pos_solute(i, :)
  END DO
  DO i = 1, n_solv
    WRITE(20, '(A, 3F15.8)') "O ", pos_solv(i, :)
  END DO

  ! Verlet Algorithm
  
  DO step = 1, nk

    ! Half-kick
    vel_solv = vel_solv + 0.5_wp * dt * (force / mass_solv)

    ! Drift
    pos_solv = pos_solv + dt * vel_solv

    ! Save Trajectory
    traj_history(:, :, step) = pos_solv

    ! Calculate new forces 
    CALL force_calculation(n_solv, pos_solv, pos_solute, force, force_solute, mu_B, omega, d0, &
                           mu_int, Omega_big, D0_big, r_cut, e_pot)
   
    ! Save Forces
    force_history(:, :, step) = force_solute

    ! Half-kick
    vel_solv = vel_solv + 0.5_wp * dt * (force / mass_solv)
  
    ! Output Results (e.g. every 100 steps)
    IF (MOD(step, 50) == 0) THEN
      
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
