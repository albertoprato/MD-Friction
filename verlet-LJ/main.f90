!========================================================================================================
!                           PROGRAM: CLASSICAL MOLECULAR DYNAMICS SIMULATION
!  
!  This is the main driver of the MD simulation. It handles:
!
!   1. Reading the simulation parameters (input.txt) and the initial configuration (system.xyz).
!
!   2. Optimization of the initial configuration via Conjugate Gradient method.
!
!   3. Initialization of the solvent molecules' velocities, given the simulation temperature.
!
!   4. An initial equilibration phase in which the system transitions from a solid crystalline 
!      structure to a liquide one.
!
!   5. Simulation of the system dynamics by time integration using the Verlet algorithm.
!
!   6. Storing the history of the forces and computing the friction tensor the friction tensor.
!
! INPUTS:
!   - input.txt: Simulation parameters.
!   - system.xyz: Initial configuration.
!
! OUTPUTS:
!   - trajectory.xyz: System trajectory file.
!   - equilibration_stats.dat: Potential energy and mean-squared displacement during equilibration phase.
!   - friction_tensor.dat: Computed friction tensor components.
!========================================================================================================

PROGRAM main
  USE kinds, ONLY: wp => dp                                           
  USE force_module
  USE friction_module
  USE minimization_module
  USE velocity_init_module
                                                        
  IMPLICIT NONE

  ! Input Variables  
  INTEGER :: nk, n_solv, i, step, n_total_atoms, atom_idx
  REAL (KIND=wp) :: dt
  REAL (KIND=wp) :: mass_solv
  REAL (KIND=wp) :: epsilon_ss, sigma_ss
  REAL (KIND=wp) :: epsilon_int, sigma_int
  REAL(KIND=wp) :: temp, k_boltz

  ! Array
  REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: pos_solv, vel_solv, force      ! (N, 3)
  REAL (KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: traj_history                 ! (N, 3, Steps)
  REAL (KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: force_history                ! (4, 3, Steps)
  REAL(KIND=wp), DIMENSION(4, 3) :: pos_solute, force_solute 
  
  REAL(KIND=wp) :: kin_en, e_pot, e_tot
  CHARACTER(LEN=2) :: atom_label  
  REAL(KIND=wp) :: inp_x, inp_y, inp_z
  REAL(KIND=wp) :: current_temp, lambda, tau_T 
  REAL(KIND=wp) :: box_L, inv_box
  REAL(KIND=wp) :: time_val
 
  ! Read Input
  OPEN(UNIT=10, FILE='input/input.txt', STATUS='old')
  
  READ(10, *) nk, dt
  READ(10, *) mass_solv, epsilon_ss, sigma_ss
  READ(10, *) epsilon_int, sigma_int
  READ(10, *) temp, k_boltz
  READ(10, *) box_L
  inv_box = 1.0_wp / box_L
  
  CLOSE(10)  
  
  PRINT *, ""
  PRINT *, "=========================================="
  PRINT *, "Molecular Dynamics with Lennard-Jones"
  PRINT *, "=========================================="
  PRINT *, ""
  PRINT '(A, I6, A, E12.5)', "Time steps: ", nk, "  dt =", dt
  PRINT '(A, G12.5)', "Solvent mass: ", mass_solv
  PRINT *, ""
  PRINT '(A, /, A, G12.5, 3X, A, G12.5)', &
        "LJ (solvent-solvent):", &
        "epsilon = ", epsilon_ss, "sigma = ", sigma_ss
  PRINT *, ""
  PRINT '(A, /, A, G12.5, 3X, A, G12.5)', &
        "LJ (solute-solvent):", &
        "epsilon = ", epsilon_int, "sigma = ", sigma_int
  PRINT *, ""
  PRINT '(A, F10.2, A)', "Simulation Temperature: ", temp, " K"
  PRINT *, ""

  OPEN(UNIT=11, FILE='input/system.xyz', STATUS='old')

  READ(11, *) n_total_atoms
  READ(11, *)
  
  n_solv = n_total_atoms - 4
  
  PRINT *, "System initialized with", n_solv, "solvent particles and", 4, "solute particles"
  PRINT *, ""
 
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
 
  PRINT *, "Starting Initial Energy Optimization..."
  
  CALL conjugate_gradient_minimize(n_solv, pos_solv, pos_solute, &
                                   epsilon_ss, sigma_ss, epsilon_int, sigma_int, &
                                   box_L, 1.0d-3, 2000)
  PRINT *, "Optimization Completed."
  PRINT *, ""  
 
  PRINT *, "Initializing Velocities..."

  CALL init_velocities(n_solv, vel_solv, mass_solv, temp, k_boltz)

  PRINT *, "Velocities assigned and scaled. Drift removed."
  PRINT *, ""  
  
  PRINT *, "Starting Equilibration..."
  
  OPEN(UNIT=40, FILE='equilibration_stats.dat', STATUS='replace')
  WRITE(40, '(A)') "# Time       E_pot"
 
  ! Calculation of the initial forces post-optimization
  CALL force_calculation(n_solv, pos_solv, pos_solute, force, force_solute, &
                         epsilon_ss, sigma_ss, epsilon_int, sigma_int, box_L, e_pot)

  ! Coupling constant
  tau_T = 1.0_wp  

  DO step = 1, 5000
    time_val = DBLE(step) * dt   
    
    ! Half-kick
    vel_solv = vel_solv + 0.5_wp * dt * (force / mass_solv)
     
    ! Drift
    pos_solv = pos_solv + dt * vel_solv
    
    DO i = 1, n_solv
      pos_solv(i, 1) = pos_solv(i, 1) - box_L * ANINT(pos_solv(i, 1) * inv_box)
      pos_solv(i, 2) = pos_solv(i, 2) - box_L * ANINT(pos_solv(i, 2) * inv_box)
      pos_solv(i, 3) = pos_solv(i, 3) - box_L * ANINT(pos_solv(i, 3) * inv_box)
    END DO

    ! New Forces
    CALL force_calculation(n_solv, pos_solv, pos_solute, force, force_solute, &
                            epsilon_ss, sigma_ss, epsilon_int, sigma_int, box_L, e_pot)
                            
    ! Half-kick
    vel_solv = vel_solv + 0.5_wp * dt * (force / mass_solv)
    
    ! Berendsen thermostat (only during equilibration)
    ! Current kinetic energy
    kin_en = 0.0_wp
    DO i = 1, n_solv
      kin_en = kin_en + SUM(vel_solv(i, :)**2)
    END DO
    kin_en = 0.5_wp * mass_solv * kin_en

    ! Current temperature
    current_temp = (2.0_wp * kin_en) / (3.0_wp * DBLE(n_solv - 1) * k_boltz)
    
    ! Berendsen scaling factor
    lambda = SQRT(1.0_wp + (dt / tau_T) * ((temp / current_temp) - 1.0_wp))
    
    ! Scale velocities
    vel_solv = vel_solv * lambda

    IF (MOD(step, 10) == 0) THEN
      WRITE(40, '(F12.4, 2X, ES14.6, 2X)') time_val, e_pot
    END IF
    
  END DO
  
  CLOSE(40)
 
  PRINT *, "Equilibration completed. Stats saved in equilibration_stats.dat"
  PRINT *, ""

  PRINT *, "Starting Molecular Dynamics..."

  OPEN(UNIT=50, FILE='energy.dat', STATUS='replace')
  WRITE(50, '(A)') "# Time(ps)       E_Kin          E_Pot          E_Tot"

  ! Calculate the forces at t = 0
  CALL force_calculation(n_solv, pos_solv, pos_solute, force, force_solute, &
                         epsilon_ss, sigma_ss, epsilon_int, sigma_int, box_L, e_pot) 

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
    
    time_val = DBLE(step) * dt

    ! Half-kick
    vel_solv = vel_solv + 0.5_wp * dt * (force / mass_solv)

    ! Drift
    pos_solv = pos_solv + dt * vel_solv

    DO i = 1, n_solv
      pos_solv(i, 1) = pos_solv(i, 1) - box_L * ANINT(pos_solv(i, 1) * inv_box)
      pos_solv(i, 2) = pos_solv(i, 2) - box_L * ANINT(pos_solv(i, 2) * inv_box)
      pos_solv(i, 3) = pos_solv(i, 3) - box_L * ANINT(pos_solv(i, 3) * inv_box)
    END DO

    ! Save Trajectory
    traj_history(:, :, step) = pos_solv

    ! Calculate new forces 
    CALL force_calculation(n_solv, pos_solv, pos_solute, force, force_solute, &
                           epsilon_ss, sigma_ss, epsilon_int, sigma_int, box_L, e_pot)
 
    ! Save Forces
    force_history(:, :, step) = force_solute

    ! Half-kick
    vel_solv = vel_solv + 0.5_wp * dt * (force / mass_solv)
  

    kin_en = 0.0_wp
    DO i = 1, n_solv
      kin_en = kin_en + SUM(vel_solv(i, :)**2)
    END DO
    kin_en = 0.5_wp * mass_solv * kin_en
     
    ! Energia Totale
    e_tot = kin_en + e_pot

    ! Salvataggio su file (ogni 10 step per non creare file enormi)
    IF (MOD(step, 10) == 0) THEN
      WRITE(50, '(F12.4, 3(2X, ES14.6))') time_val, kin_en, e_pot, e_tot
    END IF


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
  
  PRINT *, "Simulation Completed. Trajectory saved to: trajectory.xyz"
  PRINT *, ""
  PRINT *, "Computing friction tensor..."

  ! Calculate Friction Tensor
  CALL friction_tensor(nk, force_history, dt, temp, k_boltz)

  DEALLOCATE(pos_solv, vel_solv, force, traj_history, force_history)

  PRINT *, "Friction calculation completed. Results in friction_tensor.dat"
  PRINT *, ""
  PRINT *, "Program terminated successfully."

END PROGRAM main
