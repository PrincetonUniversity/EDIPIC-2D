
!===================================================================================================
! this procedure must be called before INITIATE_ION_NEUTRAL_COLLISIONS
!
SUBROUTINE INITIATE_ELECTRON_NEUTRAL_COLLISIONS

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : kB_JK, e_Cl, N_max_vel, T_e_eV, N_max_vel
  USE IonParticles, ONLY : N_spec, Ms
!  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
 
  LOGICAL exists

  CHARACTER(1) buf
  INTEGER n

  CHARACTER(23) initneutral_filename   ! init_neutral_AAAAAA.dat
                                       ! ----x----I----x----I---
  INTEGER ALLOC_ERR
  INTEGER p, i, s
  INTEGER colflag

  CHARACTER(49) initneutral_crsect_filename  ! init_neutral_AAAAAA_crsect_coll_id_NN_type_NN.dat
                                             ! ----x----I----x----I----x----I----x----I----x----

  INTEGER count

  INTEGER j

  REAL(8) energy_eV
  INTEGER indx_energy
  INTEGER indx_energy_max_prob
  REAL(8) temp

! functions
  REAL(8) frequency_of_en_collision

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  N_neutral_spec = 0
  en_collisions_turned_off = .TRUE.
  no_ionization_collisions = .TRUE.

  INQUIRE (FILE = 'init_neutrals.dat', EXIST = exists)

  IF (.NOT.exists) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### file init_neutrals.dat not found, electron-neutral collisions are turned off ###")'
     RETURN
  END IF

! specify neutral species included in simulation

  OPEN (9, FILE = 'init_neutrals.dat')
  READ (9, '(A1)') buf    !----------dd--- number of neutral species
  READ (9, '(10x,i2)') N_neutral_spec

  IF (N_neutral_spec.LT.0) N_neutral_spec = 0  ! will be used for checks later

  IF (N_neutral_spec.EQ.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### file init_neutrals.dat found, neutrals deactivated, electron-neutral collisions are turned off ###")'
     RETURN
  END IF

  ALLOCATE(neutral(1:N_neutral_spec), STAT=ALLOC_ERR)

  DO n = 1, N_neutral_spec
     READ (9, '(A1)') buf !------AAAAAA--- code/abbreviation of the neutral gas, character string 
     READ (9, '(6x,A6)') neutral(n)%name
     READ (9, '(A1)') buf !--#d.dddE#dd--- Density [m^-3]
     READ (9, '(2x,e10.3)') neutral(n)%N_m3
     READ (9, '(A1)') buf !-----ddddd.d--- Temperature [K]
     READ (9, '(5x,f7.1)') neutral(n)%T_K
  END DO
  CLOSE (9, STATUS = 'KEEP')

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
!print *, "read init_neutrals.dat"

! for each neutral species included, specify activated collisions

  DO n = 1, N_neutral_spec
     initneutral_filename = 'init_neutral_AAAAAA.dat'
     initneutral_filename(14:19) = neutral(n)%name
     INQUIRE (FILE = initneutral_filename, EXIST = exists)
     IF (.NOT.exists) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### file ",A23," NOT found, electron-neutral collisions for neutral species ",i2," (",A6,") are turned OFF ###")', &
             & initneutral_filename, n, neutral(n)%name
        neutral(n)%N_en_colproc = 0
        CYCLE
     END IF

     OPEN (9, FILE = initneutral_filename)
     READ (9, '(A1)') buf !---ddd.ddd--- mass [a.m.u.]
     READ (9, '(3x,f7.3)') neutral(n)%M_amu
     READ (9, '(A1)') buf !--------dd--- number of all possible collisional processes
     READ (9, '(8x,i2)') neutral(n)%N_en_colproc
     IF (neutral(n)%N_en_colproc.LE.0) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### file ",A23," does NOT specify electron-neutral collisions for neutral species ",i2," (",A6,"), the collisions are turned OFF ###")', &
             & initneutral_filename, n, neutral(n)%name
        CLOSE (9, STATUS = 'KEEP')
        CYCLE
     END IF

     ALLOCATE(neutral(n)%en_colproc(1:neutral(n)%N_en_colproc), STAT=ALLOC_ERR)
     DO p = 1, neutral(n)%N_en_colproc
        READ (9, '(A1)') buf !---dd--d--dd--- collision #NN :: type / activated (1/0 = Yes/No) / ion species produced (ionization collisions only) 
        READ (9, '(3x,i2,2x,i1,2x,i2)') neutral(n)%en_colproc(p)%type, colflag, neutral(n)%en_colproc(p)%ion_species_produced
        neutral(n)%en_colproc(p)%activated = .FALSE.
        IF (colflag.NE.0) neutral(n)%en_colproc(p)%activated = .TRUE.
     END DO

     READ (9, '(A1)') buf !--------dd--- number of energy segments for collision probabilities (>0)
     READ (9, '(8x,i2)') neutral(n)%N_of_energy_segments
     ALLOCATE(neutral(n)%energy_segment_boundary_value(0:neutral(n)%N_of_energy_segments), STAT=ALLOC_ERR)
     ALLOCATE(neutral(n)%energy_segment_step(1:neutral(n)%N_of_energy_segments), STAT=ALLOC_ERR)
     READ (9, '(A1)') buf !--ddddd.ddd------------- minimal energy [eV]
     READ (9, '(2x,f9.3)') neutral(n)%energy_segment_boundary_value(0)
     DO i = 1, neutral(n)%N_of_energy_segments
        READ (9, '(A1)') buf !--ddddd.ddd---ddd.ddd--- energy segment NN :: upper boundary [eV] / resolution [eV]
        READ (9, '(2x,f9.3,3x,f7.3)') neutral(n)%energy_segment_boundary_value(i), neutral(n)%energy_segment_step(i)
     END DO

     CLOSE (9, STATUS = 'KEEP')
  END DO

!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
!print *, "read ", initneutral_filename

! collision types
! 10 = elastic
! 20 = inelastic
!   21, 22, 23, etc.
! 30 = ionization, ion+ and e-
! 40 = ionization, ion++ and e- e-
! etc, we do the 3 first types only, though allow multiple inelastic collisions

! init_neutral_Xenon-.dat
! cross sections
! init_neutral_Xenon-_crsect_colltype_10.dat
! init_neutral_Xenon-_crsect_colltype_20.dat
! init_neutral_Xenon-_crsect_colltype_30.dat
! init_neutral_Xenon-_crsect_colltype_40.dat etc
! 

! for each neutral species included, for each activated collisional process, read cross sections

  DO n = 1, N_neutral_spec
     DO p = 1, neutral(n)%N_en_colproc

        IF (.NOT.neutral(n)%en_colproc(p)%activated) CYCLE

        en_collisions_turned_off = .FALSE.                                           !### flip the general collision switch

        IF (neutral(n)%en_colproc(p)%type.GE.30) no_ionization_collisions = .FALSE.  !### enable ionization collision counting

        initneutral_crsect_filename = 'init_neutral_AAAAAA_crsect_coll_id_NN_type_NN.dat'
        initneutral_crsect_filename(14:19) = neutral(n)%name
        initneutral_crsect_filename(36:37) = convert_int_to_txt_string(p, 2)
        initneutral_crsect_filename(44:45) = convert_int_to_txt_string(neutral(n)%en_colproc(p)%type, 2)
        INQUIRE( FILE = initneutral_crsect_filename, EXIST = exists)
        IF (.NOT.exists) THEN
           IF (Rank_of_process.EQ.0) PRINT '("###ERROR :: file ",A49," not found, for neutrals ",A6," while collisions ",i2," of type ",i2," are activated, program terminated")', &
                & initneutral_crsect_filename, neutral(n)%name, p, neutral(n)%en_colproc(p)%type
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        END IF

        OPEN(9, FILE = initneutral_crsect_filename)
        READ (9, '(2x,i4)') neutral(n)%en_colproc(p)%N_crsect_points
        ALLOCATE(neutral(n)%en_colproc(p)%energy_eV(1:neutral(n)%en_colproc(p)%N_crsect_points), STAT = ALLOC_ERR)
        ALLOCATE(neutral(n)%en_colproc(p)%crsect_m2(1:neutral(n)%en_colproc(p)%N_crsect_points), STAT = ALLOC_ERR)
        DO i = 1, neutral(n)%en_colproc(p)%N_crsect_points
           READ (9, '(4x,f9.3,2x,e10.3)') neutral(n)%en_colproc(p)%energy_eV(i), neutral(n)%en_colproc(p)%crsect_m2(i)
        END DO
        CLOSE (9, STATUS = 'KEEP')

! set thresholds for inelastic/ionization collisions
        neutral(n)%en_colproc(p)%threshold_energy_eV = 0.0_8   ! for elastic collisions
        IF (neutral(n)%en_colproc(p)%type.GE.20) neutral(n)%en_colproc(p)%threshold_energy_eV = neutral(n)%en_colproc(p)%energy_eV(1)

     END DO
  END DO

  ALLOCATE(collision_e_neutral(1:N_neutral_spec), STAT = ALLOC_ERR)
  DO n = 1, N_neutral_spec

! count activated collision procedures
     collision_e_neutral(n)%N_of_activated_colproc = 0
     DO p = 1, neutral(n)%N_en_colproc
        IF (neutral(n)%en_colproc(p)%activated) collision_e_neutral(n)%N_of_activated_colproc = collision_e_neutral(n)%N_of_activated_colproc + 1
     END DO

     IF (collision_e_neutral(n)%N_of_activated_colproc.EQ.0) CYCLE

! setup information about the activated collisions
     ALLOCATE(collision_e_neutral(n)%colproc_info(1:collision_e_neutral(n)%N_of_activated_colproc), STAT = ALLOC_ERR)     
     count=0
     DO p = 1, neutral(n)%N_en_colproc
        IF (.NOT.neutral(n)%en_colproc(p)%activated) CYCLE
        count=count+1
        collision_e_neutral(n)%colproc_info(count)%id_number = p
        collision_e_neutral(n)%colproc_info(count)%type = neutral(n)%en_colproc(p)%type
        collision_e_neutral(n)%colproc_info(count)%ion_species_produced = neutral(n)%en_colproc(p)%ion_species_produced
        collision_e_neutral(n)%colproc_info(count)%threshold_energy_eV = neutral(n)%en_colproc(p)%threshold_energy_eV 

        s = neutral(n)%en_colproc(p)%ion_species_produced
        IF ((s.GE.1).AND.(s.LE.N_spec)) collision_e_neutral(n)%colproc_info(count)%ion_velocity_factor = SQRT(neutral(n)%T_K * kB_JK / (T_e_eV * e_Cl)) / (N_max_vel * SQRT(Ms(s)))

     END DO

! setup general counters of collision events for minimal diagnostics
     ALLOCATE(collision_e_neutral(n)%counter(1:collision_e_neutral(n)%N_of_activated_colproc), STAT = ALLOC_ERR)     
     DO p = 1, neutral(n)%N_en_colproc
        collision_e_neutral(n)%counter(p) = 0
     END DO

! copy energy segments and steps
     collision_e_neutral(n)%N_of_energy_segments = neutral(n)%N_of_energy_segments
     ALLOCATE(collision_e_neutral(n)%energy_segment_boundary_value(0:collision_e_neutral(n)%N_of_energy_segments), STAT=ALLOC_ERR)
     collision_e_neutral(n)%energy_segment_boundary_value(0:collision_e_neutral(n)%N_of_energy_segments) = &
               & neutral(n)%energy_segment_boundary_value(0:            neutral(n)%N_of_energy_segments)
     ALLOCATE(collision_e_neutral(n)%energy_segment_step(1:collision_e_neutral(n)%N_of_energy_segments), STAT=ALLOC_ERR)
     collision_e_neutral(n)%energy_segment_step(1:collision_e_neutral(n)%N_of_energy_segments) = &
               & neutral(n)%energy_segment_step(1:            neutral(n)%N_of_energy_segments)

! identify total number of probability points / energy values
     collision_e_neutral(n)%N_of_energy_values = 0
     DO i = 1, collision_e_neutral(n)%N_of_energy_segments
        collision_e_neutral(n)%N_of_energy_values = collision_e_neutral(n)%N_of_energy_values + &
                                                  & INT( ( collision_e_neutral(n)%energy_segment_boundary_value(i) - &
                                                  &        collision_e_neutral(n)%energy_segment_boundary_value(i-1) ) / &
                                                  &      collision_e_neutral(n)%energy_segment_step(i) + 0.001_8 )
     END DO

! set index values for energy segment boundaries
     ALLOCATE(collision_e_neutral(n)%energy_segment_boundary_index(0:collision_e_neutral(n)%N_of_energy_segments), STAT=ALLOC_ERR)
     collision_e_neutral(n)%energy_segment_boundary_index(0) = 0
     DO i = 1, collision_e_neutral(n)%N_of_energy_segments
        collision_e_neutral(n)%energy_segment_boundary_index(i) = collision_e_neutral(n)%energy_segment_boundary_index(i-1) + &
                                                                & INT( ( collision_e_neutral(n)%energy_segment_boundary_value(i) - &
                                                                &        collision_e_neutral(n)%energy_segment_boundary_value(i-1) ) / &
                                                                &      collision_e_neutral(n)%energy_segment_step(i) + 0.001_8 )
     END DO

! fool proof (paranoidal)
     IF (collision_e_neutral(n)%N_of_energy_values.NE.collision_e_neutral(n)%energy_segment_boundary_index(collision_e_neutral(n)%N_of_energy_segments)) THEN
! error
        PRINT '("error")'
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF

! create array of energy values corresponding to the probability array
     ALLOCATE(collision_e_neutral(n)%energy_eV(0:collision_e_neutral(n)%N_of_energy_values), STAT=ALLOC_ERR)
     collision_e_neutral(n)%energy_eV(0) = collision_e_neutral(n)%energy_segment_boundary_value(0)
     DO i = 1, collision_e_neutral(n)%N_of_energy_segments
        DO j = collision_e_neutral(n)%energy_segment_boundary_index(i-1)+1, collision_e_neutral(n)%energy_segment_boundary_index(i)-1
           collision_e_neutral(n)%energy_eV(j) = collision_e_neutral(n)%energy_segment_boundary_value(i-1) + &
                                               & (j-collision_e_neutral(n)%energy_segment_boundary_index(i-1)) * &
                                               & collision_e_neutral(n)%energy_segment_step(i)
        END DO
        collision_e_neutral(n)%energy_eV(collision_e_neutral(n)%energy_segment_boundary_index(i)) = collision_e_neutral(n)%energy_segment_boundary_value(i)
     END DO

! set the probability arrays --------------------------------------------------------->>>>>>>>>>>>>>>

     ALLOCATE(collision_e_neutral(n)%prob_colproc_energy(1:collision_e_neutral(n)%N_of_activated_colproc, 0:collision_e_neutral(n)%N_of_energy_values), STAT = ALLOC_ERR)

     DO indx_energy = 0, collision_e_neutral(n)%N_of_energy_values

        energy_eV = collision_e_neutral(n)%energy_eV(indx_energy)

        collision_e_neutral(n)%prob_colproc_energy(1, indx_energy) = frequency_of_en_collision(energy_eV, n, collision_e_neutral(n)%colproc_info(1)%id_number)
        DO p = 2, collision_e_neutral(n)%N_of_activated_colproc
           collision_e_neutral(n)%prob_colproc_energy(p, indx_energy) = collision_e_neutral(n)%prob_colproc_energy(p-1, indx_energy) + &
                                                                   & frequency_of_en_collision(energy_eV, n, collision_e_neutral(n)%colproc_info(p)%id_number)
        END DO
     END DO

! find the maximum
     indx_energy_max_prob = 0
     DO indx_energy = 1, collision_e_neutral(n)%N_of_energy_values
        IF ( collision_e_neutral(n)%prob_colproc_energy(collision_e_neutral(n)%N_of_activated_colproc, indx_energy).GT. &
           & collision_e_neutral(n)%prob_colproc_energy(collision_e_neutral(n)%N_of_activated_colproc, indx_energy_max_prob) ) indx_energy_max_prob = indx_energy
     END DO

! find maximal fraction of colliding particles for the neutral density neutral(n)%N_m3
     collision_e_neutral(n)%max_colliding_fraction = 1.0_8 - EXP(-collision_e_neutral(n)%prob_colproc_energy(collision_e_neutral(n)%N_of_activated_colproc, indx_energy_max_prob))

! convert accumulated collision frequencies into boundaries of probability ranges for collision processes
     temp = collision_e_neutral(n)%prob_colproc_energy(collision_e_neutral(n)%N_of_activated_colproc, indx_energy_max_prob)
     DO indx_energy = 0, collision_e_neutral(n)%N_of_energy_values
        DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
           collision_e_neutral(n)%prob_colproc_energy(p, indx_energy) = MAX(0.0_8, MIN(1.0_8, collision_e_neutral(n)%prob_colproc_energy(p, indx_energy)/temp))
        END DO
     END DO
     collision_e_neutral(n)%prob_colproc_energy(collision_e_neutral(n)%N_of_activated_colproc, indx_energy_max_prob) = 1.0_8

! fool proof (paranoidal again) - make sure that boundaries of probability ranges do not decrease
     DO indx_energy = 0, collision_e_neutral(n)%N_of_energy_values
        DO p = 1, collision_e_neutral(n)%N_of_activated_colproc-1
           collision_e_neutral(n)%prob_colproc_energy(p+1, indx_energy) = MAX( collision_e_neutral(n)%prob_colproc_energy(p,   indx_energy), &
                                                                             & collision_e_neutral(n)%prob_colproc_energy(p+1, indx_energy) )
        END DO
     END DO

  END DO   !###  DO n = 1, N_neutral_spec

END SUBROUTINE INITIATE_ELECTRON_NEUTRAL_COLLISIONS

!### remember that collisions occur with the ion!!!!! time step
!### also figure out how to use dim-less v^2 instead of energy when finding collision probability

!-----------------------------------------------------------------
!
SUBROUTINE PERFORM_ELECTRON_NEUTRAL_COLLISIONS

  USE ParallelOperationValues
  USE MCCollisions
  USE ElectronParticles
  USE CurrentProblemValues, ONLY : V_scale_ms, m_e_kg, e_Cl
  USE ClusterAndItsBoundaries
  USE rng_wrapper
  USE Snapshots

!  USE ParallelOperationValues

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER ALLOC_ERR

  INTEGER n, p

  LOGICAL first_collision_not_done_yet

  REAL(8) R_collided, F_collided
  INTEGER I_collided

  INTEGER j

  REAL(8) random_r, random_n
  INTEGER random_j

  REAL(8) energy_eV
  INTEGER indx_energy
  REAL(8) a1, a0
  INTEGER i, indx_segment

  INTEGER k, indx_coll

  INTEGER n1, n2, n3, bufsize, pos
  REAL, ALLOCATABLE :: rbufer(:), rbufer2(:)

! functions
  LOGICAL Find_in_stored_list
  REAL(8) neutral_density_normalized

  INTERFACE
     RECURSIVE SUBROUTINE Node_Killer(node)
       USE MCCollisions
       TYPE (binary_tree), POINTER :: node
     END SUBROUTINE Node_Killer

     RECURSIVE SUBROUTINE Transfer_collisions_from_stored_list(node, n_neutral, indx_coll, bufsize, n1, n2, n3, rbufer)
       USE MCCollisions
       USE ElectronParticles
       USE ClusterAndItsBoundaries
       TYPE (binary_tree), POINTER :: node
       INTEGER, INTENT(IN) :: n_neutral, indx_coll, bufsize, n1, n2, n3
       REAL, DIMENSION(bufsize), INTENT(INOUT) :: rbufer
     END SUBROUTINE Transfer_collisions_from_stored_list
  END INTERFACE

  IF (en_collisions_turned_off) RETURN

! Allocate binary tree to store the numbers of particles which have collided already
  NULLIFY(Collided_particle)
  IF (.NOT.ASSOCIATED(Collided_particle)) THEN
     ALLOCATE(Collided_particle, STAT=ALLOC_ERR)
  END IF
  NULLIFY(Collided_particle%Larger)
  NULLIFY(Collided_particle%Smaller)

! clear all collision counters
  DO n = 1, N_neutral_spec
     DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
        collision_e_neutral(n)%counter(p) = 0
     END DO
  END DO

  first_collision_not_done_yet = .TRUE.
  
  DO n = 1, N_neutral_spec

     IF (collision_e_neutral(n)%N_of_activated_colproc.EQ.0) CYCLE

     R_collided = collision_e_neutral(n)%max_colliding_fraction * N_electrons
     I_collided = INT(R_collided)
     F_collided = R_collided - I_collided

     DO j = 0, I_collided

        IF (j.EQ.0) THEN
! here we process the "fractional" collisional event
           IF (well_random_number().GT.F_collided) CYCLE
        END IF

!------------- Determine the kind of collision for the selected particle
        random_r = well_random_number()
        random_n = well_random_number()

        IF (first_collision_not_done_yet) THEN
           random_j = INT(well_random_number() * N_electrons)
           random_j = MIN(MAX(random_j, 1), N_electrons)
        ELSE
           DO                        ! search will be repeated until a number will be successfully obtained
              random_j = INT(well_random_number() * N_electrons)
              random_j = MIN(MAX(random_j, 1), N_electrons)
              IF (.NOT.Find_in_stored_list(random_j)) EXIT    !#### needs some safety mechanism to avoid endless cycling
           END DO
        END IF

! account for reduced neutral density
        IF (random_n.GT.neutral_density_normalized(n, electron(random_j)%X, electron(random_j)%Y)) CYCLE      

        energy_eV = (electron(random_j)%VX**2 + electron(random_j)%VY**2 + electron(random_j)%VZ**2) * V_scale_ms**2 * m_e_kg * 0.5_8 / e_Cl

        IF (energy_eV.GE.collision_e_neutral(n)%energy_segment_boundary_value(collision_e_neutral(n)%N_of_energy_segments)) THEN
           indx_energy = collision_e_neutral(n)%N_of_energy_values-1
           a1 = 1.0_8
           a0 = 0.0_8
        ELSE IF (energy_eV.LT.collision_e_neutral(n)%energy_segment_boundary_value(0)) THEN
           indx_energy = 0
           a1 = 0.0_8
           a0 = 1.0_8
        ELSE
           DO i = collision_e_neutral(n)%N_of_energy_segments-1, 0, -1
              IF (energy_eV.GE.collision_e_neutral(n)%energy_segment_boundary_value(i)) THEN
                 indx_segment = i+1
                 EXIT
              END IF
           END DO
           indx_energy =              collision_e_neutral(n)%energy_segment_boundary_index(indx_segment-1) + &
                       & (energy_eV - collision_e_neutral(n)%energy_segment_boundary_value(indx_segment-1)) / collision_e_neutral(n)%energy_segment_step(indx_segment)
           indx_energy = MAX(0, MIN(indx_energy, collision_e_neutral(n)%N_of_energy_values-1))
! double check
           IF ((energy_eV.LT.collision_e_neutral(n)%energy_eV(indx_energy)).OR.(energy_eV.GT.collision_e_neutral(n)%energy_eV(indx_energy+1))) THEN
! error
              PRINT '("Proc ",i4," error in PERFORM_ELECTRON_NEUTRAL_COLLISIONS")', Rank_of_process
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF
           a1 = (energy_eV - collision_e_neutral(n)%energy_eV(indx_energy)) / collision_e_neutral(n)%energy_segment_step(indx_segment)
           a0 = 1.0_8 - a1
        END IF

        DO k = collision_e_neutral(n)%N_of_activated_colproc, 1, -1 
           IF (random_r.GT.(a0 * collision_e_neutral(n)%prob_colproc_energy(k, indx_energy) + &
                         &  a1 * collision_e_neutral(n)%prob_colproc_energy(k, indx_energy+1)) ) EXIT
        END DO
        indx_coll = k + 1        

        IF (indx_coll.GT.collision_e_neutral(n)%N_of_activated_colproc) CYCLE   ! the null collision

        IF (first_collision_not_done_yet) THEN
           CALL Add_to_stored_list_first_time(random_j, n, indx_coll)
           first_collision_not_done_yet = .FALSE.
        ELSE
           CALL Add_to_stored_list(random_j, n, indx_coll)
        END IF

        SELECT CASE (collision_e_neutral(n)%colproc_info(indx_coll)%type)
        CASE (10)
! scattering with the angle Ksi from Okhrimovsky
           CALL en_Collision_Elastic_10( n, random_j, energy_eV, &
                                       & collision_e_neutral(n)%counter(indx_coll))
        CASE (11)
! isotropic scattering
           CALL en_Collision_Elastic_11( n, random_j, energy_eV, &
                                       & collision_e_neutral(n)%counter(indx_coll))
        CASE (20)
           CALL en_Collision_Inelastic_20( n, random_j, energy_eV, &
                                          & collision_e_neutral(n)%colproc_info(indx_coll)%threshold_energy_eV, &
                                          & collision_e_neutral(n)%counter(indx_coll) )
        CASE (30)
           CALL en_Collision_Ionization_30( n, random_j, energy_eV, &
                                          & collision_e_neutral(n)%colproc_info(indx_coll)%threshold_energy_eV, &
                                          & collision_e_neutral(n)%colproc_info(indx_coll)%ion_species_produced, &
                                          & collision_e_neutral(n)%colproc_info(indx_coll)%ion_velocity_factor, &
                                          & collision_e_neutral(n)%counter(indx_coll) )

        END SELECT

     END DO
  END DO

!### exit if there is no ionization collisions
  IF (no_ionization_collisions) THEN
     CALL Node_Killer(Collided_particle)
     RETURN
  END IF

!### exit if the last snapshot has been saved alread, also works when no snapshots are requested (N_of_all_snaps=0)
  IF (current_snap.GT.N_of_all_snaps) THEN
     CALL Node_Killer(Collided_particle)
     RETURN
  END IF

!### exit if 2d-collision diagnostics not requested
  IF (.NOT.save_ionization_rates_2d(current_snap)) THEN
     CALL Node_Killer(Collided_particle)
     RETURN
  END IF

! account for collisions in diagnostics arrays

  n1 = c_indx_y_max - c_indx_y_min + 1
  n3 = c_indx_x_max - c_indx_x_min + 1
  n2 = -c_indx_x_min + 1 - c_indx_y_min * n3
  bufsize = n1 * n3
  ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)
  IF (Rank_cluster.EQ.0) THEN
     ALLOCATE(rbufer2(1:bufsize), STAT=ALLOC_ERR)
  ELSE
! MPI_REDUCE requires a separate receiver array for all processes, 
! but in all processes except the one where the final result (e.g. sum) is assembled 
! it is sufficient to have a valid array of minimal size. 
! The same is done in CREATE_SNAPSHOT.
     ALLOCATE(rbufer2(1), STAT=ALLOC_ERR)
  END IF

  DO n = 1, N_neutral_spec
     DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
        IF (collision_e_neutral(n)%colproc_info(p)%type.LT.30) CYCLE  ! skip non-ionizing collisions

        rbufer = 0.0
        rbufer2 = 0.0

! note, if first_collision_not_done_yet is .TRUE. then no collisions took place above and there is nothing to transfer
        IF (.NOT.first_collision_not_done_yet) THEN
           CALL Transfer_collisions_from_stored_list(Collided_particle, n, p, bufsize, n1, n2, n3, rbufer)
        END IF

        CALL MPI_REDUCE(rbufer, rbufer2, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

        IF (cluster_rank_key.EQ.0) THEN
           pos=1
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 diagnostics_neutral(n)%activated_collision(p)%counter_local(i,j) = diagnostics_neutral(n)%activated_collision(p)%counter_local(i,j) + rbufer2(pos)
                 pos=pos+1
              END DO
           END DO
        END IF
  
        CALL MPI_BARRIER(COMM_CLUSTER, ierr) !??????????

     END DO
  END DO

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT = ALLOC_ERR)

  CALL Node_Killer(Collided_particle)

END SUBROUTINE PERFORM_ELECTRON_NEUTRAL_COLLISIONS

!-------------------------------------------------------------------------------------------
!
SUBROUTINE INITIATE_en_COLL_DIAGNOSTICS

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : Start_T_cntr, N_subcycles, delta_t_s
  USE Checkpoints, ONLY : use_checkpoint

  IMPLICIT NONE

  INTEGER n
                                      ! ----x----I----x----I----x--
  CHARACTER(27) historycoll_filename  ! history_coll_e_n_AAAAAA.dat

  LOGICAL exists
  CHARACTER(1) buf
  INTEGER k, i
  INTEGER i_dummy, ios

  IF (Rank_of_process.NE.0) RETURN

  IF (en_collisions_turned_off) RETURN

  IF (use_checkpoint.EQ.1) THEN
! start from checkpoint, must trim the time dependences

     DO n = 1, N_neutral_spec

        IF (collision_e_neutral(n)%N_of_activated_colproc.EQ.0) CYCLE

        historycoll_filename = 'history_coll_e_n_AAAAAA.dat'
        historycoll_filename(18:23) = neutral(n)%name

        INQUIRE (FILE = historycoll_filename, EXIST = exists)
        IF (exists) THEN                                                       

           OPEN (21, FILE = historycoll_filename, STATUS = 'OLD')          
! skip header
           READ (21, '(A1)') buf ! WRITE (21, '("# electron time step is ",e14.7," s")'), delta_t_s
           READ (21, '(A1)') buf ! WRITE (21, '("#      ion time step is ",e14.7," s")'), delta_t_s * N_subcycles
           READ (21, '(A1)') buf ! WRITE (21, '("# column  1 is the electron step counter")')
           READ (21, '(A1)') buf ! WRITE (21, '("# column  2 is the total number of electron macroparticles in the whole system")')
           DO k = 1, collision_e_neutral(n)%N_of_activated_colproc
              READ (21, '(A1)') buf ! WRITE (21, '(# column ",i2," is the number of collision events during past ion time step for collision process with id ",i2," type ",i2)') &
                                    ! &  k+2, collision_e_neutral(n)%colproc_info(k)%id_number, collision_e_neutral(n)%colproc_info(k)%type
           END DO

           DO !i = 1, Start_T_cntr / N_subcycles            ! these files are updated at every ion timestep
              READ (21, '(2x,i8,10(2x,i8))', iostat = ios) i_dummy
              IF (ios.NE.0) EXIT
              IF (i_dummy.GE.Start_T_cntr) EXIT
           END DO
           BACKSPACE(21)
           ENDFILE 21       
           CLOSE (21, STATUS = 'KEEP')  

        ELSE

           OPEN  (21, FILE = historycoll_filename, STATUS = 'REPLACE')
! save header
           WRITE (21, '("# electron time step is ",e14.7," s")') delta_t_s
           WRITE (21, '("#      ion time step is ",e14.7," s")') delta_t_s * N_subcycles
           WRITE (21, '("# column  1 is the electron step counter")')
           WRITE (21, '("# column  2 is the total number of electron macroparticles in the whole system")')
           DO k = 1, collision_e_neutral(n)%N_of_activated_colproc
              WRITE (21, '("# column ",i2," is the number of collision events during past ion time step for collision process with id ",i2," type ",i2)') &
                   &  k+2, collision_e_neutral(n)%colproc_info(k)%id_number, collision_e_neutral(n)%colproc_info(k)%type
           END DO
           CLOSE (21, STATUS = 'KEEP')

        END IF

     END DO

  ELSE
! fresh start, empty files

     DO n = 1, N_neutral_spec

        IF (collision_e_neutral(n)%N_of_activated_colproc.EQ.0) CYCLE

        historycoll_filename = 'history_coll_e_n_AAAAAA.dat'
        historycoll_filename(18:23) = neutral(n)%name

        OPEN  (21, FILE = historycoll_filename, STATUS = 'REPLACE')
! save header
        WRITE (21, '("# electron time step is ",e14.7," s")') delta_t_s
        WRITE (21, '("#      ion time step is ",e14.7," s")') delta_t_s * N_subcycles
        WRITE (21, '("# column  1 is the electron step counter")')
        WRITE (21, '("# column  2 is the total number of electron macroparticles in the whole system")')
        DO k = 1, collision_e_neutral(n)%N_of_activated_colproc
           WRITE (21, '("# column ",i2," is the number of collision events during past ion time step for collision process with id ",i2," type ",i2)') &
                &  k+2, collision_e_neutral(n)%colproc_info(k)%id_number, collision_e_neutral(n)%colproc_info(k)%type
        END DO
        CLOSE (21, STATUS = 'KEEP')

     END DO

  END IF

END SUBROUTINE INITIATE_en_COLL_DIAGNOSTICS

!-------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_en_COLLISIONS

  USE ParallelOperationValues
  USE MCCollisions
  USE ElectronParticles, ONLY : N_electrons
  USE CurrentProblemValues, ONLY : T_cntr

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER buflen, pos, n, p
  INTEGER, ALLOCATABLE :: ibufer_send(:)
  INTEGER, ALLOCATABLE :: ibufer_receive(:)
  INTEGER ALLOC_ERR

  INTEGER global_N_of_electrons

                                      ! ----x----I----x----I----x--
  CHARACTER(27) historycoll_filename  ! history_coll_e_n_AAAAAA.dat

! report all collision counters to the process with zero global rank

  buflen=1             ! include N_electrons
  DO n = 1, N_neutral_spec
     buflen = buflen + collision_e_neutral(n)%N_of_activated_colproc
  END DO

  ALLOCATE (ibufer_send(1:buflen), STAT = ALLOC_ERR)
  ALLOCATE (ibufer_receive(1:buflen), STAT = ALLOC_ERR)

  pos=1
  ibufer_send(pos) = N_electrons
  DO n = 1, N_neutral_spec
     DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
        pos = pos+1
        ibufer_send(pos) = collision_e_neutral(n)%counter(p)
     END DO
  END DO
  
  ibufer_receive = 0

  CALL MPI_REDUCE(ibufer_send, ibufer_receive, buflen, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  IF (Rank_of_process.EQ.0) THEN
! translate the message
     pos=1
     global_N_of_electrons = ibufer_receive(pos)
     DO n = 1, N_neutral_spec
        IF (collision_e_neutral(n)%N_of_activated_colproc.LT.1) CYCLE
        DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
           pos = pos+1
           collision_e_neutral(n)%counter(p) = ibufer_receive(pos)
        END DO
        PRINT '("Total collisions with neutral species ",i2," (",A6,") :: ",10(2x,i6))', n, neutral(n)%name, collision_e_neutral(n)%counter(1:collision_e_neutral(n)%N_of_activated_colproc)
     END DO
  END IF

  DEALLOCATE(ibufer_send, STAT = ALLOC_ERR)
  DEALLOCATE(ibufer_receive, STAT = ALLOC_ERR)

  IF (Rank_of_process.NE.0) RETURN

  DO n = 1, N_neutral_spec

     IF (collision_e_neutral(n)%N_of_activated_colproc.EQ.0) CYCLE

     historycoll_filename = 'history_coll_e_n_AAAAAA.dat'
     historycoll_filename(18:23) = neutral(n)%name

     OPEN (21, FILE = historycoll_filename, POSITION = 'APPEND')
     WRITE (21, '(2x,i8,10(2x,i8))') &
          & T_cntr, &
          & global_N_of_electrons, &
          & collision_e_neutral(n)%counter(1:collision_e_neutral(n)%N_of_activated_colproc)
     CLOSE (21, STATUS = 'KEEP')
  
  END DO

END SUBROUTINE SAVE_en_COLLISIONS

!----------------------------------------
!
REAL(8) FUNCTION neutral_density_normalized(n, x, y)

!  USE CurrentProblemValues, ONLY : N_cells, delta_x_m

  IMPLICIT NONE

  INTEGER n
  REAL(8) x, y
  
  neutral_density_normalized  = 1.0_8  ! use this in case of uniform neutral density 

END FUNCTION

!-----------------------------------------------------------------
! function's value equals .TRUE. if the particle is already stored,
! otherwise function's value is .FALSE. (i.e. particle does not collide yet)
LOGICAL FUNCTION Find_in_stored_list(number)

  USE MCCollisions
  IMPLICIT NONE

  INTEGER number
  TYPE (binary_tree), POINTER :: current

  Find_in_stored_list = .FALSE.

  current => Collided_particle

  DO 

     IF (number.GT.current%number) THEN
        IF (ASSOCIATED(current%Larger)) THEN
           current => current%Larger               ! 
           CYCLE                                   ! go to the next node, with larger "number"
        ELSE
           EXIT
        END IF
     END IF

     IF (number.LT.current%number) THEN
        IF (ASSOCIATED(current%Smaller)) THEN
           current => current%Smaller              ! 
           CYCLE                                   ! go to the next node, with smaller "number"
        ELSE
           EXIT
        END IF
     END IF

     Find_in_stored_list = .TRUE.                  ! number.EQ.current%number
     EXIT                                          ! if we are here, then we found the match
     
  END DO

END FUNCTION Find_in_stored_list

!-----------------------------------------------------------------
! subroutine adds number to the binary tree
! we assume that there are no nodes in the tree with the same value yet
SUBROUTINE Add_to_stored_list_first_time(number, n_neutral, indx_coll)

  USE MCCollisions

  IMPLICIT NONE

  INTEGER number
  INTEGER n_neutral
  INTEGER indx_coll

!  TYPE (binary_tree), POINTER :: current
!  INTEGER ALLOC_ERR

!  current => Collided_particle                  ! start from the head node of the binary tree

  Collided_particle%number    = number                       ! collided electron particle number
  Collided_particle%neutral   = n_neutral                    ! neutral species number
  Collided_particle%indx_coll = indx_coll                    ! index of activated collision [not collision id]

!  NULLIFY(current%Larger) !### was done when Collided_particle was allocated
!  NULLIFY(current%Smaller)

END SUBROUTINE Add_to_stored_list_first_time

!-----------------------------------------------------------------
! subroutine adds number to the binary tree
! we assume that there are no nodes in the tree with the same value yet
SUBROUTINE Add_to_stored_list(number, n_neutral, indx_coll)

  USE MCCollisions

  IMPLICIT NONE

  INTEGER number
  INTEGER n_neutral
  INTEGER indx_coll

  TYPE (binary_tree), POINTER :: current
  INTEGER ALLOC_ERR

  current => Collided_particle                  ! start from the head node of the binary tree

  DO                                            ! go through the allocated nodes to the end of the branch

     IF (number.GT.current%number) THEN         
        IF (ASSOCIATED(current%Larger)) THEN        
           current => current%Larger               ! 
           CYCLE                                   ! go to the next node, with larger "number"
        ELSE
           ALLOCATE(current%Larger, STAT=ALLOC_ERR)
           current => current%Larger
           EXIT
        END IF
     END IF

     IF (number.LT.current%number) THEN
        IF (ASSOCIATED(current%Smaller)) THEN        
           current => current%Smaller              ! 
           CYCLE                                   ! go to the next node, with smaller "number"
        ELSE
           ALLOCATE(current%Smaller, STAT=ALLOC_ERR)
           current => current%Smaller
           EXIT
        END IF
     END IF
     
  END DO

  current%number    = number                       ! collided electron particle number
  current%neutral   = n_neutral                    ! neutral species number
  current%indx_coll = indx_coll                    ! index of activated collision [not collision id]

  NULLIFY(current%Larger)
  NULLIFY(current%Smaller)

END SUBROUTINE Add_to_stored_list

!-----------------------------------------------------------------
!
RECURSIVE SUBROUTINE Transfer_collisions_from_stored_list(node, n_neutral, indx_coll, bufsize, n1, n2, n3, rbufer)

  USE MCCollisions
  USE ElectronParticles
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  TYPE (binary_tree), POINTER :: node

  INTEGER, INTENT(IN) :: n_neutral, indx_coll, bufsize, n1, n2, n3
  REAL, DIMENSION(bufsize), INTENT(INOUT) :: rbufer

  INTEGER i, j, k
  INTEGER pos_i_j, pos_ip1_j, pos_i_jp1, pos_ip1_jp1
  REAL ax_ip1, ax_i, ay_jp1, ay_j
  REAL vij, vip1j, vijp1

  IF ((node%neutral.EQ.n_neutral).AND.(node%indx_coll.EQ.indx_coll)) THEN

     k = node%number

     i = INT(electron(k)%X)
     j = INT(electron(k)%Y)
     IF (electron(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (electron(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1
     
     pos_i_j     = i + j * n3 + n2
     pos_ip1_j   = pos_i_j + 1
     pos_i_jp1   = pos_i_j + n3
     pos_ip1_jp1 = pos_i_jp1 + 1

     ax_ip1 = REAL(electron(k)%X) - REAL(i)
     ax_i   = 1.0 - ax_ip1

     ay_jp1 = REAL(electron(k)%Y) - REAL(j)
     ay_j   = 1.0 - ay_jp1

     vij   = ax_i   * ay_j
     vip1j = ax_ip1 * ay_j
     vijp1 = ax_i   * ay_jp1

     rbufer(pos_i_j)     = rbufer(pos_i_j)     + vij                         !ax_i   * ay_j
     rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   + vip1j                       !ax_ip1 * ay_j
     rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   + vijp1                       !ax_i   * ay_jp1
     rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) + 1.0 - vij - vip1j - vijp1   !ax_ip1 * ay_jp1

  END IF

  IF (ASSOCIATED(node%Larger)) CALL Transfer_collisions_from_stored_list(node%Larger, n_neutral, indx_coll, bufsize, n1, n2, n3, rbufer)

  IF (ASSOCIATED(node%Smaller)) CALL Transfer_collisions_from_stored_list(node%Smaller, n_neutral, indx_coll, bufsize, n1, n2, n3, rbufer)

  RETURN

END SUBROUTINE Transfer_collisions_from_stored_list

!---------------------------------------------------
! this subroutine kills the nodes of the binary tree
RECURSIVE SUBROUTINE Node_Killer(node)

  USE MCCollisions
  IMPLICIT NONE
  
  TYPE (binary_tree), POINTER :: node
  INTEGER DEALLOC_ERR

  IF (ASSOCIATED(node%Larger))  CALL Node_Killer(node%Larger)
  IF (ASSOCIATED(node%Smaller)) CALL Node_Killer(node%Smaller)

  DEALLOCATE(node, STAT=DEALLOC_ERR)

  RETURN

END SUBROUTINE Node_Killer
