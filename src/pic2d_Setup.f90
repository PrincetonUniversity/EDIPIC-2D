!-----------------------------------------
!
SUBROUTINE PREPARE_SETUP_VALUES

  USE ParallelOperationValues
  USE CurrentProblemValues

  IMPLICIT NONE

  INTEGER n

  CHARACTER(14) initbo_filename
  LOGICAl exists
  
  CHARACTER(1) buf 

  REAL(8) Te_normal_constant_emit_eV
  REAL(8) Te_parallel_constant_emit_eV
  REAL(8) We_beam_constant_emit_eV

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  DO n = 1, N_of_boundary_objects

     whole_object(n)%material = '      '
     whole_object(n)%phi_const = 0.0_8
     whole_object(n)%phi_var = 0.0_8
     whole_object(n)%omega = 0.0_8
     whole_object(n)%phase = 0.0_8
     whole_object(n)%N_electron_constant_emit = 0
     whole_object(n)%model_constant_emit = 0
     whole_object(n)%factor_convert_vinj_normal_constant_emit = 0.0_8
     whole_object(n)%factor_convert_vinj_parallel_constant_emit = 0.0_8
     whole_object(n)%v_ebeam_constant_emit = 0.0_8

     initbo_filename = 'init_bo_NN.dat'
     initbo_filename(9:10) = convert_int_to_txt_string(n, 2)

     INQUIRE (FILE = initbo_filename, EXIST = exists)
     IF (.NOT.exists) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### PREPARE_SETUP_VALUES :: file ",A14," not found, use default setup values for boundary object ",i2)', initbo_filename, n
        CYCLE
     END IF

     OPEN (9, FILE = initbo_filename)

     READ (9, '(A1)') buf !------AAAAAA--- code/abbreviation of the material, character string 
     READ (9, '(6x,A6)') whole_object(n)%material
     READ (9, '(A1)') buf !---ddddd.ddd--- constant potential [V]
     READ (9, '(3x,f9.3)') whole_object(n)%phi_const
     READ (9, '(A1)') buf !---ddddd.ddd--- amplitude of potential oscillations [V]
     READ (9, '(3x,f9.3)') whole_object(n)%phi_var
     READ (9, '(A1)') buf !---ddddd.ddd--- frequency [MHz]
     READ (9, '(3x,f9.3)') whole_object(n)%omega
     READ (9, '(A1)') buf !---ddddd.ddd--- phase [deg] for sin(omega*t+phase)
     READ (9, '(3x,f9.3)') whole_object(n)%phase
     READ (9, '(A1)') buf !----dddd------- number of electron macroparticles injected each timestep, constant [dim-less]
     READ (9, '(4x,i4)') whole_object(n)%N_electron_constant_emit
     READ (9, '(A1)') buf !-------d------- emission model (0 = thermal emission, 1 = electron beam)
     READ (9, '(7x,i1)') whole_object(n)%model_constant_emit
     READ (9, '(A1)') buf !----dddd.ddd--- temperature of emitted electrons / half-energy-spread of the beam (Tb) normal to the wall [eV] (>=0)
     READ (9, '(4x,f8.3)') Te_normal_constant_emit_eV
     READ (9, '(A1)') buf !----dddd.ddd--- temperature of emitted electrons parallel to the wall [eV] (>=0)
     READ (9, '(4x,f8.3)') Te_parallel_constant_emit_eV
     READ (9, '(A1)') buf !----dddd.ddd--- energy of the electron beam [eV] (>3Tb/2)
     READ (9, '(4x,f8.3)') We_beam_constant_emit_eV

     CLOSE (9, STATUS = 'KEEP')

     IF (whole_object(n)%model_constant_emit.EQ.0) THEN
! thermal emission
        whole_object(n)%factor_convert_vinj_normal_constant_emit = SQRT(Te_normal_constant_emit_eV / T_e_eV) / N_max_vel
     ELSE
! electron beam
!
! For beam injection, we take velocity from a Maxwellian distribution with certain (see below) temperature, then add the beam velocity.
! We must ensure that particle velocity is directed away from the wall, therefore the beam velocity must exceed 
! 3 thermal velocities for the temperature used to initialize beam particles.
! When the beam is requested, the normal temperature is treated as the scale of energy spread of beam electrons in the laboratory frame.
! In order to achieve this, it is necessary to use a different (smaller) value of temperature at the first step of selecting beam velocities.
! 
! check that the beam energy is above the threshold
        IF (We_beam_constant_emit_eV.LE.(1.5_8 * Te_normal_constant_emit_eV)) THEN
           IF (Rank_of_process.EQ.0) PRINT '("Error in PREPARE_SETUP_VALUES, boundary object ",i2," :: for beam energy ",f8.3," eV the beam energy spread ",f7.3," eV exceeds threshold ",f8.3," eV")', &
                & n, We_beam_constant_emit_eV, Te_normal_constant_emit_eV, We_beam_constant_emit_eV/1.5
           STOP
        END IF
! scaling factor is calculated with a reduced temperature to achieve required energy spread in the beam frame
        whole_object(n)%factor_convert_vinj_normal_constant_emit = SQRT(0.25_8 * Te_normal_constant_emit_eV**2 / (We_beam_constant_emit_eV * T_e_eV)) / N_max_vel
        whole_object(n)%v_ebeam_constant_emit = SQRT(We_beam_constant_emit_eV / T_e_eV) / N_max_vel
     END IF

     whole_object(n)%factor_convert_vinj_parallel_constant_emit = SQRT(Te_parallel_constant_emit_eV / T_e_eV) / N_max_vel
  END DO

END SUBROUTINE PREPARE_SETUP_VALUES

!--------------------------------------------
!
SUBROUTINE PERFORM_ELECTRON_EMISSION_SETUP

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE SetupValues

  USE rng_wrapper

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER n, m, nwo
  INTEGER add_N_e_to_emit

  INTEGER k
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER, ALLOCATABLE :: ibuf_send(:)
  INTEGER, ALLOCATABLE :: ibuf_receive(:)
  INTEGER ALLOC_ERR

  IF (ht_use_e_emission_from_cathode.OR.ht_use_e_emission_from_cathode_zerogradf.OR.ht_emission_constant) RETURN   ! either PERFORM_ELECTRON_EMISSION_HT_SETUP or PERFORM_ELECTRON_EMISSION_HT_SETUP_ZERO_GRAD_F will be called instead

  IF (c_N_of_local_object_parts.LE.0) RETURN

! in a cluster, all processes have same copy of c_local_object_part (except c_local_object_part%segment_number) and c_index_of_local_object_part_*
! all processes have same copy of whole_object
! therefore each process in a cluster can calculate the number of particles to inject due to constant emission from bondary objects itself
! that is without additional communications with the master

! boundary objects along the left edge of the cluster
  DO n = 1, c_N_of_local_object_parts_left
        m = c_index_of_local_object_part_left(n)
        add_N_e_to_emit = 0
        nwo = c_local_object_part(m)%object_number
        IF (whole_object(nwo)%N_electron_constant_emit.LE.0) CYCLE
        IF (c_left_top_corner_type.EQ.FLAT_WALL_LEFT) THEN
! account for overlapping
           add_N_e_to_emit = INT( whole_object(nwo)%N_electron_constant_emit * REAL(MIN(c_local_object_part(m)%jend, c_indx_y_max-1) - c_local_object_part(m)%jstart ) / REAL(whole_object(nwo)%L) )
        ELSE
           add_N_e_to_emit = INT( whole_object(nwo)%N_electron_constant_emit * REAL(c_local_object_part(m)%jend - c_local_object_part(m)%jstart ) / REAL(whole_object(nwo)%L) )
        END IF

! account for emission split between multiple processes
     IF (Rank_cluster.EQ.N_processes_cluster-1) THEN
        add_N_e_to_emit = add_N_e_to_emit - (N_processes_cluster-1) * (add_N_e_to_emit / N_processes_cluster)
     ELSE
        add_N_e_to_emit = add_N_e_to_emit / N_processes_cluster
     END IF

! emission
     DO k = 1, add_N_e_to_emit
        IF (c_left_top_corner_type.EQ.FLAT_WALL_LEFT) THEN 
           y = DBLE(c_local_object_part(m)%jstart) + well_random_number() * DBLE( MIN(c_local_object_part(m)%jend, c_indx_y_max-1) - c_local_object_part(m)%jstart )
           y = MIN(MAX(y, DBLE(c_indx_y_min)), DBLE(c_indx_y_max-1))
        ELSE
! cluster which has no overlapping from above at the top left corner
           y = DBLE(c_local_object_part(m)%jstart) + well_random_number() * DBLE( c_local_object_part(m)%jend - c_local_object_part(m)%jstart )
           y = MIN(MAX(y, DBLE(c_indx_y_min)), DBLE(c_indx_y_max))
        END IF
        x = DBLE(c_indx_x_min) + 1.0d-6   !???

        IF (whole_object(nwo)%model_constant_emit.EQ.0) THEN
! thermal emission
           CALL GetInjMaxwellVelocity(vx)
           vx = vx * whole_object(nwo)%factor_convert_vinj_normal_constant_emit
        ELSE
! warm beam
           CALL GetMaxwellVelocity(vx)
           vx = MAX(0.0_8, whole_object(nwo)%v_ebeam_constant_emit + vx * whole_object(nwo)%factor_convert_vinj_normal_constant_emit)
        END IF
        CALL GetMaxwellVelocity(vy)
        CALL GetMaxwellVelocity(vz)
        vy = vy * whole_object(nwo)%factor_convert_vinj_parallel_constant_emit
        vz = vz * whole_object(nwo)%factor_convert_vinj_parallel_constant_emit
        tag = nwo !0

        CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
     END DO

     whole_object(nwo)%electron_emit_count = whole_object(nwo)%electron_emit_count + add_N_e_to_emit

  END DO

! boundary objects along the top edge of the cluster
  DO n = 1, c_N_of_local_object_parts_above
     m = c_index_of_local_object_part_above(n)
     add_N_e_to_emit = 0
     nwo = c_local_object_part(m)%object_number
     IF (whole_object(nwo)%N_electron_constant_emit.LE.0) CYCLE
     IF (c_right_top_corner_type.EQ.FLAT_WALL_ABOVE) THEN
! account for overlapping
        add_N_e_to_emit = INT( whole_object(nwo)%N_electron_constant_emit * REAL(MIN(c_local_object_part(m)%iend, c_indx_x_max-1) - c_local_object_part(m)%istart ) / REAL(whole_object(nwo)%L) )
     ELSE
        add_N_e_to_emit = INT( whole_object(nwo)%N_electron_constant_emit * REAL(c_local_object_part(m)%iend - c_local_object_part(m)%istart) / REAL(whole_object(nwo)%L) )
     END IF

! account for emission split between multiple processes
     IF (Rank_cluster.EQ.N_processes_cluster-1) THEN
        add_N_e_to_emit = add_N_e_to_emit - (N_processes_cluster-1) * (add_N_e_to_emit / N_processes_cluster)
     ELSE
        add_N_e_to_emit = add_N_e_to_emit / N_processes_cluster
     END IF

! emission
     DO k = 1, add_N_e_to_emit
        IF (c_right_top_corner_type.EQ.FLAT_WALL_ABOVE) THEN  !Rank_of_master_right.GE.0) THEN
           x = DBLE(c_local_object_part(m)%istart) + well_random_number() * DBLE( MIN(c_local_object_part(m)%iend, c_indx_x_max-1) - c_local_object_part(m)%istart )
           x = MIN(MAX(x, DBLE(c_indx_x_min)), DBLE(c_indx_x_max-1))
        ELSE
! cluster which has no overlapping from right at the top right corner
           x = DBLE(c_local_object_part(m)%istart) + well_random_number() * DBLE( MIN(c_local_object_part(m)%iend, c_indx_x_max) - c_local_object_part(m)%istart )
           x = MIN(MAX(x, DBLE(c_indx_x_min)), DBLE(c_indx_x_max))
        END IF
        y = DBLE(c_indx_y_max) - 1.0d-6   !???

        IF (whole_object(nwo)%model_constant_emit.EQ.0) THEN
! thermal emission
           CALL GetInjMaxwellVelocity(vy)
           vy = -vy * whole_object(nwo)%factor_convert_vinj_normal_constant_emit
        ELSE
! warm beam
           CALL GetMaxwellVelocity(vy)
           vy = -MAX(0.0_8, whole_object(nwo)%v_ebeam_constant_emit + vy * whole_object(nwo)%factor_convert_vinj_normal_constant_emit)
        END IF
        CALL GetMaxwellVelocity(vx)
        CALL GetMaxwellVelocity(vz)
        vx = vx * whole_object(nwo)%factor_convert_vinj_parallel_constant_emit
        vz = vz * whole_object(nwo)%factor_convert_vinj_parallel_constant_emit
        tag = nwo !0

        CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
     END DO

     whole_object(nwo)%electron_emit_count = whole_object(nwo)%electron_emit_count + add_N_e_to_emit

  END DO

! boundary objects along the right edge of the cluster
  DO n = 1, c_N_of_local_object_parts_right
     m = c_index_of_local_object_part_right(n)
     add_N_e_to_emit = 0
     nwo = c_local_object_part(m)%object_number
     IF (whole_object(nwo)%N_electron_constant_emit.LE.0) CYCLE
     IF (c_right_top_corner_type.EQ.FLAT_WALL_RIGHT) THEN
! account for overlapping
        add_N_e_to_emit = INT( whole_object(nwo)%N_electron_constant_emit * REAL(MIN(c_local_object_part(m)%jend, c_indx_y_max-1) - c_local_object_part(m)%jstart ) / REAL(whole_object(nwo)%L) )
     ELSE
        add_N_e_to_emit = INT( whole_object(nwo)%N_electron_constant_emit * REAL(c_local_object_part(m)%jend - c_local_object_part(m)%jstart ) / REAL(whole_object(nwo)%L) )
     END IF

! account for emission split between multiple processes
     IF (Rank_cluster.EQ.N_processes_cluster-1) THEN
        add_N_e_to_emit = add_N_e_to_emit - (N_processes_cluster-1) * (add_N_e_to_emit / N_processes_cluster)
     ELSE
        add_N_e_to_emit = add_N_e_to_emit / N_processes_cluster
     END IF

! emission
     DO k = 1, add_N_e_to_emit
        IF (c_right_top_corner_type.EQ.FLAT_WALL_RIGHT) THEN 
           y = DBLE(c_local_object_part(m)%jstart) + well_random_number() * DBLE( MIN(c_local_object_part(m)%jend, c_indx_y_max-1) - c_local_object_part(m)%jstart )
           y = MIN(MAX(y, DBLE(c_indx_y_min)), DBLE(c_indx_y_max-1))
        ELSE
! cluster which has no overlapping from above at the top left corner
           y = DBLE(c_local_object_part(m)%jstart) + well_random_number() * DBLE( c_local_object_part(m)%jend - c_local_object_part(m)%jstart )
           y = MIN(MAX(y, DBLE(c_indx_y_min)), DBLE(c_indx_y_max))
        END IF
        x = DBLE(c_indx_x_max) - 1.0d-6   !???

        IF (whole_object(nwo)%model_constant_emit.EQ.0) THEN
! thermal emission
           CALL GetInjMaxwellVelocity(vx)
           vx = -vx * whole_object(nwo)%factor_convert_vinj_normal_constant_emit
        ELSE
! warm beam
           CALL GetMaxwellVelocity(vx)
           vx = -MAX(0.0_8, whole_object(nwo)%v_ebeam_constant_emit + vx * whole_object(nwo)%factor_convert_vinj_normal_constant_emit)
        END IF
        CALL GetMaxwellVelocity(vy)
        CALL GetMaxwellVelocity(vz)
        vy = vy * whole_object(nwo)%factor_convert_vinj_parallel_constant_emit
        vz = vz * whole_object(nwo)%factor_convert_vinj_parallel_constant_emit
        tag = nwo !0

        CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
     END DO

     whole_object(nwo)%electron_emit_count = whole_object(nwo)%electron_emit_count + add_N_e_to_emit

  END DO

! boundary objects along the bottom edge of the cluster
  DO n = 1, c_N_of_local_object_parts_below
     m = c_index_of_local_object_part_below(n)
     add_N_e_to_emit = 0
     nwo = c_local_object_part(m)%object_number
     IF (whole_object(nwo)%N_electron_constant_emit.LE.0) CYCLE
     IF (c_right_bottom_corner_type.EQ.FLAT_WALL_BELOW) THEN
! account for overlapping
        add_N_e_to_emit = INT( whole_object(nwo)%N_electron_constant_emit * REAL(MIN(c_local_object_part(m)%iend, c_indx_x_max-1) - c_local_object_part(m)%istart ) / REAL(whole_object(nwo)%L) )
     ELSE
        add_N_e_to_emit = INT( whole_object(nwo)%N_electron_constant_emit * REAL(c_local_object_part(m)%iend - c_local_object_part(m)%istart ) / REAL(whole_object(nwo)%L) )
     END IF

! account for emission split between multiple processes
     IF (Rank_cluster.EQ.N_processes_cluster-1) THEN
        add_N_e_to_emit = add_N_e_to_emit - (N_processes_cluster-1) * (add_N_e_to_emit / N_processes_cluster)
     ELSE
        add_N_e_to_emit = add_N_e_to_emit / N_processes_cluster
     END IF

! emission
     DO k = 1, add_N_e_to_emit
        IF (c_right_bottom_corner_type.EQ.FLAT_WALL_BELOW) THEN 
           x = DBLE(c_local_object_part(m)%istart) + well_random_number() * DBLE( MIN(c_local_object_part(m)%iend, c_indx_x_max-1) - c_local_object_part(m)%istart )
           x = MIN(MAX(x, DBLE(c_indx_x_min)), DBLE(c_indx_x_max-1))
        ELSE
! cluster which has no overlapping from right at the top right corner
           x = DBLE(c_local_object_part(m)%istart) + well_random_number() * DBLE( c_local_object_part(m)%iend - c_local_object_part(m)%istart )
           x = MIN(MAX(x, DBLE(c_indx_x_min)), DBLE(c_indx_x_max))
        END IF
        y = DBLE(c_indx_y_min) + 1.0d-6   !???

        IF (whole_object(nwo)%model_constant_emit.EQ.0) THEN
! thermal emission
           CALL GetInjMaxwellVelocity(vy)
           vy = vy * whole_object(nwo)%factor_convert_vinj_normal_constant_emit
        ELSE
! warm beam
           CALL GetMaxwellVelocity(vy)
           vy = MAX(0.0_8, whole_object(nwo)%v_ebeam_constant_emit + vy * whole_object(nwo)%factor_convert_vinj_normal_constant_emit)
        END IF
        CALL GetMaxwellVelocity(vx)
        CALL GetMaxwellVelocity(vz)
        vx = vx * whole_object(nwo)%factor_convert_vinj_parallel_constant_emit
        vz = vz * whole_object(nwo)%factor_convert_vinj_parallel_constant_emit
        tag = nwo !0

        CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
     END DO

     whole_object(nwo)%electron_emit_count = whole_object(nwo)%electron_emit_count + add_N_e_to_emit

  END DO

! save number of emitted particles, similar to COLLECT_ELECTRON_BOUNDARY_HITS  ???? make it a separate routine???

  ALLOCATE(ibuf_send(1:N_of_boundary_objects), STAT = ALLOC_ERR)
  ALLOCATE(ibuf_receive(1:N_of_boundary_objects), STAT = ALLOC_ERR)

! each cluster adjacent to a boundary assembles electron-boundary hit counters from all cluster members in the master of the cluster

  ibuf_send(1:N_of_boundary_objects) = whole_object(1:N_of_boundary_objects)%electron_emit_count
  ibuf_receive = 0

  CALL MPI_REDUCE(ibuf_send, ibuf_receive, N_of_boundary_objects, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)

  IF ((cluster_rank_key.EQ.0).AND.(c_N_of_local_object_parts.GT.0)) THEN

     whole_object(1:N_of_boundary_objects)%electron_emit_count = ibuf_receive(1:N_of_boundary_objects)

     ibuf_send(1:N_of_boundary_objects) = whole_object(1:N_of_boundary_objects)%electron_emit_count
     ibuf_receive = 0
! 
     CALL MPI_REDUCE(ibuf_send, ibuf_receive, N_of_boundary_objects, MPI_INTEGER, MPI_SUM, 0, COMM_BOUNDARY, ierr)

     IF (Rank_of_process.EQ.0) THEN
        whole_object(1:N_of_boundary_objects)%electron_emit_count = ibuf_receive(1:N_of_boundary_objects)

print '("electrons emitted by boundaries :: ",10(2x,i8))', whole_object(1:N_of_boundary_objects)%electron_emit_count  

!! set the ion hit counters here to zero because when this subroutine is called the ions do not move
!           DO k = 1, N_of_boundary_objects
!              whole_object(k)%ion_hit_count(1:N_spec) = 0
!           END DO

     END IF
  END IF

  DEALLOCATE(ibuf_send, STAT = ALLOC_ERR)
  DEALLOCATE(ibuf_receive, STAT = ALLOC_ERR)

!print '("process ",i4," of cluster ",i4," emitted ",i5," electrons from top boundary")', Rank_of_process, particle_master, add_N_e_to_emit

END SUBROUTINE PERFORM_ELECTRON_EMISSION_SETUP
