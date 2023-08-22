
!--------------------

SUBROUTINE INITIATE_SNAPSHOTS

  USE ParallelOperationValues
  USE Snapshots
  USE Diagnostics,          ONLY : Save_probes_data_T_cntr_rff, Save_probes_data_step
  USE CurrentProblemValues, ONLY : delta_t_s, Max_T_cntr
  USE Checkpoints, ONLY : use_checkpoint, current_snap_check
  USE MCCollisions
  USE ClusterAndItsBoundaries

  use mpi

  IMPLICIT NONE

  INTEGER ierr

  INTEGER T2_old, N2_old
  LOGICAL exists
  CHARACTER (1) buf
  INTEGER flag

  INTEGER saveflagi(1:38)           ! integer flags used to set values of logical flags controlling saving of data files

  INTEGER N_of_snap_groups          ! number of sets of snapshots, read from file
  INTEGER i                         ! set of snapshots

  REAL(8) Rqst_snap_start_ns        ! requested start of current set of snapshots [ns], read from file
  REAL(8) Rqst_snap_finish_ns       ! requested finish of current set of snapshots [ns], read from file 
  INTEGER Rqst_n_of_snaps           ! requested number of snapshots in current set, read from file
  INTEGER Rqst_evdf_flag            ! requested flag defining which velocity distribution functions to save (0/1/2/3 = No/1d only/2d only/1d and 2d)
  INTEGER Rqst_pp_flag              ! requested flag defining which phase planes to save (0/1/2/3 = No/electrons only/ions only/electrons and ions)
  INTEGER Rqst_ionrate2d_flag       ! requested flag defining whether to save ionization rates (0/1 = No/Yes)
  INTEGER Rqst_part_coll_walls      ! requested flag defining whether to save particles collided with walls (0/1/2/3 = No/electrons only/ions only/electrons and ions)

  INTEGER T1, T2, N1, N2 

  INTEGER large_step 

  INTEGER Fact_n_of_snaps           ! calculated number of snapshots in one set

  INTEGER n                             ! ordering number of snapshot in the set
  INTEGER, ALLOCATABLE ::           timestep(:)   ! array for temporary storage of moments (timesteps) of snapshots
  INTEGER, ALLOCATABLE :: evdf_flag_timestep(:)   ! array for temporary storage of vdf save flags
  INTEGER, ALLOCATABLE ::   pp_flag_timestep(:)   ! array for temporary storage of phase planes save flags
  INTEGER, ALLOCATABLE :: ionrate2d_flag_timestep(:)   ! array for temporary storage of ionization rate save flags
  INTEGER, ALLOCATABLE :: part_coll_walls_flag_timestep(:)   ! array for temporary storage of particles collided with walls save flags

  INTEGER ALLOC_ERR

  INTEGER p

! default values
  N_of_all_snaps = 0
  save_data = .FALSE.
  current_snap = 1

! create structure for collision diagnostics arrays, to be stored between snapshots by cluster masters only. I need this because it can be used by averaged snapshots as well!
  IF ( .NOT.en_collisions_turned_off .AND. cluster_rank_key==0 ) THEN
   ALLOCATE(diagnostics_neutral(1:N_neutral_spec), STAT = ALLOC_ERR)
   DO n = 1, N_neutral_spec
      IF (collision_e_neutral(n)%N_of_activated_colproc.LE.0) CYCLE
      ALLOCATE(diagnostics_neutral(n)%activated_collision(1:collision_e_neutral(n)%N_of_activated_colproc), STAT = ALLOC_ERR)
   END DO   
  END IF

! read / write the data file 
  INQUIRE (FILE = 'init_snapshots.dat', EXIST = exists)

  IF (.NOT.exists) THEN
     IF (Rank_of_process.EQ.0) PRINT '(2x,"### file init_snapshots.dat not found. Snapshots will NOT be created ... ###")'
     RETURN
  END IF

  T2_old = -1
  N2_old = 0

  ALLOCATE(                     timestep(1:9999), STAT = ALLOC_ERR)
  ALLOCATE(           evdf_flag_timestep(1:9999), STAT = ALLOC_ERR)
  ALLOCATE(             pp_flag_timestep(1:9999), STAT = ALLOC_ERR)
  ALLOCATE(      ionrate2d_flag_timestep(1:9999), STAT = ALLOC_ERR)
  ALLOCATE(part_coll_walls_flag_timestep(1:9999), STAT = ALLOC_ERR)
     
  IF (Rank_of_process.EQ.0) PRINT '("### File init_snapshots is found. Reading the data file... ###")'

  OPEN (9, FILE = 'init_snapshots.dat')
     
  saveflagi = 0

  READ(9, '(A1)') buf  !--- save 2D maps of the following parameters? (1=yes, 0=no)
  READ(9, '(A1)') buf  !-----F----EX----EY--JXsum--JYsum--JZsum
  READ(9, '(A1)') buf  !----dd----dd----dd----dd----dd----dd
  READ(9, '(6(4x,i2))') saveflagi(1:6) 
  READ(9, '(A1)') buf  !----Ne----JXe---JYe---JZe---VXe---VYe---VZe---WXe---WYe---WZe---TXe---TYe---TZe---QXe---QYe---QZe
  READ(9, '(A1)') buf  !----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd
  READ(9, '(16(4x,i2))') saveflagi(7:22)
  READ(9, '(A1)') buf  !----Ni----JXi---JYi---JZi---VXi---VYi---VZi---WXi---WYi---WZi---TXi---TYi---TZi---QXi---QYi---QZi
  READ(9, '(A1)') buf  !----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd
  READ (9, '(16(4x,i2))') saveflagi(23:38)

  save_data = .TRUE.
  DO i = 1, 38
     IF (saveflagi(i).EQ.0) save_data(i) = .FALSE.
  END DO

  READ (9, '(A1)') buf !---dd--- Number of groups of snapshots ( >= 0 )
  READ (9, '(3x,i2)') N_of_snap_groups    ! below:: save VDFs (0/1/2/3 = No/1d/2d/1d+2d) / save phase planes (0/1/2/3 = No/e/i/e+i) / save ionization rates (0/1 = No/Yes) / save particles collided with walls (0/1/2/3 = No/e/i/e+i)
  READ (9, '(A1)') buf !---ddddddd.ddd---ddddddd.ddd---dddd---d---d---d---d--- group: start (ns) / finish (ns) / number of snapshots / VDFs / phase planes / ionization rates / particles collided with walls

  DO i = 1, N_of_snap_groups
! read the parameters of current set of snapshot from the data file
     READ (9, '(3x,f11.3,3x,f11.3,3x,i4,4(3x,i1))') Rqst_snap_start_ns, Rqst_snap_finish_ns, Rqst_n_of_snaps, Rqst_evdf_flag, Rqst_pp_flag, Rqst_ionrate2d_flag, Rqst_part_coll_walls
! try the next group of snapshots if the current group snapshot number is zero
     IF (Rqst_n_of_snaps.LT.1) CYCLE
! get the timestep, coinciding with the diagnostic output timestep and closest to Rqst_snap_start_ns
     T1 = Rqst_snap_start_ns / (delta_t_s * 1.0d9)
     IF (T1.LT.Save_probes_data_T_cntr_rff) T1 = Save_probes_data_T_cntr_rff
     N1 = (T1 - Save_probes_data_T_cntr_rff) / Save_probes_data_step
     T1 = Save_probes_data_T_cntr_rff + N1 * Save_probes_data_step
! get the timestep, coinciding with the diagnostic output timestep and closest to Rqst_snap_finish_ns
     T2 = Rqst_snap_finish_ns / (delta_t_s * 1.0d9)
     IF (T2.LT.Save_probes_data_T_cntr_rff) T2 = Save_probes_data_T_cntr_rff
     N2 = (T2 - Save_probes_data_T_cntr_rff) / Save_probes_data_step
     T2 = Save_probes_data_T_cntr_rff + N2 * Save_probes_data_step 
! adjust, if necessary, the start timesteps if it is before (less than) the finish timestep for the previous set of snapshots T2_old
     IF (T1.LE.T2_old) THEN 
        T1 = T2_old + Save_probes_data_step
        N1 = N2_old + 1
     END IF
! adjust, if necessary, the finish timestep if it is after (larger than) the last simulation timestep
     DO WHILE (T2.GT.Max_T_cntr)
        T2 = T2 - Save_probes_data_step
        N2 = N2 - 1
     END DO
! skip the line if the start moment is after (larger than) the finish moment
     IF (T1.GT.T2) CYCLE
! if we are here then T1 and T2 can be used for calculation of moments of snapshots
! calculate the number of snapshots which can be made in current set
     IF (Rqst_n_of_snaps.EQ.1) THEN
        large_step = 0
        Fact_n_of_snaps = 1
        T2 = T1
        N2 = N1
     ELSE
        large_step = (N2 - N1) / (Rqst_n_of_snaps - 1)
        IF (large_step.EQ.0) THEN 
           large_step = 1
           Fact_n_of_snaps = N2 - N1 + 1
        ELSE
           Fact_n_of_snaps = Rqst_n_of_snaps
           N2 = N1 + large_step * (Fact_n_of_snaps - 1)
           T2 = Save_probes_data_T_cntr_rff + N2 * Save_probes_data_step 
        END IF
     END IF
! save the final moment of the current snapshot set
     T2_old = T2
     N2_old = N2
! for all possible snapshots of the current set
     DO n = 1, Fact_n_of_snaps
! Calculate and save the snapshot moment in the temporary array
        N_of_all_snaps = N_of_all_snaps + 1
        timestep(N_of_all_snaps) =  T1 + (n - 1) * large_step * Save_probes_data_step
        evdf_flag_timestep(N_of_all_snaps) = MAX(0,MIN(3,Rqst_evdf_flag))
        pp_flag_timestep(N_of_all_snaps) = MAX(0,MIN(3,Rqst_pp_flag))
        ionrate2d_flag_timestep(N_of_all_snaps) = MAX(0,MIN(1,Rqst_ionrate2d_flag))
        part_coll_walls_flag_timestep(N_of_all_snaps) = MAX(0,MIN(3,Rqst_part_coll_walls))
     END DO        ! end of cycle over snapshots in one set
  END DO           ! end of cycle over sets of snapshots     

! for calculation of velocity distribution functions inside the plasma volume, 
! the simulation domain of each cluster will be split into sub-domains called boxes
! maximal velocity in each direction for electrons will be N_max_vel electron thermal velocities
! maximal velocity in each direction for ion species s will be N_max_vel / sqrt(Ms(s)) electron thermal velocities

  READ (9, '(A1)') buf !-------- Parameters for calculation of velocity distribution functions
  READ (9, '(A1)') buf !--d-d--- Number of spatial boxes along the X / Y direction in a cluster (>=0)
  READ (9, '(2x,i1,1x,i1)') N_vdfbox_x, N_vdfbox_y
  READ (9, '(A1)') buf !--ddd--- Electrons, maximal velocity [in units of v_Te_ms]
  READ (9, '(2x,i3)') N_max_vel_e
  READ (9, '(A1)') buf !--ddd--- Electrons, number of velocity bins per v_Te_ms
  READ (9, '(2x,i3)') N_vbins_e
  READ (9, '(A1)') buf !--ddd--- Ions, maximal velocity [in units of v_Te_ms*sqrt(me/Ms)]
  READ (9, '(2x,i3)') N_max_vel_i
  READ (9, '(A1)') buf !--ddd--- Number of velocity bins per v_Te_ms*sqrt(me/Ms) for ions
  READ (9, '(2x,i3)') N_vbins_i

  READ (9, '(A1)') buf !-------- Parameters for saving phase planes
  READ (9, '(A1)') buf !--dd---- Number of rectangular spatial boxes (>=0)
  READ (9, '(2x,i2)') N_pp_boxes
  ALLOCATE(pp_box(1:N_pp_boxes), STAT=ALLOC_ERR)
  READ (9, '(A1)') buf !--dddd--dddd----dddd--dddd--- left bottom corner X/Y (node index) / right top corner X/Y (node index)
  DO n = 1, N_pp_boxes
     READ (9, '(2x,i4,2x,i4,4x,i4,2x,i4)') pp_box(n)%imin, pp_box(n)%jmin, pp_box(n)%imax, pp_box(n)%jmax 
  END DO

  CLOSE (9, STATUS = 'KEEP')

!  current_snap = 1   ! default value

! overrite if the system is initialized using a checkpoint
  IF (use_checkpoint.EQ.1) current_snap = current_snap_check

! report about the general status of snapshot creation
  IF (N_of_all_snaps.EQ.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### Snapshots will NOT be created ... ###")'
     RETURN
  END IF 

! if we are here, snapshots will be created ...

! index limits for velocity bins array
  indx_v_min_e = -N_max_vel_e * N_vbins_e
  indx_v_max_e =  N_max_vel_e * N_vbins_e - 1

  indx_v_min_i = -N_max_vel_i * N_vbins_i
  indx_v_max_i =  N_max_vel_i * N_vbins_i - 1

! report about the general status of saving the velocity distribution functions
  N_vdfbox_all = N_vdfbox_x * N_vdfbox_y
  IF (N_vdfbox_all.EQ.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### Velocity distributions will NOT be saved ... ###")'
  ELSE
     IF (Rank_of_process.EQ.0) PRINT '("### Each cluster will create velocity distributions in ",i2," location(s) ###")', N_vdfbox_all
  END IF

! allocate the array of moments of snapshots
  ALLOCATE(Tcntr_snapshot(1:N_of_all_snaps), STAT=ALLOC_ERR)
! allocate control arrays
  ALLOCATE(        save_evdf_snapshot(1:N_of_all_snaps), STAT=ALLOC_ERR)
  ALLOCATE(          save_pp_snapshot(1:N_of_all_snaps), STAT=ALLOC_ERR)
  ALLOCATE(  save_ionization_rates_2d(1:N_of_all_snaps), STAT=ALLOC_ERR)
  ALLOCATE(save_ions_collided_with_bo(1:N_of_all_snaps), STAT=ALLOC_ERR)
  ALLOCATE(   save_e_collided_with_bo(1:N_of_all_snaps), STAT=ALLOC_ERR)

! move the calculated snapshot moments from the temporary array to the allocated array 
      Tcntr_snapshot(1:N_of_all_snaps) =           timestep(1:N_of_all_snaps)
! copy control flags
  save_evdf_snapshot(1:N_of_all_snaps) = evdf_flag_timestep(1:N_of_all_snaps)
    save_pp_snapshot(1:N_of_all_snaps) =   pp_flag_timestep(1:N_of_all_snaps)
! set logical switches
  save_ionization_rates_2d = .FALSE.
  save_ions_collided_with_bo = .FALSE.
  save_e_collided_with_bo = .FALSE.
  DO i = 1, N_of_all_snaps
     IF (ionrate2d_flag_timestep(i).GT.0) save_ionization_rates_2d(i) = .TRUE.
     SELECT CASE (part_coll_walls_flag_timestep(i))
        CASE (1)
           save_e_collided_with_bo(i) = .TRUE.
        CASE (2)
           save_ions_collided_with_bo(i) = .TRUE.
        CASE (3)
           save_e_collided_with_bo(i) = .TRUE.
           save_ions_collided_with_bo(i) = .TRUE.
     END SELECT    
  END DO
 
  DEALLOCATE(               timestep, STAT = ALLOC_ERR)
  DEALLOCATE(     evdf_flag_timestep, STAT = ALLOC_ERR)
  DEALLOCATE(       pp_flag_timestep, STAT = ALLOC_ERR)
  DEALLOCATE(ionrate2d_flag_timestep, STAT = ALLOC_ERR)
  DEALLOCATE(part_coll_walls_flag_timestep, STAT = ALLOC_ERR)
 
  IF (Rank_of_process.EQ.0) THEN 
     PRINT '("### The program will create ",i4," snapshots ###")', N_of_all_snaps

! write moments of snapshot creation into the file
     OPEN (41, FILE = '_snapmoments.dat')
!                 "--****-----*******.*****----********----*----*----*----*----*"
     WRITE (41, '(" number       time(ns)       T_cntr    vdf  pp  ioniz icbo ecbo ")')
     DO i = 1, N_of_all_snaps
        WRITE (41, '(2x,i4,5x,f13.5,4x,i8,4x,i1,4x,i1,4x,L1,4x,L1,4x,L1)') &
             & i, &
             & Tcntr_snapshot(i) * 1.0d9 * delta_t_s, &
             & Tcntr_snapshot(i), &
             & save_evdf_snapshot(i), &
             & save_pp_snapshot(i), &
             & save_ionization_rates_2d(i), &
             & save_ions_collided_with_bo(i), &
             & save_e_collided_with_bo(i)
     END DO
     CLOSE (41, STATUS = 'KEEP')
  END IF

  IF (en_collisions_turned_off) RETURN
  IF (cluster_rank_key.NE.0) RETURN

  IF (no_ionization_collisions) RETURN
  IF (.NOT.save_ionization_rates_2d(current_snap))RETURN

! create ionization rates diagnostics arrays, to be stored between snapshots by cluster masters only
  DO n = 1, N_neutral_spec
     DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
        IF (collision_e_neutral(n)%colproc_info(p)%type.LT.30) CYCLE
        ALLOCATE(diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local = 0.0
     END DO
  END DO
  
END SUBROUTINE INITIATE_SNAPSHOTS

!------------------------------------------------------------------------------------------------------------
! 
SUBROUTINE CREATE_SNAPSHOT

  USE ParallelOperationValues
  USE Snapshots
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries
  USE IonParticles, ONLY : N_spec, Qs, Ms
  USE MCCollisions

  use mpi

  IMPLICIT NONE

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER ALLOC_ERR
  INTEGER i, j, k
  INTEGER recsize, bufsize
  REAL, ALLOCATABLE :: rbufer(:)

  REAL, ALLOCATABLE :: cs_temp(:,:)

  REAL, ALLOCATABLE :: cs_JXsum(:,:)
  REAL, ALLOCATABLE :: cs_JYsum(:,:)
  REAL, ALLOCATABLE :: cs_JZsum(:,:)

  INTEGER pos1, pos2

!                                 ! ----x----I----x----I
  CHARACTER(16) filename_F        ! _NNNN_F_V_2D.bin
  CHARACTER(18) filename_E        ! _NNNN_EX_Vm_2D.bin
!                                 ! ----x----I----x----I
  CHARACTER(18) filename_Ne       ! _NNNN_Ne_m3_2D.bin
  CHARACTER(20) filename_Je       ! _NNNN_JXe_Am2_2D.bin
  CHARACTER(19) filename_Ve       ! _NNNN_VXe_ms_2D.bin
  CHARACTER(19) filename_We       ! _NNNN_WXe_eV_2D.bin
  CHARACTER(19) filename_Te       ! _NNNN_TXe_eV_2D.bin
  CHARACTER(20) filename_Qe       ! _NNNN_QXe_Wm2_2D.bin
!                                 ! ----x----I----x----I--
  CHARACTER(20) filename_Ni       ! _NNNN_Ni_s_m3_2D.bin
  CHARACTER(22) filename_Ji       ! _NNNN_JXi_s_Am2_2D.bin
  CHARACTER(21) filename_Vi       ! _NNNN_VXi_s_ms_2D.bin
  CHARACTER(21) filename_Wi       ! _NNNN_WXi_s_eV_2D.bin
  CHARACTER(21) filename_Ti       ! _NNNN_TXi_s_eV_2D.bin
  CHARACTER(22) filename_Qi       ! _NNNN_QXi_s_Wm2_2D.bin
!                                 ! ----x----I----x----I--
  CHARACTER(22) filename_Jsum     ! _NNNN_JXsum_Am2_2D.bin

  INTEGER s
  INTEGER n, p

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

! quit if all snapshots were created or if due to any reasons the snapshot counter is 
! larger than the declared number of snapshots (e.g., when no snapshots are requested) 
  IF (current_snap.GT.N_of_all_snaps) RETURN

! quit if the current moment of time is not the moment when it is necessary to create the snapshot
  IF (T_cntr.NE.Tcntr_snapshot(current_snap)) RETURN

  IF (Rank_of_process.EQ.0) PRINT '("### ^^^^^^^^^^^^^^^ Snapshot ",i4," will be created now ... ^^^^^^^^^^^^^^^^ ###")', current_snap

  IF (cluster_rank_key.EQ.0) ALLOCATE(cs_temp(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)

! potential

  IF (save_data(1)) THEN

     IF (cluster_rank_key.EQ.0) THEN

        IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN
! PETSc-based field solver
! account for own contribution of the master as a field calculator
           DO j = indx_y_min, indx_y_max
              DO i = indx_x_min, indx_x_max
                 cs_temp(i,j) = REAL(phi(i,j) * F_scale_V)
              END DO
           END DO

! account for contributions from other field calculators
           DO k = 2, cluster_N_blocks
              recsize = field_calculator(k)%indx_x_max - field_calculator(k)%indx_x_min + 1
              bufsize = recsize * (field_calculator(k)%indx_y_max - field_calculator(k)%indx_y_min + 1)
              ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)
              CALL MPI_RECV(rbufer, bufsize, MPI_REAL, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)
              pos2 = 0
              DO j = field_calculator(k)%indx_y_min, field_calculator(k)%indx_y_max
                 pos1 = pos2 + 1
                 pos2 = pos2 + recsize
                 cs_temp(field_calculator(k)%indx_x_min:field_calculator(k)%indx_x_max,j) = rbufer(pos1:pos2) 
              END DO
              DEALLOCATE(rbufer, STAT = ALLOC_ERR)
           END DO

        ELSE
! FFT-based field solver
! the cluster already knows the potential within its domain, use it -----------------
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_temp(i,j) = REAL(c_phi(i,j) * F_scale_V)
              END DO
           END DO
        END IF   !### IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

        filename_F = '_NNNN_F_V_2D.bin'
        filename_F(2:5) = convert_int_to_txt_string(current_snap, 4)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_F)

     ELSE   !### IF (cluster_rank_key.EQ.0) THEN

        IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

           recsize = indx_x_max - indx_x_min + 1
           bufsize = recsize * (indx_y_max - indx_y_min + 1)

           ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)

           pos1 = 1 - indx_x_min
           DO j = indx_y_min, indx_y_max
              DO i = indx_x_min, indx_x_max
                 rbufer(pos1 + i) = REAL(phi(i, j) * F_scale_V)
              END DO
              pos1 = pos1 + recsize
           END DO

           CALL MPI_SEND(rbufer, bufsize, MPI_REAL, field_master, Rank_of_process, MPI_COMM_WORLD, ierr) 

           DEALLOCATE(rbufer, STAT = ALLOC_ERR)
        
        END IF

     END IF   !###   IF (cluster_rank_key.EQ.0) THEN

  END IF   !###   IF (save_data(1)) THEN

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! electric field components

  IF (cluster_rank_key.EQ.0) THEN

     IF (save_data(2)) THEN
        filename_E = '_NNNN_EX_Vm_2D.bin'
        filename_E(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = REAL(EX(i,j) * E_scale_Vm)
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_E)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr) 
     END IF

     IF (save_data(3)) THEN
        filename_E = '_NNNN_EY_Vm_2D.bin'
        filename_E(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = REAL(EY(i,j) * E_scale_Vm)
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_E)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr) 
     END IF

     ALLOCATE(cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ALLOCATE(cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ALLOCATE(cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ALLOCATE(cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ALLOCATE(cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

  ELSE

! the moments of the distribution functions are calculated using MPI_REDUCE which takes sum of values from all particle calculators
! the results are stored in master processes in arrays cs_N, cs_VX, etc
! the particle calculators, in general, don't need these arrays at all, 
! so allocating these arrays in non-master processes is just a waste of memory
! however, the compiler reports an error if the code is compiled with -C flag (check everything)
! this can be avoided if at least some minimal size arrays are allocated in the non-master processes
     ALLOCATE(cs_N(1,1), STAT=ALLOC_ERR)

     ALLOCATE(cs_VX(1,1), STAT=ALLOC_ERR)
     ALLOCATE(cs_VY(1,1), STAT=ALLOC_ERR)
     ALLOCATE(cs_VZ(1,1), STAT=ALLOC_ERR)

     ALLOCATE(cs_WX(1,1), STAT=ALLOC_ERR)
     ALLOCATE(cs_WY(1,1), STAT=ALLOC_ERR)
     ALLOCATE(cs_WZ(1,1), STAT=ALLOC_ERR)

     ALLOCATE(cs_VXVY(1,1), STAT=ALLOC_ERR)
     ALLOCATE(cs_VXVZ(1,1), STAT=ALLOC_ERR)
     ALLOCATE(cs_VYVZ(1,1), STAT=ALLOC_ERR)

     ALLOCATE(cs_QX(1,1), STAT=ALLOC_ERR)
     ALLOCATE(cs_QY(1,1), STAT=ALLOC_ERR)
     ALLOCATE(cs_QZ(1,1), STAT=ALLOC_ERR)

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! collect electron and ion velocity distribution function moments (density, flows, energies)

! write electron moments to files ------------------------------------------

  IF ( save_data(4).OR. &
     & save_data(5).OR. &
     & save_data(6).OR. &
     & save_data(7).OR. &
     & save_data(8).OR. &
     & save_data(9).OR. &
     & save_data(10).OR. &
     & save_data(11).OR. &
     & save_data(12).OR. &
     & save_data(13).OR. &
     & save_data(14).OR. &
     & save_data(15).OR. &
     & save_data(16).OR. &
     & save_data(17).OR. &
     & save_data(18).OR. &
     & save_data(19).OR. &
     & save_data(20).OR. &
     & save_data(21).OR. &
     & save_data(22) )  CALL COLLECT_ELECTRON_MOMENTS

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (cluster_rank_key.EQ.0) THEN

! full electric current, electron contribution

     IF (save_data(4)) THEN
        ALLOCATE(cs_JXsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_JXsum(i,j) = -cs_N(i,j) * cs_VX(i,j)
           END DO
        END DO
     END IF

     IF (save_data(5)) THEN
        ALLOCATE(cs_JYsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_JYsum(i,j) = -cs_N(i,j) * cs_VY(i,j)
           END DO
        END DO
     END IF

     IF (save_data(6)) THEN
        ALLOCATE(cs_JZsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_JZsum(i,j) = -cs_N(i,j) * cs_VY(i,j)
           END DO
        END DO
     END IF

! number density

     IF (save_data(7)) THEN
        filename_Ne = '_NNNN_Ne_m3_2D.bin'
        filename_Ne(2:5) = convert_int_to_txt_string(current_snap, 4)
        cs_temp = cs_N * REAL(N_scale_part_m3)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ne)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

! electric current along each coordinate direction

     IF (save_data(8)) THEN
        filename_Je = '_NNNN_JXe_Am2_2D.bin'
        filename_Je(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = -cs_N(i,j) * cs_VX(i,j) * REAL(current_factor_Am2)              ! current_factor_Am2 = e_Cl * V_scale_ms * N_scale_part_m3 
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Je)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

     IF (save_data(9)) THEN
        filename_Je = '_NNNN_JYe_Am2_2D.bin'
        filename_Je(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = -cs_N(i,j) * cs_VY(i,j) * REAL(current_factor_Am2)   ! 
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Je)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

     IF (save_data(10)) THEN
        filename_Je = '_NNNN_JZe_Am2_2D.bin'
        filename_Je(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = -cs_N(i,j) * cs_VZ(i,j) * REAL(current_factor_Am2)   ! 
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Je)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

! average velocity along each coordinate direction

     IF (save_data(11)) THEN
        filename_Ve = '_NNNN_VXe_ms_2D.bin'
        filename_Ve(2:5) = convert_int_to_txt_string(current_snap, 4)
        cs_temp = cs_VX * REAL(V_scale_ms)                                                 ! V_scale_ms = N_max_vel * v_Te_ms
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ve)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

     IF (save_data(12)) THEN
        filename_Ve = '_NNNN_VYe_ms_2D.bin'
        filename_Ve(2:5) = convert_int_to_txt_string(current_snap, 4)
        cs_temp = cs_VY * REAL(V_scale_ms)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ve)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

     IF (save_data(13)) THEN
        filename_Ve = '_NNNN_VZe_ms_2D.bin'
        filename_Ve(2:5) = convert_int_to_txt_string(current_snap, 4)
        cs_temp = cs_VZ * REAL(V_scale_ms)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ve)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

! average energy of motion along each coordinate direction

     IF (save_data(14)) THEN
        filename_We = '_NNNN_WXe_eV_2D.bin'
        filename_We(2:5) = convert_int_to_txt_string(current_snap, 4)
        cs_temp = cs_WX * REAL(energy_factor_eV)                                           ! energy_factor_eV = 0.5_8 * m_e_kg * V_scale_ms**2 / e_Cl
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_We)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

     IF (save_data(15)) THEN
        filename_We = '_NNNN_WYe_eV_2D.bin'
        filename_We(2:5) = convert_int_to_txt_string(current_snap, 4)
        cs_temp = cs_WY * REAL(energy_factor_eV)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_We)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

     IF (save_data(16)) THEN
        filename_We = '_NNNN_WZe_eV_2D.bin'
        filename_We(2:5) = convert_int_to_txt_string(current_snap, 4)
        cs_temp = cs_WZ * REAL(energy_factor_eV)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_We)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

! temperature along each coordinate direction

     IF (save_data(17)) THEN
        filename_Te = '_NNNN_TXe_eV_2D.bin'
        filename_Te(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = 0.0
              IF (cs_N(i,j).LE.0.0) CYCLE
              cs_temp(i,j) = MAX(0.0, cs_WX(i,j) - cs_VX(i,j)**2) * REAL(temperature_factor_eV)   ! temperature_factor_eV = m_e_kg * V_scale_ms**2 / e_Cl
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Te)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

     IF (save_data(18)) THEN
        filename_Te = '_NNNN_TYe_eV_2D.bin'
        filename_Te(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = 0.0
              IF (cs_N(i,j).LE.0.0) CYCLE
              cs_temp(i,j) = MAX(0.0, cs_WY(i,j) - cs_VY(i,j)**2) * REAL(temperature_factor_eV)
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Te)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

     IF (save_data(19)) THEN
        filename_Te = '_NNNN_TZe_eV_2D.bin'
        filename_Te(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = 0.0
              IF (cs_N(i,j).LE.0.0) CYCLE
              cs_temp(i,j) = MAX(0.0, cs_WZ(i,j) - cs_VZ(i,j)**2) * REAL(temperature_factor_eV)
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Te)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

! heat flow along each coordinate direction

     IF (save_data(20)) THEN
        filename_Qe = '_NNNN_QXe_Wm2_2D.bin'
        filename_Qe(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = cs_N(i,j) * cs_QX(i,j) * REAL(heat_flow_factor_Wm2)   ! heat_flow_factor_Wm2 = 0.5_8 * m_e_kg * V_scale_ms**2 * N_scale_part_m3
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Qe)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

     IF (save_data(21)) THEN
        filename_Qe = '_NNNN_QYe_Wm2_2D.bin'
        filename_Qe(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = cs_N(i,j) * cs_QY(i,j) * REAL(heat_flow_factor_Wm2)
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Qe)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

     IF (save_data(22)) THEN
        filename_Qe = '_NNNN_QZe_Wm2_2D.bin'
        filename_Qe(2:5) = convert_int_to_txt_string(current_snap, 4)
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_temp(i,j) = cs_N(i,j) * cs_QZ(i,j) * REAL(heat_flow_factor_Wm2)
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Qe)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END IF

  END IF   !### IF (cluster_rank_key.EQ.0) THEN

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! write ion moments to files ------------------------------------------

  DO s = 1, N_spec

     IF ( save_data(4).OR. &
        & save_data(5).OR. &
        & save_data(6).OR. &
        & save_data(23).OR. &
        & save_data(24).OR. &
        & save_data(25).OR. &
        & save_data(26).OR. &
        & save_data(27).OR. &
        & save_data(28).OR. &
        & save_data(29).OR. &
        & save_data(30).OR. &
        & save_data(31).OR. &
        & save_data(32).OR. &
        & save_data(33).OR. &
        & save_data(34).OR. &
        & save_data(35).OR. &
        & save_data(36).OR. &
        & save_data(37).OR. &
        & save_data(38) ) CALL COLLECT_ION_MOMENTS(s)

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
     
     IF (cluster_rank_key.EQ.0) THEN

! full electric current, ion contribution

        IF (save_data(4)) THEN
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_JXsum(i,j) = cs_JXsum(i,j) + REAL(Qs(s)) * cs_N(i,j) * cs_VX(i,j)
              END DO
           END DO
        END IF

        IF (save_data(5)) THEN
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_JYsum(i,j) = cs_JYsum(i,j) + REAL(Qs(s)) * cs_N(i,j) * cs_VY(i,j)
              END DO
           END DO
        END IF

        IF (save_data(6)) THEN
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_JZsum(i,j) = cs_JZsum(i,j) + REAL(Qs(s)) * cs_N(i,j) * cs_VZ(i,j)
              END DO
           END DO
        END IF

! number density

        IF (save_data(23)) THEN
           filename_Ni = '_NNNN_Ni_s_m3_2D.bin'
           filename_Ni(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Ni(10:10) = convert_int_to_txt_string(s, 1)
           cs_temp = cs_N * REAL(N_scale_part_m3)
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ni)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

! electric current along each coordinate direction

        IF (save_data(24)) THEN
           filename_Ji = '_NNNN_JXi_s_Am2_2D.bin'
           filename_Ji(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Ji(11:11) = convert_int_to_txt_string(s, 1)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_temp(i,j) = cs_N(i,j) * cs_VX(i,j) * REAL(Qs(s) * current_factor_Am2)              ! current_factor_Am2 = e_Cl * V_scale_ms * N_scale_part_m3 
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ji)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

        IF (save_data(25)) THEN
           filename_Ji = '_NNNN_JYi_s_Am2_2D.bin'
           filename_Ji(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Ji(11:11) = convert_int_to_txt_string(s, 1)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_temp(i,j) = cs_N(i,j) * cs_VY(i,j) * REAL(Qs(s) * current_factor_Am2)
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ji)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

        IF (save_data(26)) THEN
           filename_Ji = '_NNNN_JZi_s_Am2_2D.bin'
           filename_Ji(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Ji(11:11) = convert_int_to_txt_string(s, 1)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_temp(i,j) = cs_N(i,j) * cs_VZ(i,j) * REAL(Qs(s) * current_factor_Am2)
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ji)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

! average velocity along each coordinate direction

        IF (save_data(27)) THEN
           filename_Vi = '_NNNN_VXi_s_ms_2D.bin'
           filename_Vi(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Vi(11:11) = convert_int_to_txt_string(s, 1)
           cs_temp = cs_VX * REAL(V_scale_ms)                                                 ! V_scale_ms = N_max_vel * v_Te_ms
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Vi)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

        IF (save_data(28)) THEN
           filename_Vi = '_NNNN_VYi_s_ms_2D.bin'
           filename_Vi(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Vi(11:11) = convert_int_to_txt_string(s, 1)
           cs_temp = cs_VY * REAL(V_scale_ms)
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Vi)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

        IF (save_data(29)) THEN
           filename_Vi = '_NNNN_VZi_s_ms_2D.bin'
           filename_Vi(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Vi(11:11) = convert_int_to_txt_string(s, 1)
           cs_temp = cs_VZ * REAL(V_scale_ms)
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Vi)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

! average energy of motion along each coordinate direction

        IF (save_data(30)) THEN
           filename_Wi = '_NNNN_WXi_s_eV_2D.bin'
           filename_Wi(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Wi(11:11) = convert_int_to_txt_string(s, 1)
           cs_temp = cs_WX * REAL(Ms(s) * energy_factor_eV)                                    ! energy_factor_eV = 0.5_8 * m_e_kg * V_scale_ms**2 / e_Cl
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Wi)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

        IF (save_data(31)) THEN
           filename_Wi = '_NNNN_WYi_s_eV_2D.bin'
           filename_Wi(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Wi(11:11) = convert_int_to_txt_string(s, 1)
           cs_temp = cs_WY * REAL(Ms(s) * energy_factor_eV)
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Wi)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

        IF (save_data(32)) THEN
           filename_Wi = '_NNNN_WZi_s_eV_2D.bin'
           filename_Wi(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Wi(11:11) = convert_int_to_txt_string(s, 1)
           cs_temp = cs_WZ * REAL(Ms(s) * energy_factor_eV)
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Wi)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

! temperature along each coordinate direction

        IF (save_data(33)) THEN
           filename_Ti = '_NNNN_TXi_s_eV_2D.bin'
           filename_Ti(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Ti(11:11) = convert_int_to_txt_string(s, 1)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_temp(i,j) = 0.0
                 IF (cs_N(i,j).LE.0.0) CYCLE
                 cs_temp(i,j) = MAX(0.0, cs_WX(i,j) - cs_VX(i,j)**2) * REAL(Ms(s) * temperature_factor_eV)   ! temperature_factor_eV = m_e_kg * V_scale_ms**2 / e_Cl
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ti)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

        IF (save_data(34)) THEN
           filename_Ti = '_NNNN_TYi_s_eV_2D.bin'
           filename_Ti(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Ti(11:11) = convert_int_to_txt_string(s, 1)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_temp(i,j) = 0.0
                 IF (cs_N(i,j).LE.0.0) CYCLE
                 cs_temp(i,j) = MAX(0.0, cs_WY(i,j) - cs_VY(i,j)**2) * REAL(Ms(s) * temperature_factor_eV)
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ti)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

        IF (save_data(35)) THEN
           filename_Ti = '_NNNN_TZi_s_eV_2D.bin'
           filename_Ti(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Ti(11:11) = convert_int_to_txt_string(s, 1)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_temp(i,j) = 0.0
                 IF (cs_N(i,j).LE.0.0) CYCLE
                 cs_temp(i,j) = MAX(0.0, cs_WZ(i,j) - cs_VZ(i,j)**2) * REAL(Ms(s) * temperature_factor_eV)
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Ti)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

! heat flow along each coordinate direction

        IF (save_data(36)) THEN
           filename_Qi = '_NNNN_QXi_s_Wm2_2D.bin'
           filename_Qi(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Qi(11:11) = convert_int_to_txt_string(s, 1)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_temp(i,j) = cs_N(i,j) * cs_QX(i,j) * REAL(Ms(s) * heat_flow_factor_Wm2)   ! heat_flow_factor_Wm2 = 0.5_8 * m_e_kg * V_scale_ms**2 * N_scale_part_m3
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Qi)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

        IF (save_data(37)) THEN
           filename_Qi = '_NNNN_QYi_s_Wm2_2D.bin'
           filename_Qi(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Qi(11:11) = convert_int_to_txt_string(s, 1)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_temp(i,j) = cs_N(i,j) * cs_QY(i,j) * REAL(Ms(s) * heat_flow_factor_Wm2)
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Qi)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

        IF (save_data(38)) THEN
           filename_Qi = '_NNNN_QZi_s_Wm2_2D.bin'
           filename_Qi(2:5) = convert_int_to_txt_string(current_snap, 4)
           filename_Qi(11:11) = convert_int_to_txt_string(s, 1)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_temp(i,j) = cs_N(i,j) * cs_QZ(i,j) * REAL(Ms(s) * heat_flow_factor_Wm2)
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Qi)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
        END IF

     END IF   ! ### IF (cluster_rank_key.EQ.0) THEN

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  END DO   !### DO s = 1, N_spec

! cleanup

  DEALLOCATE(cs_N, STAT = ALLOC_ERR)

  DEALLOCATE(cs_VX, STAT = ALLOC_ERR)
  DEALLOCATE(cs_VY, STAT = ALLOC_ERR)
  DEALLOCATE(cs_VZ, STAT = ALLOC_ERR)

  DEALLOCATE(cs_WX, STAT = ALLOC_ERR)
  DEALLOCATE(cs_WY, STAT = ALLOC_ERR)
  DEALLOCATE(cs_WZ, STAT = ALLOC_ERR)

  DEALLOCATE(cs_VXVY, STAT = ALLOC_ERR)
  DEALLOCATE(cs_VXVZ, STAT = ALLOC_ERR)
  DEALLOCATE(cs_VYVZ, STAT = ALLOC_ERR)

  DEALLOCATE(cs_QX, STAT = ALLOC_ERR)
  DEALLOCATE(cs_QY, STAT = ALLOC_ERR)
  DEALLOCATE(cs_QZ, STAT = ALLOC_ERR)

! FULL electric current (sum of electron and ion currents) along each coordinate direction
  IF (cluster_rank_key.EQ.0) THEN

     IF (save_data(4)) THEN
        filename_Jsum = '_NNNN_JXsum_Am2_2D.bin'
        filename_Jsum(2:5) = convert_int_to_txt_string(current_snap, 4)
        cs_temp = cs_JXsum * REAL(current_factor_Am2)                           ! current_factor_Am2 = e_Cl * V_scale_ms * N_scale_part_m3 
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Jsum)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
        DEALLOCATE(cs_JXsum, STAT = ALLOC_ERR)
     END IF

     IF (save_data(5)) THEN
        filename_Jsum = '_NNNN_JYsum_Am2_2D.bin'
        filename_Jsum(2:5) = convert_int_to_txt_string(current_snap, 4)
        cs_temp = cs_JYsum * REAL(current_factor_Am2)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Jsum)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr) 
! cleanup
        DEALLOCATE(cs_JYsum, STAT = ALLOC_ERR)
     END IF

     IF (save_data(6)) THEN
        filename_Jsum = '_NNNN_JZsum_Am2_2D.bin'
        filename_Jsum(2:5) = convert_int_to_txt_string(current_snap, 4)
        cs_temp = cs_JZsum * REAL(current_factor_Am2)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_temp, filename_Jsum)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr) 
! cleanup
        DEALLOCATE(cs_JZsum, STAT = ALLOC_ERR)
     END IF

! cleanup
     DEALLOCATE(cs_temp, STAT = ALLOC_ERR)

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! velocity distribution functions -----------------------

  IF (N_vdfbox_all.GT.0) THEN
     IF (cluster_rank_key.EQ.0) THEN
        IF ((save_evdf_snapshot(current_snap).EQ.ONLY1D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
           ALLOCATE( evxdf(indx_v_min_e:indx_v_max_e, 1:N_vdfbox_all), STAT=ALLOC_ERR)
           ALLOCATE( evydf(indx_v_min_e:indx_v_max_e, 1:N_vdfbox_all), STAT=ALLOC_ERR)
           ALLOCATE( evzdf(indx_v_min_e:indx_v_max_e, 1:N_vdfbox_all), STAT=ALLOC_ERR)
           ALLOCATE(isvxdf(indx_v_min_i:indx_v_max_i, 1:N_vdfbox_all, 1:N_spec), STAT=ALLOC_ERR)
           ALLOCATE(isvydf(indx_v_min_i:indx_v_max_i, 1:N_vdfbox_all, 1:N_spec), STAT=ALLOC_ERR)
           ALLOCATE(isvzdf(indx_v_min_i:indx_v_max_i, 1:N_vdfbox_all, 1:N_spec), STAT=ALLOC_ERR)
        END IF

        IF ((save_evdf_snapshot(current_snap).EQ.ONLY2D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
           ALLOCATE(evxvydf(indx_v_min_e:indx_v_max_e, indx_v_min_e:indx_v_max_e, 1:N_vdfbox_all), STAT=ALLOC_ERR)
        END IF
     ELSE
        IF ((save_evdf_snapshot(current_snap).EQ.ONLY1D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
           ALLOCATE( evxdf(1,1), STAT=ALLOC_ERR)
           ALLOCATE( evydf(1,1), STAT=ALLOC_ERR)
           ALLOCATE( evzdf(1,1), STAT=ALLOC_ERR)
           ALLOCATE(isvxdf(1,1,1), STAT=ALLOC_ERR)
           ALLOCATE(isvydf(1,1,1), STAT=ALLOC_ERR)
           ALLOCATE(isvzdf(1,1,1), STAT=ALLOC_ERR)
        END IF

        IF ((save_evdf_snapshot(current_snap).EQ.ONLY2D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
           ALLOCATE(evxvydf(1,1,1), STAT=ALLOC_ERR)
        END IF
     END IF
  END IF

  CALL CALCULATE_ELECTRON_VDF

  CALL MPI_BARRIER(COMM_CLUSTER, ierr)

  CALL CALCULATE_ION_VDF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL SAVE_ALL_VDF1D

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! cleanup

  IF (ALLOCATED(evxdf)) DEALLOCATE(evxdf, STAT=ALLOC_ERR)
  IF (ALLOCATED(evydf)) DEALLOCATE(evydf, STAT=ALLOC_ERR)
  IF (ALLOCATED(evzdf)) DEALLOCATE(evzdf, STAT=ALLOC_ERR)

  IF (ALLOCATED(isvxdf)) DEALLOCATE(isvxdf, STAT=ALLOC_ERR)
  IF (ALLOCATED(isvydf)) DEALLOCATE(isvydf, STAT=ALLOC_ERR)
  IF (ALLOCATED(isvzdf)) DEALLOCATE(isvzdf, STAT=ALLOC_ERR)

  CALL SAVE_ALL_VDF2D

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! cleanup
  IF (ALLOCATED(evxvydf)) DEALLOCATE(evxvydf, STAT=ALLOC_ERR)

  CALL SAVE_ELECTRON_PHASE_PLANES

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL SAVE_ION_PHASE_PLANES

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (Rank_of_process.EQ.0) PRINT '(/2x,"### ^^^^^^^^^^^^^^^^^^^^ Snapshot ",i4," completed :) ^^^^^^^^^^^^^^^^^^^ ###")', current_snap

  current_snap = current_snap + 1           ! increase the snapshots counter 

  IF (cluster_rank_key.NE.0) RETURN
  IF (en_collisions_turned_off) RETURN
  IF (no_ionization_collisions) RETURN

! the cleanup is performed below instead of the end of SAVE_en_COLLISIONS_2D because:
! (a) CREATE_SNAPSHOT is always called ###after### SAVE_en_COLLISIONS_2D (for the same snapshot number), so the data are saved already
! (b) when snapshots are taken at intervals less than N_subcycles, 
!     there may be no call of SAVE_en_COLLISIONS_2D between consecutive calls of CREATE_SNAPSHOT
!     then diagnostics_neutral may be still allocated and the code may attempt to allocate it again, causing an error
!
  IF (ALLOCATED(diagnostics_neutral)) THEN
     DO n = 1, N_neutral_spec
        IF (ALLOCATED(diagnostics_neutral(n)%activated_collision)) THEN
           DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
              IF (ALLOCATED(diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local)) DEALLOCATE(diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local, STAT=ALLOC_ERR)
           END DO
!           DEALLOCATE(diagnostics_neutral(n)%activated_collision, STAT=ALLOC_ERR)
        END IF
     END DO
!     DEALLOCATE(diagnostics_neutral, STAT=ALLOC_ERR)
  END IF

  IF (current_snap.GT.N_of_all_snaps) RETURN
  IF (.NOT.save_ionization_rates_2d(current_snap)) RETURN

! memory for ionization rate diagnostics arrays is allocated only
! (a) by cluster masters 
! (b) if the next snapshot (# current_snap) will be created 
! (c) if the next snapshot saves the ionization rate
! (d) if the ionization collisions are on

!  ALLOCATE(diagnostics_neutral(1:N_neutral_spec), STAT = ALLOC_ERR)
  DO n = 1, N_neutral_spec
     IF (collision_e_neutral(n)%N_of_activated_colproc.LE.0) CYCLE
!     ALLOCATE(diagnostics_neutral(n)%activated_collision(1:collision_e_neutral(n)%N_of_activated_colproc), STAT = ALLOC_ERR)
     DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
! presently do ionization collisions only
        IF (collision_e_neutral(n)%colproc_info(p)%type.LT.30) CYCLE
        ALLOCATE(diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local = 0.0
     END DO
  END DO

END SUBROUTINE CREATE_SNAPSHOT

!---------------------------------------------------------------------------------------------------
! this program is called only by master processes
!
SUBROUTINE SAVE_GLOBAL_2D_ARRAY(arr, filename)

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles
  USE Snapshots

  use mpi

  IMPLICIT NONE

  REAL arr(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max)

  CHARACTER*(*) filename

  INTEGER errcode,ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

!  INTEGER COMM_SAVE_CLUSTERS
!  INTEGER Rank_save_clusters
!  INTEGER N_processes_save_clusters

  INTEGER reclenx(1:N_clusters_x)
  INTEGER recleny(1:N_clusters_y)

  INTEGER col, row

  INTEGER file_handle

  INTEGER out_indx_x_min, out_indx_x_max, out_indx_y_min, out_indx_y_max

  INTEGER ALLOC_ERR
  REAL, ALLOCATABLE :: rbufer(:)

  INTEGER(kind=MPI_OFFSET_KIND) shifthead 

  INTEGER i, j

  IF (cluster_rank_key.NE.0) THEN
! just in case 
     PRINT '("Proc ",i4," Error in SAVE_GLOBAL_2D_ARRAY :: this process is not a cluster master, cluster_rank_key = ",i2," it should not call this subroutine")', Rank_of_process, cluster_rank_key
     errcode=400
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

  DO col = 1, N_clusters_x
     reclenx(col) = cluster_N_blocks_x * N_grid_block_x
  END DO
  reclenx(1)            = reclenx(1) + 2              ! include the leftmost boundary point and a y-coordinate 
  reclenx(N_clusters_x) = reclenx(N_clusters_x) + 1   ! include the rightmost bondary point

  DO row = 1, N_clusters_y
     recleny(row) = cluster_N_blocks_y * N_grid_block_y
  END DO
  recleny(1)            = recleny(1) + 1              ! include the bottom boundary point BUT doess NOT include the x-coordinate line
  recleny(N_clusters_y) = recleny(N_clusters_y) + 1   ! include the top boundary point

  CALL MPI_FILE_OPEN( COMM_HORIZONTAL, &
                    & filename,  &
                    & MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                    & MPI_INFO_NULL, &
                    & file_handle, &
                    & ierr )

  IF (c_column.EQ.1) THEN
     out_indx_x_min = c_indx_x_min-1
  ELSE
     out_indx_x_min = c_indx_x_min+1
  END IF

  IF (c_column.EQ.N_clusters_x) THEN
     out_indx_x_max = c_indx_x_max
  ELSE
     out_indx_x_max = c_indx_x_max-1
  END IF

  ALLOCATE(rbufer(out_indx_x_min:out_indx_x_max), STAT = ALLOC_ERR)

  IF (c_row.EQ.1) THEN
     out_indx_y_min = c_indx_y_min
  ELSE 
     out_indx_y_min = c_indx_y_min+1
  END IF

  IF (c_row.EQ.N_clusters_y) THEN
     out_indx_y_max = c_indx_y_max
  ELSE
     out_indx_y_max = c_indx_y_max-1
  END IF

  IF (c_row.EQ.1) THEN
! must save a line with x-coordinates, as required by gnuplot

     shifthead = 0
     DO col = 1, c_column-1
        shifthead = shifthead + reclenx(col)
     END DO
     shifthead = shifthead * 4

     CALL MPI_FILE_SEEK( file_handle, shifthead, MPI_SEEK_SET, ierr)

     IF (c_column.EQ.1) THEN
        rbufer(out_indx_x_min) = REAL(global_maximal_i + 1)    ! the very first element is the size of the x-line, as required by gnuplot
        DO i = out_indx_x_min+1, out_indx_x_max
           rbufer(i) = REAL(i * delta_x_m)
        END DO
     ELSE
        DO i = out_indx_x_min, out_indx_x_max
           rbufer(i) = REAL(i * delta_x_m)
        END DO
     END IF
     
     CALL MPI_FILE_WRITE( file_handle, rbufer(out_indx_x_min:out_indx_x_max), out_indx_x_max-out_indx_x_min+1, MPI_REAL, stattus, ierr )

  END IF

  DO j = out_indx_y_min, out_indx_y_max
        
     shifthead = global_maximal_i + 2
     DO row = 1, c_row-1
        shifthead = shifthead + recleny(row) * (global_maximal_i + 2)   ! +1 since indexing starts with 0 and another +1 since y-coordinate is included
     END DO
     shifthead = shifthead + (j - out_indx_y_min) * (global_maximal_i + 2) 
     DO col = 1, c_column-1
        shifthead = shifthead + reclenx(col)
     END DO
     shifthead = shifthead * 4

     CALL MPI_FILE_SEEK( file_handle, shifthead, MPI_SEEK_SET, ierr)
        
     IF (c_column.EQ.1) THEN
        rbufer(out_indx_x_min) = REAL(j * delta_x_m)    ! the very first element is the y-coordinate, as required by gnuplot
        rbufer(out_indx_x_min+1:out_indx_x_max) = arr(out_indx_x_min+1:out_indx_x_max, j)
     ELSE
        rbufer(out_indx_x_min:out_indx_x_max) = arr(out_indx_x_min:out_indx_x_max, j)
     END IF

     CALL MPI_FILE_WRITE( file_handle, rbufer(out_indx_x_min:out_indx_x_max), out_indx_x_max-out_indx_x_min+1, MPI_REAL, stattus, ierr )

  END DO

  CALL MPI_FILE_CLOSE(file_handle, ierr)

  DEALLOCATE(rbufer, STAT = ALLOC_ERR)

  IF (Rank_of_process.EQ.0) PRINT '("created file ",A40)', filename

END SUBROUTINE SAVE_GLOBAL_2D_ARRAY

!---------------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_ALL_VDF1D

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles
  USE Snapshots

  use mpi

  IMPLICIT NONE

  CHARACTER(15) filename      ! _NNNN_vdf1d.bin
                              ! ----x----I----x

  REAL, ALLOCATABLE :: xsplit (:)
  REAL, ALLOCATABLE :: ysplit (:)
  INTEGER ALLOC_ERR

  INTEGER i, j, s, n, pos1, pos2

  INTEGER file_handle

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER(kind=MPI_OFFSET_KIND) shifthead 

  INTEGER ibufsize
  INTEGER one_loc_rec_len

  INTEGER, ALLOCATABLE :: ibufer(:)
  REAL, ALLOCATABLE :: rbufer(:)

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (N_vdfbox_all.EQ.0) RETURN

  IF (save_evdf_snapshot(current_snap).EQ.NOANYVDF) RETURN
  IF (save_evdf_snapshot(current_snap).EQ.ONLY2D) RETURN

  IF (cluster_rank_key.NE.0) RETURN

  filename = '_NNNN_vdf1d.bin'
  filename(2:5) = convert_int_to_txt_string(current_snap, 4)

! calculate xsplit and ysplit (they differ from the ones used in CALCULATE_ELECTRON_VDF)

  ALLOCATE(xsplit(0:N_vdfbox_x), STAT = ALLOC_ERR)
  xsplit(0) = REAL(c_X_area_min * delta_x_m)
  DO i = 1, N_vdfbox_x-1
     xsplit(i) = REAL((c_X_area_min + DBLE(i) * (c_X_area_max - c_X_area_min) / N_vdfbox_x) * delta_x_m)
  END DO
  xsplit(N_vdfbox_x) = REAL(c_X_area_max * delta_x_m)

!print '(2x,i4,10(2x,f9.6))', Rank_of_process, xsplit(0:N_vdfbox_x)

  ALLOCATE(ysplit(0:N_vdfbox_y), STAT = ALLOC_ERR)
  ysplit(0) = REAL(c_Y_area_min * delta_x_m)
  DO j = 1, N_vdfbox_y-1
     ysplit(j) = (REAL(c_Y_area_min + DBLE(j) * (c_Y_area_max - c_Y_area_min) / N_vdfbox_y) * delta_x_m)
  END DO
  ysplit(N_vdfbox_y) = REAL(c_Y_area_max * delta_x_m)

! size of the integer bufer with the electron vdf (3: x,y,z) and ion vdf (3*N_spec) 
  ibufsize = 3 * (indx_v_max_e - indx_v_min_e + 1) + N_spec * 3 * (indx_v_max_i - indx_v_min_i + 1)

! length of vdf record for one sub-domain (box)
! includes 4 real numbers (boundaries of the box) and the integer bufer with all distribution functions
  one_loc_rec_len = 4 + ibufsize

  CALL MPI_FILE_OPEN( COMM_HORIZONTAL, &
                    & filename,  &
                    & MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                    & MPI_INFO_NULL, &
                    & file_handle, &
                    & ierr )

  ALLOCATE(ibufer(1:MAX(9,ibufsize)), STAT = ALLOC_ERR)
  ALLOCATE(rbufer(1:MAX(4,N_spec+1)), STAT = ALLOC_ERR)

  IF (Rank_of_process.EQ.0) THEN

     ibufer(1) = N_vdfbox_x * N_clusters_x
     ibufer(2) = N_vdfbox_y * N_clusters_y
     ibufer(3) = indx_v_min_e
     ibufer(4) = indx_v_max_e
     ibufer(5) = indx_v_min_i
     ibufer(6) = indx_v_max_i
     ibufer(7) = N_vbins_e
     ibufer(8) = N_vbins_i
     ibufer(9) = N_spec

     shifthead = 0
     CALL MPI_FILE_SEEK( file_handle, shifthead, MPI_SEEK_SET, ierr)

     CALL MPI_FILE_WRITE( file_handle, ibufer(1:9), 9, MPI_INTEGER, stattus, ierr )

     rbufer(1) = REAL(T_e_eV)
     DO s = 1, N_spec
        rbufer(s+1) = REAL(M_i_amu(s))
     END DO

     shifthead = 9*4
     CALL MPI_FILE_SEEK( file_handle, shifthead, MPI_SEEK_SET, ierr)

     CALL MPI_FILE_WRITE( file_handle, rbufer(1:(1+N_spec)), 1+N_spec, MPI_REAL, stattus, ierr )

  END IF

  DO j = 1, N_vdfbox_y

     DO i = 1, N_vdfbox_x

        rbufer(1) = xsplit(i-1)
        rbufer(2) = xsplit(i)
        rbufer(3) = ysplit(j-1)
        rbufer(4) = ysplit(j)

        n = i + (j-1) * N_vdfbox_x         ! number of the spatial box

        pos1 = 1
        pos2 = indx_v_max_e - indx_v_min_e + 1
        ibufer(pos1:pos2) = evxdf(indx_v_min_e:indx_v_max_e, n)

        pos1 = pos2 + 1
        pos2 = pos2 + (indx_v_max_e - indx_v_min_e + 1)
        ibufer(pos1:pos2) = evydf(indx_v_min_e:indx_v_max_e, n)

        pos1 = pos2 + 1
        pos2 = pos2 + (indx_v_max_e - indx_v_min_e + 1)
        ibufer(pos1:pos2) = evzdf(indx_v_min_e:indx_v_max_e, n)

        DO s = 1, N_spec
           pos1 = pos2 + 1
           pos2 = pos2 + (indx_v_max_i - indx_v_min_i + 1)
           ibufer(pos1:pos2) = isvxdf(indx_v_min_i:indx_v_max_i, n, s)

           pos1 = pos2 + 1
           pos2 = pos2 + (indx_v_max_i - indx_v_min_i + 1)
           ibufer(pos1:pos2) = isvydf(indx_v_min_i:indx_v_max_i, n, s)

           pos1 = pos2 + 1
           pos2 = pos2 + (indx_v_max_i - indx_v_min_i + 1)
           ibufer(pos1:pos2) = isvzdf(indx_v_min_i:indx_v_max_i, n, s)
        END DO

        shifthead = 9 + 1 + N_spec + &
                  & ((c_row   -1) * N_vdfbox_y + j - 1) * N_vdfbox_x * N_clusters_x * one_loc_rec_len + &
                  & ((c_column-1) * N_vdfbox_x + i - 1) * one_loc_rec_len
        shifthead = shifthead * 4

! note: (c_row-1)*N_vdfbox_y+j-1 changes from 0 at c_row=j=1 to N_clusters_y*N_vdfbox_y-1 at c_row=N_clusters_y and j=N_vdfbox_y
!
! note: ((c_column-1) * N_vdfbox_x + i - 1) changes from 0 at c_column=i=1 to N_clusters_x*N_vdfbox_x-1 at c_column=N_clusters_x and i=N_vdfbox_x

        CALL MPI_FILE_SEEK( file_handle, shifthead, MPI_SEEK_SET, ierr)

        CALL MPI_FILE_WRITE( file_handle, rbufer(1:4), 4, MPI_REAL, stattus, ierr )

        shifthead = shifthead + 4 * 4

        CALL MPI_FILE_SEEK( file_handle, shifthead, MPI_SEEK_SET, ierr)

        CALL MPI_FILE_WRITE( file_handle, ibufer(1:ibufsize), ibufsize, MPI_INTEGER, stattus, ierr )

     END DO

  END DO

  CALL MPI_FILE_CLOSE(file_handle, ierr)

  DEALLOCATE(xsplit, STAT=ALLOC_ERR)
  DEALLOCATE(ysplit, STAT=ALLOC_ERR)

  DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  DEALLOCATE(ibufer, STAT = ALLOC_ERR)

  IF (Rank_of_process.EQ.0) PRINT '("created file ",A15)', filename

END SUBROUTINE SAVE_ALL_VDF1D

!---------------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_ALL_VDF2D

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles
  USE Snapshots

  use mpi

  IMPLICIT NONE

  CHARACTER(15) filename      ! _NNNN_vdf2d.bin
                              ! ----x----I----x

  REAL, ALLOCATABLE :: xsplit (:)
  REAL, ALLOCATABLE :: ysplit (:)
  INTEGER ALLOC_ERR

  INTEGER i, j, n, pos1, pos2, jv

  INTEGER file_handle

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER(kind=MPI_OFFSET_KIND) shifthead 

  INTEGER ibufsize
  INTEGER one_loc_rec_len

  INTEGER, ALLOCATABLE :: ibufer(:)
  REAL, ALLOCATABLE :: rbufer(:)

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (N_vdfbox_all.EQ.0) RETURN

  IF (save_evdf_snapshot(current_snap).EQ.NOANYVDF) RETURN
  IF (save_evdf_snapshot(current_snap).EQ.ONLY1D) RETURN

  IF (cluster_rank_key.NE.0) RETURN

  filename = '_NNNN_vdf2d.bin'
  filename(2:5) = convert_int_to_txt_string(current_snap, 4)

! calculate xsplit and ysplit (they differ from the ones used in CALCULATE_ELECTRON_VDF)

  ALLOCATE(xsplit(0:N_vdfbox_x), STAT = ALLOC_ERR)
  xsplit(0) = REAL(c_X_area_min * delta_x_m)
  DO i = 1, N_vdfbox_x-1
     xsplit(i) = REAL((c_X_area_min + DBLE(i) * (c_X_area_max - c_X_area_min) / N_vdfbox_x) * delta_x_m)
  END DO
  xsplit(N_vdfbox_x) = REAL(c_X_area_max * delta_x_m)

  ALLOCATE(ysplit(0:N_vdfbox_y), STAT = ALLOC_ERR)
  ysplit(0) = REAL(c_Y_area_min * delta_x_m)
  DO j = 1, N_vdfbox_y-1
     ysplit(j) = (REAL(c_Y_area_min + DBLE(j) * (c_Y_area_max - c_Y_area_min) / N_vdfbox_y) * delta_x_m)
  END DO
  ysplit(N_vdfbox_y) = REAL(c_Y_area_max * delta_x_m)

! size of the integer bufer with a 2d electron vdf (Nvx*Nvy) 
  ibufsize = (indx_v_max_e - indx_v_min_e + 1) * (indx_v_max_e - indx_v_min_e + 1)

! length of vdf record for one sub-domain (box)
! includes 4 real numbers (boundaries of the box) and the integer bufer with the distribution function
  one_loc_rec_len = 4 + ibufsize

  CALL MPI_FILE_OPEN( COMM_HORIZONTAL, &
                    & filename,  &
                    & MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                    & MPI_INFO_NULL, &
                    & file_handle, &
                    & ierr )

  ALLOCATE(ibufer(1:MAX(5,ibufsize)), STAT = ALLOC_ERR)
  ALLOCATE(rbufer(1:4), STAT = ALLOC_ERR)

  IF (Rank_of_process.EQ.0) THEN

     ibufer(1) = N_vdfbox_x * N_clusters_x
     ibufer(2) = N_vdfbox_y * N_clusters_y
     ibufer(3) = indx_v_min_e
     ibufer(4) = indx_v_max_e
     ibufer(5) = N_vbins_e

     shifthead = 0
     CALL MPI_FILE_SEEK( file_handle, shifthead, MPI_SEEK_SET, ierr)

     CALL MPI_FILE_WRITE( file_handle, ibufer(1:5), 5, MPI_INTEGER, stattus, ierr )

     rbufer(1) = REAL(T_e_eV)

     shifthead = 5*4
     CALL MPI_FILE_SEEK( file_handle, shifthead, MPI_SEEK_SET, ierr)

     CALL MPI_FILE_WRITE( file_handle, rbufer(1:1), 1, MPI_REAL, stattus, ierr )

  END IF

  DO j = 1, N_vdfbox_y

     DO i = 1, N_vdfbox_x

        rbufer(1) = xsplit(i-1)
        rbufer(2) = xsplit(i)
        rbufer(3) = ysplit(j-1)
        rbufer(4) = ysplit(j)

        n = i + (j-1) * N_vdfbox_x         ! number of the spatial box

        pos2 = 0
        DO jv = indx_v_min_e, indx_v_max_e
           pos1 = pos2 + 1
           pos2 = pos2 + (indx_v_max_e - indx_v_min_e + 1)
           ibufer(pos1:pos2) = evxvydf(indx_v_min_e:indx_v_max_e, jv, n)
        END DO

        shifthead = 5 + 1 + &
                  & ((c_row   -1) * N_vdfbox_y + j - 1) * N_vdfbox_x * N_clusters_x * one_loc_rec_len + &
                  & ((c_column-1) * N_vdfbox_x + i - 1) * one_loc_rec_len
        shifthead = shifthead * 4

! note: (c_row-1)*N_vdfbox_y+j-1 changes from 0 at c_row=j=1 to N_clusters_y*N_vdfbox_y-1 at c_row=N_clusters_y and j=N_vdfbox_y
!
! note: ((c_column-1) * N_vdfbox_x + i - 1) changes from 0 at c_column=i=1 to N_clusters_x*N_vdfbox_x-1 at c_column=N_clusters_x and i=N_vdfbox_x

        CALL MPI_FILE_SEEK( file_handle, shifthead, MPI_SEEK_SET, ierr)

        CALL MPI_FILE_WRITE( file_handle, rbufer(1:4), 4, MPI_REAL, stattus, ierr )

        shifthead = shifthead + 4 * 4

        CALL MPI_FILE_SEEK( file_handle, shifthead, MPI_SEEK_SET, ierr)

        CALL MPI_FILE_WRITE( file_handle, ibufer(1:ibufsize), ibufsize, MPI_INTEGER, stattus, ierr )

     END DO

  END DO

  CALL MPI_FILE_CLOSE(file_handle, ierr)

  DEALLOCATE(xsplit, STAT=ALLOC_ERR)
  DEALLOCATE(ysplit, STAT=ALLOC_ERR)

  DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  DEALLOCATE(ibufer, STAT = ALLOC_ERR)

  IF (Rank_of_process.EQ.0) PRINT '("created file ",A15)', filename

END SUBROUTINE SAVE_ALL_VDF2D

!------------------------------------
!
SUBROUTINE SAVE_ELECTRON_PHASE_PLANES

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE Snapshots

  use mpi
  
  IMPLICIT NONE

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER count
  INTEGER i, j, k, n
  INTEGER ibufer_length, rbufer_length

  INTEGER, ALLOCATABLE :: ibufer(:)
  REAL,    ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR
  
  INTEGER pos

  INTEGER file_handle

  CHARACTER(13) filename_epp      ! _NNNN_epp.bin
                                  ! ----x----I---

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (N_pp_boxes.LE.0) RETURN

  IF (save_pp_snapshot(current_snap).EQ.NOANYPP) RETURN
  IF (save_pp_snapshot(current_snap).EQ.ONLYionPP) RETURN

! find number of particles to be saved
  count = 0
  DO k = 1, N_electrons
     i = INT(electron(k)%X)
     j = INT(electron(k)%Y)
     DO n = 1, N_pp_boxes
        IF ((i.GE.pp_box(n)%imin).AND.(i.LT.pp_box(n)%imax).AND.(j.GE.pp_box(n)%jmin).AND.(j.LT.pp_box(n)%jmax)) THEN
           count = count+1
! exit as soon as one box is found, even if the particle is in several boxes we save it only once
           EXIT
        END IF
     END DO
  END DO

  IF (Rank_of_process.EQ.0) THEN
     ibufer_length = 1+4*N_pp_boxes+1+1
     ALLOCATE(ibufer(ibufer_length), STAT=ALLOC_ERR)
     ibufer(1) = N_pp_boxes
     pos=2
     DO n = 1, N_pp_boxes
        ibufer(pos)   = pp_box(n)%imin
        ibufer(pos+1) = pp_box(n)%jmin
        ibufer(pos+2) = pp_box(n)%imax
        ibufer(pos+3) = pp_box(n)%jmax
        pos=pos+4
     END DO
     ibufer(pos) = N_of_processes
     ibufer(pos+1) = count
  ELSE
     ibufer_length = 1
     ALLOCATE(ibufer(ibufer_length), STAT=ALLOC_ERR)
     ibufer(ibufer_length) = count
  END IF

! create filename
  filename_epp = '_NNNN_epp.bin'
  filename_epp(2:5) = convert_int_to_txt_string(current_snap, 4)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL MPI_FILE_OPEN( MPI_COMM_WORLD, &
                    & filename_epp,  &
                    & MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                    & MPI_INFO_NULL, &
                    & file_handle, &
                    & ierr )
    
  CALL MPI_FILE_WRITE_ORDERED( file_handle, ibufer, ibufer_length, MPI_INTEGER, stattus, ierr )

  rbufer_length = 6*count
  IF (Rank_of_process.EQ.0) rbufer_length = rbufer_length + 1  ! to save delta_x_m

  ALLOCATE(rbufer(1:MAX(rbufer_length,1)), STAT=ALLOC_ERR)

! save particles into the buffers
  pos=1
  IF (Rank_of_process.EQ.0) THEN
     rbufer(pos) = REAL(delta_x_m)
     pos = pos+1
  END IF
  DO k = 1, N_electrons
     i = INT(electron(k)%X)
     j = INT(electron(k)%Y)
     DO n = 1, N_pp_boxes
        IF ((i.GE.pp_box(n)%imin).AND.(i.LT.pp_box(n)%imax).AND.(j.GE.pp_box(n)%jmin).AND.(j.LT.pp_box(n)%jmax)) THEN
           rbufer(pos)   = REAL(electron(k)%X)
           rbufer(pos+1) = REAL(electron(k)%Y)
           rbufer(pos+2) = REAL(electron(k)%VX * V_scale_ms)
           rbufer(pos+3) = REAL(electron(k)%VY * V_scale_ms)
           rbufer(pos+4) = REAL(electron(k)%VZ * V_scale_ms)
           rbufer(pos+5) = REAL(electron(k)%tag)
           pos = pos+6
! exit as soon as one box is found, even if the particle is in several boxes we save it only once
           EXIT
        END IF
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL MPI_FILE_WRITE_ORDERED( file_handle, rbufer, rbufer_length, MPI_REAL, stattus, ierr )

  CALL MPI_FILE_CLOSE(file_handle, ierr)

!  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  IF (ALLOCATED(ibufer)) DEALLOCATE(ibufer, STAT=ALLOC_ERR)

  IF (Rank_of_process.EQ.0) PRINT '("created file ",A13)', filename_epp

END SUBROUTINE SAVE_ELECTRON_PHASE_PLANES

!------------------------------------
!
SUBROUTINE SAVE_ION_PHASE_PLANES

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles
  USE Snapshots

  use mpi
  
  IMPLICIT NONE

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER s
  INTEGER count
  INTEGER i, j, k, n
  INTEGER ibufer_length, rbufer_length

  INTEGER, ALLOCATABLE :: ibufer(:)
  REAL,    ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR
  
  INTEGER pos

  INTEGER file_handle

  CHARACTER(13) filename_ipp      ! _NNNN_ipp.bin
                                  ! ----x----I---

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (N_pp_boxes.LE.0) RETURN

  IF (save_pp_snapshot(current_snap).EQ.NOANYPP) RETURN
  IF (save_pp_snapshot(current_snap).EQ.ONLYelectronPP) RETURN

  IF (Rank_of_process.EQ.0) THEN
     ibufer_length = 1+4*N_pp_boxes+1+1
     ALLOCATE(ibufer(ibufer_length), STAT=ALLOC_ERR)
     ibufer(1) = N_pp_boxes
     pos=2
     DO n = 1, N_pp_boxes
        ibufer(pos)   = pp_box(n)%imin
        ibufer(pos+1) = pp_box(n)%jmin
        ibufer(pos+2) = pp_box(n)%imax
        ibufer(pos+3) = pp_box(n)%jmax
        pos=pos+4
     END DO
     ibufer(pos) = N_of_processes
     ibufer(pos+1) = N_spec
  ELSE
     ibufer_length = 0
     ALLOCATE(ibufer(1), STAT=ALLOC_ERR)
     ibufer = 0  ! just to have some value
  END IF

! create filename
  filename_ipp = '_NNNN_ipp.bin'
  filename_ipp(2:5) = convert_int_to_txt_string(current_snap, 4)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL MPI_FILE_OPEN( MPI_COMM_WORLD, &
                    & filename_ipp,  &
                    & MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                    & MPI_INFO_NULL, &
                    & file_handle, &
                    & ierr )
    
  CALL MPI_FILE_WRITE_ORDERED( file_handle, ibufer, ibufer_length, MPI_INTEGER, stattus, ierr )

  DO s = 1, N_spec

! find number of particles to be saved
     count = 0
     DO k = 1, N_ions(s)
        i = INT(ion(s)%part(k)%X)
        j = INT(ion(s)%part(k)%Y)
        DO n = 1, N_pp_boxes
           IF ((i.GE.pp_box(n)%imin).AND.(i.LT.pp_box(n)%imax).AND.(j.GE.pp_box(n)%jmin).AND.(j.LT.pp_box(n)%jmax)) THEN
              count = count+1
! exit as soon as one box is found, even if the particle is in several boxes we save it only once
              EXIT
           END IF
        END DO
     END DO

     ibufer(1) = count

     CALL MPI_FILE_WRITE_ORDERED( file_handle, ibufer(1:1), 1, MPI_INTEGER, stattus, ierr )

     rbufer_length = 6*count
     IF ((Rank_of_process.EQ.0).AND.(s.EQ.1)) rbufer_length = rbufer_length + 1  ! to save delta_x_m

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
     ALLOCATE(rbufer(1:MAX(rbufer_length,1)), STAT=ALLOC_ERR)

! save particles into the buffer
     pos=1
     IF ((Rank_of_process.EQ.0).AND.(s.EQ.1)) THEN
        rbufer(pos) = REAL(delta_x_m)
        pos = pos+1
     END IF
     DO k = 1, N_ions(s)
        i = INT(ion(s)%part(k)%X)
        j = INT(ion(s)%part(k)%Y)
        DO n = 1, N_pp_boxes
           IF ((i.GE.pp_box(n)%imin).AND.(i.LT.pp_box(n)%imax).AND.(j.GE.pp_box(n)%jmin).AND.(j.LT.pp_box(n)%jmax)) THEN
              rbufer(pos)   = REAL(ion(s)%part(k)%X)
              rbufer(pos+1) = REAL(ion(s)%part(k)%Y)
              rbufer(pos+2) = REAL(ion(s)%part(k)%VX * V_scale_ms)
              rbufer(pos+3) = REAL(ion(s)%part(k)%VY * V_scale_ms)
              rbufer(pos+4) = REAL(ion(s)%part(k)%VZ * V_scale_ms)
              rbufer(pos+5) = REAL(ion(s)%part(k)%tag)
              pos = pos+6
! exit as soon as one box is found, even if the particle is in several boxes we save it only once
              EXIT
           END IF
        END DO
     END DO
     
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL MPI_FILE_WRITE_ORDERED( file_handle, rbufer, rbufer_length, MPI_REAL, stattus, ierr )

  END DO

  CALL MPI_FILE_CLOSE(file_handle, ierr)

!  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  IF (ALLOCATED(ibufer)) DEALLOCATE(ibufer, STAT=ALLOC_ERR)

  IF (Rank_of_process.EQ.0) PRINT '("created file ",A13)', filename_ipp

END SUBROUTINE SAVE_ION_PHASE_PLANES

!---------------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_en_COLLISIONS_2D

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : N_subcycles, T_cntr, N_plasma_m3, N_of_particles_cell, delta_t_s
  USE ClusterAndItsBoundaries
!  USE IonParticles
  USE MCCollisions
  USE Snapshots

  use mpi

  IMPLICIT NONE

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER T_start
  REAL conversion_factor_m3s

  CHARACTER(40) filename      ! _NNNN_neutral_AAAAAA_coll_id_NN_i_NN.bin
                              ! ----x----I----x----I----x----I----x----I

  INTEGER n, p, i, j

  INTEGER n1  ! number of nodes in the y-direction
  INTEGER n3  ! number of nodes in the x-direction

  REAL, ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (cluster_rank_key.NE.0) RETURN
  IF (en_collisions_turned_off) RETURN
  IF (no_ionization_collisions) RETURN

  IF (current_snap.GT.N_of_all_snaps) RETURN
  IF (.NOT.save_ionization_rates_2d(current_snap)) RETURN

! this procedure is called at the end of the ion advance 
! save ionization rates at the ion advance time step which is the closest to (and precede) the time step of the current snapshot
  IF ((Tcntr_snapshot(current_snap)-T_cntr).GT.N_subcycles) RETURN

! calculate conversion factor
  IF (current_snap.EQ.1) THEN
     T_start = -1
  ELSE
     T_start = N_subcycles * (Tcntr_snapshot(current_snap-1) / N_subcycles)-1  ! timestep when ions were advanced before the previous snapshot
  END IF
  conversion_factor_m3s = 0.0
  IF (T_cntr.GT.T_start) conversion_factor_m3s = REAL(N_plasma_m3 / (N_of_particles_cell * (T_cntr - T_start) * delta_t_s))

  n1 = c_indx_y_max - c_indx_y_min + 1
  n3 = c_indx_x_max - c_indx_x_min + 1

  DO n = 1, N_neutral_spec

     DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
        IF (collision_e_neutral(n)%colproc_info(p)%type.LT.30) CYCLE  ! skip non-ionizing collisions

! exchange information about collisions in overlapping nodes    

        IF (WHITE_CLUSTER) THEN  
! "white processes"

           IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

           IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right densities in the right edge
              rbufer(1:n1) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, c_indx_y_min:c_indx_y_max)
              CALL MPI_SEND(rbufer, n1, MPI_REAL, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, ierr) 
           END IF

           IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left densities in the left edge
              rbufer(1:n1) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, c_indx_y_min:c_indx_y_max)
              CALL MPI_SEND(rbufer, n1, MPI_REAL, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, ierr) 
           END IF

           IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left densities in the vertical line next to the left edge
              CALL MPI_RECV(rbufer, n1, MPI_REAL, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
              DO j = c_indx_y_min, c_indx_y_max
                 diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min+1, j) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min+1, j) + rbufer(j-c_indx_y_min+1)
              END DO
           END IF

           IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right densities in the vertical line next to the right edge
              CALL MPI_RECV(rbufer, n1, MPI_REAL, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
              DO j = c_indx_y_min, c_indx_y_max
                 diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max-1, j) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max-1, j) + rbufer(j-c_indx_y_min+1)
              END DO
           END IF

           IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

           IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up densities in the top edge
              rbufer(1:n3) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min:c_indx_x_max, c_indx_y_max)
              CALL MPI_SEND(rbufer, n3, MPI_REAL, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, ierr) 
           END IF

           IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down densities in the bottom edge
              rbufer(1:n3) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min:c_indx_x_max, c_indx_y_min)
              CALL MPI_SEND(rbufer, n3, MPI_REAL, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, ierr) 
           END IF

           IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below densities in the vertical line above the bottom line
              CALL MPI_RECV(rbufer, n3, MPI_REAL, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
              DO i = c_indx_x_min, c_indx_x_max
                 diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_min+1) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_min+1) + rbufer(i-c_indx_x_min+1)
              END DO
           END IF

           IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above densities in the vertical line under the top line
              CALL MPI_RECV(rbufer, n3, MPI_REAL, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
              DO i = c_indx_x_min, c_indx_x_max
                 diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_max-1) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_max-1) + rbufer(i-c_indx_x_min+1)
              END DO
           END IF

        ELSE
! "black" processes

           IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

           IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left densities in the vertical line next to the left edge
              CALL MPI_RECV(rbufer, n1, MPI_REAL, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
              DO j = c_indx_y_min, c_indx_y_max
                 diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min+1, j) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min+1, j) + rbufer(j-c_indx_y_min+1)
              END DO
           END IF

           IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right densities in the vertical line next to the right edge
              CALL MPI_RECV(rbufer, n1, MPI_REAL, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
              DO j = c_indx_y_min, c_indx_y_max
                 diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max-1, j) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max-1, j) + rbufer(j-c_indx_y_min+1)
              END DO
           END IF

           IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right densities in the right edge
              rbufer(1:n1) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, c_indx_y_min:c_indx_y_max)
              CALL MPI_SEND(rbufer, n1, MPI_REAL, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, ierr) 
           END IF

           IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left densities in the left edge
              rbufer(1:n1) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, c_indx_y_min:c_indx_y_max)
              CALL MPI_SEND(rbufer, n1, MPI_REAL, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, ierr) 
           END IF

           IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

           IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below densities in the vertical line above the bottom line
              CALL MPI_RECV(rbufer, n3, MPI_REAL, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
              DO i = c_indx_x_min, c_indx_x_max
                 diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_min+1) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_min+1) + rbufer(i-c_indx_x_min+1)
              END DO
           END IF

           IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above densities in the vertical line under the top line
              CALL MPI_RECV(rbufer, n3, MPI_REAL, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
              DO i = c_indx_x_min, c_indx_x_max
                 diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_max-1) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_max-1) + rbufer(i-c_indx_x_min+1)
              END DO
           END IF

           IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up densities in the top edge
              rbufer(1:n3) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min:c_indx_x_max, c_indx_y_max)
              CALL MPI_SEND(rbufer, n3, MPI_REAL, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, ierr) 
           END IF

           IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down densities in the bottom edge
              rbufer(1:n3) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min:c_indx_x_max, c_indx_y_min)
              CALL MPI_SEND(rbufer, n3, MPI_REAL, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, ierr) 
           END IF

        END IF

! adjust values at the boundaries with material walls (where the cell area may be 3/4, 1/2, or 1/4 of an inner cell)

        IF (Rank_of_master_left.LT.0) THEN
           DO j = c_indx_y_min+1, c_indx_y_max-1
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, j) = 2.0 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, j)
           END DO
        END IF

        IF (Rank_of_master_right.LT.0) THEN
           DO j = c_indx_y_min+1, c_indx_y_max-1
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, j) = 2.0 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, j)
           END DO
        END IF
  
        IF (Rank_of_master_above.LT.0) THEN
           DO i = c_indx_x_min+1, c_indx_x_max-1
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_max) = 2.0 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_max)
           END DO
        END IF
  
        IF (Rank_of_master_below.LT.0) THEN
           DO i = c_indx_x_min+1, c_indx_x_max-1
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_min) = 2.0 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i, c_indx_y_min)
           END DO
        END IF

        SELECT CASE (c_left_bottom_corner_type)
           CASE (SURROUNDED_BY_WALL)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, c_indx_y_min) = 4.0 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, c_indx_y_min) 
           CASE (EMPTY_CORNER_WALL_LEFT)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, c_indx_y_min+1) = 0.66666666666666 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_BELOW)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min+1, c_indx_y_min) = 0.66666666666666 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min+1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
        END SELECT
     
        SELECT CASE (c_left_top_corner_type)
           CASE (SURROUNDED_BY_WALL)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, c_indx_y_max) = 4.0 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, c_indx_y_max) 
           CASE (EMPTY_CORNER_WALL_LEFT)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, c_indx_y_max-1) = 0.66666666666666 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_ABOVE)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min+1, c_indx_y_max) = 0.66666666666666 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_min+1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
        END SELECT

        SELECT CASE (c_right_bottom_corner_type)
           CASE (SURROUNDED_BY_WALL)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, c_indx_y_min) = 4.0 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, c_indx_y_min) 
           CASE (EMPTY_CORNER_WALL_RIGHT)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, c_indx_y_min+1) = 0.66666666666666 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_BELOW)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max-1, c_indx_y_min) = 0.66666666666666 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max-1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
        END SELECT
     
        SELECT CASE (c_right_top_corner_type)
           CASE (SURROUNDED_BY_WALL)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, c_indx_y_max) = 4.0 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, c_indx_y_max) 
           CASE (EMPTY_CORNER_WALL_RIGHT)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, c_indx_y_max-1) = 0.66666666666666 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_ABOVE)
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max-1, c_indx_y_max) = 0.66666666666666 * diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(c_indx_x_max-1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
        END SELECT

        filename = '_NNNN_neutral_AAAAAA_coll_id_NN_i_NN.bin'
        filename(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename(15:20) = neutral(n)%name
        filename(30:31) = convert_int_to_txt_string(collision_e_neutral(n)%colproc_info(p)%id_number, 2)
        filename(35:36) = convert_int_to_txt_string(collision_e_neutral(n)%colproc_info(p)%ion_species_produced, 2)

        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i,j) = diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local(i,j) * conversion_factor_m3s
           END DO
        END DO

        CALL SAVE_GLOBAL_2D_ARRAY(diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local, filename)

! complete cleanup (deallocate + zeroing) is performed in the end of CREATE_SNAPSHOT
!        diagnostics_neutral(n)%activated_collision(p)%ionization_rate_local = 0.0

        IF (Rank_of_process.EQ.0) PRINT '("created file ",A40)', filename

     END DO  !###    DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
  END DO     !### DO n = 1, N_neutral_spec


END SUBROUTINE SAVE_en_COLLISIONS_2D

!---------------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_IONS_COLLIDED_WITH_BOUNDARY_OBJECTS

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : N_subcycles, T_cntr, N_of_boundary_and_inner_objects, delta_t_s, &
                                 & delta_x_m, V_scale_ms, N_scale_part_m3, ion_colls_with_bo, whole_object
  USE IonParticles, ONLY : N_spec, M_i_amu
  USE Snapshots

  use mpi

  IMPLICIT NONE

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  LOGICAL start_new_bo_coll_file
  LOGICAL this_was_last_record_to_bo_coll_file

  INTEGER shortibuf1(1), shortibuf2(1)

  CHARACTER(34) filename      ! _NNNN_ions_collided_with_bo_NN.bin
                              ! ----x----I----x----I----x----I----
  INTEGER file_handle

  INTEGER n, pos, bufsize, ibufsize, m, s, k

  REAL, ALLOCATABLE :: rbufer(:)
  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (current_snap.GT.N_of_all_snaps) RETURN

  IF (.NOT.save_ions_collided_with_bo(current_snap)) RETURN

  start_new_bo_coll_file = .FALSE.
  this_was_last_record_to_bo_coll_file = .FALSE.

! this procedure is called after the ion advance 
  IF (current_snap.EQ.1) THEN
! for the very first snapshot, make the first record at the first ion advance time step
     IF (T_cntr.LT.N_subcycles) start_new_bo_coll_file = .TRUE.
  ELSE
! make the first record at the ion advance time step which is the closest to (and follows) the time step of the previous snapshot
     IF ((T_cntr-Tcntr_snapshot(current_snap-1)).LT.N_subcycles) start_new_bo_coll_file = .TRUE.
  END IF

! the last record is made at the ion advance time step which is the closest to (and precede) the time step of the current  snapshot
  IF ((Tcntr_snapshot(current_snap)-T_cntr).LE.N_subcycles) this_was_last_record_to_bo_coll_file = .TRUE.

  DO n = 1, N_of_boundary_and_inner_objects

     IF (.NOT.ion_colls_with_bo(n)%must_be_saved) CYCLE

     shortibuf1(1) = ion_colls_with_bo(n)%N_of_saved_parts
     shortibuf2(1) = 0

     CALL MPI_REDUCE(shortibuf1, shortibuf2, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

     filename = '_NNNN_ions_collided_with_bo_NN.bin'
     filename(2:5) = convert_int_to_txt_string(current_snap, 4)
     filename(29:30) = convert_int_to_txt_string(n, 2)

     pos = 1
     bufsize = 5 * ion_colls_with_bo(n)%N_of_saved_parts
     ALLOCATE(rbufer(1:(bufsize+4+N_spec+2)), STAT=ALLOC_ERR)  ! 4+N_spec+2 allows same process save 
                                                               ! start time, mesh size, scale density to particle number ratio, scale_velocity, 
                                                               ! ion masses, total number of particles (from all processes), and record time

     IF (start_new_bo_coll_file) THEN
        CALL MPI_FILE_OPEN( MPI_COMM_WORLD, &
                          & filename,  &
                          & MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                          & MPI_INFO_NULL, &
                          & file_handle, &
                          & ierr )

        IF (Rank_of_process.EQ.0) THEN

           PRINT '("Started creating file ",A34)', filename

! the first process saves the following additional values
! integers: the number of segments, the segment ends coordinates, the number of ion species
           ibufsize = 2 + 4 * whole_object(n)%number_of_segments
           ALLOCATE(ibufer(1:ibufsize), STAT=ALLOC_ERR)
           ibufer(1) = whole_object(n)%number_of_segments
           DO m = 1, whole_object(n)%number_of_segments
              ibufer(1+4*(m-1)+1) = whole_object(n)%segment(m)%istart
              ibufer(1+4*(m-1)+2) = whole_object(n)%segment(m)%jstart
              ibufer(1+4*(m-1)+3) = whole_object(n)%segment(m)%iend
              ibufer(1+4*(m-1)+4) = whole_object(n)%segment(m)%jend
           END DO
           ibufer(ibufsize) = N_spec
           
! reals : time when the set begins (ns), mesh size (m), density of one macroparticle, scale velocity
           bufsize = bufsize+4+N_spec
           rbufer(1) = (T_cntr - N_subcycles) * (delta_t_s * 1.0d9) ! start time [ns]
           rbufer(2) = delta_x_m
           rbufer(3) = N_scale_part_m3
           rbufer(4) = V_scale_ms
           DO s = 1, N_spec
              rbufer(4+s) = M_i_amu(s)
           END DO
           pos = 4 + N_spec + 1

        ELSE
           ALLOCATE(ibufer(1), STAT=ALLOC_ERR)
           ibufsize=0
        END IF
! here only the zero rank process actually writes integer data, other processes do a zero-length record, that is nothing
! and this happens only when the file was opened or the very first time (created)
        CALL MPI_FILE_WRITE_ORDERED( file_handle, ibufer, ibufsize, MPI_INTEGER, stattus, ierr )
     ELSE
        CALL MPI_FILE_OPEN( MPI_COMM_WORLD, &
                          & filename,  &
                          & MPI_MODE_WRONLY + MPI_MODE_APPEND, & 
                          & MPI_INFO_NULL, &
                          & file_handle, &
                          & ierr )
     END IF

     IF (Rank_of_process.EQ.0) THEN
        rbufer(pos) = T_cntr * (delta_t_s * 1.0d9) ! record time [ns]
        rbufer(pos+1) = shortibuf2(1)
        bufsize = bufsize+2
        pos = pos+2
     END IF

     DO k = 1, ion_colls_with_bo(n)%N_of_saved_parts
        rbufer(pos)   = ion_colls_with_bo(n)%part(k)%token
        rbufer(pos+1) = ion_colls_with_bo(n)%part(k)%coll_coord
        rbufer(pos+2) = ion_colls_with_bo(n)%part(k)%VX
        rbufer(pos+3) = ion_colls_with_bo(n)%part(k)%VY
        rbufer(pos+4) = ion_colls_with_bo(n)%part(k)%VZ
        pos = pos+5
     END DO

     CALL MPI_FILE_WRITE_ORDERED( file_handle, rbufer, bufsize, MPI_REAL, stattus, ierr )

     CALL MPI_FILE_CLOSE(file_handle, ierr)

     IF (Rank_of_process.EQ.0) THEN
        IF (this_was_last_record_to_bo_coll_file) PRINT '("Finished creating file ",A34)', filename
     END IF

! cleanup
     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     IF (ALLOCATED(ibufer)) DEALLOCATE(ibufer, STAT=ALLOC_ERR)

  END DO

END SUBROUTINE SAVE_IONS_COLLIDED_WITH_BOUNDARY_OBJECTS

!---------------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_ELECTRONS_COLLIDED_WITH_BOUNDARY_OBJECTS

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : N_subcycles, T_cntr, N_of_boundary_and_inner_objects, delta_t_s, &
                                 & delta_x_m, V_scale_ms, N_scale_part_m3, e_colls_with_bo, whole_object
  USE Snapshots

  use mpi

  IMPLICIT NONE

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  LOGICAL start_new_bo_coll_file
  LOGICAL this_was_last_record_to_bo_coll_file

  INTEGER shortibuf1(1), shortibuf2(1)

  CHARACTER(31) filename      ! _NNNN_e_collided_with_bo_NN.bin
                              ! ----x----I----x----I----x----I-
  INTEGER file_handle

  INTEGER n, pos, bufsize, ibufsize, m, k

  REAL, ALLOCATABLE :: rbufer(:)
  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (current_snap.GT.N_of_all_snaps) RETURN

  IF (.NOT.save_e_collided_with_bo(current_snap)) RETURN

  start_new_bo_coll_file = .FALSE.
  this_was_last_record_to_bo_coll_file = .FALSE.

! note, in the main cycle the order of subroutines is
! CREATE_SNAPSHOT => ADVANCE_ELECTRONS => SAVE_ELECTRONS_COLLIDED_WITH_BOUNDARY_OBJECTS
!
! if Tcntr_snapshot(1)>0 then the very first call of SAVE_ELECTRONS_COLLIDED_WITH_BOUNDARY_OBJECTS @ Tcntr=0 is when current_snap == 1
! in this case the first branch of the IF below works
!
! if Tcntr_snapshot(1)==0 then the very first call of SAVE_ELECTRONS_COLLIDED_WITH_BOUNDARY_OBJECTS @ T_cntr=0 is when current_snap == 2 (or larger)
! in this case the second branch of the IF below works

! this procedure is called after the electron advance 
  IF (current_snap.EQ.1) THEN
! for the very first snapshot, make the first record at the first time step
     IF (T_cntr.EQ.0) start_new_bo_coll_file = .TRUE.
  ELSE
! make the first record at the time step when the previous snapshot was created
     IF (T_cntr.EQ.Tcntr_snapshot(current_snap-1)) start_new_bo_coll_file = .TRUE.
  END IF

! the last record is made at the time step which precedes the time step of the current snapshot
  IF (T_cntr.EQ.(Tcntr_snapshot(current_snap)-1)) this_was_last_record_to_bo_coll_file = .TRUE.

  DO n = 1, N_of_boundary_and_inner_objects

     IF (.NOT.e_colls_with_bo(n)%must_be_saved) CYCLE

     shortibuf1(1) = e_colls_with_bo(n)%N_of_saved_parts
     shortibuf2(1) = 0

     CALL MPI_REDUCE(shortibuf1, shortibuf2, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

     filename = '_NNNN_e_collided_with_bo_NN.bin'
     filename(2:5) = convert_int_to_txt_string(current_snap, 4)
     filename(26:27) = convert_int_to_txt_string(n, 2)

     pos = 1
     bufsize = 5 * e_colls_with_bo(n)%N_of_saved_parts
     ALLOCATE(rbufer(1:(bufsize+6)), STAT=ALLOC_ERR)  ! 6=4+2 allows same process save 
                                                      ! start time, mesh size, scale density to particle number ratio, scale_velocity, 
                                                      ! total number of particles (from all processes), and record time
                                                      ! (reserved fora an improbable case when one process has to save everything)

     IF (start_new_bo_coll_file) THEN
        CALL MPI_FILE_OPEN( MPI_COMM_WORLD, &
                          & filename,  &
                          & MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                          & MPI_INFO_NULL, &
                          & file_handle, &
                          & ierr )

        IF (Rank_of_process.EQ.0) THEN

           PRINT '("Started creating file ",A31)', filename

! the first process saves the following additional values
! integers: the number of segments, the segment ends coordinates, the number of ion species
           ibufsize = 1 + 4 * whole_object(n)%number_of_segments
           ALLOCATE(ibufer(1:ibufsize), STAT=ALLOC_ERR)
           ibufer(1) = whole_object(n)%number_of_segments
           DO m = 1, whole_object(n)%number_of_segments
              ibufer(1+4*(m-1)+1) = whole_object(n)%segment(m)%istart
              ibufer(1+4*(m-1)+2) = whole_object(n)%segment(m)%jstart
              ibufer(1+4*(m-1)+3) = whole_object(n)%segment(m)%iend
              ibufer(1+4*(m-1)+4) = whole_object(n)%segment(m)%jend
           END DO
           
! reals : time when the set begins (ns), mesh size (m), density of one macroparticle, scale velocity
           bufsize = bufsize+4
           rbufer(1) = (T_cntr - 1) * (delta_t_s * 1.0d9) ! start time [ns]
           rbufer(2) = delta_x_m
           rbufer(3) = N_scale_part_m3
           rbufer(4) = V_scale_ms
           pos = 5

        ELSE
           ALLOCATE(ibufer(1), STAT=ALLOC_ERR)
           ibufsize=0
        END IF
! here only the zero rank process actually writes integer data, other processes do a zero-length record, that is nothing
! and this happens only when the file was opened or the very first time (created)
        CALL MPI_FILE_WRITE_ORDERED( file_handle, ibufer, ibufsize, MPI_INTEGER, stattus, ierr )
     ELSE
        CALL MPI_FILE_OPEN( MPI_COMM_WORLD, &
                          & filename,  &
                          & MPI_MODE_WRONLY + MPI_MODE_APPEND, & 
                          & MPI_INFO_NULL, &
                          & file_handle, &
                          & ierr )
     END IF

     IF (Rank_of_process.EQ.0) THEN
        rbufer(pos) = T_cntr * (delta_t_s * 1.0d9) ! record time [ns]
        rbufer(pos+1) = shortibuf2(1)
        bufsize = bufsize+2
        pos = pos+2
     END IF

     DO k = 1, e_colls_with_bo(n)%N_of_saved_parts
        rbufer(pos)   = e_colls_with_bo(n)%part(k)%token
        rbufer(pos+1) = e_colls_with_bo(n)%part(k)%coll_coord
        rbufer(pos+2) = e_colls_with_bo(n)%part(k)%VX
        rbufer(pos+3) = e_colls_with_bo(n)%part(k)%VY
        rbufer(pos+4) = e_colls_with_bo(n)%part(k)%VZ
        pos = pos+5
     END DO

     CALL MPI_FILE_WRITE_ORDERED( file_handle, rbufer, bufsize, MPI_REAL, stattus, ierr )

     CALL MPI_FILE_CLOSE(file_handle, ierr)

     IF (Rank_of_process.EQ.0) THEN
        IF (this_was_last_record_to_bo_coll_file) PRINT '("Finished creating file ",A31)', filename
     END IF

! cleanup
     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     IF (ALLOCATED(ibufer)) DEALLOCATE(ibufer, STAT=ALLOC_ERR)

  END DO

END SUBROUTINE SAVE_ELECTRONS_COLLIDED_WITH_BOUNDARY_OBJECTS

!---------------------------------------------------------------------------------------------------
!
SUBROUTINE FINISH_SNAPSHOTS

  USE Snapshots
  IMPLICIT NONE

  INTEGER DEALLOC_ERR

  IF (ALLOCATED(    Tcntr_snapshot)) DEALLOCATE(    Tcntr_snapshot, STAT=DEALLOC_ERR)
  IF (ALLOCATED(save_evdf_snapshot)) DEALLOCATE(save_evdf_snapshot, STAT=DEALLOC_ERR)
  IF (ALLOCATED(save_pp_snapshot)) DEALLOCATE(save_pp_snapshot, STAT=DEALLOC_ERR)
  IF (ALLOCATED(save_ionization_rates_2d)) DEALLOCATE(save_ionization_rates_2d, STAT=DEALLOC_ERR)
  IF (ALLOCATED(save_ions_collided_with_bo)) DEALLOCATE(save_ions_collided_with_bo, STAT=DEALLOC_ERR)
  IF (ALLOCATED(save_e_collided_with_bo)) DEALLOCATE(save_e_collided_with_bo, STAT=DEALLOC_ERR)

END SUBROUTINE FINISH_SNAPSHOTS
