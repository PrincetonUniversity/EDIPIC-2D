
!--------------------

SUBROUTINE INITIATE_SNAPSHOTS

  USE ParallelOperationValues
  USE Snapshots
  USE Diagnostics,          ONLY : Save_probes_data_T_cntr_rff, Save_probes_data_step
  USE CurrentProblemValues, ONLY : delta_t_s, Max_T_cntr
  USE Checkpoints, ONLY : use_checkpoint, current_snap_check
  
  IMPLICIT NONE

  INTEGER T2_old, N2_old
  LOGICAL exists
  CHARACTER (1) buf
  INTEGER flag

  INTEGER saveflagi(1:32)           ! integer flags used to set values of logical flags controlling saving of data files

  INTEGER N_of_snap_groups          ! number of sets of snapshots, read from file
  INTEGER i                         ! set of snapshots

  REAL(8) Rqst_snap_start_ns        ! requested start of current set of snapshots [ns], read from file
  REAL(8) Rqst_snap_finish_ns       ! requested finish of current set of snapshots [ns], read from file 
  INTEGER Rqst_n_of_snaps           ! requested number of snapshots in current set, read from file
  INTEGER Rqst_evdf_flag            ! requested flag defining which evdfs to save (0/1/2/3 = No/1d only/2d only/1d and 2d)

  INTEGER T1, T2, N1, N2 

  INTEGER large_step 

  INTEGER Fact_n_of_snaps           ! calculated number of snapshots in one set

  INTEGER n                             ! ordering number of snapshot in the set
  INTEGER, ALLOCATABLE ::           timestep(:)   ! array for temporary storage of moments (timesteps) of snapshots
  INTEGER, ALLOCATABLE :: evdf_flag_timestep(:)   ! array for temporary storage of evdf save flags

  INTEGER ALLOC_ERR

  N_of_all_snaps = 0
  T2_old = -1
  N2_old = 0

! read / write the data file 
  INQUIRE (FILE = 'init_snapshots.dat', EXIST = exists)

  IF (exists) THEN

     ALLOCATE(          timestep(1:9999), STAT = ALLOC_ERR)
     ALLOCATE(evdf_flag_timestep(1:9999), STAT = ALLOC_ERR)
     
     IF (Rank_of_process.EQ.0) PRINT '("### File init_snapshots is found. Reading the data file... ###")'

     OPEN (9, FILE = 'init_snapshots.dat')
     
     READ(9, '(A1)') buf  !--- save 2D maps of the following parameters? (1=yes, 0=no)
     READ(9, '(A1)') buf  !-----F----EX----EY--JXsum--JYsum--JZsum
     READ(9, '(A1)') buf  !----dd----dd----dd----dd----dd----dd
     READ(9, '(6(4x,i2))') saveflagi(1:6) 
     READ(9, '(A1)') buf  !----Ne----JXe---JYe---JZe---VXe---VYe---VZe---WXe---WYe---WZe---TXe---TYe---TZe
     READ(9, '(A1)') buf  !----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd
     READ(9, '(13(4x,i2))') saveflagi(7:19)
     READ(9, '(A1)') buf  !----Ni----JXi---JYi---JZi---VXi---VYi---VZi---WXi---WYi---WZi---TXi---TYi---TZi
     READ(9, '(A1)') buf  !----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd----dd
     READ (9, '(13(4x,i2))') saveflagi(20:32)

     save_data = .TRUE.
     DO i = 1, 32
        IF (saveflagi(i).EQ.0) save_data(i) = .FALSE.
     END DO

     READ (9, '(A1)') buf !---dd--- Number of groups of snapshots ( >= 0 )
     READ (9, '(3x,i2)') N_of_snap_groups
     READ (9, '(A1)') buf !---ddddddd.ddd---ddddddd.ddd---dddd---d--- group: start (ns) / finish (ns) / number of snapshots / save evdfs (0/1/2/3 = No/1d only/2d only/1d and 2d)

     DO i = 1, N_of_snap_groups
! read the parameters of current set of snapshot from the data file
        READ (9, '(3x,f11.3,3x,f11.3,3x,i4,3x,i1)') Rqst_snap_start_ns, Rqst_snap_finish_ns, Rqst_n_of_snaps, Rqst_evdf_flag    !!!
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

     CLOSE (9, STATUS = 'KEEP')

  ELSE

     PRINT '(2x,"Process ",i4," : ERROR : file init_snapshots.dat not found. Program terminated")', Rank_of_process
     STOP

  END IF

  current_snap = 1   ! default value

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
  ALLOCATE(    Tcntr_snapshot(1:N_of_all_snaps), STAT=ALLOC_ERR)
  ALLOCATE(save_evdf_snapshot(1:N_of_all_snaps), STAT=ALLOC_ERR)

! move the calculated snapshot moments from the temporary array to the allocated array 
      Tcntr_snapshot(1:N_of_all_snaps) =           timestep(1:N_of_all_snaps)
  save_evdf_snapshot(1:N_of_all_snaps) = evdf_flag_timestep(1:N_of_all_snaps)
 
  DEALLOCATE(          timestep, STAT = ALLOC_ERR)
  DEALLOCATE(evdf_flag_timestep, STAT = ALLOC_ERR)
 
  IF (Rank_of_process.EQ.0) THEN 
     PRINT '("### The program will create ",i4," snapshots ###")', N_of_all_snaps

! write moments of snapshot creation into the file
     OPEN (41, FILE = '_snapmoments.dat')
!                 "--****-----*******.*****----********----*"
     WRITE (41, '(" number       time(ns)       T_cntr  vdf-flag")')
     DO i = 1, N_of_all_snaps
        WRITE (41, '(2x,i4,5x,f13.5,4x,i8,4x,i1)') i, Tcntr_snapshot(i) * 1.0d9 * delta_t_s, Tcntr_snapshot(i), save_evdf_snapshot(i)
     END DO
     CLOSE (41, STATUS = 'KEEP')

  END IF
  
END SUBROUTINE INITIATE_SNAPSHOTS

!------------------------------------------------------------------------------------------------------------
! 
SUBROUTINE CREATE_SNAPSHOT

  USE ParallelOperationValues
  USE Snapshots
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries
  USE IonParticles, ONLY : N_spec

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER ALLOC_ERR
  INTEGER i, j, k
  INTEGER recsize, bufsize
  REAL, ALLOCATABLE :: rbufer(:)

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
!                                 ! ----x----I----x----I--
  CHARACTER(20) filename_Ni       ! _NNNN_Ni_s_m3_2D.bin
  CHARACTER(22) filename_Ji       ! _NNNN_JXi_s_Am2_2D.bin
  CHARACTER(21) filename_Vi       ! _NNNN_VXi_s_ms_2D.bin
  CHARACTER(21) filename_Wi       ! _NNNN_WXi_s_eV_2D.bin
  CHARACTER(21) filename_Ti       ! _NNNN_TXi_s_eV_2D.bin
!                                 ! ----x----I----x----I--
  CHARACTER(22) filename_Jsum     ! _NNNN_JXsum_Am2_2D.bin

  INTEGER s

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

  IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

! collect potential from field calculators ---------------------------------------------
     IF (cluster_rank_key.EQ.0) THEN

        ALLOCATE(cs_phi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)

! account for own contribution of the master as a field calculator
        DO j = indx_y_min, indx_y_max
           DO i = indx_x_min, indx_x_max
              cs_phi(i,j) = REAL(phi(i,j) * F_scale_V)
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
              cs_phi(field_calculator(k)%indx_x_min:field_calculator(k)%indx_x_max,j) = rbufer(pos1:pos2) 
           END DO

           DEALLOCATE(rbufer, STAT = ALLOC_ERR)

        END DO

     ELSE

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

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, field_master, Rank_of_process, MPI_COMM_WORLD, request, ierr) 

        DEALLOCATE(rbufer, STAT = ALLOC_ERR)
        
     END IF

  ELSE

! the cluster already knows the potential within its domain, use it -----------------

     IF (cluster_rank_key.EQ.0) THEN

        ALLOCATE(cs_phi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)

        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_phi(i,j) = REAL(c_phi(i,j) * F_scale_V)
           END DO
        END DO

     END IF

  END IF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (cluster_rank_key.EQ.0) THEN

     ALLOCATE(cs_EX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_EY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_EX(i,j) = REAL(EX(i,j) * E_scale_Vm)
        END DO
     END DO
  
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_EY(i,j) = REAL(EY(i,j) * E_scale_Vm)
        END DO
     END DO

     ALLOCATE(cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ALLOCATE(cs_JX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_JY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_JZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ALLOCATE(cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ALLOCATE(cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ALLOCATE(cs_TX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_TY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_TZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ALLOCATE(cs_JXsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_JYsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(cs_JZsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     IF (N_vdfbox_all.GT.0) THEN

        IF ((save_evdf_snapshot(current_snap).EQ.ONLY1D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
           ALLOCATE(evxdf(indx_v_min_e:indx_v_max_e, 1:N_vdfbox_all), STAT=ALLOC_ERR)
           ALLOCATE(evydf(indx_v_min_e:indx_v_max_e, 1:N_vdfbox_all), STAT=ALLOC_ERR)
           ALLOCATE(evzdf(indx_v_min_e:indx_v_max_e, 1:N_vdfbox_all), STAT=ALLOC_ERR)

           ALLOCATE(isvxdf(indx_v_min_i:indx_v_max_i, 1:N_vdfbox_all, 1:N_spec), STAT=ALLOC_ERR)
           ALLOCATE(isvydf(indx_v_min_i:indx_v_max_i, 1:N_vdfbox_all, 1:N_spec), STAT=ALLOC_ERR)
           ALLOCATE(isvzdf(indx_v_min_i:indx_v_max_i, 1:N_vdfbox_all, 1:N_spec), STAT=ALLOC_ERR)
        END IF

        IF ((save_evdf_snapshot(current_snap).EQ.ONLY2D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
           ALLOCATE(evxvydf(indx_v_min_e:indx_v_max_e, indx_v_min_e:indx_v_max_e, 1:N_vdfbox_all), STAT=ALLOC_ERR)
        END IF

     END IF

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

     IF (N_vdfbox_all.GT.0) THEN

        IF ((save_evdf_snapshot(current_snap).EQ.ONLY1D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
           ALLOCATE(evxdf(1,1), STAT=ALLOC_ERR)
           ALLOCATE(evydf(1,1), STAT=ALLOC_ERR)
           ALLOCATE(evzdf(1,1), STAT=ALLOC_ERR)

           ALLOCATE(isvxdf(1,1,1), STAT=ALLOC_ERR)
           ALLOCATE(isvydf(1,1,1), STAT=ALLOC_ERR)
           ALLOCATE(isvzdf(1,1,1), STAT=ALLOC_ERR)
        END IF

        IF ((save_evdf_snapshot(current_snap).EQ.ONLY2D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
           ALLOCATE(evxvydf(1,1,1), STAT=ALLOC_ERR)
        END IF

     END IF

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! write electrostatic field parameters to files ----------------------------

! potential

  IF (save_data(1)) THEN
     filename_F = '_NNNN_F_V_2D.bin'
     filename_F(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_phi, filename_F)
  END IF

! electric field components

  IF (save_data(2)) THEN
     filename_E = '_NNNN_EX_Vm_2D.bin'
     filename_E(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_EX, filename_E)
  END IF

  IF (save_data(3)) THEN
     filename_E = '_NNNN_EY_Vm_2D.bin'
     filename_E(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_EY, filename_E)
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! collect electron and ion velocity distribution function moments (density, flows, energies)

! write electron moments to files ------------------------------------------

  CALL COLLECT_ELECTRON_MOMENTS

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! number density

  IF (save_data(7)) THEN
     filename_Ne = '_NNNN_Ne_m3_2D.bin'
     filename_Ne(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_N, filename_Ne)
  END IF

! electric current along each coordinate direction

  IF (save_data(8)) THEN
     filename_Je = '_NNNN_JXe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_JX, filename_Je)
  END IF

  IF (save_data(9)) THEN
     filename_Je = '_NNNN_JYe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_JY, filename_Je)
  END IF

  IF (save_data(10)) THEN
     filename_Je = '_NNNN_JZe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_JZ, filename_Je)
  END IF

! average velocity along each coordinate direction

  IF (save_data(11)) THEN
     filename_Ve = '_NNNN_VXe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_VX, filename_Ve)
  END IF

  IF (save_data(12)) THEN
     filename_Ve = '_NNNN_VYe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_VY, filename_Ve)
  END IF

  IF (save_data(13)) THEN
     filename_Ve = '_NNNN_VZe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_VZ, filename_Ve)
  END IF

! average energy of FULL motion along each coordinate direction

  IF (save_data(14)) THEN
     filename_We = '_NNNN_WXe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_WX, filename_We)
  END IF

  IF (save_data(15)) THEN
     filename_We = '_NNNN_WYe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_WY, filename_We)
  END IF

  IF (save_data(16)) THEN
     filename_We = '_NNNN_WZe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_WZ, filename_We)
  END IF

! average energy of THERMAL motion along each coordinate direction

  IF (save_data(17)) THEN
     filename_Te = '_NNNN_TXe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_TX, filename_Te)
  END IF

  IF (save_data(18)) THEN
     filename_Te = '_NNNN_TYe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_TY, filename_Te)
  END IF

  IF (save_data(19)) THEN
     filename_Te = '_NNNN_TZe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_TZ, filename_Te)
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! write ion moments to files ------------------------------------------

  DO s = 1, N_spec

     CALL COLLECT_ION_MOMENTS(s)

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! number density

     IF (save_data(20)) THEN
        filename_Ni = '_NNNN_Ni_s_m3_2D.bin'
        filename_Ni(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Ni(10:10) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_N, filename_Ni)
     END IF

! electric current along each coordinate direction

     IF (save_data(21)) THEN
        filename_Ji = '_NNNN_JXi_s_Am2_2D.bin'
        filename_Ji(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Ji(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_JX, filename_Ji)
     END IF

     IF (save_data(22)) THEN
        filename_Ji = '_NNNN_JYi_s_Am2_2D.bin'
        filename_Ji(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Ji(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_JY, filename_Ji)
     END IF

     IF (save_data(23)) THEN
        filename_Ji = '_NNNN_JZi_s_Am2_2D.bin'
        filename_Ji(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Ji(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_JZ, filename_Ji)
     END IF

! average velocity along each coordinate direction

     IF (save_data(24)) THEN
        filename_Vi = '_NNNN_VXi_s_ms_2D.bin'
        filename_Vi(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Vi(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_VX, filename_Vi)
     END IF

     IF (save_data(25)) THEN
        filename_Vi = '_NNNN_VYi_s_ms_2D.bin'
        filename_Vi(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Vi(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_VY, filename_Vi)
     END IF

     IF (save_data(26)) THEN
        filename_Vi = '_NNNN_VZi_s_ms_2D.bin'
        filename_Vi(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Vi(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_VZ, filename_Vi)
     END IF

! average energy of FULL motion along each coordinate direction

     IF (save_data(27)) THEN
        filename_Wi = '_NNNN_WXi_s_eV_2D.bin'
        filename_Wi(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Wi(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_WX, filename_Wi)
     END IF

     IF (save_data(28)) THEN
        filename_Wi = '_NNNN_WYi_s_eV_2D.bin'
        filename_Wi(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Wi(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_WY, filename_Wi)
     END IF

     IF (save_data(29)) THEN
        filename_Wi = '_NNNN_WZi_s_eV_2D.bin'
        filename_Wi(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Wi(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_WZ, filename_Wi)
     END IF

! average energy of THERMAL motion along each coordinate direction

     IF (save_data(30)) THEN
        filename_Ti = '_NNNN_TXi_s_eV_2D.bin'
        filename_Ti(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Ti(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_TX, filename_Ti)
     END IF

     IF (save_data(31)) THEN
        filename_Ti = '_NNNN_TYi_s_eV_2D.bin'
        filename_Ti(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Ti(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_TY, filename_Ti)
     END IF

     IF (save_data(32)) THEN
        filename_Ti = '_NNNN_TZi_s_eV_2D.bin'
        filename_Ti(2:5) = convert_int_to_txt_string(current_snap, 4)
        filename_Ti(11:11) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_TZ, filename_Ti)
     END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  END DO

! FULL electric current (sum of electron and ion currents) along each coordinate direction

  IF (save_data(4)) THEN
     filename_Jsum = '_NNNN_JXsum_Am2_2D.bin'
     filename_Jsum(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_JXsum, filename_Jsum)
  END IF

  IF (save_data(5)) THEN
     filename_Jsum = '_NNNN_JYsum_Am2_2D.bin'
     filename_Jsum(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_JYsum, filename_Jsum)
  END IF

  IF (save_data(6)) THEN
     filename_Jsum = '_NNNN_JZsum_Am2_2D.bin'
     filename_Jsum(2:5) = convert_int_to_txt_string(current_snap, 4)
     CALL SAVE_GLOBAL_2D_ARRAY(cs_JZsum, filename_Jsum)
  END IF

! write electron velocity distribution function -----------------------

  CALL CALCULATE_ELECTRON_VDF

  CALL MPI_BARRIER(COMM_CLUSTER, ierr)

  CALL CALCULATE_ION_VDF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL SAVE_ALL_VDF1D

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  CALL SAVE_ALL_VDF2D

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  IF (cluster_rank_key.NE.0) THEN

     DEALLOCATE(cs_N, STAT=ALLOC_ERR)
     DEALLOCATE(cs_VX, STAT=ALLOC_ERR)
     DEALLOCATE(cs_VY, STAT=ALLOC_ERR)
     DEALLOCATE(cs_VZ, STAT=ALLOC_ERR)
     DEALLOCATE(cs_WX, STAT=ALLOC_ERR)
     DEALLOCATE(cs_WY, STAT=ALLOC_ERR)
     DEALLOCATE(cs_WZ, STAT=ALLOC_ERR)

     IF (ALLOCATED(evxdf)) DEALLOCATE(evxdf, STAT=ALLOC_ERR)
     IF (ALLOCATED(evydf)) DEALLOCATE(evydf, STAT=ALLOC_ERR)
     IF (ALLOCATED(evzdf)) DEALLOCATE(evzdf, STAT=ALLOC_ERR)

     IF (ALLOCATED(evxvydf)) DEALLOCATE(evxvydf, STAT=ALLOC_ERR)

     IF (ALLOCATED(isvxdf)) DEALLOCATE(isvxdf, STAT=ALLOC_ERR)
     IF (ALLOCATED(isvydf)) DEALLOCATE(isvydf, STAT=ALLOC_ERR)
     IF (ALLOCATED(isvzdf)) DEALLOCATE(isvzdf, STAT=ALLOC_ERR)

  ELSE

     DEALLOCATE(cs_phi, STAT = ALLOC_ERR)
     DEALLOCATE(cs_EX, STAT = ALLOC_ERR)
     DEALLOCATE(cs_EY, STAT = ALLOC_ERR)

     DEALLOCATE(cs_N, STAT = ALLOC_ERR)

     DEALLOCATE(cs_JX, STAT = ALLOC_ERR)
     DEALLOCATE(cs_JY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_JZ, STAT = ALLOC_ERR)

     DEALLOCATE(cs_VX, STAT = ALLOC_ERR)
     DEALLOCATE(cs_VY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_VZ, STAT = ALLOC_ERR)

     DEALLOCATE(cs_TX, STAT = ALLOC_ERR)
     DEALLOCATE(cs_TY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_TZ, STAT = ALLOC_ERR)

     DEALLOCATE(cs_JXsum, STAT = ALLOC_ERR)
     DEALLOCATE(cs_JYsum, STAT = ALLOC_ERR)
     DEALLOCATE(cs_JZsum, STAT = ALLOC_ERR)

     IF (ALLOCATED(evxdf)) DEALLOCATE(evxdf, STAT = ALLOC_ERR)
     IF (ALLOCATED(evydf)) DEALLOCATE(evydf, STAT = ALLOC_ERR)
     IF (ALLOCATED(evzdf)) DEALLOCATE(evzdf, STAT = ALLOC_ERR)

     IF (ALLOCATED(evxvydf)) DEALLOCATE(evxvydf, STAT = ALLOC_ERR)

     IF (ALLOCATED(isvxdf)) DEALLOCATE(isvxdf, STAT = ALLOC_ERR)
     IF (ALLOCATED(isvydf)) DEALLOCATE(isvydf, STAT = ALLOC_ERR)
     IF (ALLOCATED(isvzdf)) DEALLOCATE(isvzdf, STAT = ALLOC_ERR)

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (Rank_of_process.EQ.0) PRINT '(/2x,"### ^^^^^^^^^^^^^^^^^^^^ Snapshot ",i4," completed :) ^^^^^^^^^^^^^^^^^^^ ###")', current_snap

  current_snap = current_snap + 1           ! increase the snapshots counter 

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

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL arr(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max)

  CHARACTER*(*) filename

  INTEGER ierr
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

!  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, cluster_rank_key, Rank_of_process, COMM_SAVE_CLUSTERS, ierr)
!  CALL MPI_COMM_RANK(COMM_SAVE_CLUSTERS, Rank_save_clusters, ierr)
!  CALL MPI_COMM_SIZE(COMM_SAVE_CLUSTERS, N_processes_save_clusters, ierr)

  IF (cluster_rank_key.EQ.0) THEN

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

!if (Rank_of_process.eq.0) print '("Process ",i4,"/",i4," shift ",i8," B")', Rank_of_process, Rank_horizontal,  shifthead

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

!if (Rank_of_process.eq.0) print '("Process ",i4,"/",i4," shift ",i8," B")', Rank_of_process, Rank_horizontal,  shifthead

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

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! all processes must participate in destroying the communicator, I guess...
!  CALL MPI_COMM_FREE(COMM_SAVE_CLUSTERS, ierr)

  IF (Rank_of_process.EQ.0) PRINT '("created file ",A40)', filename
!stop

END SUBROUTINE SAVE_GLOBAL_2D_ARRAY

!---------------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_ALL_VDF1D

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles
  USE Snapshots

  IMPLICIT NONE

  INCLUDE 'mpif.h'

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

  IMPLICIT NONE

  INCLUDE 'mpif.h'

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

!---------------------------------------------------------------------------------------------------
!
SUBROUTINE FINISH_SNAPSHOTS

  USE Snapshots
  IMPLICIT NONE

  INTEGER DEALLOC_ERR

  IF (ALLOCATED(    Tcntr_snapshot)) DEALLOCATE(    Tcntr_snapshot, STAT=DEALLOC_ERR)
  IF (ALLOCATED(save_evdf_snapshot)) DEALLOCATE(save_evdf_snapshot, STAT=DEALLOC_ERR)

END SUBROUTINE FINISH_SNAPSHOTS
