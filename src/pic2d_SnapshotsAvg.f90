
!--------------------

SUBROUTINE INITIATE_AVERAGED_SNAPSHOTS

  USE ParallelOperationValues
  USE AvgSnapshots
  USE CurrentProblemValues, ONLY : delta_t_s, N_subcycles, Start_T_cntr
  USE Checkpoints, ONLY : use_checkpoint
  USE MCCollisions, ONLY : N_neutral_spec, collision_e_neutral, en_collisions_turned_off

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  LOGICAL exists
  CHARACTER (1) buf

  INTEGER ios
  INTEGER saveflagi(1:39)           ! integer flags used to set values of logical flags controlling saving of data files

  INTEGER N_of_snap_groups          ! number of sets of snapshots, read from file
  INTEGER i, n, p, count

  REAL(8) Rqst_snap_start_ns        ! requested start of current set of snapshots [ns], read from file
  REAL(8) Rqst_snap_finish_ns       ! requested finish of current set of snapshots [ns], read from file 
  INTEGER Rqst_n_of_snaps           ! requested number of snapshots in current set, read from file
  
! temporary arrays
  INTEGER, ALLOCATABLE :: timestep_begin(:)   ! moments when average data collection for a snapshot begins
  INTEGER, ALLOCATABLE :: timestep_end(:)     ! moments when average data collection for a snapshot ends and the snapshot is saved
  INTEGER ALLOC_ERR

  INTEGER T1, T2, T2prev
  INTEGER large_step
  INTEGER time_begin, time_end

! default values ensure that if init_avgsnapshots.dat is not found 
! procedures COLLECT_DATA_FOR_AVERAGED_SNAPSHOT and CREATE_AVERAGED_SNAPSHOT do nothing
  N_of_all_avgsnaps = 0
  current_avgsnap = 1
  save_avg_data = .FALSE.
!  avg_data_collection_offset = -1

! read / write the data file 
  INQUIRE (FILE = 'init_avgsnapshots.dat', EXIST = exists)

  IF (.NOT.exists) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### File init_avgsnapshots.dat not found. Time-averaged snapshots will not be created. ###")'
     RETURN
  END IF

  IF (Rank_of_process.EQ.0) PRINT '("### File init_avgsnapshots.dat is found. Reading the data file... ###")'

  OPEN (9, FILE = 'init_avgsnapshots.dat')

  saveflagi = 0

  READ(9, '(A1)') buf  !--- save 2D maps of the following TIME-AVERAGED parameters? (1=yes, 0=no)
  READ(9, '(A1)') buf  !-----F----EX----EY--JXsum--JYsum--JZsum--e-n collision frequencies   (type flag values below)
  READ(9, *, iostat = ios) saveflagi(1:6) , saveflagi(39)

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### WARNING iostat ",i8," : in file init_avgsnapshots.dat first set of flags is probably incomplete, missing flag(s) set to zero ###")', ios
     BACKSPACE(9)
  END IF

  READ(9, '(A1)') buf  !----Ne----JXe---JYe---JZe---VXe---VYe---VZe---WXe---WYe---WZe---TXe---TYe---TZe---QXe---QYe---QZe  (type flag values below)
  READ(9, *, iostat = ios) saveflagi(7:22)

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### WARNING iostat ",i8," : in file init_avgsnapshots.dat second set of flags is probably incomplete, missing flag(s) set to zero ###")', ios
     BACKSPACE(9)
  END IF

  READ(9, '(A1)') buf  !----Ni----JXi---JYi---JZi---VXi---VYi---VZi---WXi---WYi---WZi---TXi---TYi---TZi---QXi---QYi---QZi  (type flag values below)
  READ (9, *, iostat = ios) saveflagi(23:38)

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### WARNING iostat ",i8," : in file init_avgsnapshots.dat third set of flags is probably incomplete, missing flag(s) set to zero ###")', ios
     BACKSPACE(9)
  END IF

  IF (Rank_of_process.EQ.0) THEN
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: saveflagi(1:6,39)     = ",7(2x,i2))', saveflagi(1:6) , saveflagi(39)
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: saveflagi(7:22)       = ",16(2x,i2))', saveflagi(7:22)
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: saveflagi(23:38)      = ",16(2x,i2))', saveflagi(7:22)
  END IF

  DO i = 1, 39
     IF (saveflagi(i).GT.0) save_avg_data(i) = .TRUE.
  END DO

  IF (en_collisions_turned_off) save_avg_data(39) = .FALSE.

! if electron-neutral collision frequencies are requested, check that at least one active collisional process has save_collfreq_2d == .TRUE.
  IF (save_avg_data(39)) THEN
     count = 0
     DO n = 1, N_neutral_spec
        DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
           IF (collision_e_neutral(n)%colproc_info(p)%save_collfreq_2d) count = count + 1
        END DO
     END DO
     IF (count.EQ.0) THEN
        save_avg_data(39) = .FALSE.
        IF (Rank_of_process.EQ.0) PRINT '("### WARNING :: Saving e-n collision frequencies was requested in init_avgsnapshots.dat but not set in init_neutral_AAAAAA.dat so it is turned off ###")'
     END IF
  END IF

  IF (Rank_of_process.EQ.0) THEN
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: save_avg_data(1:6,39) = ",7(3x,L1))', save_avg_data(1:6) , save_avg_data(39)
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: save_avg_data(7:22)   = ",16(3x,L1))', save_avg_data(7:22)
     PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: save_avg_data(23:38)  = ",16(3x,L1))', save_avg_data(7:22)
  END IF

  READ (9, '(A1)') buf !--- number of groups of snapshots (>0, no averaged snapshots if <=0), type below
  READ (9, *, iostat=ios) N_of_snap_groups

  IF (ios.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### WARNING iostat ",i8," : in file init_avgsnapshots.dat cannot read the number of snapshot groups, set it to zero ###")', ios
     N_of_snap_groups = 0
  END IF

  IF (N_of_snap_groups.LE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### ### Time-averaged snapshots will not be created. ### ###")'
     CLOSE (9, STATUS = 'KEEP')
     RETURN
  END IF

  READ (9, '(A1)') buf !--- For each group, type below start (ns), finish (ns), number of snapshots (>=0)

  T2prev = -1

  ALLOCATE(timestep_begin(1:9999), STAT = ALLOC_ERR)
  ALLOCATE(timestep_end(1:9999), STAT = ALLOC_ERR)
     
  DO i = 1, N_of_snap_groups
! read the parameters of current set of snapshot from the data file
     READ (9, *, iostat = ios) Rqst_snap_start_ns, Rqst_snap_finish_ns, Rqst_n_of_snaps

     IF (ios.NE.0) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### WARNING iostat ",i8," : in file init_avgsnapshots.dat while reading snapshot group ",i3,", skip ###")', ios, i
        CYCLE
     END IF

! try the next group of snapshots if the current group snapshot number is zero
     IF (Rqst_n_of_snaps.LT.1) CYCLE
     IF (Rqst_snap_start_ns.GE.Rqst_snap_finish_ns) CYCLE

! timestep when averaging for the first snapshot in the group begins 
     T1 = Rqst_snap_start_ns / (delta_t_s * 1.0d9)

! timestep when averaging for the last snapshot in the group ends
     T2 = Rqst_snap_finish_ns / (delta_t_s * 1.0d9)

     T1 = (T1 / N_subcycles) * N_subcycles     ! + avg_data_collection_offset
     T2 = (T2 / N_subcycles) * N_subcycles     ! + avg_data_collection_offset

     T1 = MAX(T1, T2prev)

     large_step = (T2 - T1) / Rqst_n_of_snaps
     large_step = (large_step / N_subcycles) * N_subcycles
     IF (large_step.LE.N_subcycles) CYCLE

! for all possible snapshots of the current set
     DO n = 1, Rqst_n_of_snaps  !Fact_n_of_snaps
! Calculate and save the snapshot moment in the temporary array
        time_begin = T1 + (n - 1) * large_step
        time_end = time_begin + large_step - 1 !N_subcycles
        IF (time_begin.GE.T2) EXIT
        IF (time_end.GT.T2) EXIT  !time_end = T2

        N_of_all_avgsnaps = N_of_all_avgsnaps + 1
        timestep_begin(N_of_all_avgsnaps) =  time_begin
        timestep_end(N_of_all_avgsnaps) = time_end
     END DO        ! end of cycle over snapshots in one set

     IF (N_of_all_avgsnaps.LE.0) CYCLE
     T2prev = timestep_end(N_of_all_avgsnaps)

  END DO           ! end of cycle over sets of snapshots     

  CLOSE (9, STATUS = 'KEEP')

  IF (N_of_all_avgsnaps.EQ.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### Time-averaged snapshots will NOT be created ... ###")'
! cleanup
     DEALLOCATE(timestep_begin, STAT = ALLOC_ERR)
     DEALLOCATE(timestep_end, STAT = ALLOC_ERR)
     RETURN
  END IF 

! if we are here, snapshots will be created ...

! allocate the array of moments of snapshots
  ALLOCATE(avgsnapshot(1:N_of_all_avgsnaps), STAT=ALLOC_ERR)

! move the calculated snapshot moments from the temporary array to the allocated array 
  DO n = 1, N_of_all_avgsnaps
     avgsnapshot(n)%T_cntr_begin = timestep_begin(n)
     avgsnapshot(n)%T_cntr_end = timestep_end(n)
  END DO

  DEALLOCATE(timestep_begin, STAT = ALLOC_ERR)
  DEALLOCATE(timestep_end, STAT = ALLOC_ERR)

  IF (Rank_of_process.EQ.0) THEN 
     PRINT '("### The program will create ",i4," averaged snapshots ###")', N_of_all_avgsnaps

! write moments of snapshot creation into the file
     OPEN (41, FILE = '_avg_snapmoments.dat')
     WRITE (41, '(" number   start_time(ns)   end_time(ns)   start_T_cntr   end_T_cntr   N_of_avg_points,e/i")')
     DO n = 1, N_of_all_avgsnaps
        WRITE (41, '(2x,i4,2x,2(2x,f13.5),2x,2(2x,i9),2(4x,i6))') &
             & n, &
             & avgsnapshot(n)%T_cntr_begin * 1.0d9 * delta_t_s, &
             & avgsnapshot(n)%T_cntr_end * 1.0d9 * delta_t_s, &
             & avgsnapshot(n)%T_cntr_begin, &
             & avgsnapshot(n)%T_cntr_end, &
             & avgsnapshot(n)%T_cntr_end - avgsnapshot(n)%T_cntr_begin + 1, &
             & (avgsnapshot(n)%T_cntr_end - avgsnapshot(n)%T_cntr_begin + 1) / N_subcycles
     END DO
     CLOSE (41, STATUS = 'KEEP')
  END IF

! consistency check 1
  DO n = 1, N_of_all_avgsnaps
     IF (avgsnapshot(n)%T_cntr_begin.GE.avgsnapshot(n)%T_cntr_end) THEN
        IF (Rank_of_process.EQ.0) PRINT '("ERROR-1 in INITIATE_AVERAGED_SNAPSHOTS, snapshot "i4," begins at ",i9," and ends at ",i9)', &
             & n, avgsnapshot(n)%T_cntr_begin, avgsnapshot(n)%T_cntr_end
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF
  END DO

  IF (Rank_of_process.EQ.0) PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: average snapshot timing passed consistency check 1")'

! consistency check 2
  DO n = 1, N_of_all_avgsnaps-1
     IF (avgsnapshot(n)%T_cntr_end.GE.avgsnapshot(n+1)%T_cntr_begin) THEN
        IF (Rank_of_process.EQ.0) PRINT '("ERROR-2 in INITIATE_AVERAGED_SNAPSHOTS, snapshot ",i4,2x,i9,2x,i9," overlaps with snapshot "i4,2x,i9,2x,i9)', &
             & n,   avgsnapshot(n)%T_cntr_begin,   avgsnapshot(n)%T_cntr_end, &
             & n+1, avgsnapshot(n+1)%T_cntr_begin, avgsnapshot(n+1)%T_cntr_end
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF
  END DO

  IF (Rank_of_process.EQ.0) PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: average snapshot timing passed consistency check 2")'

! note, since we are here, N_of_all_avgsnaps is not zero
! overwrite default value of the first snapshot number if the system is initialized using a checkpoint

  IF (use_checkpoint.EQ.1) THEN
     current_avgsnap = N_of_all_avgsnaps + 1  ! default assumption is that all average snapshots are already created
                                              ! this can change in the cycle below (also, see the comment after the cycle)
     DO n = 1, N_of_all_avgsnaps
        IF (avgsnapshot(n)%T_cntr_begin.GE.Start_T_cntr) THEN
           current_avgsnap = n
           IF (Rank_of_process.EQ.0) PRINT '("INITIATE_AVERAGED_SNAPSHOTS :: adjusted number of the first snapshot ",i4)', current_avgsnap
           EXIT
        END IF
     END DO
! note that when a simulation is restarted from a checkpoint, it is possible that 
! Start_T_cntr > maximal avgsnapshot%T_cntr_begin
! then the cycle above does not change current_avgsnap
! if current_avgsnap equals to 1 before the cycle then the call of COLLECT_F_EX_EY_FOR_AVERAGED_SNAPSHOT
! proceeds without allocating cs_avg_* arrays which results in memory error
! this bug was found and fixed by Willca Villafana
  END IF

END SUBROUTINE INITIATE_AVERAGED_SNAPSHOTS

!------------------------------------------------------------------------------------------------------------
! 
SUBROUTINE COLLECT_F_EX_EY_FOR_AVERAGED_SNAPSHOT

  USE ParallelOperationValues
  USE AvgSnapshots
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries
  USE IonParticles, ONLY : N_spec
  USE MCCollisions, ONLY : N_neutral_spec, collision_e_neutral
  USE Snapshots, ONLY : diagnostics_neutral


  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER n, p

  INTEGER ALLOC_ERR

  REAL, ALLOCATABLE :: cs_phi(:,:) 

  INTEGER i, j, k
  INTEGER recsize, bufsize
  REAL, ALLOCATABLE :: rbufer(:)

  INTEGER pos1, pos2

! quit if all snapshots were created or if due to any reasons the snapshot counter is 
! larger than the declared number of snapshots (e.g., when no snapshots are requested) 
  IF (current_avgsnap.GT.N_of_all_avgsnaps) RETURN

  IF (T_cntr.LT.avgsnapshot(current_avgsnap)%T_cntr_begin) RETURN

  IF ((T_cntr.EQ.avgsnapshot(current_avgsnap)%T_cntr_begin).AND.(cluster_rank_key.EQ.0)) THEN

     IF (save_avg_data(1)) THEN
        ALLOCATE(cs_avg_phi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_phi = 0.0
     END IF

     IF (save_avg_data(2)) THEN 
        ALLOCATE(cs_avg_EX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_EX = 0.0
     END IF

     IF (save_avg_data(3)) THEN
        ALLOCATE(cs_avg_EY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_EY = 0.0
     END IF

     IF (save_avg_data(7)) THEN 
        ALLOCATE(cs_avg_Ne(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_Ne = 0.0
     END IF

     IF (save_avg_data(4).OR.save_avg_data(8)) THEN
        ALLOCATE(cs_avg_JXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_JXe = 0.0
     END IF

     IF (save_avg_data(5).OR.save_avg_data(9)) THEN
        ALLOCATE(cs_avg_JYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_JYe = 0.0
     END IF

     IF (save_avg_data(6).OR.save_avg_data(10)) THEN
        ALLOCATE(cs_avg_JZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_JZe = 0.0
     END IF

     IF (save_avg_data(11)) THEN
        ALLOCATE(cs_avg_VXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_VXe = 0.0
     END IF

     IF (save_avg_data(12)) THEN
        ALLOCATE(cs_avg_VYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_VYe = 0.0
     END IF

     IF (save_avg_data(13)) THEN
        ALLOCATE(cs_avg_VZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_VZe = 0.0
     END IF

     IF (save_avg_data(14)) THEN
        ALLOCATE(cs_avg_WXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_WXe = 0.0
     END IF

     IF (save_avg_data(15)) THEN
        ALLOCATE(cs_avg_WYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_WYe = 0.0
     END IF

     IF (save_avg_data(16)) THEN
        ALLOCATE(cs_avg_WZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_WZe = 0.0
     END IF

     IF (save_avg_data(17)) THEN
        ALLOCATE(cs_avg_TXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_TXe = 0.0
     END IF

     IF (save_avg_data(18)) THEN
        ALLOCATE(cs_avg_TYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_TYe = 0.0
     END IF

     IF (save_avg_data(19)) THEN
        ALLOCATE(cs_avg_TZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_TZe = 0.0
     END IF

     IF (save_avg_data(20)) THEN
        ALLOCATE(cs_avg_QXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_QXe = 0.0
     END IF

     IF (save_avg_data(21)) THEN
        ALLOCATE(cs_avg_QYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_QYe = 0.0
     END IF

     IF (save_avg_data(22)) THEN
        ALLOCATE(cs_avg_QZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
        cs_avg_QZe = 0.0
     END IF

     IF (save_avg_data(23)) THEN
        ALLOCATE(cs_avg_Ni(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_Ni = 0.0
     END IF

     IF (save_avg_data(4).OR.save_avg_data(24)) THEN
        ALLOCATE(cs_avg_JXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_JXi = 0.0
     END IF

     IF (save_avg_data(5).OR.save_avg_data(25)) THEN
        ALLOCATE(cs_avg_JYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_JYi = 0.0
     END IF

     IF (save_avg_data(6).OR.save_avg_data(26)) THEN
        ALLOCATE(cs_avg_JZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_JZi = 0.0
     END IF

     IF (save_avg_data(27)) THEN
        ALLOCATE(cs_avg_VXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_VXi = 0.0
     END IF

     IF (save_avg_data(28)) THEN
        ALLOCATE(cs_avg_VYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_VYi = 0.0
     END IF

     IF (save_avg_data(29)) THEN
        ALLOCATE(cs_avg_VZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_VZi = 0.0
     END IF

     IF (save_avg_data(30)) THEN
        ALLOCATE(cs_avg_WXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_WXi = 0.0
     END IF

     IF (save_avg_data(31)) THEN
        ALLOCATE(cs_avg_WYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_WYi = 0.0
     END IF

     IF (save_avg_data(32)) THEN
        ALLOCATE(cs_avg_WZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_WZi = 0.0
     END IF

     IF (save_avg_data(33)) THEN
        ALLOCATE(cs_avg_TXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_TXi = 0.0
     END IF

     IF (save_avg_data(34)) THEN
        ALLOCATE(cs_avg_TYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_TYi = 0.0
     END IF

     IF (save_avg_data(35)) THEN
        ALLOCATE(cs_avg_TZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_TZi = 0.0
     END IF

     IF (save_avg_data(36)) THEN
        ALLOCATE(cs_avg_QXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_QXi = 0.0
     END IF

     IF (save_avg_data(37)) THEN
        ALLOCATE(cs_avg_QYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_QYi = 0.0
     END IF

     IF (save_avg_data(38)) THEN
        ALLOCATE(cs_avg_QZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, 1:N_spec), STAT = ALLOC_ERR)
        cs_avg_QZi = 0.0
     END IF

     IF (save_avg_data(39)) THEN
        DO n = 1, N_neutral_spec
           DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
              IF (.NOT.collision_e_neutral(n)%colproc_info(p)%save_collfreq_2d) CYCLE
              ALLOCATE(diagnostics_neutral(n)%activated_collision(p)%coll_freq_local(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
              diagnostics_neutral(n)%activated_collision(p)%coll_freq_local = 0.0
           END DO
        END DO
     END IF

  END IF   !### IF ((T_cntr.EQ.avgsnapshot(current_avgsnap)%T_cntr_begin).AND.(cluster_rank_key.EQ.0)) THEN

! potential

  IF (save_avg_data(1)) THEN

     IF (cluster_rank_key.EQ.0) THEN

        ALLOCATE(cs_phi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)

        IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN
! PETSc-based field solver
! account for own contribution of the master as a field calculator
           DO j = indx_y_min, indx_y_max
              DO i = indx_x_min, indx_x_max
                 cs_phi(i,j) = REAL(phi(i,j))
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
! FFT-based field solver
! the cluster already knows the potential within its domain, use it -----------------
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 cs_phi(i,j) = REAL(c_phi(i,j))
              END DO
           END DO
        END IF   !### IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

        cs_avg_phi = cs_avg_phi + cs_phi
! cleanup
        DEALLOCATE(cs_phi, STAT = ALLOC_ERR)

     ELSE   !### IF (cluster_rank_key.EQ.0) THEN

        IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

           recsize = indx_x_max - indx_x_min + 1
           bufsize = recsize * (indx_y_max - indx_y_min + 1)

           ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)

           pos1 = 1 - indx_x_min
           DO j = indx_y_min, indx_y_max
              DO i = indx_x_min, indx_x_max
                 rbufer(pos1 + i) = REAL(phi(i, j))
              END DO
              pos1 = pos1 + recsize
           END DO

           CALL MPI_SEND(rbufer, bufsize, MPI_REAL, field_master, Rank_of_process, MPI_COMM_WORLD, request, ierr) 

           DEALLOCATE(rbufer, STAT = ALLOC_ERR)
        
        END IF

     END IF   !###   IF (cluster_rank_key.EQ.0) THEN

  END IF   !###   IF (save_avg_data(1)) THEN

  IF (cluster_rank_key.NE.0) THEN
! non-master processes do nothing and can leave
     RETURN
  END IF

! electric field components

  IF (save_avg_data(2)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_EX(i,j) = cs_avg_EX(i,j) + REAL(EX(i,j))
        END DO
     END DO
  END IF

  IF (save_avg_data(3)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_EY(i,j) = cs_avg_EY(i,j) + REAL(EY(i,j))
        END DO
     END DO
  END IF

END SUBROUTINE COLLECT_F_EX_EY_FOR_AVERAGED_SNAPSHOT

!------------------------------------------------------------------------------------------------------------
! 
SUBROUTINE COLLECT_ELECTRON_DATA_FOR_AVERAGED_SNAPSHOT

  USE AvgSnapshots
  USE Snapshots, ONLY : cs_N, cs_VX, cs_VY, cs_VZ, cs_WX, cs_WY, cs_WZ, cs_VXVY, cs_VXVZ, cs_VYVZ, cs_QX, cs_QY, cs_QZ
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INTEGER i, j

! checks are performed before calling this procedure
! arrays cs_* are allocated before calling this procedure

! electron VDF moments are already collected in ADVANCE_ELECTRONS_AND_COLLECT_MOMENTS

! electric current along each coordinate direction

  IF (save_avg_data(4).OR.save_avg_data(8)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JXe(i,j) = cs_avg_JXe(i,j) - cs_N(i,j) * cs_VX(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(5).OR.save_avg_data(9)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JYe(i,j) = cs_avg_JYe(i,j) - cs_N(i,j) * cs_VY(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(6).OR.save_avg_data(10)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JZe(i,j) = cs_avg_JZe(i,j) - cs_N(i,j) * cs_VZ(i,j)
        END DO
     END DO
  END IF

! number density

  IF (save_avg_data(7)) THEN
     cs_avg_Ne = cs_avg_Ne + cs_N
  END IF

! average velocity along each coordinate direction

  IF (save_avg_data(11)) THEN
     cs_avg_VXe = cs_avg_VXe + cs_VX
  END IF

  IF (save_avg_data(12)) THEN
     cs_avg_VYe = cs_avg_VYe + cs_VY
  END IF

  IF (save_avg_data(13)) THEN
     cs_avg_VZe = cs_avg_VZe + cs_VZ
  END IF

! average energy of motion along each coordinate direction

  IF (save_avg_data(14)) THEN
     cs_avg_WXe = cs_avg_WXe + cs_WX
  END IF

  IF (save_avg_data(15)) THEN
     cs_avg_WYe = cs_avg_WYe + cs_WY
  END IF

  IF (save_avg_data(16)) THEN
     cs_avg_WZe = cs_avg_WZe + cs_WZ
  END IF

! temperature along each coordinate direction

  IF (save_avg_data(17)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TXe(i,j) = cs_avg_TXe(i,j) + MAX(0.0, cs_WX(i,j) - cs_VX(i,j)**2)
        END DO
     END DO
  END IF

  IF (save_avg_data(18)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TYe(i,j) = cs_avg_TYe(i,j) + MAX(0.0, cs_WY(i,j) - cs_VY(i,j)**2)
        END DO
     END DO
  END IF

  IF (save_avg_data(19)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TZe(i,j) = cs_avg_TZe(i,j) + MAX(0.0, cs_WZ(i,j) - cs_VZ(i,j)**2)
        END DO
     END DO
  END IF

! heat flow along each coordinate direction

  IF (save_avg_data(20)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QXe(i,j) = cs_avg_QXe(i,j) + cs_N(i,j) * cs_QX(i,j)
        END DO
     END DO
  END IF
     
  IF (save_avg_data(21)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QYe(i,j) = cs_avg_QYe(i,j) + cs_N(i,j) * cs_QY(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(22)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QZe(i,j) = cs_avg_QZe(i,j) + cs_N(i,j) * cs_QZ(i,j)
        END DO
     END DO
  END IF

END SUBROUTINE COLLECT_ELECTRON_DATA_FOR_AVERAGED_SNAPSHOT

!------------------------------------------------------------------------------------------------------------
! 
SUBROUTINE COLLECT_ION_DATA_FOR_AVERAGED_SNAPSHOT(s)

  USE AvgSnapshots
  USE Snapshots, ONLY : cs_N, cs_VX, cs_VY, cs_VZ, cs_WX, cs_WY, cs_WZ, cs_VXVY, cs_VXVZ, cs_VYVZ, cs_QX, cs_QY, cs_QZ
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: s

  INTEGER i, j

! checks are performed before calling this procedure
! arrays cs_* are allocated before calling this procedure

! VDF moments for ion species s are already collected in ADVANCE_IONS_AND_COLLECT_MOMENTS

! number density

  IF (save_avg_data(23)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_Ni(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_Ni(c_indx_x_min:c_indx_x_max,j,s) + cs_N(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

! electric current along each coordinate direction

  IF (save_avg_data(4).OR.save_avg_data(24)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JXi(i,j,s) = cs_avg_JXi(i,j,s) + cs_N(i,j) * cs_VX(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(5).OR.save_avg_data(25)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JYi(i,j,s) = cs_avg_JYi(i,j,s) + cs_N(i,j) * cs_VY(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(6).OR.save_avg_data(26)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_JZi(i,j,s) = cs_avg_JZi(i,j,s) + cs_N(i,j) * cs_VZ(i,j)
        END DO
     END DO
  END IF

! average velocity along each coordinate direction

  IF (save_avg_data(27)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_VXi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_VXi(c_indx_x_min:c_indx_x_max,j,s) + cs_VX(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

  IF (save_avg_data(28)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_VYi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_VYi(c_indx_x_min:c_indx_x_max,j,s) + cs_VY(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

  IF (save_avg_data(29)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_VZi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_VZi(c_indx_x_min:c_indx_x_max,j,s) + cs_VZ(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

! average energy of motion along each coordinate direction

  IF (save_avg_data(30)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_WXi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_WXi(c_indx_x_min:c_indx_x_max,j,s) + cs_WX(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

  IF (save_avg_data(31)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_WYi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_WYi(c_indx_x_min:c_indx_x_max,j,s) + cs_WY(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

  IF (save_avg_data(32)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        cs_avg_WZi(c_indx_x_min:c_indx_x_max,j,s) = cs_avg_WZi(c_indx_x_min:c_indx_x_max,j,s) + cs_WZ(c_indx_x_min:c_indx_x_max,j)
     END DO
  END IF

! temperature along each coordinate direction

  IF (save_avg_data(33)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TXi(i,j,s) = cs_avg_TXi(i,j,s) + MAX(0.0, cs_WX(i,j) - cs_VX(i,j)**2)
        END DO
     END DO
  END IF

  IF (save_avg_data(34)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TYi(i,j,s) = cs_avg_TYi(i,j,s) + MAX(0.0, cs_WY(i,j) - cs_VY(i,j)**2)
        END DO
     END DO
  END IF

  IF (save_avg_data(35)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_TZi(i,j,s) = cs_avg_TZi(i,j,s) + MAX(0.0, cs_WZ(i,j) - cs_VZ(i,j)**2)
        END DO
     END DO
  END IF

! heat flow along each coordinate direction

  IF (save_avg_data(36)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QXi(i,j,s) = cs_avg_QXi(i,j,s) + cs_N(i,j) * cs_QX(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(37)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QYi(i,j,s) = cs_avg_QYi(i,j,s) + cs_N(i,j) * cs_QY(i,j)
        END DO
     END DO
  END IF

  IF (save_avg_data(38)) THEN
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           cs_avg_QZi(i,j,s) = cs_avg_QZi(i,j,s) + cs_N(i,j) * cs_QZ(i,j)
        END DO
     END DO
  END IF

END SUBROUTINE COLLECT_ION_DATA_FOR_AVERAGED_SNAPSHOT

!------------------------------------------------------------------------------------------------------------
! 
SUBROUTINE CREATE_AVERAGED_SNAPSHOT

  USE ParallelOperationValues
  USE AvgSnapshots
  USE CurrentProblemValues, ONLY : N_subcycles, F_scale_V, E_scale_Vm, current_factor_Am2, N_scale_part_m3, V_scale_ms, &
                                 & energy_factor_eV, temperature_factor_eV, heat_flow_factor_Wm2, T_cntr, delta_t_s
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : N_spec, Qs, Ms
  USE MCCollisions, ONLY : N_neutral_spec, neutral, collision_e_neutral
  USE Snapshots, ONLY : diagnostics_neutral

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER N_averaged_timesteps, N_averaged_timesteps_i  ! different for electrons and ions
  REAL avg_factor
  INTEGER ALLOC_ERR

  REAL, ALLOCATABLE :: cs_avg_Jsum(:,:)

!                                 ! ----x----I----x----I--
  CHARACTER(20) filename_F        ! _NNNN_avg_F_V_2D.bin
  CHARACTER(22) filename_E        ! _NNNN_avg_EX_Vm_2D.bin
!                                 ! ----x----I----x----I----x-
  CHARACTER(26) filename_Jsum     ! _NNNN_avg_JXsum_Am2_2D.bin
!                                 ! ----x----I----x----I----
  CHARACTER(22) filename_Ne       ! _NNNN_avg_Ne_m3_2D.bin
  CHARACTER(24) filename_Je       ! _NNNN_avg_JXe_Am2_2D.bin
  CHARACTER(23) filename_Ve       ! _NNNN_avg_VXe_ms_2D.bin
  CHARACTER(23) filename_We       ! _NNNN_avg_WXe_eV_2D.bin
  CHARACTER(23) filename_Te       ! _NNNN_avg_TXe_eV_2D.bin
  CHARACTER(24) filename_Qe       ! _NNNN_avg_QXe_Wm2_2D.bin
!                                 ! ----x----I----x----I----x-
  CHARACTER(24) filename_Ni       ! _NNNN_avg_Ni_s_m3_2D.bin
  CHARACTER(26) filename_Ji       ! _NNNN_avg_JXi_s_Am2_2D.bin
  CHARACTER(25) filename_Vi       ! _NNNN_avg_VXi_s_ms_2D.bin
  CHARACTER(25) filename_Wi       ! _NNNN_avg_WXi_s_eV_2D.bin
  CHARACTER(25) filename_Ti       ! _NNNN_avg_TXi_s_eV_2D.bin
  CHARACTER(26) filename_Qi       ! _NNNN_avg_QXi_s_Wm2_2D.bin
!                                   ----x----I----x----I----x----I----x----I----x
  CHARACTER(45) filename_encoll   ! _NNNN_avg_frequency_e_n_AAAAAA_coll_id_NN.bin

  INTEGER s, i, j, n, p

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

! quit if all snapshots were created or if due to any reasons the snapshot counter is 
! larger than the declared number of snapshots (e.g., when no snapshots are requested) 
  IF (current_avgsnap.GT.N_of_all_avgsnaps) RETURN

! quit if the current moment of time is not the moment when it is necessary to create the snapshot
  IF (T_cntr.NE.avgsnapshot(current_avgsnap)%T_cntr_end) RETURN

  IF (cluster_rank_key.NE.0) THEN
! processes which are not cluster masters do not participate in saving data
     current_avgsnap = current_avgsnap + 1           ! increase the snapshots counter 
     RETURN
  END IF

  IF (Rank_of_process.EQ.0) PRINT '("### ^^^^^^^^^^^^^^^ Averaged Snapshot ",i4," will be created now ... ^^^^^^^^^^^^^^^^ ###")', current_avgsnap

  N_averaged_timesteps   =  avgsnapshot(current_avgsnap)%T_cntr_end - avgsnapshot(current_avgsnap)%T_cntr_begin + 1                  ! for fields and electron moments

!  N_averaged_timesteps_i = (avgsnapshot(current_avgsnap)%T_cntr_end - avgsnapshot(current_avgsnap)%T_cntr_begin) / N_subcycles + 1   ! for ion moments
  N_averaged_timesteps_i = (avgsnapshot(current_avgsnap)%T_cntr_end - avgsnapshot(current_avgsnap)%T_cntr_begin + 1) / N_subcycles   ! for ion moments

! potential

  IF (save_avg_data(1)) THEN
     filename_F = '_NNNN_avg_F_V_2D.bin'
     filename_F(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor =  REAL(F_scale_V) / N_averaged_timesteps
     cs_avg_phi = cs_avg_phi * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_phi, filename_F)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_phi, STAT = ALLOC_ERR)
  END IF

! electric field components

  IF (save_avg_data(2)) THEN
     filename_E = '_NNNN_avg_EX_Vm_2D.bin'
     filename_E(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(E_scale_Vm) / N_averaged_timesteps
     cs_avg_EX = cs_avg_EX * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_EX, filename_E)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_EX, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(3)) THEN
     filename_E = '_NNNN_avg_EY_Vm_2D.bin'
     filename_E(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(E_scale_Vm) / N_averaged_timesteps
     cs_avg_EY = cs_avg_EY * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_EY, filename_E)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_EY, STAT = ALLOC_ERR)
  END IF

! FULL electric current (sum of electron and ion currents) along each coordinate direction
! restored here from electron and ion currents

  IF (save_avg_data(4)) THEN
     filename_Jsum = '_NNNN_avg_JXsum_Am2_2D.bin'
     filename_Jsum(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     ALLOCATE(cs_avg_Jsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps 
     cs_avg_Jsum = cs_avg_JXe * avg_factor
     DO s = 1, N_spec
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_Jsum(i, j) = cs_avg_Jsum(i, j) + cs_avg_JXi(i,j,s) * avg_factor
           END DO
        END DO
     END DO
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_Jsum, filename_Jsum)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_Jsum, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(8) ) DEALLOCATE(cs_avg_JXe, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(24)) DEALLOCATE(cs_avg_JXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(5)) THEN
     filename_Jsum = '_NNNN_avg_JYsum_Am2_2D.bin'
     filename_Jsum(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     ALLOCATE(cs_avg_Jsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps 
     cs_avg_Jsum = cs_avg_JYe * avg_factor
     DO s = 1, N_spec
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_Jsum(i, j) = cs_avg_Jsum(i, j) + cs_avg_JYi(i,j,s) * avg_factor
           END DO
        END DO
     END DO
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_Jsum, filename_Jsum)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_Jsum, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(9) ) DEALLOCATE(cs_avg_JYe, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(25)) DEALLOCATE(cs_avg_JYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(6)) THEN
     filename_Jsum = '_NNNN_avg_JZsum_Am2_2D.bin'
     filename_Jsum(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     ALLOCATE(cs_avg_Jsum(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps 
     cs_avg_Jsum = cs_avg_JZe * avg_factor
     DO s = 1, N_spec
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_Jsum(i, j) = cs_avg_Jsum(i, j) + cs_avg_JZi(i,j,s) * avg_factor
           END DO
        END DO
     END DO
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_Jsum, filename_Jsum)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_Jsum, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(10)) DEALLOCATE(cs_avg_JZe, STAT = ALLOC_ERR)
     IF (.NOT.save_avg_data(26)) DEALLOCATE(cs_avg_JZi, STAT = ALLOC_ERR)
  END IF

! electron number density

  IF (save_avg_data(7)) THEN
     filename_Ne = '_NNNN_avg_Ne_m3_2D.bin'
     filename_Ne(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(N_scale_part_m3) / N_averaged_timesteps
     cs_avg_Ne = cs_avg_Ne * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_Ne, filename_Ne)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_Ne, STAT = ALLOC_ERR)
  END IF

! electron electric current along each coordinate direction

  IF (save_avg_data(8)) THEN
     filename_Je = '_NNNN_avg_JXe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps
     cs_avg_JXe = cs_avg_JXe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JXe, filename_Je)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_JXe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(9)) THEN
     filename_Je = '_NNNN_avg_JYe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps
     cs_avg_JYe = cs_avg_JYe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JYe, filename_Je)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_JYe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(10)) THEN
     filename_Je = '_NNNN_avg_JZe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(current_factor_Am2) / N_averaged_timesteps
     cs_avg_JZe = cs_avg_JZe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JZe, filename_Je)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_JZe, STAT = ALLOC_ERR)
  END IF

! electron average velocity along each coordinate direction

  IF (save_avg_data(11)) THEN
     filename_Ve = '_NNNN_avg_VXe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps
     cs_avg_VXe = cs_avg_VXe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VXe, filename_Ve)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_VXe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(12)) THEN
     filename_Ve = '_NNNN_avg_VYe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps
     cs_avg_VYe = cs_avg_VYe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VYe, filename_Ve)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_VYe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(13)) THEN
     filename_Ve = '_NNNN_avg_VZe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps
     cs_avg_VZe = cs_avg_VZe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VZe, filename_Ve)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_VZe, STAT = ALLOC_ERR)
  END IF

! electron average energy of motion along each coordinate direction

  IF (save_avg_data(14)) THEN
     filename_We = '_NNNN_avg_WXe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(energy_factor_eV) / N_averaged_timesteps
     cs_avg_WXe = cs_avg_WXe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WXe, filename_We)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_WXe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(15)) THEN
     filename_We = '_NNNN_avg_WYe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(energy_factor_eV) / N_averaged_timesteps
     cs_avg_WYe = cs_avg_WYe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WYe, filename_We)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_WYe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(16)) THEN
     filename_We = '_NNNN_avg_WZe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(energy_factor_eV) / N_averaged_timesteps
     cs_avg_WZe = cs_avg_WZe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WZe, filename_We)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_WZe, STAT = ALLOC_ERR)
  END IF

! electron temperature along each coordinate direction

  IF (save_avg_data(17)) THEN
     filename_Te = '_NNNN_avg_TXe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(temperature_factor_eV) / N_averaged_timesteps
     cs_avg_TXe = cs_avg_TXe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TXe, filename_Te)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_TXe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(18)) THEN
     filename_Te = '_NNNN_avg_TYe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(temperature_factor_eV) / N_averaged_timesteps
     cs_avg_TYe = cs_avg_TYe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TYe, filename_Te)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_TYe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(19)) THEN
     filename_Te = '_NNNN_avg_TZe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(temperature_factor_eV) / N_averaged_timesteps
     cs_avg_TZe = cs_avg_TZe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TZe, filename_Te)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_TZe, STAT = ALLOC_ERR)
  END IF

! heat flow along each coordinate direction

  IF (save_avg_data(20)) THEN
     filename_Qe = '_NNNN_avg_QXe_Wm2_2D.bin'
     filename_Qe(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(heat_flow_factor_Wm2) / N_averaged_timesteps
     cs_avg_QXe = cs_avg_QXe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QXe, filename_Qe)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_QXe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(21)) THEN
     filename_Qe = '_NNNN_avg_QYe_Wm2_2D.bin'
     filename_Qe(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(heat_flow_factor_Wm2) / N_averaged_timesteps
     cs_avg_QYe = cs_avg_QYe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QYe, filename_Qe)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_QYe, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(22)) THEN
     filename_Qe = '_NNNN_avg_QZe_Wm2_2D.bin'
     filename_Qe(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
     avg_factor = REAL(heat_flow_factor_Wm2) / N_averaged_timesteps
     cs_avg_QZe = cs_avg_QZe * avg_factor
     CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QZe, filename_Qe)
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
! cleanup
     DEALLOCATE(cs_avg_QZe, STAT = ALLOC_ERR)
  END IF

! ion number density

  IF (save_avg_data(23)) THEN
     avg_factor = REAL(N_scale_part_m3) / N_averaged_timesteps_i
     cs_avg_Ni = cs_avg_Ni * avg_factor
     DO s = 1, N_spec
        filename_Ni = '_NNNN_avg_Ni_s_m3_2D.bin'
        filename_Ni(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ni(14:14) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_Ni(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ni)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr) 
     END DO
! cleanup
     DEALLOCATE(cs_avg_Ni, STAT = ALLOC_ERR)
  END IF

! ion electric current along each coordinate direction

  IF (save_avg_data(24)) THEN
     DO s = 1, N_spec
        filename_Ji = '_NNNN_avg_JXi_s_Am2_2D.bin'
        filename_Ji(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ji(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_JXi(i,j,s) = cs_avg_JXi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ji)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_JXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(25)) THEN
     DO s = 1, N_spec
        filename_Ji = '_NNNN_avg_JYi_s_Am2_2D.bin'
        filename_Ji(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ji(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_JYi(i,j,s) = cs_avg_JYi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ji)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_JYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(26)) THEN
     DO s = 1, N_spec
        filename_Ji = '_NNNN_avg_JZi_s_Am2_2D.bin'
        filename_Ji(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ji(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Qs(s) * current_factor_Am2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_JZi(i,j,s) = cs_avg_JZi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_JZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ji)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_JZi, STAT = ALLOC_ERR)
  END IF

! ion average velocity along each coordinate direction

  IF (save_avg_data(27)) THEN
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps_i
     cs_avg_VXi = cs_avg_VXi * avg_factor
     DO s = 1, N_spec
        filename_Vi = '_NNNN_avg_VXi_s_ms_2D.bin'
        filename_Vi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Vi(15:15) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Vi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_VXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(28)) THEN
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps_i
     cs_avg_VYi = cs_avg_VYi * avg_factor
     DO s = 1, N_spec
        filename_Vi = '_NNNN_avg_VYi_s_ms_2D.bin'
        filename_Vi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Vi(15:15) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Vi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_VYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(29)) THEN
     avg_factor = REAL(V_scale_ms) / N_averaged_timesteps_i
     cs_avg_VZi = cs_avg_VZi * avg_factor
     DO s = 1, N_spec
        filename_Vi = '_NNNN_avg_VZi_s_ms_2D.bin'
        filename_Vi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Vi(15:15) = convert_int_to_txt_string(s, 1)
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_VZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Vi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr) 
     END DO
! cleanup
     DEALLOCATE(cs_avg_VZi, STAT = ALLOC_ERR)
  END IF

! ion average energy of motion along each coordinate direction

  IF (save_avg_data(30)) THEN
     DO s = 1, N_spec
        filename_Wi = '_NNNN_avg_WXi_s_eV_2D.bin'
        filename_Wi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Wi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * energy_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_WXi(i,j,s) = cs_avg_WXi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Wi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_WXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(31)) THEN
     DO s = 1, N_spec
        filename_Wi = '_NNNN_avg_WYi_s_eV_2D.bin'
        filename_Wi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Wi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * energy_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_WYi(i,j,s) = cs_avg_WYi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Wi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_WYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(32)) THEN
     DO s = 1, N_spec
        filename_Wi = '_NNNN_avg_WZi_s_eV_2D.bin'
        filename_Wi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Wi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * energy_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_WZi(i,j,s) = cs_avg_WZi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_WZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Wi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_WZi, STAT = ALLOC_ERR)
  END IF

! ion temperature along each coordinate direction

  IF (save_avg_data(33)) THEN
     DO s = 1, N_spec
        filename_Ti = '_NNNN_avg_TXi_s_eV_2D.bin'
        filename_Ti(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ti(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * temperature_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_TXi(i,j,s) = cs_avg_TXi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ti)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_TXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(34)) THEN
     DO s = 1, N_spec
        filename_Ti = '_NNNN_avg_TYi_s_eV_2D.bin'
        filename_Ti(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ti(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * temperature_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_TYi(i,j,s) = cs_avg_TYi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ti)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_TYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(35)) THEN
     DO s = 1, N_spec
        filename_Ti = '_NNNN_avg_TZi_s_eV_2D.bin'
        filename_Ti(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Ti(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * temperature_factor_eV) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_TZi(i,j,s) = cs_avg_TZi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_TZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Ti)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_TZi, STAT = ALLOC_ERR)
  END IF

! heat flow along each coordinate direction

  IF (save_avg_data(36)) THEN
     DO s = 1, N_spec
        filename_Qi = '_NNNN_avg_QXi_s_Wm2_2D.bin'
        filename_Qi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Qi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * heat_flow_factor_Wm2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_QXi(i,j,s) = cs_avg_QXi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QXi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Qi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_QXi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(37)) THEN
     DO s = 1, N_spec
        filename_Qi = '_NNNN_avg_QYi_s_Wm2_2D.bin'
        filename_Qi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Qi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * heat_flow_factor_Wm2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_QYi(i,j,s) = cs_avg_QYi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QYi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Qi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_QYi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(38)) THEN
     DO s = 1, N_spec
        filename_Qi = '_NNNN_avg_QZi_s_Wm2_2D.bin'
        filename_Qi(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
        filename_Qi(15:15) = convert_int_to_txt_string(s, 1)
        avg_factor = REAL(Ms(s) * heat_flow_factor_Wm2) / N_averaged_timesteps_i
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              cs_avg_QZi(i,j,s) = cs_avg_QZi(i,j,s) * avg_factor
           END DO
        END DO
        CALL SAVE_GLOBAL_2D_ARRAY(cs_avg_QZi(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max, s), filename_Qi)
        CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
     END DO
! cleanup
     DEALLOCATE(cs_avg_QZi, STAT = ALLOC_ERR)
  END IF

  IF (save_avg_data(39)) THEN
     DO n = 1, N_neutral_spec
        DO p = 1, collision_e_neutral(n)%N_of_activated_colproc
           IF (.NOT.collision_e_neutral(n)%colproc_info(p)%save_collfreq_2d) CYCLE       
           filename_encoll = '_NNNN_avg_frequency_e_n_AAAAAA_coll_id_NN.bin'
           filename_encoll(2:5) = convert_int_to_txt_string(current_avgsnap, 4)
           filename_encoll(25:30) = neutral(n)%name
           filename_encoll(40:41) = convert_int_to_txt_string(collision_e_neutral(n)%colproc_info(p)%id_number, 2)
           avg_factor = 1.0 / REAL(delta_t_s * N_averaged_timesteps)
           DO j = c_indx_y_min, c_indx_y_max
              DO i = c_indx_x_min, c_indx_x_max
                 diagnostics_neutral(n)%activated_collision(p)%coll_freq_local(i,j) = diagnostics_neutral(n)%activated_collision(p)%coll_freq_local(i,j) * avg_factor
              END DO
           END DO
           CALL SAVE_GLOBAL_2D_ARRAY(diagnostics_neutral(n)%activated_collision(p)%coll_freq_local, filename_encoll)
           CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)
           DEALLOCATE(diagnostics_neutral(n)%activated_collision(p)%coll_freq_local, STAT = ALLOC_ERR)
        END DO
     END DO
  END IF

  IF (Rank_of_process.EQ.0) PRINT '(/2x,"### ^^^^^^^^^^^^^^^^^^^^ Averaged Snapshot ",i4," completed :) ^^^^^^^^^^^^^^^^^^^ ###")', current_avgsnap

  current_avgsnap = current_avgsnap + 1           ! increase the snapshots counter (cluster masters only, other processes did this already)

END SUBROUTINE CREATE_AVERAGED_SNAPSHOT
