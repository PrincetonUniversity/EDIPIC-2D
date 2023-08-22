!------------------------------
!
SUBROUTINE CALCULATE_ELECTRON_VDF

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries
  USE Snapshots

  use mpi

  IMPLICIT NONE


  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER, ALLOCATABLE :: ibufer_evxdf(:)
  INTEGER, ALLOCATABLE :: ibufer_evydf(:)
  INTEGER, ALLOCATABLE :: ibufer_evzdf(:)
  INTEGER ALLOC_ERR

  INTEGER reclen
  INTEGER bufsize
  INTEGER shift

  INTEGER, ALLOCATABLE :: ibufer_evxvydf(:)
  INTEGER reclen2
  INTEGER bufsize2

  REAL(8), ALLOCATABLE :: xsplit(:)
  REAL(8), ALLOCATABLE :: ysplit(:)

  INTEGER i, j, k
  REAL(8) fvc
  INTEGER skip_count, ibufer(1), ibufer2(1)  ! bufers will be used to collect skip_count from all processes in a cluster
  INTEGER ibox, jbox, n

  INTEGER bin_vx, bin_vy, bin_vz

  INTEGER pos, pos1, pos2

  IF (N_vdfbox_all.EQ.0) RETURN

  IF (save_evdf_snapshot(current_snap).EQ.NOANYVDF) RETURN

  reclen = 2 * N_max_vel_e * N_vbins_e     ! number of points in a single 1d evdf
  bufsize = reclen * N_vdfbox_all          ! total number of values in the array which includes vdfs for one velocity component in all locations
  shift = N_max_vel_e * N_vbins_e + 1      ! if velocity bin index in evdf ranges from -N_max_vel*N_vbins_e to N_max_vel*N_vbins_e-1
                                           ! adding this "shift" modifies the range so that it is from 1 to reclen

  IF ((save_evdf_snapshot(current_snap).EQ.ONLY1D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
     ALLOCATE(ibufer_evxdf(1:bufsize), STAT=ALLOC_ERR)
     ALLOCATE(ibufer_evydf(1:bufsize), STAT=ALLOC_ERR)
     ALLOCATE(ibufer_evzdf(1:bufsize), STAT=ALLOC_ERR)
     ibufer_evxdf = 0
     ibufer_evydf = 0
     ibufer_evzdf = 0
  END IF

! a 2d-velocity distribution, f(vx,vy) is stored in an 1d array in the same manner as 
! a 2d array where vx corresponds to the first (fastest) array index and vy - to the second array index
! that is 
! f(vxmin:vxmax,vymin), f(vxmin:vxmax,vymin+1),... f(vxmin:vxmax,vymax-1), f(vxmin:vxmax,vymax)

  reclen2 = reclen * reclen                ! number of points in a single 2d evdf
  bufsize2 = reclen2 * N_vdfbox_all        ! total number of values in the array which includes vdfs for two velocity components in all locations
  
  IF ((save_evdf_snapshot(current_snap).EQ.ONLY2D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
     ALLOCATE(ibufer_evxvydf(1:bufsize2), STAT=ALLOC_ERR)
     ibufer_evxvydf = 0
  END IF

! arrays xsplit and ysplit are recalculated rather than stored
! this must be cheaper since it does not require including these arrays in the list
! of parameters a master of a cluster must share with its particle calculators in DISTRIBUTE_CLUSTER_PARAMETERS

  IF (N_vdfbox_x.GT.1) THEN 
     ALLOCATE(xsplit(1:(N_vdfbox_x-1)), STAT = ALLOC_ERR)
     DO i = 1, N_vdfbox_x-1
        xsplit(i) = c_X_area_min + DBLE(i) * (c_X_area_max - c_X_area_min) / N_vdfbox_x !????????
     END DO
  END IF

  IF (N_vdfbox_y.GT.1) THEN 
     ALLOCATE(ysplit(1:(N_vdfbox_y-1)), STAT = ALLOC_ERR)
     DO j = 1, N_vdfbox_y-1
        ysplit(j) = c_Y_area_min + DBLE(j) * (c_Y_area_max - c_Y_area_min) / N_vdfbox_y !????????
     END DO
  END IF

  fvc = DBLE(N_max_vel * N_vbins_e)   ! "fvc" stands for "factor to convert velocity"
                                      ! it is used to determine a velocity bin for a particle 
                                      ! note that we use here N_max_vel, not N_max_vel_e
  skip_count = 0

  DO k = 1, N_electrons
     
! identify a spatial box where the particle is situated
! by default assume that the particle is in the right-top corner box 
! if a cluster has just one box, this will set the box instantly 

! X-direction
     ibox = N_vdfbox_x                     
     DO i = 1, N_vdfbox_x-1
        IF (electron(k)%X.LT.xsplit(i)) THEN
           ibox = i
           EXIT
        END IF
     END DO

! Y-direction
     jbox = N_vdfbox_y
     DO j = 1, N_vdfbox_y-1
        IF (electron(k)%Y.LT.ysplit(j)) THEN
           jbox = j
           EXIT
        END IF
     END DO

! array index (in 2d arrays ev[xyz]df) corresponding to the spatial box
     n = ibox + (jbox-1) * N_vdfbox_x

! identify velocity bins and 
! skip particles which are too fast to be registered correctly by the given velocity bin set

! X-direction
     bin_vx = INT(electron(k)%VX * fvc)
     IF (electron(k)%VX.LT.0.0_8) bin_vx = bin_vx - 1

     IF ((bin_vx.LT.indx_v_min_e).OR.(bin_vx.GT.indx_v_max_e)) THEN
        skip_count = skip_count + 1
        CYCLE
     END IF

! Y-direction
     bin_vy = INT(electron(k)%VY * fvc)
     IF (electron(k)%VY.LT.0.0_8) bin_vy = bin_vy - 1

     IF ((bin_vy.LT.indx_v_min_e).OR.(bin_vy.GT.indx_v_max_e)) THEN
        skip_count = skip_count + 1
        CYCLE
     END IF 

! Z-direction
     bin_vz = INT(electron(k)%VZ * fvc)
     IF (electron(k)%VZ.LT.0.0_8) bin_vz = bin_vz - 1

     IF ((bin_vz.LT.indx_v_min_e).OR.(bin_vz.GT.indx_v_max_e))  THEN
        skip_count = skip_count + 1
        CYCLE
     END IF

! include particle in 1d-vdf arrays

     IF ((save_evdf_snapshot(current_snap).EQ.ONLY1D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
! common displacement of index accounting for location (spatial box) number
        pos1 = (n-1) * reclen + shift
! X-direction
        pos = pos1 + bin_vx
        ibufer_evxdf(pos) = ibufer_evxdf(pos) + 1
! Y-direction
        pos = pos1 + bin_vy
        ibufer_evydf(pos) = ibufer_evydf(pos) + 1
! Z-direction
        pos = pos1 + bin_vz
        ibufer_evzdf(pos) = ibufer_evzdf(pos) + 1
     END IF

! include particle in 2d-vdf arrays (vx-vy only now)

     IF ((save_evdf_snapshot(current_snap).EQ.ONLY2D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
        pos2 = (n-1) * reclen2 + (shift + bin_vy - 1) * reclen + (shift + bin_vx)
        ibufer_evxvydf(pos2) = ibufer_evxvydf(pos2) + 1
     END IF

  END DO

! collect/report (if any) numbers of skipped particles from all processes in a cluster
  ibufer(1) = skip_count
  ibufer2 = 0
  CALL MPI_REDUCE(ibufer, ibufer2, 1, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)

  IF ((cluster_rank_key.EQ.0).AND.(ibufer2(1).GT.0)) PRINT '("### Warning :: Process ",i5," :: skipped ",i9," electron particles when calculating EVDFs ###")', Rank_of_process, ibufer2(1)

  IF ((save_evdf_snapshot(current_snap).EQ.ONLY1D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
     CALL MPI_REDUCE(ibufer_evxdf, evxdf, bufsize, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)
     CALL MPI_REDUCE(ibufer_evydf, evydf, bufsize, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)
     CALL MPI_REDUCE(ibufer_evzdf, evzdf, bufsize, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)
     DEALLOCATE(ibufer_evxdf, STAT=ALLOC_ERR)
     DEALLOCATE(ibufer_evydf, STAT=ALLOC_ERR)
     DEALLOCATE(ibufer_evzdf, STAT=ALLOC_ERR)
  END IF

  IF ((save_evdf_snapshot(current_snap).EQ.ONLY2D).OR.(save_evdf_snapshot(current_snap).EQ.BOTH1DAND2D)) THEN
     CALL MPI_REDUCE(ibufer_evxvydf, evxvydf, bufsize2, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)
     DEALLOCATE(ibufer_evxvydf, STAT=ALLOC_ERR)
  END IF

  IF (N_vdfbox_all.GT.1) THEN
     DEALLOCATE(xsplit, STAT = ALLOC_ERR)
     DEALLOCATE(ysplit, STAT = ALLOC_ERR)
  END IF

END SUBROUTINE CALCULATE_ELECTRON_VDF


!------------------------------
!
SUBROUTINE CALCULATE_ION_VDF

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles
  USE ClusterAndItsBoundaries
  USE Snapshots

  use mpi

  IMPLICIT NONE


  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER, ALLOCATABLE :: ibufer_isvxdf(:)
  INTEGER, ALLOCATABLE :: ibufer_isvydf(:)
  INTEGER, ALLOCATABLE :: ibufer_isvzdf(:)
  INTEGER ALLOC_ERR

  INTEGER reclen
  INTEGER bufsize
  INTEGER shift

  REAL(8), ALLOCATABLE :: xsplit(:)
  REAL(8), ALLOCATABLE :: ysplit(:)

  INTEGER i, j, k, s
  REAL(8) fvc
  INTEGER skip_count(1:N_spec), ibufer2(1:N_spec)  ! bufers will be used to collect skip_count from all processes in a cluster
  INTEGER ibox, jbox, n

  INTEGER bin_vx, bin_vy, bin_vz

  INTEGER pos, pos1

  IF (N_vdfbox_all.EQ.0) RETURN

  IF (save_evdf_snapshot(current_snap).EQ.NOANYVDF) RETURN
  IF (save_evdf_snapshot(current_snap).EQ.ONLY2D) RETURN

  reclen = 2 * N_max_vel_i * N_vbins_i              ! number of points in a single 1d evdf for one ion species
  bufsize = reclen * N_vdfbox_all * N_spec          ! total number of values in the array which includes vdfs for one velocity component in all locations
  shift = N_max_vel_i * N_vbins_i + 1               ! if velocity bin index in evdf ranges from -N_max_vel*N_vbins_i to N_max_vel*N_vbins_i-1
                                                    ! adding this "shift" modifies the range so that it is from 1 to reclen

  ALLOCATE(ibufer_isvxdf(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(ibufer_isvydf(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(ibufer_isvzdf(1:bufsize), STAT=ALLOC_ERR)

  ibufer_isvxdf = 0
  ibufer_isvydf = 0
  ibufer_isvzdf = 0

! arrays xsplit and ysplit are recalculated rather than stored
! this must be cheaper since it does not require including these arrays in the list
! of parameters a master of a cluster must share with its particle calculators in DISTRIBUTE_CLUSTER_PARAMETERS

  IF (N_vdfbox_x.GT.1) THEN 
     ALLOCATE(xsplit(1:(N_vdfbox_x-1)), STAT = ALLOC_ERR)
     DO i = 1, N_vdfbox_x-1
        xsplit(i) = c_X_area_min + DBLE(i) * (c_X_area_max - c_X_area_min) / N_vdfbox_x
     END DO
  END IF

  IF (N_vdfbox_y.GT.1) THEN 
     ALLOCATE(ysplit(1:(N_vdfbox_y-1)), STAT = ALLOC_ERR)
     DO j = 1, N_vdfbox_y-1
        ysplit(j) = c_Y_area_min + DBLE(j) * (c_Y_area_max - c_Y_area_min) / N_vdfbox_y
     END DO
  END IF

  skip_count = 0

  DO s = 1, N_spec

     fvc = DBLE(N_max_vel * N_vbins_i) * SQRT(Ms(s))   ! "fvc" stands for "factor to convert velocity"
                                                       ! it is used to determine a velocity bin for a particle 
                                                       ! note that we use here N_max_vel, not N_max_vel_i
                                                       ! and that the factor depends on ion species
     DO k = 1, N_ions(s)
     
! identify a spatial box where the particle is situated
! by default assume that the particle is in the right-top corner box 
! if a cluster has just one box, this will set the box instantly 

! X-direction
        ibox = N_vdfbox_x                     
        DO i = 1, N_vdfbox_x-1
           IF (ion(s)%part(k)%X.LT.xsplit(i)) THEN
              ibox = i
              EXIT
           END IF
        END DO
        
! Y-direction
        jbox = N_vdfbox_y
        DO j = 1, N_vdfbox_y-1
           IF (ion(s)%part(k)%Y.LT.ysplit(j)) THEN
              jbox = j
              EXIT
           END IF
        END DO
        
! array index (in 2d arrays isv[xyz]df) corresponding to the spatial box
        n = ibox + (jbox-1) * N_vdfbox_x

! identify velocity bins and 
! skip particles which are too faster to be registered correctly by the given velocity bin set

! X-direction
        bin_vx = INT(ion(s)%part(k)%VX * fvc)
        IF (ion(s)%part(k)%VX.LT.0.0_8) bin_vx = bin_vx - 1

        IF ((bin_vx.LT.indx_v_min_i).OR.(bin_vx.GT.indx_v_max_i)) THEN
           skip_count(s) = skip_count(s) + 1
           CYCLE
        END IF

! Y-direction
        bin_vy = INT(ion(s)%part(k)%VY * fvc)
        IF (ion(s)%part(k)%VY.LT.0.0_8) bin_vy = bin_vy - 1

        IF ((bin_vy.LT.indx_v_min_i).OR.(bin_vy.GT.indx_v_max_i)) THEN
           skip_count(s) = skip_count(s) + 1
           CYCLE
        END IF

! Z-direction
        bin_vz = INT(ion(s)%part(k)%VZ * fvc)
        IF (ion(s)%part(k)%VZ.LT.0.0_8) bin_vz = bin_vz - 1

        IF ((bin_vz.LT.indx_v_min_i).OR.(bin_vz.GT.indx_v_max_i))  THEN
           skip_count(s) = skip_count(s) + 1
           CYCLE
        END IF

! common displacement of index accounting for location (spatial box) number and ion species number
        pos1 = ((s-1) * N_vdfbox_all + (n-1)) * reclen + shift

! include particle in the vdf arrays

! X-direction
        pos = pos1 + bin_vx
        ibufer_isvxdf(pos) = ibufer_isvxdf(pos) + 1

! Y-direction
        pos = pos1 + bin_vy
        ibufer_isvydf(pos) = ibufer_isvydf(pos) + 1

! Z-direction
        pos = pos1 + bin_vz
        ibufer_isvzdf(pos) = ibufer_isvzdf(pos) + 1

     END DO
  END DO

! collect moments from all processes in a cluster

  ibufer2 = 0
  CALL MPI_REDUCE(skip_count, ibufer2, N_spec, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)

  DO s = 1, N_spec
     IF ((cluster_rank_key.EQ.0).AND.(ibufer2(s).GT.0)) PRINT '("### Warning :: Process ",i5," :: skipped ",i9," particles of ion species ",i2," when calculating 1D EVDFs ###")', Rank_of_process, ibufer2(s), s
  END DO

  CALL MPI_REDUCE(ibufer_isvxdf, isvxdf, bufsize, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_REDUCE(ibufer_isvydf, isvydf, bufsize, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_REDUCE(ibufer_isvzdf, isvzdf, bufsize, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)

  DEALLOCATE(ibufer_isvxdf, STAT=ALLOC_ERR)
  DEALLOCATE(ibufer_isvydf, STAT=ALLOC_ERR)
  DEALLOCATE(ibufer_isvzdf, STAT=ALLOC_ERR)

  IF (N_vdfbox_all.GT.1) THEN
     DEALLOCATE(xsplit, STAT = ALLOC_ERR)
     DEALLOCATE(ysplit, STAT = ALLOC_ERR)
  END IF

END SUBROUTINE CALCULATE_ION_VDF
