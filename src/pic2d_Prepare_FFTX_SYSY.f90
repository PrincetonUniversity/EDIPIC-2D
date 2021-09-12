
!-------------------------------------------
!
SUBROUTINE PREPARE_FFT_X

  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues
  USE ParallelFFTX
  USE ParallelOperationValues, ONLY : Rank_of_process, cluster_rank_key

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  TYPE pair_of_comm_proc
     INTEGER proc1            ! rank of one process
     INTEGER proc2            ! rank of the other process
     LOGICAL used             ! shows whether the pair was already included in the sequence or not
  END type pair_of_comm_proc

  TYPE prepared_index_limits
     INTEGER c_indx_x_min
     INTEGER c_indx_x_max
     INTEGER fftx_band_jmin
     INTEGER fftx_band_jmax
     INTEGER invfftx_band_jmin
     INTEGER invfftx_band_jmax
  END type prepared_index_limits

  INTEGER, ALLOCATABLE :: rank_in_col(:)    ! ranks of all master processes in the same X-line

  TYPE(pair_of_comm_proc), ALLOCATABLE :: comm_pair(:,:)

  INTEGER, ALLOCATABLE :: breakj(:)

  TYPE(prepared_index_limits), ALLOCATABLE :: prepared_for_column(:)

  INTEGER ALLOC_ERR

  INTEGER col
  INTEGER i, k, m, n

  INTEGER jstart, jend
  INTEGER band_width
  INTEGER extraj
  INTEGER pos
  INTEGER otp_block_row
  INTEGER otp_block_column
  INTEGER otp_c_column

  INTEGER max_comm_steps
  INTEGER, ALLOCATABLE :: prot_comm_step_proc(:,:)
  LOGICAL empty_step, proc_not_used

  INTEGER kx, ky

  INTEGER strip_width

  IF (cluster_rank_key.NE.0) RETURN

! since we are here, this is a master process

!=== calculate communication sequence for transfer of density to (and its fft from) the bands stretching along X ===
     
  ALLOCATE(rank_in_col(1:N_clusters_x), STAT=ALLOC_ERR)
  ALLOCATE(comm_pair(2:N_clusters_x,1:N_clusters_x-1), STAT=ALLOC_ERR)
  ALLOCATE(breakj(0:N_clusters_x), STAT=ALLOC_ERR)
  ALLOCATE(prepared_for_column(1:N_clusters_x), STAT=ALLOC_ERR)

! prepare ranks of all master processes in the X-line with the same c_row as the current process
  DO col = 1, N_clusters_x
     rank_in_col(col) = N_blocks_x * (c_row - 1) * cluster_N_blocks_y + (col - 1) * cluster_N_blocks_x
  END DO

! prepare temporary matrix with all possible pairs of ranks of processes in the line communicating with each other
  DO n = 1, N_clusters_x-1
     DO m = n+1, N_clusters_x
        comm_pair(m,n)%proc1=rank_in_col(n)
        comm_pair(m,n)%proc2=rank_in_col(m)
        comm_pair(m,n)%used=.FALSE.
     END DO
  END DO

! prepare j-limits for j-band for each master process in the row for direct fft transformation
  jstart = c_indx_y_min+1
  jend   = c_indx_y_max-1

! include wall nodes for clusters adjacent to walls (both metal and dielectric)
  IF (jstart.EQ.1) jstart = 0
  IF (jend.EQ.(global_maximal_j-1)) jend = global_maximal_j

! approximate width 
  band_width = (jend - jstart + 1) / N_clusters_x
! find number of extra nodes that will be if we use the band_width everywhere 
  extraj = (jend - jstart + 1) - N_clusters_x * band_width
! in the first extraj bands, increase width by 1
  breakj(0) = jstart
  DO col = 1, extraj
     breakj(col) = breakj(col-1) + band_width + 1
  END DO
! in the other bands, use standard width
  DO col = extraj+1, N_clusters_x-1
     breakj(col) = breakj(col-1) + band_width
  END DO
  breakj(N_clusters_x) = jend

  DO col = 1, N_clusters_x-1
     prepared_for_column(col)%fftx_band_jmin = breakj(col-1)
     prepared_for_column(col)%fftx_band_jmax = breakj(col)-1
  END DO
  prepared_for_column(N_clusters_x)%fftx_band_jmin = breakj(N_clusters_x-1)
  prepared_for_column(N_clusters_x)%fftx_band_jmax = breakj(N_clusters_x)

! repeat the above for the inverse fft transformation when overlapping is permitted
  jstart = c_indx_y_min
  jend   = c_indx_y_max
! include the y-end-points
!  IF (c_row.EQ.1) jstart = c_indx_y_min+1
!  IF (c_row.EQ.N_clusters_y) jend = c_indx_y_max-1

! approximate width 
  band_width = (jend - jstart + 1) / N_clusters_x
! find number of extra nodes that will be if we use the band_width everywhere 
  extraj = (jend - jstart + 1) - N_clusters_x * band_width
! in the first extraj bands, increase width by 1
  breakj(0) = jstart
  DO col = 1, extraj
     breakj(col) = breakj(col-1) + band_width + 1
  END DO
! in the other bands, use standard width
  DO col = extraj+1, N_clusters_x-1
     breakj(col) = breakj(col-1) + band_width
  END DO
  breakj(N_clusters_x) = jend

  DO col = 1, N_clusters_x-1
     prepared_for_column(col)%invfftx_band_jmin = breakj(col-1)
     prepared_for_column(col)%invfftx_band_jmax = breakj(col)-1
  END DO
  prepared_for_column(N_clusters_x)%invfftx_band_jmin = breakj(N_clusters_x-1)
  prepared_for_column(N_clusters_x)%invfftx_band_jmax = breakj(N_clusters_x)
  
  DO col = 1, N_clusters_x
     otp_block_row    = 1 + rank_in_col(col) / N_blocks_x
     otp_block_column = 1 + rank_in_col(col) - N_blocks_x * (otp_block_row - 1)
     otp_c_column = 1 + (otp_block_column-1) / cluster_N_blocks_x
     prepared_for_column(col)%c_indx_x_min = (otp_block_column - 1) * N_grid_block_x
     prepared_for_column(col)%c_indx_x_max = (otp_block_column - 1 + cluster_N_blocks_x) * N_grid_block_x + 1 
  END DO

! use the values calculated above to set band index limits in this process
  fftx_band_jmin = prepared_for_column(c_column)%fftx_band_jmin
  fftx_band_jmax = prepared_for_column(c_column)%fftx_band_jmax

  invfftx_band_jmin = prepared_for_column(c_column)%invfftx_band_jmin
  invfftx_band_jmax = prepared_for_column(c_column)%invfftx_band_jmax

! prepare general communication protocol
  max_comm_steps = 2 * N_clusters_x  ! normally everything should be done after N_clusters_x-1 cycles
  ALLOCATE(prot_comm_step_proc(1:max_comm_steps, 1:N_clusters_x), STAT=ALLOC_ERR)
  prot_comm_step_proc = -1           ! this distinguishes legal ranks from dummy placeholders, used for steps where the set of pairs is incomplete  
  N_of_comm_steps_X = 0
  DO k = 1, max_comm_steps
! cycle over communication steps
     pos=0
     empty_step = .TRUE.
     DO n = 1, N_clusters_x-1
! cycle over columns of the matrix with rank combinations
        DO m = n+1, N_clusters_x
! cycle over rows of the matrix with rank combinations
           IF (.NOT.comm_pair(m,n)%used) THEN
! consider only combinations not used at previous communication steps
              proc_not_used = .TRUE.
              DO i = 1, pos
! check whether any of processes in the pair of processes considered are already included in this communication step
                 IF ((prot_comm_step_proc(k,i).EQ.comm_pair(m,n)%proc1).OR.(prot_comm_step_proc(k,i).EQ.comm_pair(m,n)%proc2)) THEN
                    proc_not_used = .FALSE.
                    EXIT
                 END IF
              END DO
              IF (proc_not_used) THEN
! include the pair of processes, label it as used, and show that at least one pair participates in this communication step
                 prot_comm_step_proc(k,pos+1) = comm_pair(m,n)%proc1
                 prot_comm_step_proc(k,pos+2) = comm_pair(m,n)%proc2
                 pos = pos+2
                 comm_pair(m,n)%used = .TRUE.
                 empty_step = .FALSE.
              END IF
           END IF
        END DO
     END DO
! end of cycle over communication steps
     IF (empty_step) THEN
        N_of_comm_steps_X = k-1
        EXIT
     ELSE
        IF (k.EQ.max_comm_steps) THEN
! error, too many communication steps, something is wrong
           PRINT '("Process ",i4," :: Error in PREPARE_FFT_X, too many communication steps for transfer of density X-bands : k = ", i4)', Rank_of_process, k
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        END IF
     END IF
  END DO
! extract protocol for this specific process
  ALLOCATE(comm_FFTX_step(1:N_of_comm_steps_X), STAT=ALLOC_ERR)
! note, expect that N_clusters_x is even
  DO k = 1, N_of_comm_steps_X
     comm_FFTX_step(k)%rank = -1
     DO pos = 1, N_clusters_x-1, 2

        IF (prot_comm_step_proc(k,pos).EQ.Rank_of_process) THEN
           comm_FFTX_step(k)%rank = prot_comm_step_proc(k,pos+1)
        ELSE IF (prot_comm_step_proc(k,pos+1).EQ.Rank_of_process) THEN
           comm_FFTX_step(k)%rank = prot_comm_step_proc(k,pos)
        END IF

        IF (comm_FFTX_step(k)%rank.GE.0) THEN

!              otp_block_row    = 1 + comm_FFTX_step(k)%rank / N_blocks_x
!              otp_block_column = 1 + comm_FFTX_step(k)%rank - N_blocks_x * (otp_block_row - 1)
!              otp_c_column = 1 + (otp_block_column-1) / cluster_N_blocks_x

           otp_c_column = 1 + (comm_FFTX_step(k)%rank - N_blocks_x * INT(comm_FFTX_step(k)%rank / N_blocks_x)) / cluster_N_blocks_x

           comm_FFTX_step(k)%c_indx_x_min = prepared_for_column(otp_c_column)%c_indx_x_min
           comm_FFTX_step(k)%c_indx_x_max = prepared_for_column(otp_c_column)%c_indx_x_max

           comm_FFTX_step(k)%fftx_band_jmin = prepared_for_column(otp_c_column)%fftx_band_jmin
           comm_FFTX_step(k)%fftx_band_jmax = prepared_for_column(otp_c_column)%fftx_band_jmax

           comm_FFTX_step(k)%invfftx_band_jmin = prepared_for_column(otp_c_column)%invfftx_band_jmin
           comm_FFTX_step(k)%invfftx_band_jmax = prepared_for_column(otp_c_column)%invfftx_band_jmax

        END IF

     END DO
  END DO

  IF (cluster_N_blocks.EQ.1) THEN
! a special case when no additional field calculators are used in the cluster, only the master process
     fftx_strip_jmin = fftx_band_jmin
     fftx_strip_jmax = fftx_band_jmax

     invfftx_strip_jmin = invfftx_band_jmin
     invfftx_strip_jmax = invfftx_band_jmax

  ELSE
! at least one additional field calculator accompanies the master process in the cluster

! note, ranks of field_calculator members are already set in SET_CLUSTER_STRUCTURE

     DEALLOCATE(breakj, STAT=ALLOC_ERR)
     ALLOCATE(breakj(0:cluster_N_blocks-1), STAT=ALLOC_ERR)

! prepare j-limits for j-strip for each field calculator in the cluster, including the matrix process, for direct fft transformation
     jstart = fftx_band_jmin
     jend   = fftx_band_jmax

! approximate width 
     strip_width = (jend - jstart + 1) / cluster_N_blocks
! find number of extra nodes that will be if we use the strip_width everywhere 
     extraj = (jend - jstart + 1) - cluster_N_blocks * strip_width
! in the first extraj bands, increase width by 1
     breakj(0) = jstart
     DO k = 1, extraj
        breakj(k) = breakj(k-1) + strip_width + 1
     END DO
! in the other bands, use standard width
     DO k = extraj+1, cluster_N_blocks-1
        breakj(k) = breakj(k-1) + strip_width
     END DO

     fftx_strip_jmin = jstart
     fftx_strip_jmax = breakj(1)-1
     DO k = 2, cluster_N_blocks-1
        field_calculator(k)%fftx_strip_jmin = breakj(k-1)
        field_calculator(k)%fftx_strip_jmax = breakj(k)-1
     END DO
     field_calculator(cluster_N_blocks)%fftx_strip_jmin = breakj(cluster_N_blocks-1)
     field_calculator(cluster_N_blocks)%fftx_strip_jmax = jend

! prepare j-limits for j-strip for each field calculator in the cluster, including the matrix process, for inverse fft transformation
     jstart = invfftx_band_jmin
     jend   = invfftx_band_jmax

! approximate width 
     strip_width = (jend - jstart + 1) / cluster_N_blocks
! find number of extra nodes that will be if we use the strip_width everywhere 
     extraj = (jend - jstart + 1) - cluster_N_blocks * strip_width
! in the first extraj bands, increase width by 1
     breakj(0) = jstart
     DO k = 1, extraj
        breakj(k) = breakj(k-1) + strip_width + 1
     END DO
! in the other bands, use standard width
     DO k = extraj+1, cluster_N_blocks-1
        breakj(k) = breakj(k-1) + strip_width
     END DO

     invfftx_strip_jmin = jstart
     invfftx_strip_jmax = breakj(1)-1
     DO k = 2, cluster_N_blocks-1
        field_calculator(k)%invfftx_strip_jmin = breakj(k-1)
        field_calculator(k)%invfftx_strip_jmax = breakj(k)-1
     END DO
     field_calculator(cluster_N_blocks)%invfftx_strip_jmin = breakj(cluster_N_blocks-1)
     field_calculator(cluster_N_blocks)%invfftx_strip_jmax = jend

! note, index limits will be sent to associated field calculator after this subroutine, in PREPARE_SYS_Y

  END IF

  DEALLOCATE(rank_in_col, STAT=ALLOC_ERR)
  DEALLOCATE(comm_pair, STAT=ALLOC_ERR)
  DEALLOCATE(breakj, STAT=ALLOC_ERR)
  DEALLOCATE(prepared_for_column, STAT=ALLOC_ERR)
  DEALLOCATE(prot_comm_step_proc, STAT=ALLOC_ERR)
  
END SUBROUTINE PREPARE_FFT_X

!------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>
!
SUBROUTINE PREPARE_SYS_Y

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues
!  USE BlockAndItsBoundaries, ONLY : indx_x_min, indx_x_max, indx_y_min, indx_y_max
  USE ParallelFFTX
  USE SetupValues, ONLY : ht_grid_requested, grid_j

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER ibufer(6)

  TYPE pair_of_comm_proc
     INTEGER proc1            ! rank of one process
     INTEGER proc2            ! rank of the other process
     LOGICAL used             ! shows whether the pair was already included in the sequence or not
  END type pair_of_comm_proc

  TYPE prepared_index_limits
     INTEGER c_indx_y_min
     INTEGER c_indx_y_max
     INTEGER sysy_band_nmin
     INTEGER sysy_band_nmax
  END type prepared_index_limits

  INTEGER, ALLOCATABLE :: rank_in_row(:)    ! ranks of all master processes in the same X-line

  TYPE(pair_of_comm_proc), ALLOCATABLE :: comm_pair(:,:)

  INTEGER, ALLOCATABLE :: breakn(:)

  TYPE(prepared_index_limits), ALLOCATABLE :: prepared_for_row(:)

  INTEGER ALLOC_ERR

  INTEGER row
  INTEGER i, k, m, n, j

  INTEGER nstart, nend
  INTEGER band_width
  INTEGER extran
  INTEGER pos
  INTEGER otp_block_row
  INTEGER otp_block_column
!  INTEGER otp_c_column
  INTEGER otp_c_row

  INTEGER max_comm_steps
  INTEGER, ALLOCATABLE :: prot_comm_step_proc(:,:)
  LOGICAL empty_step, proc_not_used

  INTEGER kx, ky

  INTEGER strip_width

  INTEGER indx_x_begin, indx_x_end, xsize, bufsize
  INTEGER indx_y_begin, indx_y_end, ysize

  LOGICAL wall_below_is_metal, wall_above_is_metal
  REAL(8) main_diag_term, main_diag_term_prim_upper, main_diag_term_prim_lower

  IF (cluster_rank_key.EQ.0) THEN

!=== calculate communication sequence for transfer of fft components of density to/from the bands stretching along Y ===
     
     ALLOCATE(rank_in_row(1:N_clusters_y), STAT=ALLOC_ERR)
     ALLOCATE(comm_pair(2:N_clusters_y,1:N_clusters_y-1), STAT=ALLOC_ERR)
     ALLOCATE(breakn(0:N_clusters_y), STAT=ALLOC_ERR)
     ALLOCATE(prepared_for_row(1:N_clusters_y), STAT=ALLOC_ERR)

! prepare ranks of all master processes in the Y-line with the same c_column as the current process
     DO row = 1, N_clusters_y
        rank_in_row(row) = N_blocks_x * (row - 1) * cluster_N_blocks_y + (c_column - 1) * cluster_N_blocks_x
     END DO

! prepare temporary matrix with all possible pairs of ranks of processes in the column communicating with each other
     DO n = 1, N_clusters_y-1
        DO m = n+1, N_clusters_y
           comm_pair(m,n)%proc1=rank_in_row(n)
           comm_pair(m,n)%proc2=rank_in_row(m)
           comm_pair(m,n)%used=.FALSE.
        END DO
     END DO

!     c_indx_x_min = (block_column - 1) * N_grid_block_x
!     c_indx_x_max = (block_column - 1 + cluster_N_blocks_x) * N_grid_block_x + 1 

! prepare n-limits for n-band for each master process in the column for solving the equation system
     nstart = c_indx_x_min / 2
     nend   = c_indx_x_max / 2 - 1
     IF (c_column.EQ.N_clusters_x) nend=(global_maximal_i-1)/2         !#########################?????????????????????????

! approximate width 
     band_width = (nend - nstart + 1) / N_clusters_y
! find number of extra nodes that will be if we use the band_width everywhere 
     extran = (nend - nstart + 1) - N_clusters_y * band_width
! in the first extran bands, increase width by 1 (implies that extran < N_clusters_y) 
     breakn(0) = nstart
     DO row = 1, extran
        breakn(row) = breakn(row-1) + band_width + 1
     END DO
! in the other bands, use standard width
     DO row = extran+1, N_clusters_y-1
        breakn(row) = breakn(row-1) + band_width
     END DO
     breakn(N_clusters_y) = nend

!print '(10(2x,i5),4x,5(2x,i4))', Rank_of_process, c_row, c_column, c_indx_x_min, c_indx_x_max, nstart, nend, band_width, extran, N_clusters_y, breakn(0:N_clusters_y)
!call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     DO row = 1, N_clusters_y-1
        prepared_for_row(row)%sysy_band_nmin = breakn(row-1)
        prepared_for_row(row)%sysy_band_nmax = breakn(row)-1
     END DO
     prepared_for_row(N_clusters_y)%sysy_band_nmin = breakn(N_clusters_y-1)
     prepared_for_row(N_clusters_y)%sysy_band_nmax = breakn(N_clusters_y)

     DO row = 1, N_clusters_y
        otp_block_row    = 1 + rank_in_row(row) / N_blocks_x
        prepared_for_row(row)%c_indx_y_min = (otp_block_row - 1) * N_grid_block_y                      !(otp_block_column - 1) * N_grid_block_x
        prepared_for_row(row)%c_indx_y_max = (otp_block_row - 1 + cluster_N_blocks_y) * N_grid_block_y + 1 !(otp_block_column - 1 + cluster_N_blocks_x) * N_grid_block_x + 1 
     END DO

!print '(10(2x,i5),4x,5(2x,i4))', Rank_of_process, c_row, c_column, c_indx_x_min, c_indx_x_max, nstart, nend, band_width, extran, N_clusters_y, breakn(0:N_clusters_y)
!call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! use the values calculated above to set band index limits in this process
     sysy_band_nmin = prepared_for_row(c_row)%sysy_band_nmin
     sysy_band_nmax = prepared_for_row(c_row)%sysy_band_nmax

! prepare general communication protocol
     max_comm_steps = 2 * N_clusters_y  ! normally everything should be done after N_clusters_y-1 cycles
     ALLOCATE(prot_comm_step_proc(1:max_comm_steps, 1:N_clusters_y), STAT=ALLOC_ERR)
     prot_comm_step_proc = -1           ! this distinguishes legal ranks from dummy placeholders, used for steps where the set of pairs is incomplete  
     N_of_comm_steps_Y = 0
     DO k = 1, max_comm_steps
! cycle over communication steps
        pos=0
        empty_step = .TRUE.
        DO n = 1, N_clusters_y-1
! cycle over columns of the matrix with rank combinations
           DO m = n+1, N_clusters_y
! cycle over rows of the matrix with rank combinations
              IF (.NOT.comm_pair(m,n)%used) THEN
! consider only combinations not used at previous communication steps
                 proc_not_used = .TRUE.
                 DO i = 1, pos
! check whether any of processes in the pair of processes considered are already included in this communication step
                    IF ((prot_comm_step_proc(k,i).EQ.comm_pair(m,n)%proc1).OR.(prot_comm_step_proc(k,i).EQ.comm_pair(m,n)%proc2)) THEN
                       proc_not_used = .FALSE.
                       EXIT
                    END IF
                 END DO
                 IF (proc_not_used) THEN
! include the pair of processes, label it as used, and show that at least one pair participates in this communication step
                    prot_comm_step_proc(k,pos+1) = comm_pair(m,n)%proc1
                    prot_comm_step_proc(k,pos+2) = comm_pair(m,n)%proc2
                    pos = pos+2
                    comm_pair(m,n)%used = .TRUE.
                    empty_step = .FALSE.
                 END IF
              END IF
           END DO
        END DO
! end of cycle over communication steps
        IF (empty_step) THEN
           N_of_comm_steps_Y = k-1
           EXIT
        ELSE
           IF (k.EQ.max_comm_steps) THEN
! error, too many communication steps, something is wrong
              PRINT '("Process ",i4," :: Error in PREPARE_SYS_Y, too many communication steps for transfer of density X-bands : k = ", i4)', Rank_of_process, k
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF
        END IF
     END DO
! extract protocol for this specific process
     ALLOCATE(comm_SYSY_step(1:N_of_comm_steps_Y), STAT=ALLOC_ERR)
! note, expect that N_clusters_y is even
     DO k = 1, N_of_comm_steps_Y
        comm_SYSY_step(k)%rank = -1
        DO pos = 1, N_clusters_y-1, 2

           IF (prot_comm_step_proc(k,pos).EQ.Rank_of_process) THEN
              comm_SYSY_step(k)%rank = prot_comm_step_proc(k,pos+1)
           ELSE IF (prot_comm_step_proc(k,pos+1).EQ.Rank_of_process) THEN
              comm_SYSY_step(k)%rank = prot_comm_step_proc(k,pos)
           END IF

           IF (comm_SYSY_step(k)%rank.GE.0) THEN

!              otp_block_row    = 1 + comm_SYSY_step(k)%rank / N_blocks_x
!              otp_block_column = 1 + comm_SYSY_step(k)%rank - N_blocks_x * (otp_block_row - 1)
!              otp_c_column = 1 + (otp_block_column-1) / cluster_N_blocks_x
!              otp_c_row    = 1 + (otp_block_row  - 1) / cluster_N_blocks_y

              otp_c_row = 1 + INT(comm_SYSY_step(k)%rank / N_blocks_x) / cluster_N_blocks_y

              comm_SYSY_step(k)%c_indx_y_min = prepared_for_row(otp_c_row)%c_indx_y_min
              comm_SYSY_step(k)%c_indx_y_max = prepared_for_row(otp_c_row)%c_indx_y_max

              comm_SYSY_step(k)%sysy_band_nmin = prepared_for_row(otp_c_row)%sysy_band_nmin
              comm_SYSY_step(k)%sysy_band_nmax = prepared_for_row(otp_c_row)%sysy_band_nmax

           END IF

        END DO
     END DO

     IF (cluster_N_blocks.EQ.1) THEN
! a special case when no additional field calculators are used in the cluster, only the master process
        sysy_strip_nmin = sysy_band_nmin
        sysy_strip_nmax = sysy_band_nmax

     ELSE
! at least one additional field calculator accompanies the master process in the cluster
! note, field_calculator has been already allocated and partially defined in PREPARE_FFT_X

        DEALLOCATE(breakn, STAT=ALLOC_ERR)
        ALLOCATE(breakn(0:cluster_N_blocks), STAT=ALLOC_ERR)

! prepare n-limits for n-strip for each field calculator in the cluster, including the matrix process
        nstart = sysy_band_nmin
        nend   = sysy_band_nmax

! approximate width 
        band_width = (nend - nstart + 1) / cluster_N_blocks
! find number of extra nodes that will be if we use the band_width everywhere 
        extran = (nend - nstart + 1) - cluster_N_blocks * band_width
! in the first extran bands, increase width by 1
        breakn(0) = nstart
        DO row = 1, extran
           breakn(row) = breakn(row-1) + band_width + 1
        END DO
! in the other bands, use standard width
        DO row = extran+1, cluster_N_blocks-1
           breakn(row) = breakn(row-1) + band_width
        END DO
        breakn(cluster_N_blocks) = nend

        sysy_strip_nmin = nstart
        sysy_strip_nmax = breakn(1)-1
        DO k = 2, cluster_N_blocks-1
           field_calculator(k)%sysy_strip_nmin = breakn(k-1)
           field_calculator(k)%sysy_strip_nmax = breakn(k)-1
        END DO
        field_calculator(cluster_N_blocks)%sysy_strip_nmin = breakn(cluster_N_blocks-1)
        field_calculator(cluster_N_blocks)%sysy_strip_nmax = nend

! send bufer dimensions to the associated field calculators
        DO k = 2, cluster_N_blocks

           ibufer(1) = field_calculator(k)%fftx_strip_jmin
           ibufer(2) = field_calculator(k)%fftx_strip_jmax
           ibufer(3) = field_calculator(k)%invfftx_strip_jmin
           ibufer(4) = field_calculator(k)%invfftx_strip_jmax
           ibufer(5) = field_calculator(k)%sysy_strip_nmin
           ibufer(6) = field_calculator(k)%sysy_strip_nmax

           CALL MPI_SEND(ibufer, 6, MPI_INTEGER, field_calculator(k)%rank, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
        END DO

     END IF

     DEALLOCATE(rank_in_row, STAT=ALLOC_ERR)
     DEALLOCATE(comm_pair, STAT=ALLOC_ERR)
     DEALLOCATE(breakn, STAT=ALLOC_ERR)
     DEALLOCATE(prepared_for_row, STAT=ALLOC_ERR)
     DEALLOCATE(prot_comm_step_proc, STAT=ALLOC_ERR)

! define maximal buffer size to be used in SOLVE_POISSON_FFTX_LINSYSY -----------------------
     maxbufsize = 0

     DO k = 1, N_of_comm_steps_X       !####
        IF (comm_FFTX_step(k)%rank.LT.0) CYCLE

        indx_x_begin = c_indx_x_min+1
        indx_x_end   = c_indx_x_max-1
        xsize = indx_x_end - indx_x_begin + 1
        bufsize = (comm_FFTX_step(k)%fftx_band_jmax - comm_FFTX_step(k)%fftx_band_jmin + 1) * xsize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize

        indx_x_begin = comm_FFTX_step(k)%c_indx_x_min+1
        indx_x_end   = comm_FFTX_step(k)%c_indx_x_max-1
        xsize = indx_x_end - indx_x_begin + 1
        bufsize = (fftx_band_jmax - fftx_band_jmin + 1) * xsize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize
     END DO

     DO k = 2, cluster_N_blocks       !####
        bufsize = (global_maximal_i - 1) * (field_calculator(k)%fftx_strip_jmax - field_calculator(k)%fftx_strip_jmin + 1)
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize
     END DO

     DO k = 2, cluster_N_blocks       !####
        bufsize = (global_maximal_i + 1) * (field_calculator(k)%fftx_strip_jmax - field_calculator(k)%fftx_strip_jmin + 1)
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize
     END DO

     DO k = 1, N_of_comm_steps_X      !####
        IF (comm_FFTX_step(k)%rank.LT.0) CYCLE

        indx_x_begin = comm_FFTX_step(k)%c_indx_x_min         ! to make sure that we have pairs corresponding to
        indx_x_end   = comm_FFTX_step(k)%c_indx_x_max - 2     ! real and imaginary parts of fourier components
        IF (indx_x_end.EQ.global_maximal_i-2) indx_x_end = global_maximal_i  ! to include component n=N/2 as a complex number
        xsize = indx_x_end - indx_x_begin + 1
        bufsize = (fftx_band_jmax - fftx_band_jmin + 1) * xsize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize

        indx_x_begin = c_indx_x_min
        indx_x_end   = c_indx_x_max-2
        IF (indx_x_end.EQ.global_maximal_i-2) indx_x_end = global_maximal_i
        xsize = indx_x_end - indx_x_begin + 1
        bufsize = (comm_FFTX_step(k)%fftx_band_jmax - comm_FFTX_step(k)%fftx_band_jmin + 1) * xsize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize
     END DO

     DO k = 1, N_of_comm_steps_Y       !####
        IF (comm_SYSY_step(k)%rank.LT.0) CYCLE           ! added on 21-Feb-2017

        IF (c_indx_y_min.GT.0) THEN
           indx_y_begin = c_indx_y_min+1
        ELSE
           indx_y_begin = 0
        END IF
        IF (c_indx_y_max.LT.global_maximal_j) THEN
           indx_y_end   = c_indx_y_max-1
        ELSE
           indx_y_end = global_maximal_j
        END IF
        ysize = 2 * (indx_y_end - indx_y_begin + 1)
        bufsize = (comm_SYSY_step(k)%sysy_band_nmax - comm_SYSY_step(k)%sysy_band_nmin + 1) * ysize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize

        IF (comm_SYSY_step(k)%c_indx_y_min.GT.0) THEN
           indx_y_begin = comm_SYSY_step(k)%c_indx_y_min+1
        ELSE
           indx_y_begin = 0
        END IF
        IF (comm_SYSY_step(k)%c_indx_y_max.LT.global_maximal_j) THEN
           indx_y_end   = comm_SYSY_step(k)%c_indx_y_max-1
        ELSE
           indx_y_end = global_maximal_j
        END IF
        bufsize = (sysy_band_nmax - sysy_band_nmin + 1) * 2 * (indx_y_end - indx_y_begin + 1)
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize
     END DO

     DO k = 2, cluster_N_blocks       !####
        ysize = 2 * (global_maximal_j + 1)
        bufsize = ysize * (field_calculator(k)%sysy_strip_nmax - field_calculator(k)%sysy_strip_nmin + 1)        ! try to use %indx_x_min/max instead
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize
     END DO

     DO k = 1, N_of_comm_steps_Y       !####
        IF (comm_SYSY_step(k)%rank.LT.0) CYCLE         ! changed on 21-Feb-2017, was comm_FFTX_step(k)

        indx_y_begin = 2 * comm_SYSY_step(k)%c_indx_y_min  !MAX(2, 2 * comm_SYSY_step(k)%c_indx_y_min)                           ! even, real part of the first complex amplitude
        indx_y_end = 2 * comm_SYSY_step(k)%c_indx_y_max + 1 !MIN(2 * comm_SYSY_step(k)%c_indx_y_max + 1, 2 * global_maximal_j - 1)  ! odd, includes the imaginary part of the last complex amplitude
        ysize = indx_y_end - indx_y_begin + 1
        bufsize = (sysy_band_nmax - sysy_band_nmin + 1) * ysize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize

        indx_y_begin = c_indx_y_min !MAX(c_indx_y_min, 1) 
        indx_y_end   = c_indx_y_max !MIN(c_indx_y_max, global_maximal_j-1)
        ysize = 2 * (indx_y_end - indx_y_begin + 1)
        bufsize = (comm_SYSY_step(k)%sysy_band_nmax - comm_SYSY_step(k)%sysy_band_nmin + 1) * ysize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize
     END DO

     DO k = 1, N_of_comm_steps_X       !####
        IF (comm_FFTX_step(k)%rank.LT.0) CYCLE

        indx_x_begin = c_indx_x_min
        indx_x_end   = c_indx_x_max-2
        IF (c_column.EQ.N_clusters_x) indx_x_end = c_indx_x_max
        xsize = indx_x_end - indx_x_begin + 1
        bufsize = (comm_FFTX_step(k)%invfftx_band_jmax - comm_FFTX_step(k)%invfftx_band_jmin + 1) * xsize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize

        indx_x_begin = comm_FFTX_step(k)%c_indx_x_min
        indx_x_end   = comm_FFTX_step(k)%c_indx_x_max-2
        IF (indx_x_end.EQ.global_maximal_i-2) indx_x_end = global_maximal_i
        xsize = indx_x_end - indx_x_begin + 1
        bufsize = (invfftx_band_jmax - invfftx_band_jmin + 1) * xsize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize
     END DO

     DO k = 2, cluster_N_blocks       !####
        indx_y_begin = field_calculator(k)%invfftx_strip_jmin
        indx_y_end   = field_calculator(k)%invfftx_strip_jmax
        bufsize = (global_maximal_i + 1) * (indx_y_end - indx_y_begin + 1)
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize
     END DO

     DO k = 1, N_of_comm_steps_X       !####

        IF (comm_FFTX_step(k)%rank.LT.0) CYCLE
        
        indx_x_begin = comm_FFTX_step(k)%c_indx_x_min
        indx_x_end   = comm_FFTX_step(k)%c_indx_x_max
        xsize = indx_x_end - indx_x_begin + 1
        bufsize = (invfftx_band_jmax - invfftx_band_jmin + 1) * xsize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize

        indx_x_begin = c_indx_x_min
        indx_x_end   = c_indx_x_max
        xsize = indx_x_end - indx_x_begin + 1
        bufsize = (comm_FFTX_step(k)%invfftx_band_jmax - comm_FFTX_step(k)%invfftx_band_jmin + 1) * xsize
        IF (bufsize.GT.maxbufsize) maxbufsize = bufsize
     END DO

  ELSE

     CALL MPI_RECV(ibufer, 6, MPI_INTEGER, field_master, field_master, MPI_COMM_WORLD, stattus, ierr)
     fftx_strip_jmin = ibufer(1)
     fftx_strip_jmax = ibufer(2)
     invfftx_strip_jmin = ibufer(3)
     invfftx_strip_jmax = ibufer(4)
     sysy_strip_nmin = ibufer(5)
     sysy_strip_nmax = ibufer(6)
 
! define maximal buffer size to be used in SOLVE_POISSON_FFTX_LINSYSY -----------------------
     maxbufsize = 0

!####
     bufsize = (global_maximal_i - 1) * (fftx_strip_jmax - fftx_strip_jmin + 1)
     IF (bufsize.GT.maxbufsize) maxbufsize = bufsize

!####
     bufsize = (global_maximal_i + 1) * (fftx_strip_jmax - fftx_strip_jmin + 1)
     IF (bufsize.GT.maxbufsize) maxbufsize = bufsize

!####
     ysize = 2 * (global_maximal_j + 1)
     bufsize = ysize * (sysy_strip_nmax - sysy_strip_nmin + 1)
     IF (bufsize.GT.maxbufsize) maxbufsize = bufsize

!####
     bufsize = (global_maximal_i + 1) * (invfftx_strip_jmax - invfftx_strip_jmin + 1)
     IF (bufsize.GT.maxbufsize) maxbufsize = bufsize

  END IF

! prepare boundary flags
  wall_below_is_metal = .TRUE.
  DO n = 1, N_of_boundary_objects
     IF ( (whole_object(n)%segment(1)%jstart.EQ.0).AND. &
        & (whole_object(n)%segment(1)%jend.EQ.0).AND. &
        & (whole_object(n)%object_type.EQ.DIELECTRIC) ) THEN
        wall_below_is_metal = .FALSE.
     END IF
  END DO

  wall_above_is_metal = .TRUE.
  DO n = 1, N_of_boundary_objects
     IF ( (whole_object(n)%segment(1)%jstart.EQ.global_maximal_j).AND. &
        & (whole_object(n)%segment(1)%jend.EQ.global_maximal_j).AND. &
        & (whole_object(n)%object_type.EQ.DIELECTRIC) ) THEN
        wall_above_is_metal = .FALSE.
     END IF
  END DO

  two_dielectric_walls = .TRUE.
  IF (wall_below_is_metal.OR.wall_above_is_metal) two_dielectric_walls = .FALSE. 

! flags wall_below/above_is_metal work when the bottom/top boundary is a solid metal wall or a segmented wall 
! (combination of METAL_WALL and VACUUM_GAP boundary objects) with potential profile that is a known function of time

! precalculate coefficients of the linear system --------------------------------------------

  ALLOCATE(a_eq(0:global_maximal_j, sysy_strip_nmin:sysy_strip_nmax), STAT = ALLOC_ERR)  ! diagonal main+1 (above)
  
!  DO n = sysy_strip_nmin, sysy_strip_nmax !0, N_grid_x/2 !-1
!     main_diag_term = -2.0_8 - 4.0_8 * pi * pi * dble(n) * dble(n) / (dble(global_maximal_i-1) * dble(global_maximal_i-1))
!     main_diag_term_prim_upper = -1.0_8 - 2.0_8 * pi * (pi * dble(n) + whole_object(2)%eps_diel) * dble(n) / (dble(global_maximal_i-1) * dble(global_maximal_i-1))         !####
!     main_diag_term_prim_lower = -1.0_8 - 2.0_8 * pi * (pi * dble(n) + whole_object(4)%eps_diel) * dble(n) / (dble(global_maximal_i-1) * dble(global_maximal_i-1))         !####
!! diagonalization . . .
!     IF (wall_below_is_metal.OR.(n.EQ.0)) THEN
!!#### if the wall y=0 is metal or n=0
!        a_eq(0, n) = 0.0_8
!     ELSE
!!#### if the wall y=0 is dielectric and n>0
!        a_eq(0, n) = 1.0_8 / main_diag_term_prim_lower
!     END IF
!     DO j = 1, global_maximal_j - 1
!        a_eq(j, n) = 1.0_8 / (main_diag_term - a_eq(j-1, n))
!     END DO
!     IF (wall_above_is_metal.OR.(n.EQ.0)) THEN
!!#### if the wall y=ymax is metal or n=0
!        a_eq(global_maximal_j, n) = 0.0_8
!     ELSE
!!#### if the wall y=ymax is dielectric and n>0
!        a_eq(global_maximal_j, n) = 1.0_8 / (main_diag_term_prim_upper - a_eq(global_maximal_j-1, n))
!     END IF
!  END DO

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  DO n = sysy_strip_nmin, sysy_strip_nmax !0, N_grid_x/2 !-1

     main_diag_term = -2.0_8 - 4.0_8 * pi * pi * dble(n) * dble(n) / (dble(global_maximal_i-1) * dble(global_maximal_i-1))
     
     main_diag_term_prim_upper = -1.0_8 - 2.0_8 * pi * (pi * dble(n) + whole_object(2)%eps_diel) * dble(n) / (dble(global_maximal_i-1) * dble(global_maximal_i-1))         !####
     main_diag_term_prim_lower = -1.0_8 - 2.0_8 * pi * (pi * dble(n) + whole_object(4)%eps_diel) * dble(n) / (dble(global_maximal_i-1) * dble(global_maximal_i-1))         !####

! if n=0
!     main_diag_term            = -2.0_8
!     main_diag_term_prim_upper = -1.0_8         !####
!     main_diag_term_prim_lower = -1.0_8         !####

! diagonalization . . .

!was     IF (wall_below_is_metal.OR.(n.EQ.0)) THEN
     IF (wall_below_is_metal) THEN
!#### if the wall y=0 is metal
        a_eq(0, n) = 0.0_8
     ELSE
!#### if the wall y=0 is dielectric

        if (wall_above_is_metal) then
           a_eq(0, n) = 1.0_8 / main_diag_term_prim_lower
        else
!#### both walls are dielectric
           if (n.eq.0) then
!#### set potential at -infinity to zero (note that potential at +infinity may be nonzero in this case)
              a_eq(0, n) = 0.0_8
           else
              a_eq(0, n) = 1.0_8 / main_diag_term_prim_lower
           end if
        end if

!was        a_eq(0, n) = 1.0_8 / main_diag_term_prim_lower
     END IF

     IF (ht_grid_requested) THEN

        DO j = 1, grid_j - 1
           a_eq(j, n) = 1.0_8 / (main_diag_term - a_eq(j-1, n))
        END DO

        a_eq(grid_j, n) = 0.0_8

        DO j = grid_j + 1, global_maximal_j - 1
           a_eq(j, n) = 1.0_8 / (main_diag_term - a_eq(j-1, n))
        END DO

     ELSE

        DO j = 1, global_maximal_j - 1
           a_eq(j, n) = 1.0_8 / (main_diag_term - a_eq(j-1, n))
        END DO

     END IF

!was     IF (wall_above_is_metal.OR.(n.EQ.0)) THEN
     IF (wall_above_is_metal) THEN
!#### if the wall y=ymax is metal
        a_eq(global_maximal_j, n) = 0.0_8
     ELSE
!#### if the wall y=ymax is dielectric
        a_eq(global_maximal_j, n) = 1.0_8 / (main_diag_term_prim_upper - a_eq(global_maximal_j-1, n))
     END IF

  END DO

END SUBROUTINE PREPARE_SYS_Y
