!------------------------------------------
!
! identify neighbours of each block at left/right/below/above
! these neighbors will be used by the field solver
!
! this procedure is not called when periodic boundaries are involved
!
SUBROUTINE IDENTIFY_BLOCK_NEIGHBOURS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries

  IMPLICIT NONE

  INTEGER ALLOC_ERR

  character(17) boxproc_filename ! box_proc_NNNN.dat
  integer devid  ! id of device for writing the file
  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface

! moved to INITIATE_PARAMETERS
!  block_row    = 1 + Rank_of_process / N_blocks_x
!  block_column = 1 + Rank_of_process - N_blocks_x * (block_row - 1)

  indx_x_min = (block_column-1) * N_grid_block_x
  indx_x_max =  block_column    * N_grid_block_x + 1 
  indx_y_min = (block_row-1) * N_grid_block_y
  indx_y_max =  block_row    * N_grid_block_y + 1 

  IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN
! these arrays are used in SOR
     ALLOCATE(rho_i(indx_x_min:indx_x_max, indx_y_min:indx_y_max), STAT=ALLOC_ERR)  !
     ALLOCATE(rho_e(indx_x_min:indx_x_max, indx_y_min:indx_y_max), STAT=ALLOC_ERR)  ! these arrays will never be resized or removed
     ALLOCATE(  phi(indx_x_min:indx_x_max, indx_y_min:indx_y_max), STAT=ALLOC_ERR)  !
  END IF

  X_area_min = DBLE(indx_x_min)
  X_area_max = DBLE(indx_x_max)
  Y_area_min = DBLE(indx_y_min)
  Y_area_max = DBLE(indx_y_max)

!           shift_x = (m-1) * (N_grid_proc_x-1)
!           shift_y = (n-1) * (N_grid_proc_y-1)

  Rank_of_process_left  = Rank_of_process - 1
  Rank_of_process_right = Rank_of_process + 1
  Rank_of_process_below = N_blocks_x * (block_row-2) + block_column - 1
  Rank_of_process_above = N_blocks_x *  block_row    + block_column - 1

  IF (block_column.EQ.1)          Rank_of_process_left  = -1
  IF (block_column.EQ.N_blocks_x) Rank_of_process_right = -1
  IF (block_row.EQ.1)             Rank_of_process_below = -1
  IF (block_row.EQ.N_blocks_y)    Rank_of_process_above = -1

  IF (MOD(block_row,2).EQ.1) THEN
     IF (MOD(block_column,2).EQ.1) THEN
        WHITE = .TRUE.
     ELSE
        WHITE = .FALSE.
     END IF
  ELSE
     IF (MOD(block_column,2).EQ.0) THEN
        WHITE = .TRUE.
     ELSE
        WHITE = .FALSE.
     END IF
  END IF
        
  IF (WHITE) THEN
     PRINT '("Process ",i4," : WHITE : grid column number= ",i4," row number= ",i4," ; neighbours : left ",i4," right ",i4," above ",i4," below ",i4)', &
          & Rank_of_process, block_column, block_row, Rank_of_process_left, Rank_of_process_right, Rank_of_process_above, Rank_of_process_below
  ELSE
     PRINT '("Process ",i4," : BLACK : grid column number= ",i4," row number= ",i4," ; neighbours : left ",i4," right ",i4," above ",i4," below ",i4)', &
          & Rank_of_process, block_column, block_row, Rank_of_process_left, Rank_of_process_right, Rank_of_process_above, Rank_of_process_below
  END IF

  boxproc_filename = 'box_proc_NNNN.dat'
  boxproc_filename(10:13) = convert_int_to_txt_string(Rank_of_process, 4)

  devid = 9+Rank_of_process

  open (devid, file=boxproc_filename)
  write (devid, '(3(2x,i5),2(2x,e14.7))') indx_x_min, indx_y_min, Rank_of_process, X_area_min * delta_x_m, Y_area_min * delta_x_m
  write (devid, '(3(2x,i5),2(2x,e14.7))') indx_x_min, indx_y_max, Rank_of_process, X_area_min * delta_x_m, Y_area_max * delta_x_m
  write (devid, '(3(2x,i5),2(2x,e14.7))') indx_x_max, indx_y_max, Rank_of_process, X_area_max * delta_x_m, Y_area_max * delta_x_m
  write (devid, '(3(2x,i5),2(2x,e14.7))') indx_x_max, indx_y_min, Rank_of_process, X_area_max * delta_x_m, Y_area_min * delta_x_m
  write (devid, '(3(2x,i5),2(2x,e14.7))') indx_x_min, indx_y_min, Rank_of_process, X_area_min * delta_x_m, Y_area_min * delta_x_m
  close (devid, status = 'keep')

  print '("Process ",i4," : IDENTIFY_BLOCK_NEIGHBOURS : file ",A17," is ready")', Rank_of_process, boxproc_filename

END SUBROUTINE IDENTIFY_BLOCK_NEIGHBOURS

!-------------------------------------------------
! this subroutine works only with master processes
!
! this procedure is not called when periodic boundaries are involved
!!
SUBROUTINE IDENTIFY_BLOCK_BOUNDARIES
  
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries

  use mpi

  IMPLICIT NONE

  INTEGER errcode,ierr

  INTEGER n, m
  INTEGER jstart, jend, istart, iend
  INTEGER minimal_index, maximal_index

! most probable default assumption - no boundary objects

  N_of_local_object_parts = 0

  N_of_local_object_parts_left  = 0
  N_of_local_object_parts_right = 0
  N_of_local_object_parts_above = 0
  N_of_local_object_parts_below = 0

  block_has_symmetry_plane_X_left = .FALSE.

  IF ( (Rank_of_process_left.GE.0).AND. &
       (Rank_of_process_right.GE.0).AND. &
       (Rank_of_process_below.GE.0).AND. &
       (Rank_of_process_above.GE.0) ) THEN
     PRINT '(" Block ",i4," is not connected to any boundary object")', Rank_of_process
     RETURN
  END IF

! since we are here, there is at least one boundary object

! check the left edge of the block
  IF (Rank_of_process_left.EQ.-1) THEN
     DO n = 1, N_of_boundary_objects
        DO m = 1, whole_object(n)%number_of_segments
           IF ( (whole_object(n)%segment(m)%istart.EQ.indx_x_min) .AND. &
              & (whole_object(n)%segment(m)%iend.EQ.indx_x_min) ) THEN
              
              jstart = MAX(whole_object(n)%segment(m)%jstart, indx_y_min)
              jend   = MIN(whole_object(n)%segment(m)%jend, indx_y_max)

              IF (jstart.LT.jend) THEN
                 N_of_local_object_parts = N_of_local_object_parts + 1
                 IF (N_of_local_object_parts.GT.max_N_of_local_object_parts) THEN
                    PRINT '("Process ",i4," : ERROR-1 in IDENTIFY_BLOCK_BOUNDARIES : maximal number of boundary parts exceeded : ",i4)', Rank_of_process, N_of_local_object_parts
                    errcode=100
                    CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
                 END IF
                 local_object_part(N_of_local_object_parts)%object_number = n
                 local_object_part(N_of_local_object_parts)%istart = indx_x_min
                 local_object_part(N_of_local_object_parts)%jstart = jstart 
                 local_object_part(N_of_local_object_parts)%iend   = indx_x_min
                 local_object_part(N_of_local_object_parts)%jend   = jend

                 N_of_local_object_parts_left = N_of_local_object_parts_left+1
                 index_of_local_object_part_left(N_of_local_object_parts_left) = N_of_local_object_parts

                 IF (whole_object(n)%object_type.EQ.SYMMETRY_PLANE) block_has_symmetry_plane_X_left = .TRUE.

              END IF

           END IF
        END DO
     END DO
  END IF

! check the top edge of the block
  IF (Rank_of_process_above.EQ.-1) THEN
     DO n = 1, N_of_boundary_objects
        DO m = 1, whole_object(n)%number_of_segments
           IF ( (whole_object(n)%segment(m)%jstart.EQ.indx_y_max) .AND. &
              & (whole_object(n)%segment(m)%jend.EQ.indx_y_max) ) THEN
              
              istart = MAX(whole_object(n)%segment(m)%istart, indx_x_min)
              iend   = MIN(whole_object(n)%segment(m)%iend, indx_x_max)

              IF (istart.LT.iend) THEN
                 N_of_local_object_parts = N_of_local_object_parts + 1
                 IF (N_of_local_object_parts.GT.max_N_of_local_object_parts) THEN
                    PRINT '("Process ",i4," : ERROR-2 in IDENTIFY_BLOCK_BOUNDARIES : maximal number of boundary parts exceeded : ",i4)', Rank_of_process, N_of_local_object_parts
                    errcode=101
                    CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
                 END IF
                 local_object_part(N_of_local_object_parts)%object_number = n
                 local_object_part(N_of_local_object_parts)%istart = istart
                 local_object_part(N_of_local_object_parts)%jstart = indx_y_max
                 local_object_part(N_of_local_object_parts)%iend   = iend
                 local_object_part(N_of_local_object_parts)%jend   = indx_y_max

                 N_of_local_object_parts_above = N_of_local_object_parts_above+1
                 index_of_local_object_part_above(N_of_local_object_parts_above) = N_of_local_object_parts
              END IF

           END IF
        END DO
     END DO
  END IF

! check the right edge of the block
  IF (Rank_of_process_right.EQ.-1) THEN
     DO n = 1, N_of_boundary_objects
        DO m = 1, whole_object(n)%number_of_segments
           IF ( (whole_object(n)%segment(m)%istart.EQ.indx_x_max) .AND. &
              & (whole_object(n)%segment(m)%iend.EQ.indx_x_max) ) THEN
              
              jstart = MAX(whole_object(n)%segment(m)%jstart, indx_y_min)
              jend   = MIN(whole_object(n)%segment(m)%jend, indx_y_max)

              IF (jstart.LT.jend) THEN
                 N_of_local_object_parts = N_of_local_object_parts + 1
                 IF (N_of_local_object_parts.GT.max_N_of_local_object_parts) THEN
                    PRINT '("Process ",i4," : ERROR-3 in IDENTIFY_BLOCK_BOUNDARIES : maximal number of boundary parts exceeded : ",i4)', Rank_of_process, N_of_local_object_parts
                    errcode=102
                    CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
                 END IF
                 local_object_part(N_of_local_object_parts)%object_number = n
                 local_object_part(N_of_local_object_parts)%istart = indx_x_max
                 local_object_part(N_of_local_object_parts)%jstart = jstart 
                 local_object_part(N_of_local_object_parts)%iend   = indx_x_max
                 local_object_part(N_of_local_object_parts)%jend   = jend

                 N_of_local_object_parts_right = N_of_local_object_parts_right+1
                 index_of_local_object_part_right(N_of_local_object_parts_right) = N_of_local_object_parts
              END IF

           END IF
        END DO
     END DO
  END IF

! check the bottom edge of the block
  IF (Rank_of_process_below.EQ.-1) THEN
     DO n = 1, N_of_boundary_objects
        DO m = 1, whole_object(n)%number_of_segments
           IF ( (whole_object(n)%segment(m)%jstart.EQ.indx_y_min) .AND. &
              & (whole_object(n)%segment(m)%jend.EQ.indx_y_min) ) THEN
              
              istart = MAX(whole_object(n)%segment(m)%istart, indx_x_min)
              iend   = MIN(whole_object(n)%segment(m)%iend, indx_x_max)

              IF (istart.LT.iend) THEN
                 N_of_local_object_parts = N_of_local_object_parts + 1
                 IF (N_of_local_object_parts.GT.max_N_of_local_object_parts) THEN
                    PRINT '("Process ",i4," : ERROR-4 in IDENTIFY_BLOCK_BOUNDARIES : maximal number of boundary parts exceeded : ",i4)', Rank_of_process, N_of_local_object_parts
                    errcode=103
                    CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
                 END IF
                 local_object_part(N_of_local_object_parts)%object_number = n
                 local_object_part(N_of_local_object_parts)%istart = istart
                 local_object_part(N_of_local_object_parts)%jstart = indx_y_min
                 local_object_part(N_of_local_object_parts)%iend   = iend
                 local_object_part(N_of_local_object_parts)%jend   = indx_y_min

                 N_of_local_object_parts_below = N_of_local_object_parts_below+1
                 index_of_local_object_part_below(N_of_local_object_parts_below) = N_of_local_object_parts
              END IF

           END IF
        END DO
     END DO
  END IF

END SUBROUTINE IDENTIFY_BLOCK_BOUNDARIES

!-------------------------------------------------
! this subroutine works mostly with master processes
! but other processes are involved as well
! in order to create new temporary communicators
!
SUBROUTINE INCLUDE_BLOCK_PERIODICITY
  
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries

  use mpi

  IMPLICIT NONE

  INTEGER errcode,ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER periodic_message(1:12)
  INTEGER n, m, nwo

  INTEGER, ALLOCATABLE :: all_periodic_messages(:,:)
  INTEGER ALLOC_ERR
  LOGICAL pair_found
  INTEGER ibufer(1:4)

  ! check the presence of periodic boundaries

  periodic_message(1:12)=-1

  block_periodic_boundary_X_left  = .FALSE.
  block_periodic_boundary_X_right = .FALSE.
  block_periodic_boundary_Y_below = .FALSE.
  block_periodic_boundary_Y_above = .FALSE.

! left boundary (periodicity in X)
  DO n = 1, N_of_local_object_parts_left
     m = index_of_local_object_part_left(n)
     nwo = local_object_part(m)%object_number
     IF (whole_object(nwo)%object_type.EQ.PERIODIC_PIPELINE_X) THEN
        periodic_message(1) = nwo   ! 1<=nwo<=N_of_boundary_objects
        periodic_message(2) = local_object_part(m)%jstart
        periodic_message(3) = local_object_part(m)%jend
     END IF
  END DO

! right boundary (periodicity in X)
  DO n = 1, N_of_local_object_parts_right
     m = index_of_local_object_part_right(n)
     nwo = local_object_part(m)%object_number
     IF (whole_object(nwo)%object_type.EQ.PERIODIC_PIPELINE_X) THEN
        periodic_message(4) = nwo   ! 1<=nwo<=N_of_boundary_objects
        periodic_message(5) = local_object_part(m)%jstart
        periodic_message(6) = local_object_part(m)%jend
     END IF
  END DO

! bottom boundary (periodicity in Y)
  DO n = 1, N_of_local_object_parts_below
     m = index_of_local_object_part_below(n)
     nwo = local_object_part(m)%object_number
     IF (whole_object(nwo)%object_type.EQ.PERIODIC_PIPELINE_Y) THEN
        periodic_message(7) = nwo   ! 1<=nwo<=N_of_boundary_objects
        periodic_message(8) = local_object_part(m)%istart
        periodic_message(9) = local_object_part(m)%iend
     END IF
  END DO
  
! top boundary (periodicity in Y)
  DO n = 1, N_of_local_object_parts_above
     m = index_of_local_object_part_above(n)
     nwo = local_object_part(m)%object_number
     IF (whole_object(nwo)%object_type.EQ.PERIODIC_PIPELINE_Y) THEN
        periodic_message(10) = nwo   ! 1<=nwo<=N_of_boundary_objects
        periodic_message(11) = local_object_part(m)%istart
        periodic_message(12) = local_object_part(m)%iend
     END IF
  END DO

! send/receive-process periodic_message to the process with rank zero

  IF (Rank_of_process.EQ.0) THEN

     ALLOCATE(all_periodic_messages(1:14,0:N_of_processes-1), STAT=ALLOC_ERR)
     all_periodic_messages = -1
     all_periodic_messages(1:12,0) = periodic_message(1:12)
! receive periodic_message from all other clusters
     DO n = 1, N_of_processes-1
        CALL MPI_RECV(all_periodic_messages(1:12,n), 12, MPI_INTEGER, n, n, MPI_COMM_WORLD, stattus, ierr)        
     END DO

! find pairs of blocks with periodic boundary in X
     DO n = 0, N_of_processes-1
! find one such a cluster
        IF (all_periodic_messages(1,n).GT.0) THEN
! find a matching pair
           pair_found=.FALSE.
           DO m = 0, N_of_processes-1
              IF ((all_periodic_messages(4,m).GT.0).AND.(all_periodic_messages(5,m).EQ.all_periodic_messages(2,n)).AND.(all_periodic_messages(6,m).EQ.all_periodic_messages(3,n))) THEN
                 pair_found=.TRUE.
                 all_periodic_messages(13,n) = m  ! MPI_COMM_WORLD rank
                 all_periodic_messages(13,m) = n  ! MPI_COMM_WORLD rank
                 EXIT
              END IF
           END DO
           IF (.NOT.pair_found) THEN
! error message if the pair is not found
              PRINT '("Process ",i4," :: Error in INCLUDE_BLOCK_PERIODICITY :: cannot find a matching pair for X-periodic boundary of process ",i4)', Rank_of_process, n
              errcode=104
              CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
           END IF
        END IF
     END DO

! find pairs of blocks with periodic boundary in Y
     DO n = 0, N_of_processes-1
! find one such a cluster
        IF (all_periodic_messages(7,n).GT.0) THEN
! find a matching pair
           pair_found=.FALSE.
           DO m = 0, N_of_processes-1
              IF ((all_periodic_messages(10,m).GT.0).AND.(all_periodic_messages(11,m).EQ.all_periodic_messages(8,n)).AND.(all_periodic_messages(12,m).EQ.all_periodic_messages(9,n))) THEN
                 pair_found=.TRUE.
                 all_periodic_messages(14,n) = m  ! MPI_COMM_WORLD rank
                 all_periodic_messages(14,m) = n  ! MPI_COMM_WORLD rank
                 EXIT
              END IF
           END DO
           IF (.NOT.pair_found) THEN
! error message if the pair is not found
              PRINT '("Process ",i4," :: Error in INCLUDE_BLOCK_PERIODICITY :: cannot find a matching pair for Y-periodic boundary of process ",i4)', Rank_of_process, n
              errcode=105
              CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
           END IF
        END IF
     END DO

     ! send  information about periodicity pairs to each process
     DO n = 1, N_of_processes-1
        CALL MPI_SEND(all_periodic_messages(13:14,n), 2, MPI_INTEGER, n, Rank_of_process, MPI_COMM_WORLD, ierr)        
     END DO
     
! process own periodicity pair info
     IF (all_periodic_messages(13,0).GT.0) THEN
        IF (Rank_of_process_left.LT.0) THEN
           Rank_of_process_left  = all_periodic_messages(13,0)
           block_periodic_boundary_X_left = .TRUE.
        END IF
        IF (Rank_of_process_right.LT.0) THEN
           Rank_of_process_right = all_periodic_messages(13,0)
           block_periodic_boundary_X_right = .TRUE.
        END IF
     END IF
     IF (all_periodic_messages(14,0).GT.0) THEN
        IF (Rank_of_process_below.LT.0) THEN
           Rank_of_process_below = all_periodic_messages(14,0)
           block_periodic_boundary_Y_below = .TRUE.
        END IF
        IF (Rank_of_process_above.LT.0) THEN
           Rank_of_process_above = all_periodic_messages(14,0)
           block_periodic_boundary_Y_above = .TRUE.
        END IF
     END IF

     DEALLOCATE(all_periodic_messages, STAT=ALLOC_ERR)
 
  ELSE
     
! send periodic_message to zero-rank process
     CALL MPI_SEND(periodic_message, 12, MPI_INTEGER, 0, Rank_of_process, MPI_COMM_WORLD, ierr)
     
! receive information about periodicity pairs
     CALL MPI_RECV(ibufer(1:2), 2, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, stattus, ierr) 
     IF (ibufer(1).GE.0) THEN
        IF (Rank_of_process_left.LT.0) THEN
           Rank_of_process_left  = ibufer(1)
           block_periodic_boundary_X_left = .TRUE.
        END IF
        IF (Rank_of_process_right.LT.0) THEN
           Rank_of_process_right = ibufer(1)
           block_periodic_boundary_X_right = .TRUE.
        END IF
     END IF
     IF (ibufer(2).GE.0) THEN
        IF (Rank_of_process_below.LT.0) THEN
           Rank_of_process_below = ibufer(2)
           block_periodic_boundary_Y_below = .TRUE.
        END IF
        IF (Rank_of_process_above.LT.0) THEN
           Rank_of_process_above = ibufer(2)
           block_periodic_boundary_Y_above = .TRUE.
        END IF
     END IF

  END IF

print '("Block ",i4," :: P.B. L/R/B/A :: ",4(1x,L1))', Rank_of_process, block_periodic_boundary_X_left, block_periodic_boundary_X_right, block_periodic_boundary_Y_below, block_periodic_boundary_Y_above

END SUBROUTINE INCLUDE_BLOCK_PERIODICITY

!---------------------------------------------
! this procedure must be called after INCLUDE_BLOCK_PERIODICITY
! because there new neighbors may be assigned due to periodicity
!
SUBROUTINE CALCULATE_BLOCK_OFFSET
  
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries

  use mpi

  IMPLICIT NONE

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER N_of_nodes_to_solve_x, N_of_nodes_to_solve_y

  INTEGER, ALLOCATABLE :: offset_of_block(:)
  INTEGER ALLOC_ERR

  INTEGER i, itemp
  INTEGER N_to_solve_in_block_im1

  INTEGER bottom_right_inner_node
  INTEGER bottom_left_inner_node
  INTEGER left_top_inner_node
  INTEGER left_bottom_inner_node

  INTEGER ibufer(2)

  N_of_nodes_to_solve_x = indx_x_max - indx_x_min - 1
  IF (Rank_of_process_left.LT.0)  N_of_nodes_to_solve_x = N_of_nodes_to_solve_x + 1
  IF (Rank_of_process_right.LT.0) N_of_nodes_to_solve_x = N_of_nodes_to_solve_x + 1
  
  N_of_nodes_to_solve_y = indx_y_max - indx_y_min - 1
  IF (Rank_of_process_below.LT.0) N_of_nodes_to_solve_y = N_of_nodes_to_solve_y + 1
  IF (Rank_of_process_above.LT.0) N_of_nodes_to_solve_y = N_of_nodes_to_solve_y + 1
  
  block_N_of_nodes_to_solve = N_of_nodes_to_solve_x * N_of_nodes_to_solve_y

  ALLOCATE(offset_of_block(0:N_of_processes-1), STAT=ALLOC_ERR)

  CALL MPI_GATHER(block_N_of_nodes_to_solve, 1, MPI_INTEGER, offset_of_block, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  
  IF (Rank_of_process.EQ.0) THEN
     N_to_solve_in_block_im1 = block_N_of_nodes_to_solve
     offset_of_block(0) = 0
     DO i = 1, N_of_processes-1
        itemp = offset_of_block(i)
        offset_of_block(i) = offset_of_block(i-1) + N_to_solve_in_block_im1
        N_to_solve_in_block_im1 = itemp
     END DO
     N_to_solve_total = offset_of_block(N_of_processes-1) + N_to_solve_in_block_im1
  END IF

  CALL MPI_BCAST(offset_of_block, N_of_processes, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(N_to_solve_total, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  global_offset = offset_of_block(Rank_of_process)-1 !######### check, should work if petsc numbering of rows/columns begins with zero

  DEALLOCATE(offset_of_block, STAT=ALLOC_ERR)

! talk to left/right neighbor to get process_left_bottom_right_inner_node    process_right_bottom_left_inner_node 
!                                    process_left_solved_nodes_row_length    process_right_solved_nodes_row_length

! talk to neighbors above/below to get process_above_left_bottom_inner_node and process_below_left_top_inner_node

  IF (WHITE) THEN  
! "white processes"

     IF (Rank_of_process_right.GE.0) THEN
! ## 1 ## send right
        bottom_right_inner_node = global_offset + N_of_nodes_to_solve_x
        IF (Rank_of_process_below.LT.0) bottom_right_inner_node = bottom_right_inner_node + N_of_nodes_to_solve_x  ! skip line of nodes along the bottom boundary
        ibufer(1) = bottom_right_inner_node
        ibufer(2) = N_of_nodes_to_solve_x
        CALL MPI_SEND(ibufer, 2, MPI_INTEGER, Rank_of_process_right, Rank_of_process, MPI_COMM_WORLD, ierr) 
     END IF

     IF (Rank_of_process_left.GE.0) THEN
! ## 2 ## send left
        bottom_left_inner_node = global_offset + 1
        IF (Rank_of_process_below.LT.0) bottom_left_inner_node = bottom_left_inner_node + N_of_nodes_to_solve_x  ! skip line of nodes along the bottom boundary
        ibufer(1) = bottom_left_inner_node
        ibufer(2) = N_of_nodes_to_solve_x
        CALL MPI_SEND(ibufer, 2, MPI_INTEGER, Rank_of_process_left, Rank_of_process, MPI_COMM_WORLD, ierr)
     END IF

     IF (Rank_of_process_left.GE.0) THEN
! ## 3 ## receive from left
        CALL MPI_RECV(ibufer, 2, MPI_INTEGER, Rank_of_process_left, Rank_of_process_left, MPI_COMM_WORLD, stattus, ierr)
        process_left_bottom_right_inner_node = ibufer(1)
        process_left_solved_nodes_row_length = ibufer(2)
     END IF
 
     IF (Rank_of_process_right.GE.0) THEN
! ## 4 ## receive from right
        CALL MPI_RECV(ibufer, 2, MPI_INTEGER, Rank_of_process_right, Rank_of_process_right, MPI_COMM_WORLD, stattus, ierr)
        process_right_bottom_left_inner_node = ibufer(1)
        process_right_solved_nodes_row_length = ibufer(2)
     END IF

     IF (Rank_of_process_above.GE.0) THEN
! ## 5 ## send up
        left_top_inner_node = global_offset + N_of_nodes_to_solve_x * (N_of_nodes_to_solve_y-1) + 1
        IF (Rank_of_process_left.LT.0) left_top_inner_node = left_top_inner_node + 1  ! skip the node at the left boundary
        ibufer(1) = left_top_inner_node
        CALL MPI_SEND(ibufer(1), 1, MPI_INTEGER, Rank_of_process_above, Rank_of_process, MPI_COMM_WORLD, ierr) 
     END IF

     IF (Rank_of_process_below.GE.0) THEN
! ## 6 ## send down
        left_bottom_inner_node = global_offset + 1
        IF (Rank_of_process_left.LT.0) left_bottom_inner_node = left_bottom_inner_node + 1  ! skip the node at the left boundary
        ibufer(1) = left_bottom_inner_node
        CALL MPI_SEND(ibufer(1), 1, MPI_INTEGER, Rank_of_process_below, Rank_of_process, MPI_COMM_WORLD, ierr) 
     END IF

     IF (Rank_of_process_below.GE.0) THEN
! ## 7 ## receive from below
        CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, Rank_of_process_below, Rank_of_process_below, MPI_COMM_WORLD, stattus, ierr)
        process_below_left_top_inner_node = ibufer(1)
     END IF
 
     IF (Rank_of_process_above.GE.0) THEN
! ## 8 ## receive from above
        CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, Rank_of_process_above, Rank_of_process_above, MPI_COMM_WORLD, stattus, ierr)
        process_above_left_bottom_inner_node = ibufer(1)
     END IF
 
  ELSE

     IF (Rank_of_process_left.GE.0) THEN
! ## 1 ## receive from left
        CALL MPI_RECV(ibufer, 2, MPI_INTEGER, Rank_of_process_left, Rank_of_process_left, MPI_COMM_WORLD, stattus, ierr)
        process_left_bottom_right_inner_node = ibufer(1)
        process_left_solved_nodes_row_length = ibufer(2)
     END IF
 
     IF (Rank_of_process_right.GE.0) THEN
! ## 2 ## receive from right
        CALL MPI_RECV(ibufer, 2, MPI_INTEGER, Rank_of_process_right, Rank_of_process_right, MPI_COMM_WORLD, stattus, ierr)
        process_right_bottom_left_inner_node = ibufer(1)
        process_right_solved_nodes_row_length = ibufer(2)
     END IF

     IF (Rank_of_process_right.GE.0) THEN
! ## 3 ## send right
        bottom_right_inner_node = global_offset + N_of_nodes_to_solve_x
        IF (Rank_of_process_below.LT.0) bottom_right_inner_node = bottom_right_inner_node + N_of_nodes_to_solve_x  ! skip line of nodes along the bottom boundary
        ibufer(1) = bottom_right_inner_node
        ibufer(2) = N_of_nodes_to_solve_x
        CALL MPI_SEND(ibufer, 2, MPI_INTEGER, Rank_of_process_right, Rank_of_process, MPI_COMM_WORLD, ierr) 
     END IF

     IF (Rank_of_process_left.GE.0) THEN
! ## 4 ## send left
        bottom_left_inner_node = global_offset + 1
        IF (Rank_of_process_below.LT.0) bottom_left_inner_node = bottom_left_inner_node + N_of_nodes_to_solve_x  ! skip line of nodes along the bottom boundary
        ibufer(1) = bottom_left_inner_node
        ibufer(2) = N_of_nodes_to_solve_x
        CALL MPI_SEND(ibufer, 2, MPI_INTEGER, Rank_of_process_left, Rank_of_process, MPI_COMM_WORLD, ierr)
     END IF

     IF (Rank_of_process_below.GE.0) THEN
! ## 5 ## receive from below
        CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, Rank_of_process_below, Rank_of_process_below, MPI_COMM_WORLD, stattus, ierr)
        process_below_left_top_inner_node = ibufer(1)
     END IF
 
     IF (Rank_of_process_above.GE.0) THEN
! ## 6 ## receive from above
        CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, Rank_of_process_above, Rank_of_process_above, MPI_COMM_WORLD, stattus, ierr)
        process_above_left_bottom_inner_node = ibufer(1)
     END IF
 
     IF (Rank_of_process_above.GE.0) THEN
! ## 7 ## send up
        left_top_inner_node = global_offset + N_of_nodes_to_solve_x * (N_of_nodes_to_solve_y-1) + 1
        IF (Rank_of_process_left.LT.0) left_top_inner_node = left_top_inner_node + 1  ! skip the node at the left boundary
        ibufer(1) = left_top_inner_node
        CALL MPI_SEND(ibufer(1), 1, MPI_INTEGER, Rank_of_process_above, Rank_of_process, MPI_COMM_WORLD, ierr) 
     END IF

     IF (Rank_of_process_below.GE.0) THEN
! ## 8 ## send down
        left_bottom_inner_node = global_offset + 1
        IF (Rank_of_process_left.LT.0) left_bottom_inner_node = left_bottom_inner_node + 1  ! skip the node at the left boundary
        ibufer(1) = left_bottom_inner_node
        CALL MPI_SEND(ibufer(1), 1, MPI_INTEGER, Rank_of_process_below, Rank_of_process, MPI_COMM_WORLD, ierr) 
     END IF

  END IF

END SUBROUTINE CALCULATE_BLOCK_OFFSET

!------------------------------------------------
!
SUBROUTINE ANALYZE_BOUNDARY_OBJECTS

  USE CurrentProblemValues

  use mpi

  IMPLICIT NONE

  INTEGER errcode,ierr

  INTEGER nobj, nseg
  LOGICAL connection_found
  INTEGER n, m

  INTEGER connected_object
  INTEGER connected_segment

  INTEGER iv1x, iv1y
  INTEGER iv2x, iv2y

  DO nobj = 1, N_of_boundary_objects
     DO nseg = 1, whole_object(nobj)%number_of_segments

! find where the ends of this segment are connected to

! first process the start point ---------------------------------------------------------------------------------

        connection_found = .FALSE.

        DO n = 1, N_of_boundary_objects
           DO m = 1, whole_object(n)%number_of_segments
! do not compare segment to itself
              IF ((n.EQ.nobj).AND.(m.EQ.nseg)) CYCLE
! try start of the other segment
              IF ( (whole_object(nobj)%segment(nseg)%istart.EQ.whole_object(n)%segment(m)%istart).AND. &
                 & (whole_object(nobj)%segment(nseg)%jstart.EQ.whole_object(n)%segment(m)%jstart) ) THEN
                 connection_found = .TRUE.
                 connected_object = n
                 connected_segment = m
                 EXIT
              END IF
! try end of the other segment
              IF ( (whole_object(nobj)%segment(nseg)%istart.EQ.whole_object(n)%segment(m)%iend).AND. &
                 & (whole_object(nobj)%segment(nseg)%jstart.EQ.whole_object(n)%segment(m)%jend) ) THEN
                 connection_found = .TRUE.
                 connected_object = n
                 connected_segment = m
                 EXIT
              END IF
           END DO
           IF (connection_found) EXIT
        END DO

        IF (.NOT.connection_found) THEN
           PRINT '("ERROR :: cannot find connection for the start point of object ",i4," segment ",i4," with i/j= ",2(2x,i6))', &
                & nobj, nseg, whole_object(nobj)%segment(nseg)%istart, whole_object(nobj)%segment(nseg)%jstart
           errcode=106
           CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
        END IF

        n = connected_object
        m = connected_segment
! vector of object/segment nobj/nseg
        iv1x = whole_object(nobj)%segment(nseg)%iend - whole_object(nobj)%segment(nseg)%istart
        iv1y = whole_object(nobj)%segment(nseg)%jend - whole_object(nobj)%segment(nseg)%jstart
! vector of object/segment n/m
        iv2x = whole_object(n)%segment(m)%iend - whole_object(n)%segment(m)%istart
        iv2y = whole_object(n)%segment(m)%jend - whole_object(n)%segment(m)%jstart

        IF (n.EQ.nobj) THEN
! start end of the segment is connected to another segment of the same object
           IF ((iv1x * iv2y - iv1y * iv2x).EQ.0) THEN
! vectors iv1 and iv2 (i.e. the two segments) are parallel, this should not happen
              PRINT '("ERROR :: boundary object ",i4," has segments ",i4," and ",i4," connected and parallel")', &
                & nobj, nseg, m
              errcode=107
              CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
           ELSE
! vectors iv1 and iv2 (i.e. the two segments) are orthogonal
!### presently, it is assumed that this is a CONCAVE_CORNER, which works in a RECTANGULAR DOMAIN ONLY
!### must be fixed if the whole domain consists of multiple rectangles, so that CONVEX_CORNER 
              whole_object(nobj)%segment(nseg)%start_type = CONCAVE_CORNER
           END IF
        ELSE  !### IF (n.EQ.nobj) THEN
! start end of the segment is connected to a segment of a different object
           IF ((iv1x * iv2y - iv1y * iv2x).EQ.0) THEN
! vectors iv1 and iv2 (i.e. the two segments) are parallel
              whole_object(nobj)%segment(nseg)%start_type = END_FLAT
           ELSE
! vectors iv1 and iv2 (i.e. the two segments) are orthogonal
!### presently, it is assumed that this is an END_CORNER_CONCAVE, which works in a RECTANGULAR DOMAIN ONLY
!### must be fixed if the whole domain consists of multiple rectangles, where END_CORNER_CONVEX is possible
              whole_object(nobj)%segment(nseg)%start_type = END_CORNER_CONCAVE
           END IF
        END IF  !### IF (n.EQ.nobj) THEN

! now process the end point ---------------------------------------------------------------------------------

        connection_found = .FALSE.

        DO n = 1, N_of_boundary_objects
           DO m = 1, whole_object(n)%number_of_segments
! do not compare segment to itself
              IF ((n.EQ.nobj).AND.(m.EQ.nseg)) CYCLE
! try start of the other segment
              IF ( (whole_object(nobj)%segment(nseg)%iend.EQ.whole_object(n)%segment(m)%istart).AND. &
                 & (whole_object(nobj)%segment(nseg)%jend.EQ.whole_object(n)%segment(m)%jstart) ) THEN
                 connection_found = .TRUE.
                 connected_object = n
                 connected_segment = m
                 EXIT
              END IF
! try end of the other segment
              IF ( (whole_object(nobj)%segment(nseg)%iend.EQ.whole_object(n)%segment(m)%iend).AND. &
                 & (whole_object(nobj)%segment(nseg)%jend.EQ.whole_object(n)%segment(m)%jend) ) THEN
                 connection_found = .TRUE.
                 connected_object = n
                 connected_segment = m
                 EXIT
              END IF
           END DO
           IF (connection_found) EXIT
        END DO

        IF (.NOT.connection_found) THEN
           PRINT '("ERROR :: cannot find connection for the end point of object ",i4," segment ",i4," with i/j= ",2(2x,i6))', &
                & nobj, nseg, whole_object(nobj)%segment(nseg)%iend, whole_object(nobj)%segment(nseg)%jend
           errcode=108
           CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
        END IF

        n = connected_object
        m = connected_segment
! vector of object/segment nobj/nseg (already known)
!        iv1x = whole_object(nobj)%segment(nseg)%iend - whole_object(nobj)%segment(nseg)%istart
!        iv1y = whole_object(nobj)%segment(nseg)%jend - whole_object(nobj)%segment(nseg)%jstart
! vector of object/segment n/m
        iv2x = whole_object(n)%segment(m)%iend - whole_object(n)%segment(m)%istart
        iv2y = whole_object(n)%segment(m)%jend - whole_object(n)%segment(m)%jstart

        IF (n.EQ.nobj) THEN
! end of the segment is connected to another segment of the same object
           IF ((iv1x * iv2y - iv1y * iv2x).EQ.0) THEN
! vectors iv1 and iv2 (i.e. the two segments) are parallel, this should not happen
              PRINT '("ERROR :: boundary object ",i4," has segments ",i4," and ",i4," connected and parallel")', &
                & nobj, nseg, m
              errcode=109
              CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
           ELSE
! vectors iv1 and iv2 (i.e. the two segments) are orthogonal
!### presently, it is assumed that this is a CONCAVE_CORNER, which works in a RECTANGULAR DOMAIN ONLY
!### must be fixed if the whole domain consists of multiple rectangles, so that CONVEX_CORNER 
              whole_object(nobj)%segment(nseg)%end_type = CONCAVE_CORNER
           END IF
        ELSE  !### IF (n.EQ.nobj) THEN
! end of the segment is connected to a segment of a different object
           IF ((iv1x * iv2y - iv1y * iv2x).EQ.0) THEN
! vectors iv1 and iv2 (i.e. the two segments) are parallel
              whole_object(nobj)%segment(nseg)%end_type = END_FLAT
           ELSE
! vectors iv1 and iv2 (i.e. the two segments) are orthogonal
!### presently, it is assumed that this is an END_CORNER_CONCAVE, which works in a RECTANGULAR DOMAIN ONLY
!### must be fixed if the whole domain consists of multiple rectangles, where END_CORNER_CONVEX is possible
              whole_object(nobj)%segment(nseg)%end_type = END_CORNER_CONCAVE
           END IF
        END IF  !### IF (n.EQ.nobj) THEN

     END DO   !###  DO nseg = 1, whole_object(n)%number_of_segments
  END DO   !###  DO nobj = 1, N_of_boundary_objects

END SUBROUTINE ANALYZE_BOUNDARY_OBJECTS

