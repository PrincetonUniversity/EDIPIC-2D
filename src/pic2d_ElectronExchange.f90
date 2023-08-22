!-------------------------------------------------------
!
SUBROUTINE EXCHANGE_ELECTRONS_WITH_LEFT_RIGHT_NEIGHBOURS

  USE ParallelOperationValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries, ONLY : c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max

  use mpi

  IMPLICIT NONE

  INTEGER total_part_number_to_receive
  INTEGER n
  INTEGER current_N

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER N_e_to_receive_cluster(1:N_processes_cluster)

  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR

  INTEGER N_e_to_receive
  INTEGER pos, i, k
  INTEGER pos_begin, pos_end

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

!print '("Process ",i4," entered EXCHANGE_ELECTRONS_WITH_LEFT_RIGHT_NEIGHBOURS")', Rank_of_process

! before sending ### RIGHT ###, move particles from upper levels to lower levels if the upper levels have no horizontal neighbors
  IF ((N_processes_cluster_right.GT.0).AND.(N_processes_cluster_right.LT.N_processes_cluster)) THEN
     IF (Rank_cluster.LT.N_processes_cluster_right) THEN
! receive the number of particles from all participating processes to get the total number of particles
        total_part_number_to_receive = 0
        DO n = Rank_cluster + N_processes_cluster_right, N_processes_cluster-1, N_processes_cluster_right
           CALL MPI_RECV(N_e_to_receive_cluster(n), 1, MPI_INTEGER, n, 0, COMM_CLUSTER, stattus, ierr)
           total_part_number_to_receive = total_part_number_to_receive + N_e_to_receive_cluster(n)
        END DO
        current_N = N_e_to_send_right
        N_e_to_send_right = N_e_to_send_right + total_part_number_to_receive   !### horizontal communications use N_e_to_send_right
! if too many particles, have to increase the size of array electron_to_send_right
        IF (N_e_to_send_right.GT.max_N_e_to_send_right) THEN
! save my own particles first
           ALLOCATE(rbufer(1:current_N*6), STAT=ALLOC_ERR)
           pos = 1
           DO i = 1, current_N
              rbufer(pos)   = electron_to_send_right(i)%X
              rbufer(pos+1) = electron_to_send_right(i)%Y
              rbufer(pos+2) = electron_to_send_right(i)%VX
              rbufer(pos+3) = electron_to_send_right(i)%VY
              rbufer(pos+4) = electron_to_send_right(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_right(i)%tag)
              pos = pos+6
           END DO
! re-size and re-allocate array
           IF (ALLOCATED(electron_to_send_right)) DEALLOCATE(electron_to_send_right, STAT = ALLOC_ERR)
           max_N_e_to_send_right = N_e_to_send_right + MAX(50, N_e_to_send_right/10)
           ALLOCATE(electron_to_send_right(1:max_N_e_to_send_right), STAT=ALLOC_ERR)
! restore my own particles
           pos = 1
           DO k = 1, current_N
              electron_to_send_right(k)%X   = rbufer(pos)
              electron_to_send_right(k)%Y   = rbufer(pos+1)
              electron_to_send_right(k)%VX  = rbufer(pos+2) 
              electron_to_send_right(k)%VY  = rbufer(pos+3)
              electron_to_send_right(k)%VZ  = rbufer(pos+4)
              electron_to_send_right(k)%tag = INT(rbufer(pos+5))
              pos = pos+6
           END DO
           DEALLOCATE(rbufer, STAT = ALLOC_ERR)
        END IF
! receive the particles
        pos_end = current_N
        DO n = Rank_cluster + N_processes_cluster_right, N_processes_cluster-1, N_processes_cluster_right
           IF (N_e_to_receive_cluster(n).GT.0) THEN
              ALLOCATE(rbufer(1:N_e_to_receive_cluster(n)*6), STAT=ALLOC_ERR)
              CALL MPI_RECV(rbufer, N_e_to_receive_cluster(n)*6, MPI_DOUBLE_PRECISION, n, n, COMM_CLUSTER, stattus, ierr)               
              pos_begin = pos_end + 1
              pos_end   = pos_end + N_e_to_receive_cluster(n)
              pos = 1
              DO k = pos_begin, pos_end
                 electron_to_send_right(k)%X   = rbufer(pos)
                 electron_to_send_right(k)%Y   = rbufer(pos+1)
                 electron_to_send_right(k)%VX  = rbufer(pos+2) 
                 electron_to_send_right(k)%VY  = rbufer(pos+3)
                 electron_to_send_right(k)%VZ  = rbufer(pos+4)
                 electron_to_send_right(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
              DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           END IF
        END DO
     ELSE
! send the number of particles
        CALL MPI_SEND(N_e_to_send_right, 1, MPI_INTEGER, MOD(Rank_cluster,N_processes_cluster_right), 0, COMM_CLUSTER, ierr) 
        IF (N_e_to_send_right.GT.0) THEN
! prepare array to send
           ALLOCATE(rbufer(1:N_e_to_send_right*6), STAT=ALLOC_ERR)
           pos=1
           DO i = 1, N_e_to_send_right
              rbufer(pos)   = electron_to_send_right(i)%X
              rbufer(pos+1) = electron_to_send_right(i)%Y
              rbufer(pos+2) = electron_to_send_right(i)%VX
              rbufer(pos+3) = electron_to_send_right(i)%VY
              rbufer(pos+4) = electron_to_send_right(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_right(i)%tag)
              pos = pos+6
           END DO
! send particles
           CALL MPI_SEND(rbufer, N_e_to_send_right*6, MPI_DOUBLE_PRECISION, MOD(Rank_cluster,N_processes_cluster_right), Rank_cluster, COMM_CLUSTER, ierr)     
! clear the counter
           N_e_to_send_right = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF
!print '("White process ",i4," done with vertical transport to send right")', Rank_of_process
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! before sending ### LEFT ###, move particles from upper levels to lower levels if the upper levels have no horizontal neighbors
  IF ((N_processes_cluster_left.GT.0).AND.(N_processes_cluster_left.LT.N_processes_cluster)) THEN
     IF (Rank_cluster.LT.N_processes_cluster_left) THEN
! receive the number of particles from all participating processes to get the total number of particles
        total_part_number_to_receive = 0
        DO n = Rank_cluster + N_processes_cluster_left, N_processes_cluster-1, N_processes_cluster_left
           CALL MPI_RECV(N_e_to_receive_cluster(n), 1, MPI_INTEGER, n, 0, COMM_CLUSTER, stattus, ierr)
           total_part_number_to_receive = total_part_number_to_receive + N_e_to_receive_cluster(n)
        END DO
        current_N = N_e_to_send_left
        N_e_to_send_left = N_e_to_send_left + total_part_number_to_receive   !### horizontal communications use N_e_to_send_left
! if too many particles, have to increase the size of array electron_to_send_left
        IF (N_e_to_send_left.GT.max_N_e_to_send_left) THEN
! save my own particles first
           ALLOCATE(rbufer(1:current_N*6), STAT=ALLOC_ERR)
           pos = 1
           DO i = 1, current_N
              rbufer(pos)   = electron_to_send_left(i)%X
              rbufer(pos+1) = electron_to_send_left(i)%Y
              rbufer(pos+2) = electron_to_send_left(i)%VX
              rbufer(pos+3) = electron_to_send_left(i)%VY
              rbufer(pos+4) = electron_to_send_left(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_left(i)%tag)
              pos = pos+6
           END DO
! re-size and re-allocate array
           IF (ALLOCATED(electron_to_send_left)) DEALLOCATE(electron_to_send_left, STAT = ALLOC_ERR)
           max_N_e_to_send_left = N_e_to_send_left + MAX(50, N_e_to_send_left/10)
           ALLOCATE(electron_to_send_left(1:max_N_e_to_send_left), STAT=ALLOC_ERR)
! restore my own particles
           pos = 1
           DO k = 1, current_N
              electron_to_send_left(k)%X   = rbufer(pos)
              electron_to_send_left(k)%Y   = rbufer(pos+1)
              electron_to_send_left(k)%VX  = rbufer(pos+2) 
              electron_to_send_left(k)%VY  = rbufer(pos+3)
              electron_to_send_left(k)%VZ  = rbufer(pos+4)
              electron_to_send_left(k)%tag = INT(rbufer(pos+5))
              pos = pos+6
           END DO
           DEALLOCATE(rbufer, STAT = ALLOC_ERR)
        END IF
! receive the particles
        pos_end = current_N
        DO n = Rank_cluster + N_processes_cluster_left, N_processes_cluster-1, N_processes_cluster_left
           IF (N_e_to_receive_cluster(n).GT.0) THEN
              ALLOCATE(rbufer(1:N_e_to_receive_cluster(n)*6), STAT=ALLOC_ERR)
              CALL MPI_RECV(rbufer, N_e_to_receive_cluster(n)*6, MPI_DOUBLE_PRECISION, n, n, COMM_CLUSTER, stattus, ierr)               
              pos_begin = pos_end + 1
              pos_end   = pos_end + N_e_to_receive_cluster(n)
              pos = 1
              DO k = pos_begin, pos_end
                 electron_to_send_left(k)%X   = rbufer(pos)
                 electron_to_send_left(k)%Y   = rbufer(pos+1)
                 electron_to_send_left(k)%VX  = rbufer(pos+2) 
                 electron_to_send_left(k)%VY  = rbufer(pos+3)
                 electron_to_send_left(k)%VZ  = rbufer(pos+4)
                 electron_to_send_left(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
              DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           END IF
        END DO
     ELSE
! send the number of particles
        CALL MPI_SEND(N_e_to_send_left, 1, MPI_INTEGER, MOD(Rank_cluster,N_processes_cluster_left), 0, COMM_CLUSTER, ierr) 
        IF (N_e_to_send_left.GT.0) THEN
! prepare array to send
           ALLOCATE(rbufer(1:N_e_to_send_left*6), STAT=ALLOC_ERR)
           pos=1
           DO i = 1, N_e_to_send_left
              rbufer(pos)   = electron_to_send_left(i)%X
              rbufer(pos+1) = electron_to_send_left(i)%Y
              rbufer(pos+2) = electron_to_send_left(i)%VX
              rbufer(pos+3) = electron_to_send_left(i)%VY
              rbufer(pos+4) = electron_to_send_left(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_left(i)%tag)
              pos = pos+6
           END DO
! send particles
           CALL MPI_SEND(rbufer, N_e_to_send_left*6, MPI_DOUBLE_PRECISION, MOD(Rank_cluster,N_processes_cluster_left), Rank_cluster, COMM_CLUSTER, ierr)     
! clear the counter
           N_e_to_send_left = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF
!print '("White process ",i4," done with vertical transport to send left")', Rank_of_process
  END IF
 
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  IF (WHITE_CLUSTER) THEN  ! WHITE='TRUE' ############################################################################################
! "white processes"

! ##  1 ## send right the number of particles ----------------------------------
     IF (Rank_horizontal_right.GE.0) THEN
        CALL MPI_SEND(N_e_to_send_right, 1, MPI_INTEGER, Rank_horizontal_right, 0, COMM_HORIZONTAL, ierr) 

! ##  2 ## send right the particles
        IF (N_e_to_send_right.GT.0) THEN
           ALLOCATE(rbufer(1:N_e_to_send_right*6), STAT=ALLOC_ERR)
           pos=1
           DO i = 1, N_e_to_send_right
              rbufer(pos)   = electron_to_send_right(i)%X
              rbufer(pos+1) = electron_to_send_right(i)%Y
              rbufer(pos+2) = electron_to_send_right(i)%VX
              rbufer(pos+3) = electron_to_send_right(i)%VY
              rbufer(pos+4) = electron_to_send_right(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_right(i)%tag)
              pos = pos+6
           END DO
! send particles to right neighbor
           CALL MPI_SEND(rbufer, N_e_to_send_right*6, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, ierr)     
! clear the counter
           N_e_to_send_right = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("White process ",i4," done with sending right")', Rank_of_process

! ##  3 ## send left the number of particles
     IF (Rank_horizontal_left.GE.0) THEN
        CALL MPI_SEND(N_e_to_send_left, 1, MPI_INTEGER, Rank_horizontal_left, 0, COMM_HORIZONTAL, ierr) 

! ##  4 ## send left the particles
        IF (N_e_to_send_left.GT.0) THEN
           ALLOCATE(rbufer(1:N_e_to_send_left*6), STAT=ALLOC_ERR)
           pos=1
           DO i = 1, N_e_to_send_left
              rbufer(pos)   = electron_to_send_left(i)%X
              rbufer(pos+1) = electron_to_send_left(i)%Y
              rbufer(pos+2) = electron_to_send_left(i)%VX
              rbufer(pos+3) = electron_to_send_left(i)%VY
              rbufer(pos+4) = electron_to_send_left(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_left(i)%tag)
              pos = pos+6
           END DO
! send particles to left neighbor
           CALL MPI_SEND(rbufer, N_e_to_send_left*6, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, ierr)     
! clear the counter
           N_e_to_send_left = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("White process ",i4," done with sending left")', Rank_of_process

! ##  5 ## receive from left the number of particles ===========================
     IF (Rank_horizontal_left.GE.0) THEN
        CALL MPI_RECV(N_e_to_receive, 1, MPI_INTEGER, Rank_horizontal_left, 0, COMM_HORIZONTAL, stattus, ierr)

! ##  6 ## receive from left the particles
        IF (N_e_to_receive.GT.0) THEN
! receive particles from left neighbor
           ALLOCATE(rbufer(1:N_e_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, N_e_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
           pos=1
           DO i = 1, N_e_to_receive
              x   = rbufer(pos)
              y   = rbufer(pos+1)
              vx  = rbufer(pos+2)
              vy  = rbufer(pos+3)
              vz  = rbufer(pos+4)
              tag = INT(rbufer(pos+5))
              pos = pos+6
              IF ((y.GE.c_Y_area_min).AND.(y.LE.c_Y_area_max)) THEN
                 CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
              ELSE IF (y.LT.c_Y_area_min) THEN
                 IF (Rank_of_master_below.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_BELOW(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(x, y, vx, vy, vz, tag)
                 END IF
              ELSE IF (y.GT.c_Y_area_max) THEN
                 IF (Rank_of_master_above.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(x, y, vx, vy, vz, tag)
                 END IF
              END IF
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("White process ",i4," done with receiving from left")', Rank_of_process

! ##  7 ## receive from right the number of particles
     IF (Rank_horizontal_right.GE.0) THEN
        CALL MPI_RECV(N_e_to_receive, 1, MPI_INTEGER, Rank_horizontal_right, 0, COMM_HORIZONTAL, stattus, ierr)

! ##  8 ## receive from right the particles
        IF (N_e_to_receive.GT.0) THEN
! receive particles from right neighbor
           ALLOCATE(rbufer(1:N_e_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, N_e_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
           pos = 1
           DO i = 1, N_e_to_receive
              x   = rbufer(pos)
              y   = rbufer(pos+1)
              vx  = rbufer(pos+2)
              vy  = rbufer(pos+3)
              vz  = rbufer(pos+4)
              tag = INT(rbufer(pos+5))
              pos = pos+6
              IF ((y.GE.c_Y_area_min).AND.(y.LE.c_Y_area_max)) THEN
                 CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
              ELSE IF (y.LT.c_Y_area_min) THEN
                 IF (Rank_of_master_below.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_BELOW(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(x, y, vx, vy, vz, tag)
                 END IF
              ELSE IF (y.GT.c_Y_area_max) THEN
                 IF (Rank_of_master_above.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(x, y, vx, vy, vz, tag)
                 END IF
              END IF
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF
     
!print '("White process ",i4," done with receiving from right")', Rank_of_process

  ELSE             ! WHITE_CLUSTER='FALSE' ############################################################################################
! "black" processes

! ##  1 ## receive from left the number of particles ===========================
     IF (Rank_horizontal_left.GE.0) THEN
        CALL MPI_RECV(N_e_to_receive, 1, MPI_INTEGER, Rank_horizontal_left, 0, COMM_HORIZONTAL, stattus, ierr)

! ##  2 ## receive from left the particles
        IF (N_e_to_receive.GT.0) THEN
! receive particles from left neighbor
           ALLOCATE(rbufer(1:N_e_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, N_e_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
           pos = 1
           DO i = 1, N_e_to_receive
              x   = rbufer(pos)
              y   = rbufer(pos+1)
              vx  = rbufer(pos+2)
              vy  = rbufer(pos+3)
              vz  = rbufer(pos+4)
              tag = INT(rbufer(pos+5))
              pos = pos+6
              IF ((y.GE.c_Y_area_min).AND.(y.LE.c_Y_area_max)) THEN
                 CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
              ELSE IF (y.LT.c_Y_area_min) THEN
                 IF (Rank_of_master_below.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_BELOW(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(x, y, vx, vy, vz, tag)
                 END IF
              ELSE IF (y.GT.c_Y_area_max) THEN
                 IF (Rank_of_master_above.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(x, y, vx, vy, vz, tag)
                 END IF
              END IF
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("Black process ",i4," done with receiving from left")', Rank_of_process

! ##  3 ## receive from right the number of particles
     IF (Rank_horizontal_right.GE.0) THEN
        CALL MPI_RECV(N_e_to_receive, 1, MPI_INTEGER, Rank_horizontal_right, 0, COMM_HORIZONTAL, stattus, ierr)

! ##  4 ## receive from right the particles
        IF (N_e_to_receive.GT.0) THEN
! receive particles from right neighbor
           ALLOCATE(rbufer(1:N_e_to_receive*6), STAT=ALLOC_ERR)
!           CALL MPI_PROBE(Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           CALL MPI_RECV(rbufer, N_e_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
           pos=1
           DO i = 1, N_e_to_receive
              x   = rbufer(pos)
              y   = rbufer(pos+1)
              vx  = rbufer(pos+2)
              vy  = rbufer(pos+3)
              vz  = rbufer(pos+4)
              tag = INT(rbufer(pos+5))
              pos = pos+6
              IF ((y.GE.c_Y_area_min).AND.(y.LE.c_Y_area_max)) THEN
                 CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
              ELSE IF (y.LT.c_Y_area_min) THEN
                 IF (Rank_of_master_below.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_BELOW(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(x, y, vx, vy, vz, tag)
                 END IF
              ELSE IF (y.GT.c_Y_area_max) THEN
                 IF (Rank_of_master_above.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(x, y, vx, vy, vz, tag)
                 END IF
              END IF
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("Black process ",i4," done with receiving from right")', Rank_of_process

! ##  5 ## send right the number of particles ----------------------------------
     IF (Rank_horizontal_right.GE.0) THEN
        CALL MPI_SEND(N_e_to_send_right, 1, MPI_INTEGER, Rank_horizontal_right, 0, COMM_HORIZONTAL, ierr) 

! ##  6 ## send right the particles
        IF (N_e_to_send_right.GT.0) THEN
! it is assumed that if aray rbufer is allocated already, it happened when particles from upper levels were
! transferred to lower levels and the array already contains all particles to be sent
! otherwise
! prepare array to send to right neighbor
           IF (.NOT.ALLOCATED(rbufer)) THEN
              ALLOCATE(rbufer(1:N_e_to_send_right*6), STAT=ALLOC_ERR)
              pos=1
              DO i = 1, N_e_to_send_right
                 rbufer(pos)   = electron_to_send_right(i)%X
                 rbufer(pos+1) = electron_to_send_right(i)%Y
                 rbufer(pos+2) = electron_to_send_right(i)%VX
                 rbufer(pos+3) = electron_to_send_right(i)%VY
                 rbufer(pos+4) = electron_to_send_right(i)%VZ
                 rbufer(pos+5) = DBLE(electron_to_send_right(i)%tag)
                 pos = pos+6
              END DO
           END IF
! send particles to right neighbor
           CALL MPI_SEND(rbufer, N_e_to_send_right*6, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, ierr)     
! clear the counter
           N_e_to_send_right = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("Black process ",i4," done with sending right")', Rank_of_process

! ##  7 ## send left the number of particles
     IF (Rank_horizontal_left.GE.0) THEN
        CALL MPI_SEND(N_e_to_send_left, 1, MPI_INTEGER, Rank_horizontal_left, 0, COMM_HORIZONTAL, ierr) 

! ##  8 ## send left the particles
        IF (N_e_to_send_left.GT.0) THEN
! it is assumed that if array rbufer is allocated already, it happened when particles from upper levels were
! transferred to lower levels and the array already contains all particles to be sent
! otherwise
! prepare array to send to left neighbor
           IF (.NOT.ALLOCATED(rbufer)) THEN
              ALLOCATE(rbufer(1:N_e_to_send_left*6), STAT=ALLOC_ERR)
              pos=1
              DO i = 1, N_e_to_send_left
                 rbufer(pos)   = electron_to_send_left(i)%X
                 rbufer(pos+1) = electron_to_send_left(i)%Y
                 rbufer(pos+2) = electron_to_send_left(i)%VX
                 rbufer(pos+3) = electron_to_send_left(i)%VY
                 rbufer(pos+4) = electron_to_send_left(i)%VZ
                 rbufer(pos+5) = DBLE(electron_to_send_left(i)%tag)
                 pos = pos+6
              END DO
           END IF
! send particles to left neighbor
           CALL MPI_SEND(rbufer, N_e_to_send_left*6, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, ierr)     
! clear the counter
           N_e_to_send_left = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("Black process ",i4," done with sending left")', Rank_of_process

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," left EXCHANGE_ELECTRONS_WITH_LEFT_RIGHT_NEIGHBOURS, N_e_to_send l/r/a/b :: ",4(1x,i8))', Rank_of_process, N_e_to_send_left, N_e_to_send_right, N_e_to_send_above, N_e_to_send_below

END SUBROUTINE EXCHANGE_ELECTRONS_WITH_LEFT_RIGHT_NEIGHBOURS

!------------------------------------------------
!
SUBROUTINE EXCHANGE_ELECTRONS_WITH_ABOVE_BELOW_NEIGHBOURS

  USE ParallelOperationValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries, ONLY : c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max

  use mpi

  IMPLICIT NONE

  INTEGER total_part_number_to_receive
  INTEGER n
  INTEGER current_N

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER N_e_to_receive_cluster(1:N_processes_cluster)

  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR

  INTEGER N_e_to_receive
  INTEGER pos, i, k
  INTEGER pos_begin, pos_end

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

!print '("Process ",i4," entered EXCHANGE_ELECTRONS_WITH_NEIGHBOURS")', Rank_of_process
! before sending ### ABOVE ###, move particles from upper levels to lower levels if the upper levels have no horizontal neighbors
  IF ((N_processes_cluster_above.GT.0).AND.(N_processes_cluster_above.LT.N_processes_cluster)) THEN
     IF (Rank_cluster.LT.N_processes_cluster_above) THEN
! receive the number of particles from all participating processes to get the total number of particles
        total_part_number_to_receive = 0
        DO n = Rank_cluster + N_processes_cluster_above, N_processes_cluster-1, N_processes_cluster_above
           CALL MPI_RECV(N_e_to_receive_cluster(n), 1, MPI_INTEGER, n, 0, COMM_CLUSTER, stattus, ierr)
           total_part_number_to_receive = total_part_number_to_receive + N_e_to_receive_cluster(n)
        END DO
        current_N = N_e_to_send_above
        N_e_to_send_above = N_e_to_send_above + total_part_number_to_receive   !### horizontal communications use N_e_to_send_above
! if too many particles, have to increase the size of array electron_to_send_above
        IF (N_e_to_send_above.GT.max_N_e_to_send_above) THEN
! save my own particles first
           ALLOCATE(rbufer(1:current_N*6), STAT=ALLOC_ERR)
           pos = 1
           DO i = 1, current_N
              rbufer(pos)   = electron_to_send_above(i)%X
              rbufer(pos+1) = electron_to_send_above(i)%Y
              rbufer(pos+2) = electron_to_send_above(i)%VX
              rbufer(pos+3) = electron_to_send_above(i)%VY
              rbufer(pos+4) = electron_to_send_above(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_above(i)%tag)
              pos = pos+6
           END DO
! re-size and re-allocate array
           IF (ALLOCATED(electron_to_send_above)) DEALLOCATE(electron_to_send_above, STAT = ALLOC_ERR)
           max_N_e_to_send_above = N_e_to_send_above + MAX(50, N_e_to_send_above/10)
           ALLOCATE(electron_to_send_above(1:max_N_e_to_send_above), STAT=ALLOC_ERR)
! restore my own particles
           pos = 1
           DO k = 1, current_N
              electron_to_send_above(k)%X   = rbufer(pos)
              electron_to_send_above(k)%Y   = rbufer(pos+1)
              electron_to_send_above(k)%VX  = rbufer(pos+2) 
              electron_to_send_above(k)%VY  = rbufer(pos+3)
              electron_to_send_above(k)%VZ  = rbufer(pos+4)
              electron_to_send_above(k)%tag = INT(rbufer(pos+5))
              pos = pos+6
           END DO
           DEALLOCATE(rbufer, STAT = ALLOC_ERR)
        END IF
! receive the particles
        pos_end = current_N
        DO n = Rank_cluster + N_processes_cluster_above, N_processes_cluster-1, N_processes_cluster_above
           IF (N_e_to_receive_cluster(n).GT.0) THEN
              ALLOCATE(rbufer(1:N_e_to_receive_cluster(n)*6), STAT=ALLOC_ERR)
              CALL MPI_RECV(rbufer, N_e_to_receive_cluster(n)*6, MPI_DOUBLE_PRECISION, n, n, COMM_CLUSTER, stattus, ierr)               
              pos_begin = pos_end + 1
              pos_end   = pos_end + N_e_to_receive_cluster(n)
              pos = 1
              DO k = pos_begin, pos_end
                 electron_to_send_above(k)%X   = rbufer(pos)
                 electron_to_send_above(k)%Y   = rbufer(pos+1)
                 electron_to_send_above(k)%VX  = rbufer(pos+2) 
                 electron_to_send_above(k)%VY  = rbufer(pos+3)
                 electron_to_send_above(k)%VZ  = rbufer(pos+4)
                 electron_to_send_above(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
              DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           END IF
        END DO
     ELSE
! send the number of particles
        CALL MPI_SEND(N_e_to_send_above, 1, MPI_INTEGER, MOD(Rank_cluster,N_processes_cluster_above), 0, COMM_CLUSTER, ierr) 
        IF (N_e_to_send_above.GT.0) THEN
! prepare array to send
           ALLOCATE(rbufer(1:N_e_to_send_above*6), STAT=ALLOC_ERR)
           pos=1
           DO i = 1, N_e_to_send_above
              rbufer(pos)   = electron_to_send_above(i)%X
              rbufer(pos+1) = electron_to_send_above(i)%Y
              rbufer(pos+2) = electron_to_send_above(i)%VX
              rbufer(pos+3) = electron_to_send_above(i)%VY
              rbufer(pos+4) = electron_to_send_above(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_above(i)%tag)
              pos = pos+6
           END DO
! send particles
           CALL MPI_SEND(rbufer, N_e_to_send_above*6, MPI_DOUBLE_PRECISION, MOD(Rank_cluster,N_processes_cluster_above), Rank_cluster, COMM_CLUSTER, ierr)     
! clear the counter
           N_e_to_send_above = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF
!print '("White process ",i4," done with vertical transport to send above")', Rank_of_process
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! before sending ### BELOW ###, move particles from upper levels to lower levels if the upper levels have no horizontal neighbors
  IF ((N_processes_cluster_below.GT.0).AND.(N_processes_cluster_below.LT.N_processes_cluster)) THEN
     IF (Rank_cluster.LT.N_processes_cluster_below) THEN
! receive the number of particles from all participating processes to get the total number of particles
        total_part_number_to_receive = 0
        DO n = Rank_cluster + N_processes_cluster_below, N_processes_cluster-1, N_processes_cluster_below
           CALL MPI_RECV(N_e_to_receive_cluster(n), 1, MPI_INTEGER, n, 0, COMM_CLUSTER, stattus, ierr)
           total_part_number_to_receive = total_part_number_to_receive + N_e_to_receive_cluster(n)
        END DO
        current_N = N_e_to_send_below
        N_e_to_send_below = N_e_to_send_below + total_part_number_to_receive   !### horizontal communications use N_e_to_send_below
! if too many particles, have to increase the size of array electron_to_send_below
        IF (N_e_to_send_below.GT.max_N_e_to_send_below) THEN
! save my own particles first
           ALLOCATE(rbufer(1:current_N*6), STAT=ALLOC_ERR)
           pos = 1
           DO i = 1, current_N
              rbufer(pos)   = electron_to_send_below(i)%X
              rbufer(pos+1) = electron_to_send_below(i)%Y
              rbufer(pos+2) = electron_to_send_below(i)%VX
              rbufer(pos+3) = electron_to_send_below(i)%VY
              rbufer(pos+4) = electron_to_send_below(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_below(i)%tag)
              pos = pos+6
           END DO
! re-size and re-allocate array
           IF (ALLOCATED(electron_to_send_below)) DEALLOCATE(electron_to_send_below, STAT = ALLOC_ERR)
           max_N_e_to_send_below = N_e_to_send_below + MAX(50, N_e_to_send_below/10)
           ALLOCATE(electron_to_send_below(1:max_N_e_to_send_below), STAT=ALLOC_ERR)
! restore my own particles
           pos = 1
           DO k = 1, current_N
              electron_to_send_below(k)%X   = rbufer(pos)
              electron_to_send_below(k)%Y   = rbufer(pos+1)
              electron_to_send_below(k)%VX  = rbufer(pos+2) 
              electron_to_send_below(k)%VY  = rbufer(pos+3)
              electron_to_send_below(k)%VZ  = rbufer(pos+4)
              electron_to_send_below(k)%tag = INT(rbufer(pos+5))
              pos = pos+6
           END DO
           DEALLOCATE(rbufer, STAT = ALLOC_ERR)
        END IF
! receive the particles
        pos_end = current_N
        DO n = Rank_cluster + N_processes_cluster_below, N_processes_cluster-1, N_processes_cluster_below
           IF (N_e_to_receive_cluster(n).GT.0) THEN
              ALLOCATE(rbufer(1:N_e_to_receive_cluster(n)*6), STAT=ALLOC_ERR)
              CALL MPI_RECV(rbufer, N_e_to_receive_cluster(n)*6, MPI_DOUBLE_PRECISION, n, n, COMM_CLUSTER, stattus, ierr)               
              pos_begin = pos_end + 1
              pos_end   = pos_end + N_e_to_receive_cluster(n)
              pos = 1
              DO k = pos_begin, pos_end
                 electron_to_send_below(k)%X   = rbufer(pos)
                 electron_to_send_below(k)%Y   = rbufer(pos+1)
                 electron_to_send_below(k)%VX  = rbufer(pos+2) 
                 electron_to_send_below(k)%VY  = rbufer(pos+3)
                 electron_to_send_below(k)%VZ  = rbufer(pos+4)
                 electron_to_send_below(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
              DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           END IF
        END DO
     ELSE
! send the number of particles
        CALL MPI_SEND(N_e_to_send_below, 1, MPI_INTEGER, MOD(Rank_cluster,N_processes_cluster_below), 0, COMM_CLUSTER, ierr) 
        IF (N_e_to_send_below.GT.0) THEN
! prepare array to send
           ALLOCATE(rbufer(1:N_e_to_send_below*6), STAT=ALLOC_ERR)
           pos=1
           DO i = 1, N_e_to_send_below
              rbufer(pos)   = electron_to_send_below(i)%X
              rbufer(pos+1) = electron_to_send_below(i)%Y
              rbufer(pos+2) = electron_to_send_below(i)%VX
              rbufer(pos+3) = electron_to_send_below(i)%VY
              rbufer(pos+4) = electron_to_send_below(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_below(i)%tag)
              pos = pos+6
           END DO
! send particles
           CALL MPI_SEND(rbufer, N_e_to_send_below*6, MPI_DOUBLE_PRECISION, MOD(Rank_cluster,N_processes_cluster_below), Rank_cluster, COMM_CLUSTER, ierr)     
! clear the counter
           N_e_to_send_below = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF
!print '("White process ",i4," done with vertical transport to send below")', Rank_of_process
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
 
  IF (WHITE_CLUSTER) THEN  ! WHITE='TRUE' ############################################################################################
! "white processes"

! ##  9 ## send up the number of particles -------------------------------------
     IF (Rank_horizontal_above.GE.0) THEN
        CALL MPI_SEND(N_e_to_send_above, 1, MPI_INTEGER, Rank_horizontal_above, 0, COMM_HORIZONTAL, ierr) 

! ## 10 ## send up the particles
        IF (N_e_to_send_above.GT.0) THEN
           ALLOCATE(rbufer(1:N_e_to_send_above*6), STAT=ALLOC_ERR)
           pos=1
           DO i = 1, N_e_to_send_above
              rbufer(pos)   = electron_to_send_above(i)%X
              rbufer(pos+1) = electron_to_send_above(i)%Y
              rbufer(pos+2) = electron_to_send_above(i)%VX
              rbufer(pos+3) = electron_to_send_above(i)%VY
              rbufer(pos+4) = electron_to_send_above(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_above(i)%tag)
              pos = pos+6
           END DO
! send particles to neighbor above
           CALL MPI_SEND(rbufer, N_e_to_send_above*6, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, ierr)     
! clear the counter
           N_e_to_send_above = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF
     
!print '("White process ",i4," done with sending up")', Rank_of_process

! ## 11 ## send down the number of particles
     IF (Rank_horizontal_below.GE.0) THEN
        CALL MPI_SEND(N_e_to_send_below, 1, MPI_INTEGER, Rank_horizontal_below, 0, COMM_HORIZONTAL, ierr) 

! ## 12 ## send down the particles
        IF (N_e_to_send_below.GT.0) THEN
           ALLOCATE(rbufer(1:N_e_to_send_below*6), STAT=ALLOC_ERR)
           pos=1
           DO i = 1, N_e_to_send_below
              rbufer(pos)   = electron_to_send_below(i)%X
              rbufer(pos+1) = electron_to_send_below(i)%Y
              rbufer(pos+2) = electron_to_send_below(i)%VX
              rbufer(pos+3) = electron_to_send_below(i)%VY
              rbufer(pos+4) = electron_to_send_below(i)%VZ
              rbufer(pos+5) = DBLE(electron_to_send_below(i)%tag)
              pos = pos+6
           END DO
! send particles to neighbor below
           CALL MPI_SEND(rbufer, N_e_to_send_below*6, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, ierr)     
! clear the counter
           N_e_to_send_below = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("White process ",i4," done with sending down")', Rank_of_process

! ## 13 ## receive from below the number of particles ==========================
     IF (Rank_horizontal_below.GE.0) THEN
        CALL MPI_RECV(N_e_to_receive, 1, MPI_INTEGER, Rank_horizontal_below, 0, COMM_HORIZONTAL, stattus, ierr)

! ## 14 ## receive from below the particles
        IF (N_e_to_receive.GT.0) THEN
! receive particles from neighbor below
           ALLOCATE(rbufer(1:N_e_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, N_e_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
           pos=1
           DO i = 1, N_e_to_receive
              x   = rbufer(pos)
              y   = rbufer(pos+1)
              vx  = rbufer(pos+2)
              vy  = rbufer(pos+3)
              vz  = rbufer(pos+4)
              tag = INT(rbufer(pos+5))
              pos = pos+6
              IF ((x.GE.c_X_area_min).AND.(x.LE.c_X_area_max)) THEN
                 CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
              ELSE IF (x.LT.c_X_area_min) THEN
                 IF (Rank_of_master_left.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_LEFT(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(x, y, vx, vy, vz, tag)
                 END IF
              ELSE IF (x.GT.c_X_area_max) THEN
                 IF (Rank_of_master_right.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_RIGHT(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(x, y, vx, vy, vz, tag)
                 END IF
              END IF
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("White process ",i4," done with receiving from below")', Rank_of_process

! ## 15 ## receive from above the number of particles
     IF (Rank_horizontal_above.GE.0) THEN
        CALL MPI_RECV(N_e_to_receive, 1, MPI_INTEGER, Rank_horizontal_above, 0, COMM_HORIZONTAL, stattus, ierr)

! ## 16 ## receive from above the particles
        IF (N_e_to_receive.GT.0) THEN
!print '("Process ",i4," >> D-4")', Rank_horizontal
! receive particles from neighbor above
           ALLOCATE(rbufer(1:N_e_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, N_e_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)     
!print '("Process ",i4," >> D-5")', Rank_horizontal
! process the received particles
           pos=1
           DO i = 1, N_e_to_receive
              x   = rbufer(pos)
              y   = rbufer(pos+1)
              vx  = rbufer(pos+2)
              vy  = rbufer(pos+3)
              vz  = rbufer(pos+4)
              tag = INT(rbufer(pos+5))
              pos = pos+6
              IF ((x.GE.c_X_area_min).AND.(x.LE.c_X_area_max)) THEN
                 CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
              ELSE IF (x.LT.c_X_area_min) THEN
                 IF (Rank_of_master_left.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_LEFT(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(x, y, vx, vy, vz, tag)
                 END IF
              ELSE IF (x.GT.c_X_area_max) THEN
                 IF (Rank_of_master_right.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_RIGHT(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(x, y, vx, vy, vz, tag)
                 END IF
              END IF
           END DO
!print '("Process ",i4," >> D-6")', Rank_horizontal
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("White process ",i4," done with receiving from above")', Rank_of_process

  ELSE             ! WHITE_CLUSTER='FALSE' ############################################################################################
! "black" processes

! ##  9 ## receive from below the number of particles ==========================
     IF (Rank_horizontal_below.GE.0) THEN
        CALL MPI_RECV(N_e_to_receive, 1, MPI_INTEGER, Rank_horizontal_below, 0, COMM_HORIZONTAL, stattus, ierr)

! ## 10 ## receive from below the particles
        IF (N_e_to_receive.GT.0) THEN
! receive particles from neighbor below
           ALLOCATE(rbufer(1:N_e_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, N_e_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
           pos=1
           DO i = 1, N_e_to_receive
              x   = rbufer(pos)
              y   = rbufer(pos+1)
              vx  = rbufer(pos+2)
              vy  = rbufer(pos+3)
              vz  = rbufer(pos+4)
              tag = INT(rbufer(pos+5))
              pos = pos+6
              IF ((x.GE.c_X_area_min).AND.(x.LE.c_X_area_max)) THEN
                 CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
              ELSE IF (x.LT.c_X_area_min) THEN
                 IF (Rank_of_master_left.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_LEFT(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(x, y, vx, vy, vz, tag)
                 END IF
              ELSE IF (x.GT.c_X_area_max) THEN
                 IF (Rank_of_master_right.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_RIGHT(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(x, y, vx, vy, vz, tag)
                 END IF
              END IF
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("Black process ",i4," done with receiving from below")', Rank_of_process

! ## 11 ## receive from above the  number of particles
     IF (Rank_horizontal_above.GE.0) THEN
        CALL MPI_RECV(N_e_to_receive, 1, MPI_INTEGER, Rank_horizontal_above, 0, COMM_HORIZONTAL, stattus, ierr)

! ## 12 ## receive from above the particles
        IF (N_e_to_receive.GT.0) THEN
!print '("Process ",i4," >> D-4")', Rank_horizontal
! receive particles from neighbor above
           ALLOCATE(rbufer(1:N_e_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, N_e_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)     
!print '("Process ",i4," >> D-5")', Rank_horizontal
! process the received particles
           pos=1
           DO i = 1, N_e_to_receive
              x   = rbufer(pos)
              y   = rbufer(pos+1)
              vx  = rbufer(pos+2)
              vy  = rbufer(pos+3)
              vz  = rbufer(pos+4)
              tag = INT(rbufer(pos+5))
              pos = pos+6
              IF ((x.GE.c_X_area_min).AND.(x.LE.c_X_area_max)) THEN
                 CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
              ELSE IF (x.LT.c_X_area_min) THEN
                 IF (Rank_of_master_left.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_LEFT(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(x, y, vx, vy, vz, tag)
                 END IF
              ELSE IF (x.GT.c_X_area_max) THEN
                 IF (Rank_of_master_right.GE.0) THEN
                    CALL ADD_ELECTRON_TO_SEND_RIGHT(x, y, vx, vy, vz, tag)              
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(x, y, vx, vy, vz, tag)
                 END IF
              END IF
           END DO
!print '("Process ",i4," >> D-6")', Rank_horizontal
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF
       
!print '("Black process ",i4," done with receiving from above")', Rank_of_process

! ## 13 ## send up the number of particles -------------------------------------
     IF (Rank_horizontal_above.GE.0) THEN
        CALL MPI_SEND(N_e_to_send_above, 1, MPI_INTEGER, Rank_horizontal_above, 0, COMM_HORIZONTAL, ierr) 

! ## 14 ## send up the particles
        IF (N_e_to_send_above.GT.0) THEN
! it is assumed that if aray rbufer is allocated already, it happened when particles from upper levels were
! transferred to lower levels and the array already contains all particles to be sent
! otherwise
! prepare array to send to neighbor above
           IF (.NOT.ALLOCATED(rbufer)) THEN
              ALLOCATE(rbufer(1:N_e_to_send_above*6), STAT=ALLOC_ERR)
              pos=1
              DO i = 1, N_e_to_send_above
                 rbufer(pos)   = electron_to_send_above(i)%X
                 rbufer(pos+1) = electron_to_send_above(i)%Y
                 rbufer(pos+2) = electron_to_send_above(i)%VX
                 rbufer(pos+3) = electron_to_send_above(i)%VY
                 rbufer(pos+4) = electron_to_send_above(i)%VZ
                 rbufer(pos+5) = DBLE(electron_to_send_above(i)%tag)
                 pos = pos+6
              END DO
           END IF
! send particles to neighbor above
           CALL MPI_SEND(rbufer, N_e_to_send_above*6, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, ierr)     
! clear the counter
           N_e_to_send_above = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF
            
!print '("Black process ",i4," done with sending up")', Rank_of_process

! ## 15 ## send down the number of particles
     IF (Rank_horizontal_below.GE.0) THEN
        CALL MPI_SEND(N_e_to_send_below, 1, MPI_INTEGER, Rank_horizontal_below, 0, COMM_HORIZONTAL, ierr) 

! ## 16 ## send down the particles
        IF (N_e_to_send_below.GT.0) THEN
! it is assumed that if aray rbufer is allocated already, it happened when particles from upper levels were
! transferred to lower levels and the array already contains all particles to be sent
! otherwise
! prepare array to send to neighbor below
           IF (.NOT.ALLOCATED(rbufer)) THEN
              ALLOCATE(rbufer(1:N_e_to_send_below*6), STAT=ALLOC_ERR)
              pos=1
              DO i = 1, N_e_to_send_below
                 rbufer(pos)   = electron_to_send_below(i)%X
                 rbufer(pos+1) = electron_to_send_below(i)%Y
                 rbufer(pos+2) = electron_to_send_below(i)%VX
                 rbufer(pos+3) = electron_to_send_below(i)%VY
                 rbufer(pos+4) = electron_to_send_below(i)%VZ
                 rbufer(pos+5) = DBLE(electron_to_send_below(i)%tag)
                 pos = pos+6
              END DO
           END IF
! send particles to neighbor below
           CALL MPI_SEND(rbufer, N_e_to_send_below*6, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, ierr)     
! clear the counter
           N_e_to_send_below = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF
            
!print '("Black process ",i4," done with sending down")', Rank_of_process

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," left EXCHANGE_ELECTRONS_WITH_LEFT_RIGHT_NEIGHBOURS, N_e_to_send l/r/a/b :: ",4(1x,i8))', Rank_of_process, N_e_to_send_left, N_e_to_send_right, N_e_to_send_above, N_e_to_send_below

END SUBROUTINE EXCHANGE_ELECTRONS_WITH_ABOVE_BELOW_NEIGHBOURS
