!-------------------------------------------------------
!
SUBROUTINE EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS

  USE ParallelOperationValues
  USE ElectronParticles
  USE IonParticles
  USE ClusterAndItsBoundaries, ONLY : c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER total_part_number_to_receive(0:N_spec)
  INTEGER n
  INTEGER current_N

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER N_part_to_receive_cluster(0:N_spec, 1:N_processes_cluster)
  INTEGER sum_N_part_to_receive
  INTEGER sum_N_part_to_send

  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR

  INTEGER ibufer(0:N_spec)
  INTEGER pos, i, k, s
  INTEGER pos_begin(0:N_spec), pos_end(0:N_spec)

  INTEGER N_part_to_receive(0:N_spec)
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

!print '("Process ",i4," entered EXCHANGE_ELECTRONS_WITH_LEFT_RIGHT_NEIGHBOURS")', Rank_of_process

! before sending ### RIGHT ###, move particles from upper levels to lower levels if the upper levels have no horizontal neighbors
  IF ((N_processes_cluster_right.GT.0).AND.(N_processes_cluster_right.LT.N_processes_cluster)) THEN
     IF (Rank_cluster.LT.N_processes_cluster_right) THEN
! receive the number of particles from all participating processes to get the total number of particles
        total_part_number_to_receive(0:N_spec) = 0
        DO n = Rank_cluster + N_processes_cluster_right, N_processes_cluster-1, N_processes_cluster_right
           CALL MPI_RECV(N_part_to_receive_cluster(0:N_spec,n), N_spec+1, MPI_INTEGER, n, 0, COMM_CLUSTER, stattus, ierr)
           total_part_number_to_receive(0:N_spec) = total_part_number_to_receive(0:N_spec) + N_part_to_receive_cluster(0:N_spec, n)
        END DO

! if necessary, resize electron arrays
        current_N = N_e_to_send_right
        N_e_to_send_right = N_e_to_send_right + total_part_number_to_receive(0)   !### horizontal communications use N_e_to_send_right
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
        pos_end(0) = current_N

! if necessary, resize ion arrays
        DO s = 1, N_spec
           current_N = N_ions_to_send_right(s)
           N_ions_to_send_right(s) = N_ions_to_send_right(s) + total_part_number_to_receive(s)   !### horizontal communications use N_ions_to_send_right
! if too many particles, have to increase the size of array ion_to_send_right
           IF (N_ions_to_send_right(s).GT.max_N_ions_to_send_right(s)) THEN
! save my own particles first
              ALLOCATE(rbufer(1:current_N*6), STAT=ALLOC_ERR)
              pos = 1
              DO i = 1, current_N
                 rbufer(pos)   = ion_to_send_right(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_right(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_right(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_right(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_right(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_right(s)%part(i)%tag)
                 pos = pos+6
              END DO
! re-size and re-allocate array
              IF (ALLOCATED(ion_to_send_right(s)%part)) THEN
                 DEALLOCATE(ion_to_send_right(s)%part, STAT = ALLOC_ERR)
                 !NULLIFY(ion_to_send_right(s)%part)
              END IF
              max_N_ions_to_send_right(s) = N_ions_to_send_right(s) + MAX(50, N_ions_to_send_right(s)/10)
              ALLOCATE(ion_to_send_right(s)%part(1:max_N_ions_to_send_right(s)), STAT=ALLOC_ERR)
! restore my own particles
              pos = 1
              DO k = 1, current_N
                 ion_to_send_right(s)%part(k)%X   = rbufer(pos)
                 ion_to_send_right(s)%part(k)%Y   = rbufer(pos+1)
                 ion_to_send_right(s)%part(k)%VX  = rbufer(pos+2) 
                 ion_to_send_right(s)%part(k)%VY  = rbufer(pos+3)
                 ion_to_send_right(s)%part(k)%VZ  = rbufer(pos+4)
                 ion_to_send_right(s)%part(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
              DEALLOCATE(rbufer, STAT = ALLOC_ERR)
           END IF
           pos_end(s) = current_N
        END DO

! receive the particles
        DO n = Rank_cluster + N_processes_cluster_right, N_processes_cluster-1, N_processes_cluster_right

           sum_N_part_to_receive = N_part_to_receive_cluster(0,n)
           DO s = 1, N_spec
              sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive_cluster(s,n)
           END DO
           
           IF (sum_N_part_to_receive.GT.0) THEN
              ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
              CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, n, n, COMM_CLUSTER, stattus, ierr) 
! extract electrons              
              pos_begin(0) = pos_end(0) + 1
              pos_end(0)   = pos_end(0) + N_part_to_receive_cluster(0, n)
              pos = 1
              DO k = pos_begin(0), pos_end(0)
                 electron_to_send_right(k)%X   = rbufer(pos)
                 electron_to_send_right(k)%Y   = rbufer(pos+1)
                 electron_to_send_right(k)%VX  = rbufer(pos+2) 
                 electron_to_send_right(k)%VY  = rbufer(pos+3)
                 electron_to_send_right(k)%VZ  = rbufer(pos+4)
                 electron_to_send_right(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
! extract ions
              DO s = 1, N_spec
                 pos_begin(s) = pos_end(s) + 1
                 pos_end(s)   = pos_end(s) + N_part_to_receive_cluster(s, n)
                 DO k = pos_begin(s), pos_end(s)
                    ion_to_send_right(s)%part(k)%X   = rbufer(pos)
                    ion_to_send_right(s)%part(k)%Y   = rbufer(pos+1)
                    ion_to_send_right(s)%part(k)%VX  = rbufer(pos+2)
                    ion_to_send_right(s)%part(k)%VY  = rbufer(pos+3)
                    ion_to_send_right(s)%part(k)%VZ  = rbufer(pos+4)
                    ion_to_send_right(s)%part(k)%tag = INT(rbufer(pos+5))
                    pos = pos+6
                 END DO
              END DO

              DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           END IF
        END DO
     ELSE
! send the number of particles
        ibufer(0) = N_e_to_send_right
        ibufer(1:N_spec) = N_ions_to_send_right(1:N_spec)
        CALL MPI_SEND(ibufer(0:N_spec), N_spec+1, MPI_INTEGER, MOD(Rank_cluster,N_processes_cluster_right), 0, COMM_CLUSTER, request, ierr) 
        sum_N_part_to_send = N_e_to_send_right
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_right(s)
        END DO
        IF (sum_N_part_to_send.GT.0) THEN
! prepare array to send
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_right(s)
                 rbufer(pos)   = ion_to_send_right(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_right(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_right(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_right(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_right(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_right(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, MOD(Rank_cluster,N_processes_cluster_right), Rank_cluster, COMM_CLUSTER, request, ierr)     
! clear the counters
           N_e_to_send_right = 0
           N_ions_to_send_right = 0
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
        total_part_number_to_receive(0:N_spec) = 0
        DO n = Rank_cluster + N_processes_cluster_left, N_processes_cluster-1, N_processes_cluster_left
           CALL MPI_RECV(N_part_to_receive_cluster(0:N_spec,n), N_spec+1, MPI_INTEGER, n, 0, COMM_CLUSTER, stattus, ierr)
           total_part_number_to_receive(0:N_spec) = total_part_number_to_receive(0:N_spec) + N_part_to_receive_cluster(0:N_spec, n)
        END DO

! if necessary, resize electron arrays
        current_N = N_e_to_send_left
        N_e_to_send_left = N_e_to_send_left + total_part_number_to_receive(0)   !### horizontal communications use N_e_to_send_left
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
        pos_end(0) = current_N

! if necessary, resize ion arrays
        DO s = 1, N_spec
           current_N = N_ions_to_send_left(s)
           N_ions_to_send_left(s) = N_ions_to_send_left(s) + total_part_number_to_receive(s)   !### horizontal communications use N_ions_to_send_left
! if too many particles, have to increase the size of array ion_to_send_left
           IF (N_ions_to_send_left(s).GT.max_N_ions_to_send_left(s)) THEN
! save my own particles first
              ALLOCATE(rbufer(1:current_N*6), STAT=ALLOC_ERR)
              pos = 1
              DO i = 1, current_N
                 rbufer(pos)   = ion_to_send_left(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_left(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_left(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_left(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_left(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_left(s)%part(i)%tag)
                 pos = pos+6
              END DO
! re-size and re-allocate array
              IF (ALLOCATED(ion_to_send_left(s)%part)) THEN
                 DEALLOCATE(ion_to_send_left(s)%part, STAT = ALLOC_ERR)
                 !NULLIFY(ion_to_send_left(s)%part)
              END IF
              max_N_ions_to_send_left(s) = N_ions_to_send_left(s) + MAX(50, N_ions_to_send_left(s)/10)
              ALLOCATE(ion_to_send_left(s)%part(1:max_N_ions_to_send_left(s)), STAT=ALLOC_ERR)
! restore my own particles
              pos = 1
              DO k = 1, current_N
                 ion_to_send_left(s)%part(k)%X   = rbufer(pos)
                 ion_to_send_left(s)%part(k)%Y   = rbufer(pos+1)
                 ion_to_send_left(s)%part(k)%VX  = rbufer(pos+2) 
                 ion_to_send_left(s)%part(k)%VY  = rbufer(pos+3)
                 ion_to_send_left(s)%part(k)%VZ  = rbufer(pos+4)
                 ion_to_send_left(s)%part(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
              DEALLOCATE(rbufer, STAT = ALLOC_ERR)
           END IF
           pos_end(s) = current_N
        END DO

! receive the particles
        DO n = Rank_cluster + N_processes_cluster_left, N_processes_cluster-1, N_processes_cluster_left

           sum_N_part_to_receive = N_part_to_receive_cluster(0,n)
           DO s = 1, N_spec
              sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive_cluster(s,n)
           END DO
           
           IF (sum_N_part_to_receive.GT.0) THEN
              ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
              CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, n, n, COMM_CLUSTER, stattus, ierr) 
! extract electrons              
              pos_begin(0) = pos_end(0) + 1
              pos_end(0)   = pos_end(0) + N_part_to_receive_cluster(0, n)
              pos = 1
              DO k = pos_begin(0), pos_end(0)
                 electron_to_send_left(k)%X   = rbufer(pos)
                 electron_to_send_left(k)%Y   = rbufer(pos+1)
                 electron_to_send_left(k)%VX  = rbufer(pos+2) 
                 electron_to_send_left(k)%VY  = rbufer(pos+3)
                 electron_to_send_left(k)%VZ  = rbufer(pos+4)
                 electron_to_send_left(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
! extract ions
              DO s = 1, N_spec
                 pos_begin(s) = pos_end(s) + 1
                 pos_end(s)   = pos_end(s) + N_part_to_receive_cluster(s, n)
                 DO k = pos_begin(s), pos_end(s)
                    ion_to_send_left(s)%part(k)%X   = rbufer(pos)
                    ion_to_send_left(s)%part(k)%Y   = rbufer(pos+1)
                    ion_to_send_left(s)%part(k)%VX  = rbufer(pos+2)
                    ion_to_send_left(s)%part(k)%VY  = rbufer(pos+3)
                    ion_to_send_left(s)%part(k)%VZ  = rbufer(pos+4)
                    ion_to_send_left(s)%part(k)%tag = INT(rbufer(pos+5))
                    pos = pos+6
                 END DO
              END DO

              DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           END IF
        END DO
     ELSE
! send the number of particles
        ibufer(0) = N_e_to_send_left
        ibufer(1:N_spec) = N_ions_to_send_left(1:N_spec)
        CALL MPI_SEND(ibufer(0:N_spec), N_spec+1, MPI_INTEGER, MOD(Rank_cluster,N_processes_cluster_left), 0, COMM_CLUSTER, request, ierr) 
        sum_N_part_to_send = N_e_to_send_left
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_left(s)
        END DO
        IF (sum_N_part_to_send.GT.0) THEN
! prepare array to send
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_left(s)
                 rbufer(pos)   = ion_to_send_left(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_left(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_left(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_left(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_left(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_left(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, MOD(Rank_cluster,N_processes_cluster_left), Rank_cluster, COMM_CLUSTER, request, ierr)     
! clear the counter
           N_e_to_send_left = 0
           N_ions_to_send_left = 0
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

        ibufer(0) = N_e_to_send_right
        ibufer(1:N_spec) = N_ions_to_send_right(1:N_spec)
        sum_N_part_to_send = N_e_to_send_right
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_right(s)
        END DO
        
        CALL MPI_SEND(ibufer, N_spec+1, MPI_INTEGER, Rank_horizontal_right, 0, COMM_HORIZONTAL, request, ierr) 

! ##  2 ## send right the particles
        IF (sum_N_part_to_send.GT.0) THEN
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_right(s)
                 rbufer(pos)   = ion_to_send_right(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_right(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_right(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_right(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_right(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_right(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles to right neighbor
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr)     
! clear the counters
           N_e_to_send_right = 0
           N_ions_to_send_right = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("White process ",i4," done with sending right")', Rank_of_process

! ##  3 ## send left the number of particles
     IF (Rank_horizontal_left.GE.0) THEN

        ibufer(0) = N_e_to_send_left
        ibufer(1:N_spec) = N_ions_to_send_left(1:N_spec)
        sum_N_part_to_send = N_e_to_send_left
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_left(s)
        END DO
        
        CALL MPI_SEND(ibufer, N_spec+1, MPI_INTEGER, Rank_horizontal_left, 0, COMM_HORIZONTAL, request, ierr) 

! ##  4 ## send left the particles
        IF (sum_N_part_to_send.GT.0) THEN
! prepare array to send
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_left(s)
                 rbufer(pos)   = ion_to_send_left(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_left(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_left(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_left(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_left(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_left(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles to left neighbor
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr)     
! clear the counters
           N_e_to_send_left = 0
           N_ions_to_send_left = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("White process ",i4," done with sending left")', Rank_of_process

! ##  5 ## receive from left the number of particles ===========================
     IF (Rank_horizontal_left.GE.0) THEN
        CALL MPI_RECV(N_part_to_receive(0:N_spec), N_spec+1, MPI_INTEGER, Rank_horizontal_left, 0, COMM_HORIZONTAL, stattus, ierr)
        sum_N_part_to_receive = N_part_to_receive(0)
        DO s = 1, N_spec
           sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive(s)
        END DO

! ##  6 ## receive from left the particles
        IF (sum_N_part_to_receive.GT.0) THEN
! receive particles from left neighbor
           ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
! electrons first
           pos=1
           DO i = 1, N_part_to_receive(0)
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
! then ions
           DO s = 1, N_spec
              DO i = 1, N_part_to_receive(s)
                 x   = rbufer(pos)
                 y   = rbufer(pos+1)
                 vx  = rbufer(pos+2)
                 vy  = rbufer(pos+3)
                 vz  = rbufer(pos+4)
                 tag = INT(rbufer(pos+5))
                 pos = pos+6
                 IF ((y.GE.c_Y_area_min).AND.(y.LE.c_Y_area_max)) THEN
                    CALL ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, tag)
                 ELSE IF (y.LT.c_Y_area_min) THEN
                    IF (Rank_of_master_below.GE.0) THEN
                       CALL ADD_ION_TO_SEND_BELOW(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, x, y, vx, vy, vz, tag)
                    END IF
                 ELSE IF (y.GT.c_Y_area_max) THEN
                    IF (Rank_of_master_above.GE.0) THEN
                       CALL ADD_ION_TO_SEND_ABOVE(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, x, y, vx, vy, vz, tag)
                    END IF
                 END IF
              END DO   
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

!print '("White process ",i4," done with receiving from left")', Rank_of_process

! ##  7 ## receive from right the number of particles
     IF (Rank_horizontal_right.GE.0) THEN
        CALL MPI_RECV(N_part_to_receive, N_spec+1, MPI_INTEGER, Rank_horizontal_right, 0, COMM_HORIZONTAL, stattus, ierr)
        sum_N_part_to_receive = N_part_to_receive(0)
        DO s = 1, N_spec
           sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive(s)
        END DO

! ##  8 ## receive from right the particles
        IF (sum_N_part_to_receive.GT.0) THEN
! receive particles from right neighbor
           ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
! electrons first
           pos = 1
           DO i = 1, N_part_to_receive(0)
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
! then ions
           DO s = 1, N_spec
              DO i = 1, N_part_to_receive(s)
                 x   = rbufer(pos)
                 y   = rbufer(pos+1)
                 vx  = rbufer(pos+2)
                 vy  = rbufer(pos+3)
                 vz  = rbufer(pos+4)
                 tag = INT(rbufer(pos+5))
                 pos = pos+6
                 IF ((y.GE.c_Y_area_min).AND.(y.LE.c_Y_area_max)) THEN
                    CALL ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, tag)
                 ELSE IF (y.LT.c_Y_area_min) THEN
                    IF (Rank_of_master_below.GE.0) THEN
                       CALL ADD_ION_TO_SEND_BELOW(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, x, y, vx, vy, vz, tag)
                    END IF
                 ELSE IF (y.GT.c_Y_area_max) THEN
                    IF (Rank_of_master_above.GE.0) THEN
                       CALL ADD_ION_TO_SEND_ABOVE(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, x, y, vx, vy, vz, tag)
                    END IF
                 END IF
              END DO
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF
     
!print '("White process ",i4," done with receiving from right")', Rank_of_process

  ELSE             ! WHITE_CLUSTER='FALSE' ############################################################################################
! "black" processes

! ##  1 ## receive from left the number of particles ===========================
     IF (Rank_horizontal_left.GE.0) THEN
        CALL MPI_RECV(N_part_to_receive(0:N_spec), N_spec+1, MPI_INTEGER, Rank_horizontal_left, 0, COMM_HORIZONTAL, stattus, ierr)
        sum_N_part_to_receive = N_part_to_receive(0)
        DO s = 1, N_spec
           sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive(s)
        END DO

! ##  2 ## receive from left the particles
        IF (sum_N_part_to_receive.GT.0) THEN
! receive particles from left neighbor
           ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
! electrons first
           pos=1
           DO i = 1, N_part_to_receive(0)
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
! then ions
           DO s = 1, N_spec
              DO i = 1, N_part_to_receive(s)
                 x   = rbufer(pos)
                 y   = rbufer(pos+1)
                 vx  = rbufer(pos+2)
                 vy  = rbufer(pos+3)
                 vz  = rbufer(pos+4)
                 tag = INT(rbufer(pos+5))
                 pos = pos+6
                 IF ((y.GE.c_Y_area_min).AND.(y.LE.c_Y_area_max)) THEN
                    CALL ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, tag)
                 ELSE IF (y.LT.c_Y_area_min) THEN
                    IF (Rank_of_master_below.GE.0) THEN
                       CALL ADD_ION_TO_SEND_BELOW(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, x, y, vx, vy, vz, tag)
                    END IF
                 ELSE IF (y.GT.c_Y_area_max) THEN
                    IF (Rank_of_master_above.GE.0) THEN
                       CALL ADD_ION_TO_SEND_ABOVE(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, x, y, vx, vy, vz, tag)
                    END IF
                 END IF
              END DO   
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

! ##  3 ## receive from right the number of particles
     IF (Rank_horizontal_right.GE.0) THEN
        CALL MPI_RECV(N_part_to_receive, N_spec+1, MPI_INTEGER, Rank_horizontal_right, 0, COMM_HORIZONTAL, stattus, ierr)
        sum_N_part_to_receive = N_part_to_receive(0)
        DO s = 1, N_spec
           sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive(s)
        END DO

! ##  4 ## receive from right the particles
        IF (sum_N_part_to_receive.GT.0) THEN
! receive particles from right neighbor
           ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
! electrons first
           pos = 1
           DO i = 1, N_part_to_receive(0)
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
! then ions
           DO s = 1, N_spec
              DO i = 1, N_part_to_receive(s)
                 x   = rbufer(pos)
                 y   = rbufer(pos+1)
                 vx  = rbufer(pos+2)
                 vy  = rbufer(pos+3)
                 vz  = rbufer(pos+4)
                 tag = INT(rbufer(pos+5))
                 pos = pos+6
                 IF ((y.GE.c_Y_area_min).AND.(y.LE.c_Y_area_max)) THEN
                    CALL ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, tag)
                 ELSE IF (y.LT.c_Y_area_min) THEN
                    IF (Rank_of_master_below.GE.0) THEN
                       CALL ADD_ION_TO_SEND_BELOW(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, x, y, vx, vy, vz, tag)
                    END IF
                 ELSE IF (y.GT.c_Y_area_max) THEN
                    IF (Rank_of_master_above.GE.0) THEN
                       CALL ADD_ION_TO_SEND_ABOVE(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, x, y, vx, vy, vz, tag)
                    END IF
                 END IF
              END DO
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

! ##  5 ## send right the number of particles ----------------------------------
     IF (Rank_horizontal_right.GE.0) THEN

        ibufer(0) = N_e_to_send_right
        ibufer(1:N_spec) = N_ions_to_send_right(1:N_spec)
        sum_N_part_to_send = N_e_to_send_right
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_right(s)
        END DO
        
        CALL MPI_SEND(ibufer, N_spec+1, MPI_INTEGER, Rank_horizontal_right, 0, COMM_HORIZONTAL, request, ierr) 

! ##  6 ## send right the particles
        IF (sum_N_part_to_send.GT.0) THEN
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_right(s)
                 rbufer(pos)   = ion_to_send_right(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_right(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_right(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_right(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_right(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_right(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles to right neighbor
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr)     
! clear the counters
           N_e_to_send_right = 0
           N_ions_to_send_right = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

! ##  7 ## send left the number of particles
     IF (Rank_horizontal_left.GE.0) THEN

        ibufer(0) = N_e_to_send_left
        ibufer(1:N_spec) = N_ions_to_send_left(1:N_spec)
        sum_N_part_to_send = N_e_to_send_left
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_left(s)
        END DO
        
        CALL MPI_SEND(ibufer, N_spec+1, MPI_INTEGER, Rank_horizontal_left, 0, COMM_HORIZONTAL, request, ierr) 

! ##  8 ## send left the particles
        IF (sum_N_part_to_send.GT.0) THEN
! prepare array to send
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_left(s)
                 rbufer(pos)   = ion_to_send_left(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_left(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_left(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_left(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_left(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_left(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles to left neighbor
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr)     
! clear the counters
           N_e_to_send_left = 0
           N_ions_to_send_left = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," left EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS, N_e_to_send l/r/a/b :: ",4(1x,i8))', Rank_of_process, N_e_to_send_left, N_e_to_send_right, N_e_to_send_above, N_e_to_send_below

END SUBROUTINE EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS

!------------------------------------------------
!
SUBROUTINE EXCHANGE_PARTICLES_WITH_ABOVE_BELOW_NEIGHBOURS

  USE ParallelOperationValues
  USE ElectronParticles
  USE IonParticles
  USE ClusterAndItsBoundaries, ONLY : c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER total_part_number_to_receive(0:N_spec)
  INTEGER n
  INTEGER current_N

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER N_part_to_receive_cluster(0:N_spec, 1:N_processes_cluster)
  INTEGER sum_N_part_to_receive
  INTEGER sum_N_part_to_send

  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR

  INTEGER ibufer(0:N_spec)
  INTEGER pos, i, k, s
  INTEGER pos_begin(0:N_spec), pos_end(0:N_spec)

  INTEGER N_part_to_receive(0:N_spec)
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

!print '("Process ",i4," entered EXCHANGE_ELECTRONS_WITH_NEIGHBOURS")', Rank_of_process
! before sending ### ABOVE ###, move particles from upper levels to lower levels if the upper levels have no horizontal neighbors
  IF ((N_processes_cluster_above.GT.0).AND.(N_processes_cluster_above.LT.N_processes_cluster)) THEN
     IF (Rank_cluster.LT.N_processes_cluster_above) THEN
! receive the number of particles from all participating processes to get the total number of particles
        total_part_number_to_receive(0:N_spec) = 0
        DO n = Rank_cluster + N_processes_cluster_above, N_processes_cluster-1, N_processes_cluster_above
           CALL MPI_RECV(N_part_to_receive_cluster(0:N_spec,n), N_spec+1, MPI_INTEGER, n, 0, COMM_CLUSTER, stattus, ierr)
           total_part_number_to_receive(0:N_spec) = total_part_number_to_receive(0:N_spec) + N_part_to_receive_cluster(0:N_spec, n)
        END DO

! if necessary, resize electron arrays
        current_N = N_e_to_send_above
        N_e_to_send_above = N_e_to_send_above + total_part_number_to_receive(0)   !### horizontal communications use N_e_to_send_above
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
        pos_end(0) = current_N

! if necessary, resize ion arrays
        DO s = 1, N_spec
           current_N = N_ions_to_send_above(s)
           N_ions_to_send_above(s) = N_ions_to_send_above(s) + total_part_number_to_receive(s)   !### horizontal communications use N_ions_to_send_above
! if too many particles, have to increase the size of array ion_to_send_above
           IF (N_ions_to_send_above(s).GT.max_N_ions_to_send_above(s)) THEN
! save my own particles first
              ALLOCATE(rbufer(1:current_N*6), STAT=ALLOC_ERR)
              pos = 1
              DO i = 1, current_N
                 rbufer(pos)   = ion_to_send_above(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_above(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_above(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_above(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_above(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_above(s)%part(i)%tag)
                 pos = pos+6
              END DO
! re-size and re-allocate array
              IF (ALLOCATED(ion_to_send_above(s)%part)) THEN
                 DEALLOCATE(ion_to_send_above(s)%part, STAT = ALLOC_ERR)
                 !NULLIFY(ion_to_send_above(s)%part)
              END IF
              max_N_ions_to_send_above(s) = N_ions_to_send_above(s) + MAX(50, N_ions_to_send_above(s)/10)
              ALLOCATE(ion_to_send_above(s)%part(1:max_N_ions_to_send_above(s)), STAT=ALLOC_ERR)
! restore my own particles
              pos = 1
              DO k = 1, current_N
                 ion_to_send_above(s)%part(k)%X   = rbufer(pos)
                 ion_to_send_above(s)%part(k)%Y   = rbufer(pos+1)
                 ion_to_send_above(s)%part(k)%VX  = rbufer(pos+2) 
                 ion_to_send_above(s)%part(k)%VY  = rbufer(pos+3)
                 ion_to_send_above(s)%part(k)%VZ  = rbufer(pos+4)
                 ion_to_send_above(s)%part(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
              DEALLOCATE(rbufer, STAT = ALLOC_ERR)
           END IF
           pos_end(s) = current_N
        END DO

! receive the particles
        DO n = Rank_cluster + N_processes_cluster_above, N_processes_cluster-1, N_processes_cluster_above

           sum_N_part_to_receive = N_part_to_receive_cluster(0,n)
           DO s = 1, N_spec
              sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive_cluster(s,n)
           END DO
           
           IF (sum_N_part_to_receive.GT.0) THEN
              ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
              CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, n, n, COMM_CLUSTER, stattus, ierr) 
! extract electrons              
              pos_begin(0) = pos_end(0) + 1
              pos_end(0)   = pos_end(0) + N_part_to_receive_cluster(0, n)
              pos = 1
              DO k = pos_begin(0), pos_end(0)
                 electron_to_send_above(k)%X   = rbufer(pos)
                 electron_to_send_above(k)%Y   = rbufer(pos+1)
                 electron_to_send_above(k)%VX  = rbufer(pos+2) 
                 electron_to_send_above(k)%VY  = rbufer(pos+3)
                 electron_to_send_above(k)%VZ  = rbufer(pos+4)
                 electron_to_send_above(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
! extract ions
              DO s = 1, N_spec
                 pos_begin(s) = pos_end(s) + 1
                 pos_end(s)   = pos_end(s) + N_part_to_receive_cluster(s, n)
                 DO k = pos_begin(s), pos_end(s)
                    ion_to_send_above(s)%part(k)%X   = rbufer(pos)
                    ion_to_send_above(s)%part(k)%Y   = rbufer(pos+1)
                    ion_to_send_above(s)%part(k)%VX  = rbufer(pos+2)
                    ion_to_send_above(s)%part(k)%VY  = rbufer(pos+3)
                    ion_to_send_above(s)%part(k)%VZ  = rbufer(pos+4)
                    ion_to_send_above(s)%part(k)%tag = INT(rbufer(pos+5))
                    pos = pos+6
                 END DO
              END DO

              DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           END IF
        END DO
     ELSE
! send the number of particles
        ibufer(0) = N_e_to_send_above
        ibufer(1:N_spec) = N_ions_to_send_above(1:N_spec)
        CALL MPI_SEND(ibufer(0:N_spec), N_spec+1, MPI_INTEGER, MOD(Rank_cluster,N_processes_cluster_above), 0, COMM_CLUSTER, request, ierr) 
        sum_N_part_to_send = N_e_to_send_above
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_above(s)
        END DO
        IF (sum_N_part_to_send.GT.0) THEN
! prepare array to send
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_above(s)
                 rbufer(pos)   = ion_to_send_above(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_above(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_above(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_above(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_above(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_above(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, MOD(Rank_cluster,N_processes_cluster_above), Rank_cluster, COMM_CLUSTER, request, ierr)     
! clear the counters
           N_e_to_send_above = 0
           N_ions_to_send_above = 0
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
        total_part_number_to_receive(0:N_spec) = 0
        DO n = Rank_cluster + N_processes_cluster_below, N_processes_cluster-1, N_processes_cluster_below
           CALL MPI_RECV(N_part_to_receive_cluster(0:N_spec,n), N_spec+1, MPI_INTEGER, n, 0, COMM_CLUSTER, stattus, ierr)
           total_part_number_to_receive(0:N_spec) = total_part_number_to_receive(0:N_spec) + N_part_to_receive_cluster(0:N_spec, n)
        END DO

! if necessary, resize electron arrays
        current_N = N_e_to_send_below
        N_e_to_send_below = N_e_to_send_below + total_part_number_to_receive(0)   !### horizontal communications use N_e_to_send_below
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
        pos_end(0) = current_N

! if necessary, resize ion arrays
        DO s = 1, N_spec
           current_N = N_ions_to_send_below(s)
           N_ions_to_send_below(s) = N_ions_to_send_below(s) + total_part_number_to_receive(s)   !### horizontal communications use N_ions_to_send_below
! if too many particles, have to increase the size of array ion_to_send_below
           IF (N_ions_to_send_below(s).GT.max_N_ions_to_send_below(s)) THEN
! save my own particles first
              ALLOCATE(rbufer(1:current_N*6), STAT=ALLOC_ERR)
              pos = 1
              DO i = 1, current_N
                 rbufer(pos)   = ion_to_send_below(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_below(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_below(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_below(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_below(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_below(s)%part(i)%tag)
                 pos = pos+6
              END DO
! re-size and re-allocate array
              IF (ALLOCATED(ion_to_send_below(s)%part)) THEN
                 DEALLOCATE(ion_to_send_below(s)%part, STAT = ALLOC_ERR)
                 !NULLIFY(ion_to_send_below(s)%part)
              END IF
              max_N_ions_to_send_below(s) = N_ions_to_send_below(s) + MAX(50, N_ions_to_send_below(s)/10)
              ALLOCATE(ion_to_send_below(s)%part(1:max_N_ions_to_send_below(s)), STAT=ALLOC_ERR)
! restore my own particles
              pos = 1
              DO k = 1, current_N
                 ion_to_send_below(s)%part(k)%X   = rbufer(pos)
                 ion_to_send_below(s)%part(k)%Y   = rbufer(pos+1)
                 ion_to_send_below(s)%part(k)%VX  = rbufer(pos+2) 
                 ion_to_send_below(s)%part(k)%VY  = rbufer(pos+3)
                 ion_to_send_below(s)%part(k)%VZ  = rbufer(pos+4)
                 ion_to_send_below(s)%part(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
              DEALLOCATE(rbufer, STAT = ALLOC_ERR)
           END IF
           pos_end(s) = current_N
        END DO

! receive the particles
        DO n = Rank_cluster + N_processes_cluster_below, N_processes_cluster-1, N_processes_cluster_below

           sum_N_part_to_receive = N_part_to_receive_cluster(0,n)
           DO s = 1, N_spec
              sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive_cluster(s,n)
           END DO
           
           IF (sum_N_part_to_receive.GT.0) THEN
              ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
              CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, n, n, COMM_CLUSTER, stattus, ierr) 
! extract electrons              
              pos_begin(0) = pos_end(0) + 1
              pos_end(0)   = pos_end(0) + N_part_to_receive_cluster(0, n)
              pos = 1
              DO k = pos_begin(0), pos_end(0)
                 electron_to_send_below(k)%X   = rbufer(pos)
                 electron_to_send_below(k)%Y   = rbufer(pos+1)
                 electron_to_send_below(k)%VX  = rbufer(pos+2) 
                 electron_to_send_below(k)%VY  = rbufer(pos+3)
                 electron_to_send_below(k)%VZ  = rbufer(pos+4)
                 electron_to_send_below(k)%tag = INT(rbufer(pos+5))
                 pos = pos+6
              END DO
! extract ions
              DO s = 1, N_spec
                 pos_begin(s) = pos_end(s) + 1
                 pos_end(s)   = pos_end(s) + N_part_to_receive_cluster(s, n)
                 DO k = pos_begin(s), pos_end(s)
                    ion_to_send_below(s)%part(k)%X   = rbufer(pos)
                    ion_to_send_below(s)%part(k)%Y   = rbufer(pos+1)
                    ion_to_send_below(s)%part(k)%VX  = rbufer(pos+2)
                    ion_to_send_below(s)%part(k)%VY  = rbufer(pos+3)
                    ion_to_send_below(s)%part(k)%VZ  = rbufer(pos+4)
                    ion_to_send_below(s)%part(k)%tag = INT(rbufer(pos+5))
                    pos = pos+6
                 END DO
              END DO

              DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           END IF
        END DO
     ELSE
! send the number of particles
        ibufer(0) = N_e_to_send_below
        ibufer(1:N_spec) = N_ions_to_send_below(1:N_spec)
        CALL MPI_SEND(ibufer(0:N_spec), N_spec+1, MPI_INTEGER, MOD(Rank_cluster,N_processes_cluster_below), 0, COMM_CLUSTER, request, ierr) 
        sum_N_part_to_send = N_e_to_send_below
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_below(s)
        END DO
        IF (sum_N_part_to_send.GT.0) THEN
! prepare array to send
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_below(s)
                 rbufer(pos)   = ion_to_send_below(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_below(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_below(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_below(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_below(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_below(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, MOD(Rank_cluster,N_processes_cluster_below), Rank_cluster, COMM_CLUSTER, request, ierr)     
! clear the counters
           N_e_to_send_below = 0
           N_ions_to_send_below = 0
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

        ibufer(0) = N_e_to_send_above
        ibufer(1:N_spec) = N_ions_to_send_above(1:N_spec)
        sum_N_part_to_send = N_e_to_send_above
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_above(s)
        END DO
        
        CALL MPI_SEND(ibufer, N_spec+1, MPI_INTEGER, Rank_horizontal_above, 0, COMM_HORIZONTAL, request, ierr) 

! ## 10 ## send up the particles
        IF (sum_N_part_to_send.GT.0) THEN
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_above(s)
                 rbufer(pos)   = ion_to_send_above(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_above(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_above(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_above(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_above(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_above(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles to neighbor above
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr)     
! clear the counters
           N_e_to_send_above = 0
           N_ions_to_send_above = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

! ## 11 ## send down the number of particles
     IF (Rank_horizontal_below.GE.0) THEN

        ibufer(0) = N_e_to_send_below
        ibufer(1:N_spec) = N_ions_to_send_below(1:N_spec)
        sum_N_part_to_send = N_e_to_send_below
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_below(s)
        END DO
        
        CALL MPI_SEND(ibufer, N_spec+1, MPI_INTEGER, Rank_horizontal_below, 0, COMM_HORIZONTAL, request, ierr) 

! ## 12 ## send down the particles
        IF (sum_N_part_to_send.GT.0) THEN
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_below(s)
                 rbufer(pos)   = ion_to_send_below(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_below(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_below(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_below(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_below(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_below(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles to neighbor below
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr)     
! clear the counters
           N_e_to_send_below = 0
           N_ions_to_send_below = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

! ## 13 ## receive from below the number of particles ==========================
     IF (Rank_horizontal_below.GE.0) THEN
        CALL MPI_RECV(N_part_to_receive(0:N_spec), N_spec+1, MPI_INTEGER, Rank_horizontal_below, 0, COMM_HORIZONTAL, stattus, ierr)
        sum_N_part_to_receive = N_part_to_receive(0)
        DO s = 1, N_spec
           sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive(s)
        END DO

! ## 14 ## receive from below the particles
        IF (sum_N_part_to_receive.GT.0) THEN
! receive particles from neighbor below
           ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
! electrons first
           pos=1
           DO i = 1, N_part_to_receive(0)
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
! then ions
           DO s = 1, N_spec
              DO i = 1, N_part_to_receive(s)
                 x   = rbufer(pos)
                 y   = rbufer(pos+1)
                 vx  = rbufer(pos+2)
                 vy  = rbufer(pos+3)
                 vz  = rbufer(pos+4)
                 tag = INT(rbufer(pos+5))
                 pos = pos+6
                 IF ((x.GE.c_X_area_min).AND.(x.LE.c_X_area_max)) THEN
                    CALL ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, tag)
                 ELSE IF (x.LT.c_X_area_min) THEN
                    IF (Rank_of_master_left.GE.0) THEN
                       CALL ADD_ION_TO_SEND_LEFT(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, x, y, vx, vy, vz, tag)
                    END IF
                 ELSE IF (x.GT.c_X_area_max) THEN
                    IF (Rank_of_master_right.GE.0) THEN
                       CALL ADD_ION_TO_SEND_RIGHT(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, x, y, vx, vy, vz, tag)
                    END IF
                 END IF
              END DO
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

! ## 15 ## receive from above the number of particles
     IF (Rank_horizontal_above.GE.0) THEN
        CALL MPI_RECV(N_part_to_receive(0:N_spec), N_spec+1, MPI_INTEGER, Rank_horizontal_above, 0, COMM_HORIZONTAL, stattus, ierr)
        sum_N_part_to_receive = N_part_to_receive(0)
        DO s = 1, N_spec
           sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive(s)
        END DO

! ## 16 ## receive from above the particles
        IF (sum_N_part_to_receive.GT.0) THEN
! receive particles from neighbor above
           ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
! electrons first
           pos=1
           DO i = 1, N_part_to_receive(0)
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
! then ions
           DO s = 1, N_spec
              DO i = 1, N_part_to_receive(s)
                 x   = rbufer(pos)
                 y   = rbufer(pos+1)
                 vx  = rbufer(pos+2)
                 vy  = rbufer(pos+3)
                 vz  = rbufer(pos+4)
                 tag = INT(rbufer(pos+5))
                 pos = pos+6
                 IF ((x.GE.c_X_area_min).AND.(x.LE.c_X_area_max)) THEN
                    CALL ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, tag)
                 ELSE IF (x.LT.c_X_area_min) THEN
                    IF (Rank_of_master_left.GE.0) THEN
                       CALL ADD_ION_TO_SEND_LEFT(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, x, y, vx, vy, vz, tag)
                    END IF
                 ELSE IF (x.GT.c_X_area_max) THEN
                    IF (Rank_of_master_right.GE.0) THEN
                       CALL ADD_ION_TO_SEND_RIGHT(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, x, y, vx, vy, vz, tag)
                    END IF
                 END IF
              END DO
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

  ELSE             ! WHITE_CLUSTER='FALSE' ############################################################################################
! "black" processes

! ##  9 ## receive from below the number of particles ==========================
     IF (Rank_horizontal_below.GE.0) THEN
        CALL MPI_RECV(N_part_to_receive(0:N_spec), N_spec+1, MPI_INTEGER, Rank_horizontal_below, 0, COMM_HORIZONTAL, stattus, ierr)
        sum_N_part_to_receive = N_part_to_receive(0)
        DO s = 1, N_spec
           sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive(s)
        END DO

! ## 10 ## receive from below the particles
        IF (sum_N_part_to_receive.GT.0) THEN
! receive particles from neighbor below
           ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
! electrons first
           pos=1
           DO i = 1, N_part_to_receive(0)
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
! then ions
           DO s = 1, N_spec
              DO i = 1, N_part_to_receive(s)
                 x   = rbufer(pos)
                 y   = rbufer(pos+1)
                 vx  = rbufer(pos+2)
                 vy  = rbufer(pos+3)
                 vz  = rbufer(pos+4)
                 tag = INT(rbufer(pos+5))
                 pos = pos+6
                 IF ((x.GE.c_X_area_min).AND.(x.LE.c_X_area_max)) THEN
                    CALL ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, tag)
                 ELSE IF (x.LT.c_X_area_min) THEN
                    IF (Rank_of_master_left.GE.0) THEN
                       CALL ADD_ION_TO_SEND_LEFT(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, x, y, vx, vy, vz, tag)
                    END IF
                 ELSE IF (x.GT.c_X_area_max) THEN
                    IF (Rank_of_master_right.GE.0) THEN
                       CALL ADD_ION_TO_SEND_RIGHT(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, x, y, vx, vy, vz, tag)
                    END IF
                 END IF
              END DO
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

! ## 11 ## receive from above the  number of particles
     IF (Rank_horizontal_above.GE.0) THEN
        CALL MPI_RECV(N_part_to_receive(0:N_spec), N_spec+1, MPI_INTEGER, Rank_horizontal_above, 0, COMM_HORIZONTAL, stattus, ierr)
        sum_N_part_to_receive = N_part_to_receive(0)
        DO s = 1, N_spec
           sum_N_part_to_receive = sum_N_part_to_receive + N_part_to_receive(s)
        END DO

! ## 12 ## receive from above the particles
        IF (sum_N_part_to_receive.GT.0) THEN
! receive particles from neighbor above
           ALLOCATE(rbufer(1:sum_N_part_to_receive*6), STAT=ALLOC_ERR)
           CALL MPI_RECV(rbufer, sum_N_part_to_receive*6, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)     
! process the received particles
! electrons first
           pos=1
           DO i = 1, N_part_to_receive(0)
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
! then ions
           DO s = 1, N_spec
              DO i = 1, N_part_to_receive(s)
                 x   = rbufer(pos)
                 y   = rbufer(pos+1)
                 vx  = rbufer(pos+2)
                 vy  = rbufer(pos+3)
                 vz  = rbufer(pos+4)
                 tag = INT(rbufer(pos+5))
                 pos = pos+6
                 IF ((x.GE.c_X_area_min).AND.(x.LE.c_X_area_max)) THEN
                    CALL ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, tag)
                 ELSE IF (x.LT.c_X_area_min) THEN
                    IF (Rank_of_master_left.GE.0) THEN
                       CALL ADD_ION_TO_SEND_LEFT(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, x, y, vx, vy, vz, tag)
                    END IF
                 ELSE IF (x.GT.c_X_area_max) THEN
                    IF (Rank_of_master_right.GE.0) THEN
                       CALL ADD_ION_TO_SEND_RIGHT(s, x, y, vx, vy, vz, tag)              
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, x, y, vx, vy, vz, tag)
                    END IF
                 END IF
              END DO
           END DO
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

! ## 13 ## send up the number of particles -------------------------------------
     IF (Rank_horizontal_above.GE.0) THEN

        ibufer(0) = N_e_to_send_above
        ibufer(1:N_spec) = N_ions_to_send_above(1:N_spec)
        sum_N_part_to_send = N_e_to_send_above
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_above(s)
        END DO
        
        CALL MPI_SEND(ibufer, N_spec+1, MPI_INTEGER, Rank_horizontal_above, 0, COMM_HORIZONTAL, request, ierr) 

! ## 14 ## send up the particles
        IF (sum_N_part_to_send.GT.0) THEN
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_above(s)
                 rbufer(pos)   = ion_to_send_above(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_above(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_above(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_above(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_above(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_above(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles to neighbor above
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr)     
! clear the counters
           N_e_to_send_above = 0
           N_ions_to_send_above = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

! ## 15 ## send down the number of particles
     IF (Rank_horizontal_below.GE.0) THEN

        ibufer(0) = N_e_to_send_below
        ibufer(1:N_spec) = N_ions_to_send_below(1:N_spec)
        sum_N_part_to_send = N_e_to_send_below
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions_to_send_below(s)
        END DO
        
        CALL MPI_SEND(ibufer, N_spec+1, MPI_INTEGER, Rank_horizontal_below, 0, COMM_HORIZONTAL, request, ierr) 

! ## 16 ## send down the particles
        IF (sum_N_part_to_send.GT.0) THEN
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
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
! place ions
           DO s = 1, N_spec
              DO i = 1, N_ions_to_send_below(s)
                 rbufer(pos)   = ion_to_send_below(s)%part(i)%X
                 rbufer(pos+1) = ion_to_send_below(s)%part(i)%Y
                 rbufer(pos+2) = ion_to_send_below(s)%part(i)%VX
                 rbufer(pos+3) = ion_to_send_below(s)%part(i)%VY
                 rbufer(pos+4) = ion_to_send_below(s)%part(i)%VZ
                 rbufer(pos+5) = DBLE(ion_to_send_below(s)%part(i)%tag)
                 pos = pos+6
              END DO
           END DO
! send particles to neighbor below
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr)     
! clear the counters
           N_e_to_send_below = 0
           N_ions_to_send_below = 0
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END IF

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

END SUBROUTINE EXCHANGE_PARTICLES_WITH_ABOVE_BELOW_NEIGHBOURS
