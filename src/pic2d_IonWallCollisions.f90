!-----------------------------------------
!
SUBROUTINE PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : Qs
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER jbelow, jabove
  REAL(8) dqbelow, dqabove

  particle_not_processed = .TRUE.

  y = MIN(MAX(y,c_Y_area_min),c_Y_area_max)

  DO n = 1, c_N_of_local_object_parts_left

     m = c_index_of_local_object_part_left(n)

     IF ( (y.GE.c_local_object_part(m)%jstart).AND. &
        & (y.LE.c_local_object_part(m)%jend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%ion_hit_count(s) = whole_object(nwo)%ion_hit_count(s) + 1

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = (y - jbelow) * Qs(s)
              dqbelow = Qs(s) - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) + dqbelow
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) + dqabove
        END SELECT


!>>>> processing begins >>>>
! the simplest case - specular reflection everywhere
!###        CALL ADD_ELECTRON_TO_ADD_LIST(c_X_area_min, y, -vx, vy, vz, tag)
! here the long min/max operator is necessary to process a particle hitting a corner
!<<<< processing ends <<<<
        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ION_COLL_WITH_BOUNDARY_LEFT")', Rank_of_process
     PRINT '("particle s= ",i4," x= ",e14.7," y= ",e14.7)', s, x, y
     STOP
  END IF

END SUBROUTINE PROCESS_ION_COLL_WITH_BOUNDARY_LEFT

!-----------------------------------------
!
SUBROUTINE PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : Qs
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER jbelow, jabove
  REAL(8) dqbelow, dqabove

  particle_not_processed = .TRUE.

  y = MIN(MAX(y,c_Y_area_min),c_Y_area_max)

  DO n = 1, c_N_of_local_object_parts_right

     m = c_index_of_local_object_part_right(n)

     IF ( (y.GE.c_local_object_part(m)%jstart).AND. &
        & (y.LE.c_local_object_part(m)%jend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%ion_hit_count(s) = whole_object(nwo)%ion_hit_count(s) + 1

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = (y - jbelow) * Qs(s)
              dqbelow = Qs(s) - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) + dqbelow
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) + dqabove
        END SELECT

!>>>> processing begins >>>>
! the simplest case - specular reflection everywhere
!###         CALL ADD_ELECTRON_TO_ADD_LIST(c_X_area_max, y, -vx, vy, vz, tag)
!<<<< processing ends <<<<
        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT")', Rank_of_process
     PRINT '("particle s= ",i4," x= ",e14.7," y= ",e14.7)', s, x, y
     STOP
  END IF

END SUBROUTINE PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT

!-----------------------------------------
!
SUBROUTINE PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : Qs
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER ileft, iright
  REAL(8) dqleft, dqright

  particle_not_processed = .TRUE.

  x = MIN(MAX(x,c_X_area_min),c_X_area_max)

  DO n = 1, c_N_of_local_object_parts_below

     m = c_index_of_local_object_part_below(n)

     IF ( (x.GE.c_local_object_part(m)%istart).AND. &
        & (x.LE.c_local_object_part(m)%iend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%ion_hit_count(s) = whole_object(nwo)%ion_hit_count(s) + 1

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = Qs(s) * (x - ileft)
              dqleft = Qs(s) - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft) + dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) + dqright
        END SELECT

!>>>> processing begins >>>>
! the simplest case - specular reflection everywhere
!###         CALL ADD_ION_TO_ADD_LIST(x, c_Y_area_min, vx, -vy, vz, tag)
!<<<< processing ends <<<<
        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ION_COLL_WITH_BOUNDARY_BELOW")', Rank_of_process
     PRINT '("particle s= ",i4," x= ",e14.7," y= ",e14.7)', s, x, y
     STOP
  END IF

END SUBROUTINE PROCESS_ION_COLL_WITH_BOUNDARY_BELOW

!-----------------------------------------
!
SUBROUTINE PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : Qs
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER ileft, iright
  REAL(8) dqleft, dqright

  particle_not_processed = .TRUE.

  x = MIN(MAX(x,c_X_area_min),c_X_area_max)

  DO n = 1, c_N_of_local_object_parts_above

     m = c_index_of_local_object_part_above(n)

     IF ( (x.GE.c_local_object_part(m)%istart).AND. &
        & (x.LE.c_local_object_part(m)%iend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%ion_hit_count(s) = whole_object(nwo)%ion_hit_count(s) + 1

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = Qs(s) * (x - ileft)
              dqleft = Qs(s) - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft) + dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) + dqright
        END SELECT

!>>>> processing begins >>>>
! the simplest case - specular reflection everywhere
!###        CALL ADD_ION_TO_ADD_LIST(x, c_Y_area_max, vx, -vy, vz, tag)
!<<<< processing ends <<<<
        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE")', Rank_of_process
     PRINT '("particle s= ",i4," x= ",e14.7," y= ",e14.7)', s, x, y
     STOP
  END IF

END SUBROUTINE PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE

!------------------------------------------------------
!
SUBROUTINE COLLECT_PARTICLE_BOUNDARY_HITS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : N_spec

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
!  INTEGER stattus(MPI_STATUS_SIZE)
!  INTEGER request

  INTEGER bufsize
  INTEGER, ALLOCATABLE :: ibuf_send(:)
  INTEGER, ALLOCATABLE :: ibuf_receive(:)
  INTEGER ALLOC_ERR
  
  INTEGER k, pos1, pos2, s

  IF (c_N_of_local_object_parts.GT.0) THEN

     bufsize = N_of_boundary_objects * (1 + N_spec)

     ALLOCATE(ibuf_send(1:bufsize), STAT = ALLOC_ERR)
     ALLOCATE(ibuf_receive(1:bufsize), STAT = ALLOC_ERR)

! each cluster adjacent to a boundary assembles electron-boundary hit counters from all cluster members in the master of the cluster

! include electron counters first
     ibuf_send(1:N_of_boundary_objects) = whole_object(1:N_of_boundary_objects)%electron_hit_count
     pos2 = N_of_boundary_objects
     DO k = 1, N_of_boundary_objects
! include counters for the ion species
        pos1 = pos2 + 1
        pos2 = pos2 + N_spec
        ibuf_send(pos1:pos2) = whole_object(k)%ion_hit_count(1:N_spec)
     END DO
     ibuf_receive = 0

     CALL MPI_REDUCE(ibuf_send, ibuf_receive, bufsize, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)

     IF ((cluster_rank_key.EQ.0).AND.(c_N_of_local_object_parts.GT.0)) THEN
! now counters from all boundary cluster masters are assembled in the process with global rank zero

        ibuf_send = ibuf_receive
        ibuf_receive = 0
! 
        CALL MPI_REDUCE(ibuf_send, ibuf_receive, bufsize, MPI_INTEGER, MPI_SUM, 0, COMM_BOUNDARY, ierr)

        IF (Rank_of_process.EQ.0) THEN
           whole_object(1:N_of_boundary_objects)%electron_hit_count = ibuf_receive(1:N_of_boundary_objects)
           pos2 = N_of_boundary_objects
           DO k = 1, N_of_boundary_objects
              pos1 = pos2 + 1
              pos2 = pos2 + N_spec
              whole_object(k)%ion_hit_count(1:N_spec) = ibuf_receive(pos1:pos2)
           END DO

print '("electrons hit boundaries :: ",10(2x,i8))', whole_object(1:N_of_boundary_objects)%electron_hit_count           
do s = 1, N_spec
   do k = 1, N_of_boundary_objects
      ibuf_send(k) = whole_object(k)%ion_hit_count(s)
   end do
   print '("ions (",i2,") hit boundaries :: ",10(2x,i8))', s, ibuf_send(1:N_of_boundary_objects)
end do

        END IF

     END IF

     DEALLOCATE(ibuf_send, STAT = ALLOC_ERR)
     DEALLOCATE(ibuf_receive, STAT = ALLOC_ERR)

  END IF

END SUBROUTINE COLLECT_PARTICLE_BOUNDARY_HITS

!-----------------------------------------------
!
SUBROUTINE GATHER_SURFACE_CHARGE_DENSITY

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries  

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  LOGICAL dielectric_found
  INTEGER n
  INTEGER nwo

  INTEGER bufsize
  REAL(8), ALLOCATABLE :: rbufer(:), rbufer2(:)
  INTEGER ALLOC_ERR

  INTEGER pos1, pos2

! if no dielectric walls, exit
  dielectric_found = .FALSE.
  DO n = 1, c_N_of_local_object_parts
     nwo = c_local_object_part(n)%object_number
     IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
        dielectric_found = .TRUE.
        EXIT
     END IF
  END DO
  IF (.NOT.dielectric_found) RETURN

! collect contributions from all particle calculators and add them to the density stored in the master process 
  bufsize = 0
  DO n = 1, c_N_of_local_object_parts
     nwo = c_local_object_part(n)%object_number
     IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
        IF (c_local_object_part(n)%jend.GT.c_local_object_part(n)%jstart) THEN
! vertical boundary
           bufsize = bufsize + c_local_object_part(n)%jend - c_local_object_part(n)%jstart + 1
        ELSE IF (c_local_object_part(n)%iend.GT.c_local_object_part(n)%istart) THEN
! horizontal boundary
           bufsize = bufsize + c_local_object_part(n)%iend - c_local_object_part(n)%istart + 1
        ELSE
! error
        END IF
     END IF
  END DO

  ALLOCATE(rbufer(bufsize), STAT = ALLOC_ERR)
  ALLOCATE(rbufer2(bufsize), STAT = ALLOC_ERR)

  pos2 = 0
  DO n = 1, c_N_of_local_object_parts
     nwo = c_local_object_part(n)%object_number
     IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
        IF (c_local_object_part(n)%jend.GT.c_local_object_part(n)%jstart) THEN
! vertical boundary
           pos1 = pos2 + 1
           pos2 = pos2 + c_local_object_part(n)%jend - c_local_object_part(n)%jstart + 1
           rbufer(pos1:pos2) = c_local_object_part(n)%surface_charge(c_local_object_part(n)%jstart:c_local_object_part(n)%jend)
        ELSE IF (c_local_object_part(n)%iend.GT.c_local_object_part(n)%istart) THEN
! horizontal boundary
           pos1 = pos2 + 1
           pos2 = pos2 + c_local_object_part(n)%iend - c_local_object_part(n)%istart + 1
           rbufer(pos1:pos2) = c_local_object_part(n)%surface_charge(c_local_object_part(n)%istart:c_local_object_part(n)%iend)
        ELSE
! error
           PRINT '("Process ",i5," :: ERROR-1 in GATHER_SURFACE_CHARGE_DENSITY :: inconsistent limits of local part ",i3," object ",i3," istart/iend ",2(2x,i6)," jstart/jend ",2(2x,i6))', &
                & Rank_of_process, &
                & n, &
                & c_local_object_part(n)%object_number, &
                & c_local_object_part(n)%istart, &
                & c_local_object_part(n)%iend, &
                & c_local_object_part(n)%jstart, &
                & c_local_object_part(n)%jend
           STOP
        END IF
     END IF
  END DO
  
  CALL MPI_REDUCE(rbufer, rbufer2, bufsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_CLUSTER, ierr)

  IF (cluster_rank_key.EQ.0) THEN
! cluster master saves the collected data in a permanent array
     pos2 = 0
     DO n = 1, c_N_of_local_object_parts
        nwo = c_local_object_part(n)%object_number
        IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
           IF (c_local_object_part(n)%jend.GT.c_local_object_part(n)%jstart) THEN
! vertical boundary
              pos1 = pos2 + 1
              pos2 = pos2 + c_local_object_part(n)%jend - c_local_object_part(n)%jstart + 1
              c_local_object_part(n)%surface_charge(c_local_object_part(n)%jstart:c_local_object_part(n)%jend) = rbufer2(pos1:pos2)
           ELSE IF (c_local_object_part(n)%iend.GT.c_local_object_part(n)%istart) THEN
! horizontal boundary
              pos1 = pos2 + 1
              pos2 = pos2 + c_local_object_part(n)%iend - c_local_object_part(n)%istart + 1
              c_local_object_part(n)%surface_charge(c_local_object_part(n)%istart:c_local_object_part(n)%iend) = rbufer2(pos1:pos2)
           ELSE
! error
              PRINT '("Process ",i5," :: ERROR-2 in GATHER_SURFACE_CHARGE_DENSITY :: inconsistent limits of local part ",i3," object ",i3," istart/iend ",2(2x,i6)," jstart/jend ",2(2x,i6))', &
                   & Rank_of_process, &
                   & n, &
                   & c_local_object_part(n)%object_number, &
                   & c_local_object_part(n)%istart, &
                   & c_local_object_part(n)%iend, &
                   & c_local_object_part(n)%jstart, &
                   & c_local_object_part(n)%jend
              STOP
           END IF
        END IF
     END DO
  ELSE   
! non-master processes clear surface charge density arrays
     DO n = 1, c_N_of_local_object_parts
        nwo = c_local_object_part(n)%object_number
        IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
           c_local_object_part(n)%surface_charge = 0.0_8
        END IF
     END DO
! deallocate arrays and quit
     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
     IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT = ALLOC_ERR)
     RETURN
  END IF

! if the boundary object is larger than the cluster, communicate with neighbor clusters to get their surface charge contributions 
! note that two neighbor clusters may share up to two boundary objects which is why clusters will exchange two values
! the surface charge density values in this pair will be  ordered "left-to-right" for vertical boundaries and exchanges below/above
! and "bottom-to-top" for horizontal boundaries and exchanges left/right
! the situation when two clusters share two boundary objects appears  when, for example, there is only one row of clusters in the system
! if two clusters share only one boundary object, which should be the most common situation, one of the two values is just a dummy

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT = ALLOC_ERR)

  ALLOCATE(rbufer(1:2), STAT = ALLOC_ERR)

  IF (WHITE_CLUSTER) THEN  
! "white processes"

     IF (connect_right) THEN
! ## 1 ## send right charge density in the right end 
        DO n = 1, 2
           rbufer(n) = 0.0_8
           IF (n_right(n).GT.0) THEN
              rbufer(n) = c_local_object_part(n_right(n))%surface_charge(c_indx_x_max)
! clear the density in the endpoint because it will be accumulated in the neighbor cluster
              c_local_object_part(n_right(n))%surface_charge(c_indx_x_max) = 0.0_8
           END IF
        END DO
        CALL MPI_SEND(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_left) THEN
! ## 2 ## send left charge density in the left end
        DO n = 1, 2
           rbufer(n) = 0.0_8
           IF (n_left(n).GT.0) THEN
              rbufer(n) = c_local_object_part(n_left(n))%surface_charge(c_indx_x_min)
              c_local_object_part(n_left(n))%surface_charge(c_indx_x_min) = 0.0_8
           END IF
        END DO
        CALL MPI_SEND(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_left) THEN
! ## 3 ## receive from left charge density in the node next to the left end
        CALL MPI_RECV(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_left, Rank_of_master_left, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
! if a charge density in the endpoint of an object part was sent to the neighbor cluster, 
! that cluster will return its contribution to the same object part (only to a different point of it, of course)
           IF (n_left(n).GT.0) THEN
              c_local_object_part(n_left(n))%surface_charge(c_indx_x_min+1) = c_local_object_part(n_left(n))%surface_charge(c_indx_x_min+1) + rbufer(n)
           END IF
        END DO
     END IF

     IF (connect_right) THEN
! ## 4 ## receive from right charge density in the node next to the right end
        CALL MPI_RECV(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_right, Rank_of_master_right, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_right(n).GT.0) THEN
              c_local_object_part(n_right(n))%surface_charge(c_indx_x_max-1) = c_local_object_part(n_right(n))%surface_charge(c_indx_x_max-1) + rbufer(n)
           END IF
        END DO
     END IF

     IF (connect_above) THEN
! ## 5 ## send above charge density in the top end 
        DO n = 1, 2
           rbufer(n) = 0.0_8
           IF (n_above(n).GT.0) THEN
              rbufer(n) = c_local_object_part(n_above(n))%surface_charge(c_indx_y_max)
              c_local_object_part(n_above(n))%surface_charge(c_indx_y_max) = 0.0_8
           END IF
        END DO
        CALL MPI_SEND(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_below) THEN
! ## 6 ## send below charge density in the bottom end
        DO n = 1, 2
           rbufer(n) = 0.0_8
           IF (n_below(n).GT.0) THEN
              rbufer(n) = c_local_object_part(n_below(n))%surface_charge(c_indx_y_min)
              c_local_object_part(n_below(n))%surface_charge(c_indx_y_min) = 0.0_8
           END IF
        END DO
        CALL MPI_SEND(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_below) THEN
! ## 7 ## receive from below charge density in the node next to the bottom end
        CALL MPI_RECV(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_below, Rank_of_master_below, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_below(n).GT.0) THEN
              c_local_object_part(n_below(n))%surface_charge(c_indx_y_min+1) = c_local_object_part(n_below(n))%surface_charge(c_indx_y_min+1) + rbufer(n)
           END IF
        END DO
     END IF

     IF (connect_above) THEN
! ## 8 ## receive from above charge density in the node next to the top end
        CALL MPI_RECV(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_above, Rank_of_master_above, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_above(n).GT.0) THEN
              c_local_object_part(n_above(n))%surface_charge(c_indx_y_max-1) = c_local_object_part(n_above(n))%surface_charge(c_indx_y_max-1) + rbufer(n)
           END IF
        END DO
     END IF

  ELSE
! "black" processes

     IF (connect_left) THEN
! ## 1 ## receive from left charge density in the node next to the left end
        CALL MPI_RECV(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_left, Rank_of_master_left, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
! if a charge density in the endpoint of an object part was sent to the neighbor cluster, 
! that cluster will return its contribution to the same object part (only to a different point of it, of course)
           IF (n_left(n).GT.0) THEN
              c_local_object_part(n_left(n))%surface_charge(c_indx_x_min+1) = c_local_object_part(n_left(n))%surface_charge(c_indx_x_min+1) + rbufer(n)
           END IF
        END DO
     END IF

     IF (connect_right) THEN
! ## 2 ## receive from right charge density in the node next to the right end
        CALL MPI_RECV(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_right, Rank_of_master_right, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_right(n).GT.0) THEN
              c_local_object_part(n_right(n))%surface_charge(c_indx_x_max-1) = c_local_object_part(n_right(n))%surface_charge(c_indx_x_max-1) + rbufer(n)
           END IF
        END DO
     END IF

     IF (connect_right) THEN
! ## 3 ## send right charge density in the right end 
        DO n = 1, 2
           rbufer(n) = 0.0_8
           IF (n_right(n).GT.0) THEN
              rbufer(n) = c_local_object_part(n_right(n))%surface_charge(c_indx_x_max)
              c_local_object_part(n_right(n))%surface_charge(c_indx_x_max) = 0.0_8
           END IF
        END DO
        CALL MPI_SEND(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_left) THEN
! ## 4 ## send left charge density in the left end
        DO n = 1, 2
           rbufer(n) = 0.0_8
           IF (n_left(n).GT.0) THEN
              rbufer(n) = c_local_object_part(n_left(n))%surface_charge(c_indx_x_min)
              c_local_object_part(n_left(n))%surface_charge(c_indx_x_min) = 0.0_8
           END IF
        END DO
        CALL MPI_SEND(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_below) THEN
! ## 5 ## receive from below charge density in the node next to the bottom end
        CALL MPI_RECV(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_below, Rank_of_master_below, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_below(n).GT.0) THEN
              c_local_object_part(n_below(n))%surface_charge(c_indx_y_min+1) = c_local_object_part(n_below(n))%surface_charge(c_indx_y_min+1) + rbufer(n)
           END IF
        END DO
     END IF

     IF (connect_above) THEN
! ## 6 ## receive from above charge density in the node next to the top end
        CALL MPI_RECV(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_above, Rank_of_master_above, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_above(n).GT.0) THEN
              c_local_object_part(n_above(n))%surface_charge(c_indx_y_max-1) = c_local_object_part(n_above(n))%surface_charge(c_indx_y_max-1) + rbufer(n)
           END IF
        END DO
     END IF

     IF (connect_above) THEN
! ## 7 ## send above charge densitiy in the top end 
        DO n = 1, 2
           rbufer(n) = 0.0_8
           IF (n_above(n).GT.0) THEN
              rbufer(n) = c_local_object_part(n_above(n))%surface_charge(c_indx_y_max)
              c_local_object_part(n_above(n))%surface_charge(c_indx_y_max) = 0.0_8
           END IF
        END DO
        CALL MPI_SEND(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_below) THEN
! ## 8 ## send below charge density in the bottom end
        DO n = 1, 2
           rbufer(n) = 0.0_8
           IF (n_below(n).GT.0) THEN
              rbufer(n) = c_local_object_part(n_below(n))%surface_charge(c_indx_y_min)
              c_local_object_part(n_below(n))%surface_charge(c_indx_y_min) = 0.0_8
           END IF
        END DO
        CALL MPI_SEND(rbufer(1:2), 2, MPI_DOUBLE_PRECISION, Rank_of_master_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

  END IF

! deallocate arrays and quit
  DEALLOCATE(rbufer, STAT = ALLOC_ERR)

END SUBROUTINE GATHER_SURFACE_CHARGE_DENSITY
