!-----------------------------------------
!
SUBROUTINE PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : Qs
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

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

        CALL ADD_ION_TO_BO_COLLS_LIST(s, REAL(y), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        IF (whole_object(nwo)%reflects_all_ions) THEN
           CALL INJECT_REFLECTED_ION(s, x, y, vx, vy, vz, tag, whole_object(nwo), m, 1)   ! "1" is for a left wall 
           particle_not_processed = .FALSE.
           EXIT
        END IF

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

        IF (whole_object(nwo)%ion_induced_EE_enabled) THEN
           CALL PROCESS_ION_INDUCED_ELECTRON_EMISSION(s, x, y, vx, vy, vz, tag, whole_object(nwo), m, 1)   ! "1" is for a left wall 
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ION_COLL_WITH_BOUNDARY_LEFT")', Rank_of_process
     PRINT '("particle s= ",i4," x= ",e14.7," y= ",e14.7)', s, x, y
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
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

  INCLUDE 'mpif.h'

  INTEGER ierr

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

        CALL ADD_ION_TO_BO_COLLS_LIST(s, REAL(y), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        IF (whole_object(nwo)%reflects_all_ions) THEN
           CALL INJECT_REFLECTED_ION(s, x, y, vx, vy, vz, tag, whole_object(nwo), m, 3)   ! "3" is for a right wall 
           particle_not_processed = .FALSE.
           EXIT
        END IF

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

        IF (whole_object(nwo)%ion_induced_EE_enabled) THEN
           CALL PROCESS_ION_INDUCED_ELECTRON_EMISSION(s, x, y, vx, vy, vz, tag, whole_object(nwo), m, 3)   ! "3" is for a right wall 
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT")', Rank_of_process
     PRINT '("particle s= ",i4," x= ",e14.7," y= ",e14.7)', s, x, y
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
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

  INCLUDE 'mpif.h'

  INTEGER ierr

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

        CALL ADD_ION_TO_BO_COLLS_LIST(s, REAL(x), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        IF (whole_object(nwo)%reflects_all_ions) THEN
           CALL INJECT_REFLECTED_ION(s, x, y, vx, vy, vz, tag, whole_object(nwo), m, 4)   ! "4" is for a wall below
           particle_not_processed = .FALSE.
           EXIT
        END IF

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

        IF (whole_object(nwo)%ion_induced_EE_enabled) THEN
           CALL PROCESS_ION_INDUCED_ELECTRON_EMISSION(s, x, y, vx, vy, vz, tag, whole_object(nwo), m, 4)   ! "4" is for a wall below
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ION_COLL_WITH_BOUNDARY_BELOW")', Rank_of_process
     PRINT '("particle s= ",i4," x= ",e14.7," y= ",e14.7)', s, x, y
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
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

  INCLUDE 'mpif.h'

  INTEGER ierr

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

        CALL ADD_ION_TO_BO_COLLS_LIST(s, REAL(x), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        IF (whole_object(nwo)%reflects_all_ions) THEN
           CALL INJECT_REFLECTED_ION(s, x, y, vx, vy, vz, tag, whole_object(nwo), m, 2)   ! "2" is for a wall above
           particle_not_processed = .FALSE.
           EXIT
        END IF

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

        IF (whole_object(nwo)%ion_induced_EE_enabled) THEN
           CALL PROCESS_ION_INDUCED_ELECTRON_EMISSION(s, x, y, vx, vy, vz, tag, whole_object(nwo), m, 2)   ! "2" is for a wall above
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE")', Rank_of_process
     PRINT '("particle s= ",i4," x= ",e14.7," y= ",e14.7)', s, x, y
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
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

  bufsize = N_of_boundary_and_inner_objects * (1 + N_spec)

  ALLOCATE(ibuf_send(1:bufsize), STAT = ALLOC_ERR)
  ALLOCATE(ibuf_receive(1:bufsize), STAT = ALLOC_ERR)

! include electron counters first
  ibuf_send(1:N_of_boundary_and_inner_objects) = whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count
  pos2 = N_of_boundary_and_inner_objects
  DO k = 1, N_of_boundary_and_inner_objects
! include counters for the ion species
     pos1 = pos2 + 1
     pos2 = pos2 + N_spec
     ibuf_send(pos1:pos2) = whole_object(k)%ion_hit_count(1:N_spec)
  END DO
  ibuf_receive = 0

  CALL MPI_REDUCE(ibuf_send, ibuf_receive, bufsize, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr) !??? use Rank_of_bottom_left_cluster_master ???

  IF (Rank_of_process.EQ.0) THEN
! now counters from all processes are assembled in the process with global rank zero
    
     whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count = ibuf_receive(1:N_of_boundary_and_inner_objects)
     pos2 = N_of_boundary_and_inner_objects
     DO k = 1, N_of_boundary_and_inner_objects
        pos1 = pos2 + 1
        pos2 = pos2 + N_spec
        whole_object(k)%ion_hit_count(1:N_spec) = ibuf_receive(pos1:pos2)
     END DO
     
     print '("electrons hit boundaries :: ",20(2x,i8))', whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count           
     do s = 1, N_spec
        do k = 1, N_of_boundary_and_inner_objects
           ibuf_send(k) = whole_object(k)%ion_hit_count(s)
        end do
        print '("ions (",i2,") hit boundaries :: ",20(2x,i8))', s, ibuf_send(1:N_of_boundary_and_inner_objects)
     end do

  END IF

     DEALLOCATE(ibuf_send, STAT = ALLOC_ERR)
     DEALLOCATE(ibuf_receive, STAT = ALLOC_ERR)

END SUBROUTINE COLLECT_PARTICLE_BOUNDARY_HITS

!-----------------------------------------
!
SUBROUTINE TRY_ION_COLL_WITH_INNER_OBJECT(s, x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : Qs
  USE CurrentProblemValues !, ONLY : inner_object, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

!  INTEGER nio  ! number of the inner object
  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag
!  TYPE(boundary_object) myobject

  REAL(8) xorg, yorg

  INTEGER n_do     ! number of inner object that the particle collided with
  INTEGER mcross   ! number of the segment of inner object number n_do that the particle collided with
  REAL(8) xcross, ycross  ! coordinates of the crossing
  REAL(8) distorg         ! distance from the origin to the crossing (we keep the crossing with the smallest distance from the origin)

  INTEGER n_try
  INTEGER mcross_try
  REAL(8) xcross_try, ycross_try, distorg_try

  INTEGER coll_direction_flag

  REAL coll_coord   ! coordinate of collision point, y/x for collisions with vertical/horizontal segments, respectively

  xorg = x - vx*N_subcycles
  yorg = y - vy*N_subcycles

  n_do = -1
  mcross = -1

  DO n_try = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     CALL FIND_CLOSEST_INTERSECTION_WITH_OBJECT(xorg, yorg, x, y, n_try, mcross_try, xcross_try, ycross_try, distorg_try)

     IF (mcross_try.LT.0) CYCLE  ! no crossing found

     IF (mcross.EQ.-1) THEN
! the very first crossing was found
        n_do = n_try
        mcross = mcross_try
        xcross = xcross_try
        ycross = ycross_try
        distorg = distorg_try
     ELSE
        IF (distorg_try.GE.distorg) CYCLE
! the new crossing is closer to the origin than the previously found one
        n_do = n_try
        mcross = mcross_try
        xcross = xcross_try
        ycross = ycross_try
        distorg = distorg_try
     END IF

  END DO   !### DO n_try = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
   
  IF (mcross.EQ.-1) THEN
     PRINT '("Error-1 in TRY_ION_COLL_WITH_INNER_OBJECT")'
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  SELECT CASE (mcross)
     CASE (1)
        coll_direction_flag = 3
        coll_coord = REAL(ycross)
     CASE (2)
        coll_direction_flag = 4
        coll_coord = REAL(xcross)
     CASE (3)
        coll_direction_flag = 1
        coll_coord = REAL(ycross)
     CASE (4)
        coll_direction_flag = 2
        coll_coord = REAL(xcross)
  END SELECT

  CALL ADD_ION_TO_BO_COLLS_LIST(s, coll_coord, REAL(vx), REAL(vy), REAL(vz), tag, n_do, mcross)

  CALL DO_ION_COLL_WITH_INNER_OBJECT(s, xcross, ycross, vx, vy, vz, tag, whole_object(n_do), coll_direction_flag)

END SUBROUTINE TRY_ION_COLL_WITH_INNER_OBJECT

!-----------------------------------------
!
SUBROUTINE DO_ION_COLL_WITH_INNER_OBJECT(s, x, y, vx, vy, vz, tag, myobject, coll_direction_flag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : Qs
  USE CurrentProblemValues !, ONLY : inner_object, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

!  INTEGER nio  ! number of the inner object
  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag
  TYPE(boundary_object) myobject
  INTEGER coll_direction_flag

  REAL(8) xmin, xmax, ymin, ymax

  INTEGER i_left_top, i_right_top, i_right_bottom, i_left_bottom_bis, i, ip1

  REAL(8) dqip1, dqi

! identify side of the inner object hit by the particle

  xmin = myobject%xmin
  xmax = myobject%xmax
  ymin = myobject%ymin
  ymax = myobject%ymax

  myobject%ion_hit_count(s) = myobject%ion_hit_count(s) + 1

  IF (myobject%reflects_all_ions) THEN
     CALL INJECT_REFLECTED_ION(s, x, y, vx, vy, vz, tag, myobject, -1, coll_direction_flag)
     RETURN
  END IF

  IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge

! ilt   -  ---- irt
!  |             |
!  |             |
!  1 ilbb ---- irb
!
     i_left_top        =                  myobject%jtop   - myobject%jbottom + 1
     i_right_top       = i_left_top     + myobject%iright - myobject%ileft
     i_right_bottom    = i_right_top    + myobject%jtop   - myobject%jbottom
     i_left_bottom_bis = i_right_bottom + myobject%iright - myobject%ileft - 1

     SELECT CASE (coll_direction_flag)
        CASE (3)
! left wall of the object
           i = MIN(INT(y - ymin) + 1, i_left_top - 1)
           dqip1 = (y - INT(y)) * Qs(s)
           dqi = Qs(s) - dqip1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1

        CASE (4)
! top wall of the object
           i = MIN(INT(x - xmin) + i_left_top, i_right_top - 1)
           dqip1 = (x - INT(x)) * Qs(s)
           dqi = Qs(s) - dqip1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1

        CASE (1)
! right wall of the object
           i = MIN(INT(ymax - y) + i_right_top, i_right_bottom - 1)
           dqi = (y - INT(y)) * Qs(s)
           dqip1 = Qs(s) - dqi
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1

        CASE (2)
! bottom wall of the object
           i = MIN(INT(xmax - x) + i_right_bottom, i_left_bottom_bis)
           dqi = (x - INT(x)) * Qs(s)
           dqip1 = Qs(s) - dqi
           ip1 = i+1
           IF (i.EQ.i_left_bottom_bis) ip1 = 1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
           myobject%surface_charge_variation(ip1) = myobject%surface_charge_variation(ip1) + dqip1

     END SELECT
  END IF   !### IF (myobject%object_type.EQ.DIELECTRIC) THEN

  IF (myobject%ion_induced_EE_enabled) THEN
     CALL PROCESS_ION_INDUCED_ELECTRON_EMISSION(s, x, y, vx, vy, vz, tag, myobject, -1, coll_direction_flag)
  END IF

END SUBROUTINE DO_ION_COLL_WITH_INNER_OBJECT

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
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
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
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
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

!-----------------------------------------------
!
SUBROUTINE GATHER_SURFACE_CHARGE_DENSITY_INNER_OBJECTS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries  

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER bufsize
  INTEGER n, k

  REAL(8), ALLOCATABLE :: rbufer(:), rbufer2(:)
  INTEGER ALLOC_ERR

  INTEGER pos1, pos2

  INTEGER nofcopy

  INTEGER i_left_top
  INTEGER i_right_top
  INTEGER i_right_bottom
  INTEGER i_left_bottom_bis

  INTEGER i0bot, i0top

  INTEGER k1, k2

  IF (N_of_inner_objects.EQ.0) RETURN

! define size of the buffer containing charge density along edges of all dielectric inner objects
  bufsize = 0
  DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
     IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE 
     bufsize = bufsize + whole_object(n)%N_boundary_nodes
  END DO
  
! if no dielectric inner objects, exit
  IF (bufsize.EQ.0) RETURN

  ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)
  ALLOCATE(rbufer2(1:bufsize), STAT = ALLOC_ERR)

  pos2=0
  DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
     IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE
     pos1 = pos2+1
     pos2 = pos2+whole_object(n)%N_boundary_nodes
     rbufer(pos1:pos2) = whole_object(n)%surface_charge_variation(1:whole_object(n)%N_boundary_nodes)
  END DO

! add everything together and put it into process master of the bottom left cluster
  CALL MPI_REDUCE(rbufer, rbufer2, bufsize, MPI_DOUBLE_PRECISION, MPI_SUM, Rank_of_bottom_left_cluster_master, MPI_COMM_WORLD, ierr)

  IF (cluster_rank_key.EQ.0) THEN
! note that process with global rank Rank_of_left_bottom_cluster_master
! must be a cluster master and have rank zero in COMM_HORIZONTAL
     CALL MPI_BCAST(rbufer2, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_HORIZONTAL, ierr)

! translate the message / update surface charge
     pos2=0
     DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
        IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE
        pos1 = pos2+1
        pos2 = pos2+whole_object(n)%N_boundary_nodes
        whole_object(n)%surface_charge_variation(1:whole_object(n)%N_boundary_nodes) = rbufer2(pos1:pos2)
     END DO

! a dielectric object crossing one periodic X-boundary has its copy crossing the opposite periodic X-boundary
     DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
        IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE
        IF (whole_object(n)%object_copy_periodic_X_right.GT.0) THEN
! there is a periodic copy of this object placed across the opposite periodic boundary along the X direction (on the right)
           nofcopy = whole_object(n)%object_copy_periodic_X_right
           DO k = 1, whole_object(n)%N_boundary_nodes
              whole_object(n)%surface_charge_variation(k) = whole_object(n)%surface_charge_variation(k) + whole_object(nofcopy)%surface_charge_variation(k)
              whole_object(nofcopy)%surface_charge_variation(k) = whole_object(n)%surface_charge_variation(k)
           END DO
        END IF
     END DO

! dielectric object crossing one periodic Y-boundary has its copy crossing the opposite periodic Y-boundary
     DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
        IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE
        IF (whole_object(n)%object_copy_periodic_Y_above.GT.0) THEN
! there is a periodic copy of this object placed across the opposite periodic boundary along the Y direction (above)
           nofcopy = whole_object(n)%object_copy_periodic_Y_above
           DO k = 1, whole_object(n)%N_boundary_nodes
              whole_object(n)%surface_charge_variation(k) = whole_object(n)%surface_charge_variation(k) + whole_object(nofcopy)%surface_charge_variation(k)
              whole_object(nofcopy)%surface_charge_variation(k) = whole_object(n)%surface_charge_variation(k)
           END DO
        END IF
     END DO

! in a dielectric object crossing symmetry plane x=0, variation of space charge at nodes x=0 must be doubled
! and possible contribution from nodes x<0 must be accounted for
     DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

        IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE
        IF (whole_object(n)%object_does_NOT_cross_symmetry_plane_X) CYCLE

! ilt -------i0top--- irt    i0top/bot is at x=0
!  |                   |
!  |                   |
!  1 ilbb ---i0bot--- irb
!
        i_left_top        =                  whole_object(n)%jtop   - whole_object(n)%jbottom + 1
        i_right_top       = i_left_top     + whole_object(n)%iright - whole_object(n)%ileft
        i_right_bottom    = i_right_top    + whole_object(n)%jtop   - whole_object(n)%jbottom
        i_left_bottom_bis = i_right_bottom + whole_object(n)%iright - whole_object(n)%ileft - 1

! it is expected that such an object has ileft = -iright

        i0top = (i_right_top + i_left_top)/2
        IF ((2*i0top).NE.(i_right_top + i_left_top)) THEN
! this should not happen
           PRINT '("Error-A in GATHER_SURFACE_CHARGE_DENSITY_INNER_OBJECTS, object ",i3," is not symmetric relative x=0")', n
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        END IF

        i0bot = (i_left_bottom_bis + 1 + i_right_bottom)/2
        IF ((2*i0bot).NE.(i_left_bottom_bis + 1 + i_right_bottom)) THEN
! this should not happen
           PRINT '("Error-B in GATHER_SURFACE_CHARGE_DENSITY_INNER_OBJECTS, object ",i3," is not symmetric relative x=0")', n
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        END IF

! add the left wall of the object to the right wall
        DO k1 = 1, i_left_top
           k2 = i_right_bottom - k1 + 1  ! when k1=1 , k2=i_right_bottom ; when k1=i_left_top, k2=i_right_top
           whole_object(n)%surface_charge_variation(k2) = whole_object(n)%surface_charge_variation(k2) + whole_object(n)%surface_charge_variation(k1)
        END DO

! add the left half of the top wall to the right half of the top wall, this includes the upper node at the symmetry axis
        DO k1 = i_left_top+1, i0top
           k2 = i_right_top-k1+i_left_top  ! when k1=i_left_top+1, k2=i_right_top-1; when k1=i0top, k2=i0top
           whole_object(n)%surface_charge_variation(k2) = whole_object(n)%surface_charge_variation(k2) + whole_object(n)%surface_charge_variation(k1)
        END DO

! add the left half of the bottom wall to the right half of the bottom wall, this includes the lower node at the symmetry axis
        DO k1 = i0bot, i_left_bottom_bis
           k2 = i0bot + i0bot - k1        ! when k1 = i0bot, k2 = i0bot ; when k1 = i_left_bottom_bis, k2 = i_right_bottom+1
           whole_object(n)%surface_charge_variation(k2) = whole_object(n)%surface_charge_variation(k2) + whole_object(n)%surface_charge_variation(k1)
        END DO

! now make the perturbation zero for x<0

! left wall of the object
        DO k1 = 1, i_left_top
           whole_object(n)%surface_charge_variation(k1) = 0.0_8
        END DO
! left half of the top wall of the object
        DO k1 = i_left_top+1, i0top-1
           whole_object(n)%surface_charge_variation(k1) = 0.0_8
        END DO
! left half of the bottom wall of the object
        DO k1 = i0bot+1, i_left_bottom_bis
           whole_object(n)%surface_charge_variation(k1) = 0.0_8
        END DO

     END DO   !### DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
        IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE
        DO k = 1, whole_object(n)%N_boundary_nodes
           whole_object(n)%surface_charge(k) = whole_object(n)%surface_charge(k) +  whole_object(n)%surface_charge_variation(k)
        END DO
     END DO

  END IF   !### IF (cluster_rank_key.EQ.0) THEN

! clear surface charge variation arrays
  DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
     IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE
     whole_object(n)%surface_charge_variation = 0.0_8
  END DO

! cleanup
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)
  RETURN

END SUBROUTINE GATHER_SURFACE_CHARGE_DENSITY_INNER_OBJECTS

!---------------------------------------------------------
!
SUBROUTINE ADD_ION_TO_BO_COLLS_LIST(s, coll_coord, vx, vy, vz, tag, nwo, nseg)

  USE CurrentProblemValues, ONLY : ion_colls_with_bo
  USE Snapshots

  IMPLICIT NONE

  INTEGER s
  REAL coll_coord, vx, vy, vz
  INTEGER tag
  
  INTEGER nwo   ! number of the boundary object
  INTEGER nseg  ! number of the segment of the boundary object
  
  TYPE collided_particle
     INTEGER token
     REAL coll_coord
     REAL VX
     REAL VY
     REAL VZ
  END TYPE collided_particle

  TYPE(collided_particle), ALLOCATABLE :: bufer(:)
  INTEGER ALLOC_ERR

  INTEGER k
  INTEGER current_N

  IF (.NOT.ion_colls_with_bo(nwo)%must_be_saved) RETURN

  IF (current_snap.GT.N_of_all_snaps) RETURN

  IF (.NOT.save_ions_collided_with_bo(current_snap)) RETURN

  ion_colls_with_bo(nwo)%N_of_saved_parts = ion_colls_with_bo(nwo)%N_of_saved_parts+1

  IF (ion_colls_with_bo(nwo)%N_of_saved_parts.GT.ion_colls_with_bo(nwo)%max_N_of_saved_parts) THEN
! increase the size of the array
     current_N = ion_colls_with_bo(nwo)%max_N_of_saved_parts
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%token      = ion_colls_with_bo(nwo)%part(k)%token
        bufer(k)%coll_coord = ion_colls_with_bo(nwo)%part(k)%coll_coord
        bufer(k)%VX         = ion_colls_with_bo(nwo)%part(k)%VX
        bufer(k)%VY         = ion_colls_with_bo(nwo)%part(k)%VY
        bufer(k)%VZ         = ion_colls_with_bo(nwo)%part(k)%VZ
     END DO
     IF (ALLOCATED(ion_colls_with_bo(nwo)%part)) DEALLOCATE(ion_colls_with_bo(nwo)%part, STAT=ALLOC_ERR)
     ion_colls_with_bo(nwo)%max_N_of_saved_parts = ion_colls_with_bo(nwo)%max_N_of_saved_parts + MAX(50, ion_colls_with_bo(nwo)%max_N_of_saved_parts/10)
     ALLOCATE(ion_colls_with_bo(nwo)%part(1:ion_colls_with_bo(nwo)%max_N_of_saved_parts), STAT=ALLOC_ERR)
     DO k = 1, current_N
        ion_colls_with_bo(nwo)%part(k)%token      = bufer(k)%token
        ion_colls_with_bo(nwo)%part(k)%coll_coord = bufer(k)%coll_coord
        ion_colls_with_bo(nwo)%part(k)%VX         = bufer(k)%VX
        ion_colls_with_bo(nwo)%part(k)%VY         = bufer(k)%VY
        ion_colls_with_bo(nwo)%part(k)%VZ         = bufer(k)%VZ
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=ALLOC_ERR)
  END IF

  k = ion_colls_with_bo(nwo)%N_of_saved_parts

  ion_colls_with_bo(nwo)%part(k)%token = tag + 100 * nseg + 10000 * s

  ion_colls_with_bo(nwo)%part(k)%coll_coord = coll_coord
  ion_colls_with_bo(nwo)%part(k)%VX = vx
  ion_colls_with_bo(nwo)%part(k)%VY = vy
  ion_colls_with_bo(nwo)%part(k)%VZ = vz

END SUBROUTINE ADD_ION_TO_BO_COLLS_LIST
