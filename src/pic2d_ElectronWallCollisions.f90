!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr

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

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        CALL ADD_ELECTRON_TO_BO_COLLS_LIST(REAL(y), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = y - jbelow
              dqbelow = 1.0_8 - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) - dqbelow
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) - dqabove
        END SELECT

        IF (whole_object(nwo)%SEE_enabled) THEN
           CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, whole_object(nwo), m, 1)   ! "1" is for a left wall 
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     errcode=240
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT

!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr

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

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        CALL ADD_ELECTRON_TO_BO_COLLS_LIST(REAL(y), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = y - jbelow
              dqbelow = 1.0_8 - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) - dqbelow
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) - dqabove
        END SELECT

        IF (whole_object(nwo)%SEE_enabled) THEN
           CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, whole_object(nwo), m, 3)   ! "3" is for a right wall 
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     errcode=241
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT

!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr

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

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        CALL ADD_ELECTRON_TO_BO_COLLS_LIST(REAL(x), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = x - ileft
              dqleft = 1.0_8 - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft) - dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) - dqright
        END SELECT

        IF (whole_object(nwo)%SEE_enabled) THEN
           CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, whole_object(nwo), m, 4)   ! "4" is for a wall below
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     errcode=242
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW

!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr

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

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        CALL ADD_ELECTRON_TO_BO_COLLS_LIST(REAL(x), REAL(vx), REAL(vy), REAL(vz), tag, nwo, c_local_object_part(m)%segment_number)

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = x - ileft
              dqleft = 1.0_8 - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft) - dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) - dqright
        END SELECT

        IF (whole_object(nwo)%SEE_enabled) THEN
           CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, whole_object(nwo), m, 2)   ! "2" is for a wall above
        END IF

        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     errcode=243
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE

!------------------------------------------------------
!
SUBROUTINE COLLECT_ELECTRON_BOUNDARY_HITS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : N_spec

  use mpi

  IMPLICIT NONE


  INTEGER ierr
!  INTEGER stattus(MPI_STATUS_SIZE)
!  INTEGER request

  INTEGER, ALLOCATABLE :: ibuf_send(:)
  INTEGER, ALLOCATABLE :: ibuf_receive(:)
  INTEGER ALLOC_ERR

  INTEGER k

  ALLOCATE(ibuf_send(1:N_of_boundary_and_inner_objects), STAT = ALLOC_ERR)
  ALLOCATE(ibuf_receive(1:N_of_boundary_and_inner_objects), STAT = ALLOC_ERR)

! each cluster adjacent to a boundary assembles electron-boundary hit counters from all cluster members in the master of the cluster

  ibuf_send(1:N_of_boundary_and_inner_objects) = whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count
  ibuf_receive = 0

  CALL MPI_REDUCE(ibuf_send, ibuf_receive, N_of_boundary_and_inner_objects, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)  !??? use Rank_of_bottom_left_cluster_master ???

  IF (Rank_of_process.EQ.0) THEN
! now counters from all processes are assembled in the process with global rank zero
    
     whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count = ibuf_receive(1:N_of_boundary_and_inner_objects)
     print '("electrons hit boundaries :: ",10(2x,i8))', whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count  

     DO k = 1, N_of_boundary_and_inner_objects
        whole_object(k)%ion_hit_count(1:N_spec) = 0
     END DO
     
  END IF

  DEALLOCATE(ibuf_send, STAT = ALLOC_ERR)
  DEALLOCATE(ibuf_receive, STAT = ALLOC_ERR)

END SUBROUTINE COLLECT_ELECTRON_BOUNDARY_HITS

!-------------------------------------------------------------------------------------------
!
SUBROUTINE INITIATE_WALL_DIAGNOSTICS

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : N_of_boundary_and_inner_objects, Start_T_cntr
  USE Checkpoints, ONLY : use_checkpoint
!  USE Diagnostics, ONLY : N_of_saved_records
  USE SetupValues, ONLY : ht_use_e_emission_from_cathode, ht_use_e_emission_from_cathode_zerogradf, ht_emission_constant

  use mpi

  IMPLICIT NONE

                                    ! ----x----I----x--
  CHARACTER(17) historybo_filename  ! history_bo_NN.dat

  LOGICAL exists
  INTEGER i, k
  INTEGER i_dummy, ios

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (Rank_of_process.NE.0) RETURN

  IF (ht_use_e_emission_from_cathode.OR.ht_use_e_emission_from_cathode_zerogradf.OR.ht_emission_constant) RETURN

! hardwired for objects #2 (cathode) and #4 (anode)

  IF (use_checkpoint.EQ.1) THEN
! start from checkpoint, must trim the time dependences

     DO k = 1, N_of_boundary_and_inner_objects

        historybo_filename = 'history_bo_NN.dat'
        historybo_filename(12:13) = convert_int_to_txt_string(k, 2)

        INQUIRE (FILE = historybo_filename, EXIST = exists)
        IF (exists) THEN                                                       
           OPEN (21, FILE = historybo_filename, STATUS = 'OLD')          
           DO !i = 1, Start_T_cntr   !N_of_saved_records             ! these files are updated at every electron timestep
              READ (21, '(2x,i8,10(2x,i8))', iostat = ios) i_dummy
              IF (ios.NE.0) EXIT
              IF (i_dummy.GE.Start_T_cntr) EXIT
           END DO
           BACKSPACE(21)
           ENDFILE 21       
           CLOSE (21, STATUS = 'KEEP')        
        END IF

     END DO

  ELSE
! fresh start, empty files, clean up whatever garbage there might be

     DO k = 1, N_of_boundary_and_inner_objects

        historybo_filename = 'history_bo_NN.dat'
        historybo_filename(12:13) = convert_int_to_txt_string(k, 2)

        OPEN  (21, FILE = historybo_filename, STATUS = 'REPLACE')          
        CLOSE (21, STATUS = 'KEEP')

     END DO

  END IF

END SUBROUTINE INITIATE_WALL_DIAGNOSTICS

!-------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE SetupValues, ONLY : ht_use_e_emission_from_cathode, ht_use_e_emission_from_cathode_zerogradf, ht_emission_constant
  USE IonParticles, ONLY : N_spec, Qs
  USE ExternalCircuit

  use mpi

  IMPLICIT NONE

  INTEGER nn, noi, s
  INTEGER k
                                    ! ----x----I----x--
  CHARACTER(17) historybo_filename  ! history_bo_NN.dat

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (Rank_of_process.NE.0) RETURN

  IF (ht_use_e_emission_from_cathode.OR.ht_use_e_emission_from_cathode_zerogradf.OR.ht_emission_constant) RETURN

  DO nn = 1, N_of_object_potentials_to_solve
     noi = object_charge_calculation(1)%noi
     dQ_plasma_of_object(nn) = -whole_object(noi)%electron_hit_count + &
                             &  whole_object(noi)%electron_emit_count                      ! include electron emission
     DO s = 1, N_spec
        dQ_plasma_of_object(nn) = dQ_plasma_of_object(nn) + Qs(s) * whole_object(noi)%ion_hit_count(s)
     END DO
  END DO

  DO k = 1, N_of_boundary_and_inner_objects

     historybo_filename = 'history_bo_NN.dat'
     historybo_filename(12:13) = convert_int_to_txt_string(k, 2)

     OPEN (21, FILE = historybo_filename, POSITION = 'APPEND')
     WRITE (21, '(2x,i8,10(2x,i8))') &
          & T_cntr, &
          & whole_object(k)%electron_hit_count , &
          & whole_object(k)%ion_hit_count(1:N_spec), &
          & whole_object(k)%electron_emit_count

     CLOSE (21, STATUS = 'KEEP')
  
  END DO

END SUBROUTINE SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS

!-----------------------------------------
!
SUBROUTINE TRY_ELECTRON_COLL_WITH_INNER_OBJECT(x, y, vx, vy, vz, tag) !, myobject)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues !, ONLY : inner_object, METAL_WALL, DIELECTRIC

  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr

!  INTEGER nio  ! number of the inner object
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

  xorg = x - vx
  yorg = y - vy

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
     PRINT '("Error-1 in TRY_ELECTRON_COLL_WITH_INNER_OBJECT ",4(2x,f10.4))', xorg, yorg, x, y
     errcode=244
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
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

  CALL ADD_ELECTRON_TO_BO_COLLS_LIST(coll_coord, REAL(vx), REAL(vy), REAL(vz), tag, n_do, mcross)

  CALL DO_ELECTRON_COLL_WITH_INNER_OBJECT(xcross, ycross, vx, vy, vz, tag, whole_object(n_do), coll_direction_flag)

END SUBROUTINE TRY_ELECTRON_COLL_WITH_INNER_OBJECT

!-----------------------------------------
!
SUBROUTINE DO_ELECTRON_COLL_WITH_INNER_OBJECT(x, y, vx, vy, vz, tag, myobject, coll_direction_flag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues !, ONLY : inner_object, METAL_WALL, DIELECTRIC

  use mpi

  IMPLICIT NONE

!  INTEGER nio  ! number of the inner object
  REAL(8) x, y, vx, vy, vz
  INTEGER tag
  TYPE(boundary_object) myobject
  INTEGER coll_direction_flag

  REAL(8) xmin, xmax, ymin, ymax

  INTEGER i_left_top, i_right_top, i_right_bottom, i_left_bottom_bis, i, ip1

  REAL(8) dqip1, dqi

! identify side of the inner object hit by the particle

  xmin = myobject%Xmin
  xmax = myobject%Xmax
  ymin = myobject%Ymin
  ymax = myobject%Ymax

  myobject%electron_hit_count = myobject%electron_hit_count + 1

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
           dqip1 = y - INT(y)
           dqi = 1.0_8 - dqip1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   - dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) - dqip1
           
        CASE (4)
! top wall of the object
           i = MIN(INT(x - xmin) + i_left_top, i_right_top - 1)
           dqip1 = x - INT(x)
           dqi = 1.0_8 - dqip1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   - dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) - dqip1

        CASE (1)
! right wall of the object
           i = MIN(INT(ymax - y) + i_right_top, i_right_bottom - 1)
           dqi = y - INT(y)
           dqip1 = 1.0_8 - dqi
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   - dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) - dqip1

        CASE (2)
! bottom wall of the object
           i = MIN(INT(xmax - x) + i_right_bottom, i_left_bottom_bis)
           dqi = x - INT(x)
           dqip1 = 1.0_8 - dqi
           ip1 = i+1
           IF (i.EQ.i_left_bottom_bis) ip1 = 1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   - dqi
           myobject%surface_charge_variation(ip1) = myobject%surface_charge_variation(ip1) - dqip1

     END SELECT
  END IF   !### IF (myobject%object_type.EQ.DIELECTRIC) THEN

  IF (myobject%SEE_enabled) THEN
     CALL PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, myobject, -1, coll_direction_flag)
  END IF

END SUBROUTINE DO_ELECTRON_COLL_WITH_INNER_OBJECT

!--------------------------
!
SUBROUTINE FIND_CLOSEST_INTERSECTION_WITH_OBJECT(xorg, yorg, x, y, n, mcross, xcross, ycross, distorg)

  USE CurrentProblemValues !, ONLY : inner_object, METAL_WALL, DIELECTRIC

  use mpi

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: xorg, yorg, x, y
  INTEGER, INTENT(IN) :: n

  INTEGER, INTENT(OUT) :: mcross
  REAL(8), INTENT(OUT) :: xcross, ycross
  REAL(8), INTENT(OUT) :: distorg

  INTEGER m

  INTEGER jbot, jtop, jcross
  INTEGER ileft, iright, icross
  REAL(8) myxcross, myycross
  INTEGER mystatus
  REAL(8) mydistorg

  mcross = -1

  DO m = 1, whole_object(n)%number_of_segments

     IF (whole_object(n)%segment(m)%istart.EQ.whole_object(n)%segment(m)%iend) THEN
! vertical segment

        jbot = MIN(whole_object(n)%segment(m)%jstart, whole_object(n)%segment(m)%jend)
        jtop = MAX(whole_object(n)%segment(m)%jstart, whole_object(n)%segment(m)%jend)
        myxcross = DBLE(whole_object(n)%segment(m)%istart)

        CALL CHECK_INTERSECTION_WITH_VERTICAL_SEGMENT( xorg, yorg, x, y, myxcross, DBLE(jbot), DBLE(jtop), mystatus, myycross )
        IF (mystatus.NE.0) CYCLE

        jcross = MAX(jbot, MIN(jtop-1, INT(myycross)))

!print *,"aa ", jcross, n, m, whole_object(n)%segment(m)%cell_is_covered(jcross)

! check that the intersection point is not in the prohibited (covered) part of the segment
        IF (whole_object(n)%segment(m)%cell_is_covered(jcross)) CYCLE

        mydistorg = (myxcross-xorg)**2 + (myycross-yorg)**2

        IF (mcross.EQ.-1) THEN
           mcross = m
           xcross = myxcross
           ycross = myycross
           distorg = mydistorg
        ELSE
           IF (mydistorg.GE.distorg) CYCLE
           mcross = m
           xcross = myxcross
           ycross = myycross
           distorg = mydistorg
        END IF

     ELSE IF (whole_object(n)%segment(m)%jstart.EQ.whole_object(n)%segment(m)%jend) THEN
! horizontal segment

        ileft  = MIN(whole_object(n)%segment(m)%istart, whole_object(n)%segment(m)%iend)
        iright = MAX(whole_object(n)%segment(m)%istart, whole_object(n)%segment(m)%iend)
        myycross = DBLE(whole_object(n)%segment(m)%jstart)

        CALL CHECK_INTERSECTION_WITH_HORIZONTAL_SEGMENT( xorg, yorg, x, y, myycross, DBLE(ileft), DBLE(iright), mystatus, myxcross )
        IF (mystatus.NE.0) CYCLE

        icross = MAX(ileft, MIN(iright-1, INT(myxcross)))

!print *,"bb ", icross, n, m, whole_object(n)%segment(m)%cell_is_covered(icross)

! check that the intersection point is not in the prohibited (covered) part of the segment
        IF (whole_object(n)%segment(m)%cell_is_covered(icross)) CYCLE

        mydistorg = (myxcross-xorg)**2 + (myycross-yorg)**2

        IF (mcross.EQ.-1) THEN
           mcross = m
           xcross = myxcross
           ycross = myycross
           distorg = mydistorg
        ELSE
           IF (mydistorg.GE.distorg) CYCLE
           mcross = m
           xcross = myxcross
           ycross = myycross
           distorg = mydistorg
        END IF

     END IF  !### ELSE IF (whole_object(n)%segment(m)%jstart.EQ.whole_object(n)%segment(m)%jend) THEN

  END DO   !### DO m = 1, whole_object(n)%number_of_segments

END SUBROUTINE FIND_CLOSEST_INTERSECTION_WITH_OBJECT

!----------------------------------------------------------
!### assume that yminseg < ymaxseg
!
SUBROUTINE CHECK_INTERSECTION_WITH_VERTICAL_SEGMENT( xorg, yorg, x, y, xseg, yminseg, ymaxseg, mystatus, ycross )

  use, intrinsic :: ieee_arithmetic

  use mpi

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: xorg, yorg, x, y          ! coordinates of the ends of particle trajectory segment
  REAL(8), INTENT(IN) :: xseg, yminseg, ymaxseg    ! coordinates of the ends of vertical boundary segment
  INTEGER, INTENT(OUT) :: mystatus                 ! zero if crossing found, nonzero otherwise
  REAL(8), INTENT(OUT) :: ycross                   ! y-coordinate of the crossing

  mystatus = -1

! check obvious things first
  IF (MAX(xorg, x).LT.xseg) RETURN
  IF (MIN(xorg, x).GT.xseg) RETURN
  IF (MAX(yorg, y).LT.yminseg) RETURN
  IF (MIN(yorg, y).GT.ymaxseg) RETURN

! extremely unlikely situation, particle goes exactly along the surface of the object
  IF (xorg.EQ.x) RETURN

! since we are here, ends of segment {xorg,yorg}-{x,y} are on different sides of segment {xseg,yminseg}-{xseg,ymaxseg}

  ycross = yorg + (y-yorg) * (xseg-xorg) / (x-xorg)

  IF (.NOT.ieee_is_finite(ycross)) THEN
     mystatus = 1
     RETURN
  END IF

  IF (ycross.LT.yminseg) RETURN
  IF (ycross.GT.ymaxseg) RETURN
  IF (ycross.LT.MIN(y,yorg)) RETURN  ! paranoidal failsafe check?
  IF (ycross.GT.MAX(y,yorg)) RETURN  !

  mystatus = 0
  RETURN

END SUBROUTINE CHECK_INTERSECTION_WITH_VERTICAL_SEGMENT

!----------------------------------------------------------
!### assume that xminseg < xmagseg
!
SUBROUTINE CHECK_INTERSECTION_WITH_HORIZONTAL_SEGMENT( xorg, yorg, x, y, yseg, xminseg, xmaxseg, mystatus, xcross )

  use, intrinsic :: ieee_arithmetic

  use mpi

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: xorg, yorg, x, y          ! coordinates of the ends of particle trajectory segment
  REAL(8), INTENT(IN) :: yseg, xminseg, xmaxseg    ! coordinates of the ends of horizontal boundary segment
  INTEGER, INTENT(OUT) :: mystatus                 ! zero if crossing found, nonzero otherwise
  REAL(8), INTENT(OUT) :: xcross                   ! x-coordinate of the crossing

  mystatus = -1

! check obvious things first
  IF (MAX(yorg, y).LT.yseg) RETURN
  IF (MIN(yorg, y).GT.yseg) RETURN
  IF (MAX(xorg, x).LT.xminseg) RETURN
  IF (MIN(xorg, x).GT.xmaxseg) RETURN

! extremely unlikely situation, particle goes exactly along the surface of the object
  IF (yorg.EQ.y) RETURN

! since we are here, ends of segment {xorg,yorg}-{x,y} are on different sides of segment {xseg,yminseg}-{xseg,ymaxseg}

  xcross = xorg + (x-xorg) * (yseg-yorg) / (y-yorg)

  IF (.NOT.ieee_is_finite(xcross)) THEN
     mystatus = 1
     RETURN
  END IF

  IF (xcross.LT.xminseg) RETURN
  IF (xcross.GT.xmaxseg) RETURN
  IF (xcross.LT.MIN(x,xorg)) RETURN  ! paranoidal failsafe check?
  IF (xcross.GT.MAX(x,xorg)) RETURN  !

  mystatus = 0
  RETURN

END SUBROUTINE CHECK_INTERSECTION_WITH_HORIZONTAL_SEGMENT

!-------------------------------------------------------------------------------------------
! Prepares the tabulated values of integral of the maxwell distribution function
!  
SUBROUTINE PrepareMaxwellDistribIntegral

  USE ParallelOperationValues
  USE MaxwellVelocity
!  USE CurrentProblemValues, ONLY : N_box_vel
  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr

  INTEGER i
  INTEGER N_pnts
  INTEGER count
  REAL(8) V_min, V_max
  REAL(8) F(0:180003)      ! to be sure that we overcome V_max
  REAL(8) temp
  REAL(8) dV

  LOGICAL check1, check2

  check1 = .FALSE.
  check2 = .FALSE.

! ------- for symmetrical maxwellian
  N_pnts = 180000  !30000
  V_min  = -U_max
  V_max  =  U_max
  dV = (V_max - V_min) / N_pnts

  F = 0.0_8
  DO i = 1, N_pnts + 3
     F(i) = F(i-1) + EXP( - (V_min + (DBLE(i)-0.5_8) * dV)**2 )
  END DO

  temp = F(N_pnts)
  F = F * R_max / temp   ! normalize integral such that F(N_pnts) = R_max

  v(0) = V_min
  count = 0
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        v(count) = V_min + i * dV
        IF (count.EQ.R_max) THEN
           check1 = .TRUE.
           EXIT
        END IF
     END IF
  END DO

!  v = v * N_box_vel

!--------- for asymmetrical maxwellian * v (used for injection, v > 0)

  N_pnts = 90000   !15000
  V_min  = 0.0_8
  V_max  = U_max
  dV = (V_max - V_min) / N_pnts

  F = 0.0_8
  DO i = 1, N_pnts + 3
     temp = V_min + (REAL(i)-0.5_8) * dV
     F(i) = F(i-1) + EXP( - temp**2 ) * temp
  END DO

  temp = F(N_pnts)
  F(1:(N_pnts+3)) = F(1:(N_pnts+3)) * R_max_inj / temp   ! normalize integral such that F(N_pnts) = R_max_inj

  v_inj(0) = V_min
  count = 0
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        v_inj(count) = V_min + i * dV
        IF (count.EQ.R_max_inj) THEN
           check2 = .TRUE.
           EXIT
        END IF
     END IF
  END DO

!  v_inj = v_inj * N_box_vel

  IF (check1.AND.check2) THEN
!     PRINT '(2x,"Process ",i3," : Integrals for producing maxwell distributions are successfully obtained ...")', &
!                                                                                                  & Rank_of_process
  ELSE
     PRINT '(2x,"Process ",i3," : ERROR in PrepareMaxwellDistribIntegral !!!")', Rank_of_process
     PRINT '(2x,"The initialization in PrepareMaxwellDistribIntegral is not performed !!!")'
     PRINT '(2x,"The program will be terminated now :(")'
     errcode=245
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

END SUBROUTINE PrepareMaxwellDistribIntegral

!-------------------------------------------------------------------------------------------
!  
SUBROUTINE GetInjMaxwellVelocity(U) 

  USE MaxwellVelocity
 
  USE rng_wrapper

  use mpi

  IMPLICIT NONE

  REAL(8) U

  REAL(8) R
  INTEGER indx
  
  R = R_max_inj * well_random_number()

  indx = INT(R)

  IF (indx.LT.R_max_inj) THEN
     U = v_inj(indx) + (R - indx) * (v_inj(indx+1) - v_inj(indx))
  ELSE
     U = v_inj(R_max_inj)
  END IF
  RETURN
  
END SUBROUTINE GetInjMaxwellVelocity

!-------------------------------------------------------------------------------------------
!  
SUBROUTINE GetMaxwellVelocity(U) 

  USE MaxwellVelocity

  USE rng_wrapper

  use mpi

  IMPLICIT NONE

  REAL(8) U

  REAL(8) R
  INTEGER indx
  
  R = R_max * well_random_number()

  indx = INT(R)

  IF (indx.LT.R_max) THEN
     U = v(indx) + (R - indx) * (v(indx+1) - v(indx))
  ELSE
     U = v(R_max)
  END IF
  RETURN
  
END SUBROUTINE GetMaxwellVelocity

!-------------------------
!
REAL(8) FUNCTION vector_product_z(ax, ay, bx, by)

  use mpi

  IMPLICIT NONE
  REAL(8) ax, ay, bx, by

  vector_product_z = ax * by - ay * bx

END FUNCTION vector_product_z

!---------------------------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_BO_COLLS_LIST(coll_coord, vx, vy, vz, tag, nwo, nseg)

  USE CurrentProblemValues, ONLY : e_colls_with_bo
  USE Snapshots

  use mpi

  IMPLICIT NONE

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

  IF (.NOT.e_colls_with_bo(nwo)%must_be_saved) RETURN

  IF (current_snap.GT.N_of_all_snaps) RETURN

  IF (.NOT.save_e_collided_with_bo(current_snap)) RETURN

  e_colls_with_bo(nwo)%N_of_saved_parts = e_colls_with_bo(nwo)%N_of_saved_parts+1

  IF (e_colls_with_bo(nwo)%N_of_saved_parts.GT.e_colls_with_bo(nwo)%max_N_of_saved_parts) THEN
! increase the size of the array
     current_N = e_colls_with_bo(nwo)%max_N_of_saved_parts
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%token      = e_colls_with_bo(nwo)%part(k)%token
        bufer(k)%coll_coord = e_colls_with_bo(nwo)%part(k)%coll_coord
        bufer(k)%VX         = e_colls_with_bo(nwo)%part(k)%VX
        bufer(k)%VY         = e_colls_with_bo(nwo)%part(k)%VY
        bufer(k)%VZ         = e_colls_with_bo(nwo)%part(k)%VZ
     END DO
     IF (ALLOCATED(e_colls_with_bo(nwo)%part)) DEALLOCATE(e_colls_with_bo(nwo)%part, STAT=ALLOC_ERR)
     e_colls_with_bo(nwo)%max_N_of_saved_parts = e_colls_with_bo(nwo)%max_N_of_saved_parts + MAX(50, e_colls_with_bo(nwo)%max_N_of_saved_parts/10)
     ALLOCATE(e_colls_with_bo(nwo)%part(1:e_colls_with_bo(nwo)%max_N_of_saved_parts), STAT=ALLOC_ERR)
     DO k = 1, current_N
        e_colls_with_bo(nwo)%part(k)%token      = bufer(k)%token
        e_colls_with_bo(nwo)%part(k)%coll_coord = bufer(k)%coll_coord
        e_colls_with_bo(nwo)%part(k)%VX         = bufer(k)%VX
        e_colls_with_bo(nwo)%part(k)%VY         = bufer(k)%VY
        e_colls_with_bo(nwo)%part(k)%VZ         = bufer(k)%VZ
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=ALLOC_ERR)
  END IF

  k = e_colls_with_bo(nwo)%N_of_saved_parts

  e_colls_with_bo(nwo)%part(k)%token = tag + 100 * nseg

  e_colls_with_bo(nwo)%part(k)%coll_coord = coll_coord
  e_colls_with_bo(nwo)%part(k)%VX = vx
  e_colls_with_bo(nwo)%part(k)%VY = vy
  e_colls_with_bo(nwo)%part(k)%VZ = vz

END SUBROUTINE ADD_ELECTRON_TO_BO_COLLS_LIST
