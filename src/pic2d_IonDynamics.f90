!---------------------------------------
!
SUBROUTINE ADVANCE_IONS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles
  USE ClusterAndItsBoundaries
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER s, k
  INTEGER i, j
  REAL(8) ax_ip1, ax_i, ay_jp1, ay_j
  REAL(8) E_X, E_Y, E_Z
  REAL(8) alfa_x, alfa_y, alfa_z
  REAL(8) alfa_x2, alfa_y2, alfa_z2
  REAL(8) theta2, invtheta
  REAL(8) K11, K12, K13, K21, K22, K23, K31, K32, K33
  REAL(8) VX_minus, VY_minus, VZ_minus
  REAL(8) VX_plus, VY_plus, VZ_plus

  INTEGER n
  LOGICAL collision_with_inner_object_occurred

! functions
  REAL(8) Bx, By, Bz, Ez

! clear counters of particles to be sent to neighbor processes
  N_ions_to_send_left = 0
  N_ions_to_send_right = 0
  N_ions_to_send_above = 0
  N_ions_to_send_below = 0

! clear counters of particles that hit the boundary objects
  DO k = 1, N_of_boundary_and_inner_objects
     whole_object(k)%ion_hit_count(1:N_spec) = 0
     ion_colls_with_bo(k)%N_of_saved_parts = 0
  END DO

! cycle over ion species
  DO s = 1, N_spec

! cycle over particles of the ion species
     k=0
     DO WHILE (k.LT.N_ions(s))

        k = k + 1

        if ( (ion(s)%part(k)%X.lt.c_X_area_min).or. &
           & (ion(s)%part(k)%X.gt.c_X_area_max).or. &
           & (ion(s)%part(k)%Y.lt.c_Y_area_min).or. &
           & (ion(s)%part(k)%Y.gt.c_Y_area_max) ) then
           print '("Process ",i4," : Error-1 in ADVANCE_IONS : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
           print '("Process ",i4," : s/k/N_ions(s) : ",i3,2x,i8,2x,i8)', Rank_of_process, s, k, N_ions(s)
           print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        end if

! interpolate electric field

        i = INT(ion(s)%part(k)%X)
        j = INT(ion(s)%part(k)%Y)

        IF (ion(s)%part(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
        IF (ion(s)%part(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

        if ((i.lt.c_indx_x_min).or.(i.gt.(c_indx_x_max-1)).or.(j.lt.c_indx_y_min).or.(j.gt.(c_indx_y_max-1))) then
           print '("Process ",i4," : Error-2 in ADVANCE_IONS : index out of bounds")', Rank_of_process
           print '("Process ",i4," : s/k/N_ions(s) : ",i3,2x,i8,2x,i8)', Rank_of_process, s, k, N_ions(s)
           print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag
           print '("Process ",i4," : minx/maxx/miny/maxy : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        end if

        ax_ip1 = ion(s)%part(k)%X - DBLE(i)
        ax_i   = 1.0_8 - ax_ip1

        ay_jp1 = ion(s)%part(k)%Y - DBLE(j)
        ay_j = 1.0_8 - ay_jp1

        E_X = acc_EX(i,j) * ax_i * ay_j + acc_EX(i+1,j) * ax_ip1 * ay_j + acc_EX(i,j+1) * ax_i * ay_jp1 + acc_EX(i+1,j+1) * ax_ip1 * ay_jp1   ! use accumulaed (averaged) electric fields
        E_Y = acc_EY(i,j) * ax_i * ay_j + acc_EY(i+1,j) * ax_ip1 * ay_j + acc_EY(i,j+1) * ax_i * ay_jp1 + acc_EY(i+1,j+1) * ax_ip1 * ay_jp1   !

        IF (ions_sense_Ez) THEN
           E_Z = Ez(ion(s)%part(k)%X, ion(s)%part(k)%Y) * N_subcycles         ! Aug-3-2017 found a bug here, N_subcycles was missing
        ELSE
           E_Z = 0.0_8
        END IF

        IF (ions_sense_magnetic_field) THEN
! magnetic field accounted for
! calculate magnetic field factors   !##### MODIFY FOR IONS ########

           alfa_x = QM2sNsub(s) * Bx(ion(s)%part(k)%X, ion(s)%part(k)%Y)
           alfa_y = QM2sNsub(s) * By(ion(s)%part(k)%X, ion(s)%part(k)%Y)
           alfa_z = QM2sNsub(s) * Bz(ion(s)%part(k)%X, ion(s)%part(k)%Y)

           alfa_x2 = alfa_x**2
           alfa_y2 = alfa_y**2
           alfa_z2 = alfa_z**2

           theta2 = alfa_x2 + alfa_y2 + alfa_z2
           invtheta = 1.0_8 / (1.0_8 + theta2)

           K11 = (1.0_8 - theta2 + 2.0_8 * alfa_x2) * invtheta
           K12 =  2.0_8 * (alfa_x * alfa_y + alfa_z) * invtheta
           K13 =  2.0_8 * (alfa_x * alfa_z - alfa_y) * invtheta

           K21 =  2.0_8 * (alfa_x * alfa_y - alfa_z) * invtheta
           K22 = (1.0_8 - theta2 + 2.0_8 * alfa_y2) * invtheta
           K23 =  2.0_8 * (alfa_y * alfa_z + alfa_x) * invtheta

           K31 =  2.0_8 * (alfa_x * alfa_z + alfa_y) * invtheta
           K32 =  2.0_8 * (alfa_y * alfa_z - alfa_x) * invtheta
           K33 = (1.0_8 - theta2 + 2.0_8 * alfa_z2) * invtheta

!print *, K11*(K22*K33-K23*K32) - K12*(K21*K33-K23*K31) + K13*(K21*K32-K22*K31), E_X*E_scale_Vm

! velocity advance: first half-acceleration due to electric field

           VX_minus = ion(s)%part(k)%VX + QM2s(s) * E_X
           VY_minus = ion(s)%part(k)%VY + QM2s(s) * E_Y
           VZ_minus = ion(s)%part(k)%VZ + QM2s(s) * E_Z

! velocity advance: rotation in the magnetic field

           VX_plus = K11 * VX_minus + K12 * VY_minus + K13 * VZ_minus
           VY_plus = K21 * VX_minus + K22 * VY_minus + K23 * VZ_minus
           VZ_plus = K31 * VX_minus + K32 * VY_minus + K33 * VZ_minus

! velocity advance: second half-acceleration due to electric field

           ion(s)%part(k)%VX = VX_plus + QM2s(s) * E_X
           ion(s)%part(k)%VY = VY_plus + QM2s(s) * E_Y
           ion(s)%part(k)%VZ = VZ_plus + QM2s(s) * E_Z

        ELSE
! magnetic field effects omitted
! velocity advance:  combine both half-accelerations, no magnetic field

           ion(s)%part(k)%VX = ion(s)%part(k)%VX + (QM2s(s) + QM2s(s)) * E_X
           ion(s)%part(k)%VY = ion(s)%part(k)%VY + (QM2s(s) + QM2s(s)) * E_Y
           ion(s)%part(k)%VZ = ion(s)%part(k)%VZ + (QM2s(s) + QM2s(s)) * E_Z

        END IF

! coordinate advance

        ion(s)%part(k)%X = ion(s)%part(k)%X + ion(s)%part(k)%VX * N_subcycles
        ion(s)%part(k)%Y = ion(s)%part(k)%Y + ion(s)%part(k)%VY * N_subcycles

! a particle crossed symmetry plane, reflect it
        IF (symmetry_plane_X_left) THEN
           IF (ion(s)%part(k)%X.LT.c_X_area_min) THEN
              ion(s)%part(k)%X = MAX(c_X_area_min, c_X_area_min + c_X_area_min - ion(s)%part(k)%X)
              ion(s)%part(k)%VX = -ion(s)%part(k)%VX
           END IF
        END IF

! check whether a collision with an inner object occurred
        collision_with_inner_object_occurred = .FALSE.
        DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (ion(s)%part(k)%X.LE.whole_object(n)%Xmin) CYCLE
           IF (ion(s)%part(k)%X.GE.whole_object(n)%Xmax) CYCLE
           IF (ion(s)%part(k)%Y.LE.whole_object(n)%Ymin) CYCLE
           IF (ion(s)%part(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! collision detected
           CALL TRY_ION_COLL_WITH_INNER_OBJECT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag) !, whole_object(n))
           CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           collision_with_inner_object_occurred = .TRUE.
           EXIT
        END DO

        IF (collision_with_inner_object_occurred) CYCLE

! most probable situation when the particle remains inside the area
        IF ( (ion(s)%part(k)%X.GE.c_X_area_min) .AND. &
           & (ion(s)%part(k)%X.LE.c_X_area_max) .AND. &
           & (ion(s)%part(k)%Y.GE.c_Y_area_min) .AND. &
           & (ion(s)%part(k)%Y.LE.c_Y_area_max) ) CYCLE

! since we are here, a particle did not collide with an inner object but crossed an area boundary
! note, in a periodic system the particle still may collide with an inner object after transfer to the other domain
! therefore it is still necessary to check add list for collisions after process receives particles from neighbors

        IF (ion(s)%part(k)%X.LT.c_X_area_min) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the ion by the period length
              ion(s)%part(k)%X = ion(s)%part(k)%X + L_period_X
              IF (ion(s)%part(k)%Y.LT.c_Y_area_min) THEN
! particle is somewhere near the left bottom cell of the area
                 IF (Rank_of_master_below.LT.0) THEN
                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 ELSE
                    CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)                       
                 END IF
                 CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              ELSE IF (ion(s)%part(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
                 IF (Rank_of_master_above.LT.0) THEN
                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 ELSE
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 END IF
                 CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              END IF
              CYCLE
           END IF

           IF ( (ion(s)%part(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
              & (ion(s)%part(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

              IF (Rank_of_master_left.GE.0) THEN
! left neighbor cluster exists
                 CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
! left neighbor cluster does not exist
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX,  ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)   ! left
              END IF

           ELSE IF (ion(s)%part(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the left bottom cell of the area

              SELECT CASE (c_left_bottom_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((c_X_area_min-ion(s)%part(k)%X).LT.(c_Y_area_min-ion(s)%part(k)%Y)) THEN
                       CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)

                 CASE (FLAT_WALL_LEFT)
                    IF (ion(s)%part(k)%Y.GE.c_Y_area_min) THEN                 
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)                       
                    END IF

                 CASE (SURROUNDED_BY_WALL)
                    IF ((c_X_area_min-ion(s)%part(k)%X).LT.(c_Y_area_min-ion(s)%part(k)%Y)) THEN
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    END IF

                 CASE (EMPTY_CORNER_WALL_LEFT)
                    CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
               
                 CASE (EMPTY_CORNER_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)

              END SELECT

           ELSE IF (ion(s)%part(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top left cell of the area

              SELECT CASE (c_left_top_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((c_X_area_min-ion(s)%part(k)%X).LT.(ion(s)%part(k)%Y-c_Y_area_max)) THEN
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)

                 CASE (FLAT_WALL_LEFT)
                    IF (ion(s)%part(k)%Y.LE.c_Y_area_max) THEN                 
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    END IF

                 CASE (SURROUNDED_BY_WALL)
                    IF ((c_X_area_min-ion(s)%part(k)%X).LT.(ion(s)%part(k)%Y-c_Y_area_max)) THEN
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    END IF

                 CASE (EMPTY_CORNER_WALL_LEFT)
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)

                 CASE (EMPTY_CORNER_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)

              END SELECT

           ELSE
! ERROR, we shouldn't be here
              PRINT '("ERROR-1 in ADVANCE_IONS: we should not be here")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF

           CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE

        END IF

        IF (ion(s)%part(k)%X.GT.c_X_area_max) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the ion by the period length
              ion(s)%part(k)%X = ion(s)%part(k)%X - L_period_X
              IF (ion(s)%part(k)%Y.LT.c_Y_area_min) THEN
! particle is somewhere near the left bottom cell of the area
                 IF (Rank_of_master_below.LT.0) THEN
                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 ELSE
                    CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)                       
                 END IF
                 CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              ELSE IF (ion(s)%part(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
                 IF (Rank_of_master_above.LT.0) THEN
                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 ELSE
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 END IF
                 CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              END IF
              CYCLE
           END IF
           
           IF ( (ion(s)%part(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
              & (ion(s)%part(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

              IF (Rank_of_master_right.GE.0) THEN
! right neighbor cluster exists
                 CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
! right neighbor cluster does not exist
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              END IF

           ELSE IF (ion(s)%part(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the right bottom cell of the area

              SELECT CASE (c_right_bottom_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((ion(s)%part(k)%X-c_X_area_max).LT.(c_Y_area_min-ion(s)%part(k)%Y)) THEN
                       CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                
                 CASE (FLAT_WALL_RIGHT)
                    IF (ion(s)%part(k)%Y.GE.c_Y_area_min) THEN
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)                       
                    END IF

                 CASE (SURROUNDED_BY_WALL)
                    IF ((ion(s)%part(k)%X-c_X_area_max).LT.(c_Y_area_min-ion(s)%part(k)%Y)) THEN
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    END IF

                 CASE (EMPTY_CORNER_WALL_RIGHT)
                    CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)

                 CASE (EMPTY_CORNER_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)

              END SELECT

           ELSE IF (ion(s)%part(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top right cell of the area

              SELECT CASE (c_right_top_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((ion(s)%part(k)%X-c_X_area_max).LT.(ion(s)%part(k)%Y-c_Y_area_max)) THEN
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)

                 CASE (FLAT_WALL_RIGHT)
                    IF (ion(s)%part(k)%Y.LE.c_Y_area_max) THEN
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)                    
                    END IF

                 CASE (SURROUNDED_BY_WALL)
                    IF ((ion(s)%part(k)%X-c_X_area_max).LT.(ion(s)%part(k)%Y-c_Y_area_max)) THEN
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    ELSE
                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                    END IF

                 CASE (EMPTY_CORNER_WALL_RIGHT)
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)

                 CASE (EMPTY_CORNER_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)

              END SELECT

           ELSE
! ERROR, we shouldn't be here
              PRINT '("ERROR-2 in ADVANCE_IONS: we should not be here")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF

           CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE

        END IF

! since we are here, c_X_area_min <= ion(s)&part(k)%Y <= c_X_area_max

        IF (ion(s)%part(k)%Y.GT.c_Y_area_max) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
              IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
                 CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              END IF
              CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              CYCLE
           END IF

           IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
              CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           ELSE
! neighbor cluster above does not exist
              IF ((ion(s)%part(k)%X.GE.(c_X_area_min+1.0_8)).AND.(ion(s)%part(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE IF (ion(s)%part(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left top corner
                 IF (c_left_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                    CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 ELSE
                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 END IF
              ELSE IF (ion(s)%part(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right top corner
                 IF (c_right_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX,  ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 ELSE
                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 END IF
              END IF
           END IF
           CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE
        END IF

        IF (ion(s)%part(k)%Y.LT.c_Y_area_min) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
              IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster above exists
                 CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              END IF
              CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              CYCLE
           END IF

           IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster below exists, remove particle and prepare to send it to the neighbor below
              CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           ELSE
! neighbor cluster below does not exist
              IF ((ion(s)%part(k)%X.GE.(c_X_area_min+1.0_8)).AND.(ion(s)%part(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE IF (ion(s)%part(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left bottom corner
                 IF (c_left_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                    CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 ELSE
                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 END IF
              ELSE IF (ion(s)%part(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right bottom corner
                 IF (c_right_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 ELSE
                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
                 END IF
              END IF
           END IF
           CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE
        END IF
! 

     END DO  ! end of cycle over particles of ion species

  END DO ! end of cycle over ion species

!print '("Proc ",i4," will add ",i8," send l/r/a/b ",4(2x,i8))', Rank_of_process, N_e_to_add, N_ions_to_send_left, N_ions_to_send_right, N_ions_to_send_above, N_ions_to_send_below

END SUBROUTINE ADVANCE_IONS

!----------------------------------------
!
SUBROUTINE REMOVE_ION(s, k)

  USE ParallelOperationValues
  USE IonParticles

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(IN) :: s
  INTEGER, INTENT(INOUT) :: k

  IF ((s.LT.1).OR.(s.GT.N_spec)) THEN
     PRINT '("Process ",i6," : ERROR-1 in REMOVE_ION : index s invalid")', Rank_of_process
     PRINT '("Process ",i6," : s = ",i3," k= ", i7," N_ions(s)= ",i7)', Rank_of_process, s, k, N_ions(s)
     PRINT '("Process ",i6," : PROGRAM TERMINATED")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  IF ((k.LT.1).OR.(k.GT.N_ions(s))) THEN
     PRINT '("Process ",i6," : ERROR-2 in REMOVE_ION : index k invalid")', Rank_of_process
     PRINT '("Process ",i6," : s = ",i3," k= ", i7," N_ions(s)= ",i7)', Rank_of_process, s, k, N_ions(s)
     PRINT '("Process ",i6," : PROGRAM TERMINATED")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  IF (k.LT.N_ions(s)) THEN

     ion(s)%part(k)%X   = ion(s)%part(N_ions(s))%X
     ion(s)%part(k)%Y   = ion(s)%part(N_ions(s))%Y
     ion(s)%part(k)%VX  = ion(s)%part(N_ions(s))%VX
     ion(s)%part(k)%VY  = ion(s)%part(N_ions(s))%VY
     ion(s)%part(k)%VZ  = ion(s)%part(N_ions(s))%VZ
     ion(s)%part(k)%tag = ion(s)%part(N_ions(s))%tag

  END IF

  N_ions(s) = N_ions(s) - 1
  k = k - 1                  ! to ensure that the new k-particle is processed

END SUBROUTINE REMOVE_ION

!------------------------------------------
!
SUBROUTINE ADD_ION_TO_SEND_LEFT(s, x, y, vx, vy, vz, tag)

  USE IonParticles , ONLY : N_ions_to_send_left, max_N_ions_to_send_left, ion_to_send_left
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_X_left, L_period_x

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_ions_to_send_left(s) = N_ions_to_send_left(s) + 1

  IF (N_ions_to_send_left(s).GT.max_N_ions_to_send_left(s)) THEN
! increase the size of the list array
     current_N = max_N_ions_to_send_left(s)
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = ion_to_send_left(s)%part(k)%X
        bufer(k)%Y   = ion_to_send_left(s)%part(k)%Y
        bufer(k)%VX  = ion_to_send_left(s)%part(k)%VX
        bufer(k)%VY  = ion_to_send_left(s)%part(k)%VY
        bufer(k)%VZ  = ion_to_send_left(s)%part(k)%VZ
        bufer(k)%tag = ion_to_send_left(s)%part(k)%tag
     END DO
     IF (ALLOCATED(ion_to_send_left(s)%part)) THEN 
        DEALLOCATE(ion_to_send_left(s)%part, STAT=DEALLOC_ERR)
        !NULLIFY(ion_to_send_left(s)%part)
     END IF
     max_N_ions_to_send_left(s) = max_N_ions_to_send_left(s) + MAX(50, max_N_ions_to_send_left(s)/10)
     ALLOCATE(ion_to_send_left(s)%part(1:max_N_ions_to_send_left(s)), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        ion_to_send_left(s)%part(k)%X   = bufer(k)%X
        ion_to_send_left(s)%part(k)%Y   = bufer(k)%Y
        ion_to_send_left(s)%part(k)%VX  = bufer(k)%VX
        ion_to_send_left(s)%part(k)%VY  = bufer(k)%VY
        ion_to_send_left(s)%part(k)%VZ  = bufer(k)%VZ
        ion_to_send_left(s)%part(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  IF (periodic_boundary_X_left) THEN
     ion_to_send_left(s)%part(N_ions_to_send_left(s))%X = x + L_period_X
  ELSE
     ion_to_send_left(s)%part(N_ions_to_send_left(s))%X = x
  END IF
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%Y   = y
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%VX  = vx
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%VY  = vy
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%VZ  = vz
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%tag = tag

!print '("Process ",i4," called ADD_ION_TO_SEND_LEFT, T_cntr= ",i7)', Rank_of_process, T_cntr

END SUBROUTINE ADD_ION_TO_SEND_LEFT

!------------------------------------------
!
SUBROUTINE ADD_ION_TO_SEND_RIGHT(s, x, y, vx, vy, vz, tag)

  USE IonParticles, ONLY : N_ions_to_send_right, max_N_ions_to_send_right, ion_to_send_right
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_X_right, L_period_x

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_ions_to_send_right(s) = N_ions_to_send_right(s) + 1

  IF (N_ions_to_send_right(s).GT.max_N_ions_to_send_right(s)) THEN     !######## allocate somewhere ion_to_send_right first ######
! increase the size of the list array
     current_N = max_N_ions_to_send_right(s)
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = ion_to_send_right(s)%part(k)%X
        bufer(k)%Y   = ion_to_send_right(s)%part(k)%Y
        bufer(k)%VX  = ion_to_send_right(s)%part(k)%VX
        bufer(k)%VY  = ion_to_send_right(s)%part(k)%VY
        bufer(k)%VZ  = ion_to_send_right(s)%part(k)%VZ
        bufer(k)%tag = ion_to_send_right(s)%part(k)%tag
     END DO
     IF (ALLOCATED(ion_to_send_right(s)%part)) THEN 
        DEALLOCATE(ion_to_send_right(s)%part, STAT=DEALLOC_ERR)
        !NULLIFY(ion_to_send_right(s)%part)
     END IF
     max_N_ions_to_send_right(s) = max_N_ions_to_send_right(s) + MAX(50, max_N_ions_to_send_right(s)/10)
     ALLOCATE(ion_to_send_right(s)%part(1:max_N_ions_to_send_right(s)), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        ion_to_send_right(s)%part(k)%X   = bufer(k)%X
        ion_to_send_right(s)%part(k)%Y   = bufer(k)%Y
        ion_to_send_right(s)%part(k)%VX  = bufer(k)%VX
        ion_to_send_right(s)%part(k)%VY  = bufer(k)%VY
        ion_to_send_right(s)%part(k)%VZ  = bufer(k)%VZ
        ion_to_send_right(s)%part(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  IF (periodic_boundary_X_right) THEN
     ion_to_send_right(s)%part(N_ions_to_send_right(s))%X = x - L_period_X
  ELSE
     ion_to_send_right(s)%part(N_ions_to_send_right(s))%X = x
  END IF
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%Y   = y
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%VX  = vx
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%VY  = vy
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%VZ  = vz
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%tag = tag

!print '("Process ",i4," called ADD_ION_TO_SEND_RIGHT, T_cntr= ",i7)', Rank_of_process, T_cntr

END SUBROUTINE ADD_ION_TO_SEND_RIGHT

!------------------------------------------
!
SUBROUTINE ADD_ION_TO_SEND_ABOVE(s, x, y, vx, vy, vz, tag)

  USE IonParticles, ONLY : N_ions_to_send_above, max_N_ions_to_send_above, ion_to_send_above
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_Y_above, L_period_y

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_ions_to_send_above(s) = N_ions_to_send_above(s) + 1

  IF (N_ions_to_send_above(s).GT.max_N_ions_to_send_above(s)) THEN
! increase the size of the list array
     current_N = max_N_ions_to_send_above(s)
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = ion_to_send_above(s)%part(k)%X
        bufer(k)%Y   = ion_to_send_above(s)%part(k)%Y
        bufer(k)%VX  = ion_to_send_above(s)%part(k)%VX
        bufer(k)%VY  = ion_to_send_above(s)%part(k)%VY
        bufer(k)%VZ  = ion_to_send_above(s)%part(k)%VZ
        bufer(k)%tag = ion_to_send_above(s)%part(k)%tag
     END DO
     IF (ALLOCATED(ion_to_send_above(s)%part)) THEN 
        DEALLOCATE(ion_to_send_above(s)%part, STAT=DEALLOC_ERR)
        !NULLIFY(ion_to_send_above(s)%part)
     END IF
     max_N_ions_to_send_above(s) = max_N_ions_to_send_above(s) + MAX(50, max_N_ions_to_send_above(s)/10)
     ALLOCATE(ion_to_send_above(s)%part(1:max_N_ions_to_send_above(s)), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        ion_to_send_above(s)%part(k)%X   = bufer(k)%X
        ion_to_send_above(s)%part(k)%Y   = bufer(k)%Y
        ion_to_send_above(s)%part(k)%VX  = bufer(k)%VX
        ion_to_send_above(s)%part(k)%VY  = bufer(k)%VY
        ion_to_send_above(s)%part(k)%VZ  = bufer(k)%VZ
        ion_to_send_above(s)%part(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%X   = x
  IF (periodic_boundary_Y_above) THEN
     ion_to_send_above(s)%part(N_ions_to_send_above(s))%Y = y - L_period_y
  ELSE
     ion_to_send_above(s)%part(N_ions_to_send_above(s))%Y = y
  END IF
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%VX  = vx
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%VY  = vy
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%VZ  = vz
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%tag = tag

!print '("Process ",i4," called ADD_ION_TO_SEND_ABOVE, T_cntr= ",i7)', Rank_of_process, T_cntr

END SUBROUTINE ADD_ION_TO_SEND_ABOVE

!------------------------------------------
!
SUBROUTINE ADD_ION_TO_SEND_BELOW(s, x, y, vx, vy, vz, tag)

  USE IonParticles, ONLY : N_ions_to_send_below, max_N_ions_to_send_below, ion_to_send_below
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_Y_below, L_period_y

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_ions_to_send_below(s) = N_ions_to_send_below(s) + 1

  IF (N_ions_to_send_below(s).GT.max_N_ions_to_send_below(s)) THEN
! increase the size of the list array
     current_N = max_N_ions_to_send_below(s)
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = ion_to_send_below(s)%part(k)%X
        bufer(k)%Y   = ion_to_send_below(s)%part(k)%Y
        bufer(k)%VX  = ion_to_send_below(s)%part(k)%VX
        bufer(k)%VY  = ion_to_send_below(s)%part(k)%VY
        bufer(k)%VZ  = ion_to_send_below(s)%part(k)%VZ
        bufer(k)%tag = ion_to_send_below(s)%part(k)%tag
     END DO
     IF (ALLOCATED(ion_to_send_below(s)%part)) THEN 
        DEALLOCATE(ion_to_send_below(s)%part, STAT=DEALLOC_ERR)
        !NULLIFY(ion_to_send_below(s)%part)
     END IF
     max_N_ions_to_send_below(s) = max_N_ions_to_send_below(s) + MAX(50, max_N_ions_to_send_below(s)/10)
     ALLOCATE(ion_to_send_below(s)%part(1:max_N_ions_to_send_below(s)), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        ion_to_send_below(s)%part(k)%X   = bufer(k)%X
        ion_to_send_below(s)%part(k)%Y   = bufer(k)%Y
        ion_to_send_below(s)%part(k)%VX  = bufer(k)%VX
        ion_to_send_below(s)%part(k)%VY  = bufer(k)%VY
        ion_to_send_below(s)%part(k)%VZ  = bufer(k)%VZ
        ion_to_send_below(s)%part(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%X   = x
  IF (periodic_boundary_Y_below) THEN
     ion_to_send_below(s)%part(N_ions_to_send_below(s))%Y = y + L_period_y
  ELSE
     ion_to_send_below(s)%part(N_ions_to_send_below(s))%Y = y
  END IF
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%VX  = vx
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%VY  = vy
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%VZ  = vz
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%tag = tag

!print '("Process ",i4," called ADD_ION_TO_SEND_BELOW, T_cntr= ",i7)', Rank_of_process, T_cntr

END SUBROUTINE ADD_ION_TO_SEND_BELOW

!----------------------------------------
!
SUBROUTINE ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, tag)

  USE IonParticles, ONLY : N_ions_to_add, max_N_ions_to_add, ion_to_add

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_ions_to_add(s) = N_ions_to_add(s) + 1

  IF (N_ions_to_add(s).GT.max_N_ions_to_add(s)) THEN
! increase the size of the list array
     current_N = max_N_ions_to_add(s)
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = ion_to_add(s)%part(k)%X
        bufer(k)%Y   = ion_to_add(s)%part(k)%Y
        bufer(k)%VX  = ion_to_add(s)%part(k)%VX
        bufer(k)%VY  = ion_to_add(s)%part(k)%VY
        bufer(k)%VZ  = ion_to_add(s)%part(k)%VZ
        bufer(k)%tag = ion_to_add(s)%part(k)%tag
     END DO
     IF (ALLOCATED(ion_to_add(s)%part)) DEALLOCATE(ion_to_add(s)%part, STAT=DEALLOC_ERR)
     max_N_ions_to_add(s) = max_N_ions_to_add(s) + MAX(50, max_N_ions_to_add(s)/10)
     ALLOCATE(ion_to_add(s)%part(1:max_N_ions_to_add(s)), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        ion_to_add(s)%part(k)%X   = bufer(k)%X
        ion_to_add(s)%part(k)%Y   = bufer(k)%Y
        ion_to_add(s)%part(k)%VX  = bufer(k)%VX
        ion_to_add(s)%part(k)%VY  = bufer(k)%VY
        ion_to_add(s)%part(k)%VZ  = bufer(k)%VZ
        ion_to_add(s)%part(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  ion_to_add(s)%part(N_ions_to_add(s))%X   = x
  ion_to_add(s)%part(N_ions_to_add(s))%Y   = y
  ion_to_add(s)%part(N_ions_to_add(s))%VX  = vx
  ion_to_add(s)%part(N_ions_to_add(s))%VY  = vy
  ion_to_add(s)%part(N_ions_to_add(s))%VZ  = vz
  ion_to_add(s)%part(N_ions_to_add(s))%tag = tag

!print '("Process ",i4," called ADD_ION_TO_ADD_LIST, T_cntr= ",i7)', Rank_of_process, T_cntr
!print '("Process ",i4," : ADD_ION_TO_ADD_LIST : x/y/vx/vy/vz/tag ",5(2x,e14.7),2x,i4)', Rank_of_process, x, y, vx, vy, vz,tag

END SUBROUTINE ADD_ION_TO_ADD_LIST

!----------------------------------------
!
SUBROUTINE REMOVE_ION_FROM_ADD_LIST(s, k)

  USE ParallelOperationValues
  USE IonParticles

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(IN)    :: s
  INTEGER, INTENT(INOUT) :: k

  IF ((k.LT.1).OR.(k.GT.N_ions_to_add(s))) THEN
     PRINT '("Process ",i6," : ERROR in REMOVE_ION_FROM_ADD_LIST : index k invalid")', Rank_of_process
     PRINT '("Process ",i6," : k= ", i7," N_ions_to_add(",i2,")= ",i7)', Rank_of_process, k, s, N_ions_to_add(s)
     PRINT '("Process ",i6," : PROGRAM TERMINATED")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  IF (k.LT.N_ions_to_add(s)) THEN

     ion_to_add(s)%part(k)%X   = ion_to_add(s)%part(N_ions_to_add(s))%X
     ion_to_add(s)%part(k)%Y   = ion_to_add(s)%part(N_ions_to_add(s))%Y
     ion_to_add(s)%part(k)%VX  = ion_to_add(s)%part(N_ions_to_add(s))%VX
     ion_to_add(s)%part(k)%VY  = ion_to_add(s)%part(N_ions_to_add(s))%VY
     ion_to_add(s)%part(k)%VZ  = ion_to_add(s)%part(N_ions_to_add(s))%VZ
     ion_to_add(s)%part(k)%tag = ion_to_add(s)%part(N_ions_to_add(s))%tag

  END IF

  N_ions_to_add(s) = N_ions_to_add(s) - 1
  k = k-1                  ! to ensure that the new k-particle is processed

!print '("Process ",i4," called REMOVE_ION_FROM_ADD_LIST(k), T_cntr= ",i7," k= ",i4)', Rank_of_process, T_cntr, k

END SUBROUTINE REMOVE_ION_FROM_ADD_LIST

!----------------------------------------
! This subroutine is called after ADVANCE_IONS, before exchange of ions takes place.
! It removes particles which do not belong to the cluster domain.
! Such particles may appear after emission from surface of inner objects.
! It is expected that the emission is never directed into boundary objects aligned along the main simulation domain's boundary.
! The algorithm below, which places alien particles into proper SEND* lists is the same as in ADVANCE_IONS.
! However, collision with a boundary object at this stage is not expected and is considered an error.
!
SUBROUTINE FIND_ALIENS_IN_ION_ADD_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : N_ions_to_add, ion_to_add, N_spec

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER s, k, n 

  IF (N_of_inner_objects.EQ.0) RETURN

! cycle over ion species
  DO s = 1, N_spec
! cycle over particles of the ion species
     k=0
     DO WHILE (k.LT.N_ions_to_add(s))
        k = k + 1

        IF (symmetry_plane_X_left) THEN
           IF (ion_to_add(s)%part(k)%X.LT.c_X_area_min) THEN
! a particle crossed symmetry plane... this is an ion!!! YIKES !!!!
! since we are here, a particle to be added is beyond the symmetry plane at X=0
! the reason may be injection of ions from surfaces of inner material objects, which is not implemented
! so, presently, this should not happen, and we report it as an error
! note that for electrons this is possible (due to SEE), and such particles are moved symmetrically relative to plane x=0
!              ion_to_add(s)%part(k)%X = MAX(c_X_area_min, c_X_area_min + c_X_area_min - ion_to_add(s)%part(k)%X)
!              ion(s)%part(k)%VX = -ion(s)%part(k)%VX
              PRINT '("Proc ",i4," Error-00 in FIND_ALIENS_IN_ION_ADD_LIST, particle ",i8," of species ",i2," is beyond symmetry plane ",5(2x,e12.5),2x,i2)', Rank_of_process, k, s, &
                   & ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF
        END IF

! most probable situation when the particle is inside the area
        IF ( (ion_to_add(s)%part(k)%X.GE.c_X_area_min) .AND. &
           & (ion_to_add(s)%part(k)%X.LE.c_X_area_max) .AND. &
           & (ion_to_add(s)%part(k)%Y.GE.c_Y_area_min) .AND. &
           & (ion_to_add(s)%part(k)%Y.LE.c_Y_area_max) ) CYCLE

! since we are here, particle s,k is outside the domain of this cluster

        IF (ion_to_add(s)%part(k)%X.LT.c_X_area_min) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the ion by the period length
              ion_to_add(s)%part(k)%X = ion_to_add(s)%part(k)%X + L_period_X
              IF (ion_to_add(s)%part(k)%Y.LT.c_Y_area_min) THEN
! particle is below the bottom side of the area
                 IF (Rank_of_master_below.LT.0) THEN
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-1 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)                       
                 END IF
                 CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              ELSE IF (ion_to_add(s)%part(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
                 IF (Rank_of_master_above.LT.0) THEN
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-2 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                 END IF
                 CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              END IF
              CYCLE
           END IF

           IF ( (ion_to_add(s)%part(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
              & (ion_to_add(s)%part(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

              IF (Rank_of_master_left.GE.0) THEN
! left neighbor cluster exists
                 CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
              ELSE
! left neighbor cluster does not exist
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX,  ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)   ! left
! error
                 print '("error-3 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF

           ELSE IF (ion_to_add(s)%part(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the left bottom cell of the area

              SELECT CASE (c_left_bottom_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((c_X_area_min-ion_to_add(s)%part(k)%X).LT.(c_Y_area_min-ion_to_add(s)%part(k)%Y)) THEN
                       CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)

                 CASE (FLAT_WALL_LEFT)
                    IF (ion_to_add(s)%part(k)%Y.GE.c_Y_area_min) THEN                 
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                       print '("error-4 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                    ELSE
                       CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)                       
                    END IF

                 CASE (SURROUNDED_BY_WALL)
!                    IF ((c_X_area_min-ion_to_add(s)%part(k)%X).LT.(c_Y_area_min-ion_to_add(s)%part(k)%Y)) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    ELSE
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    END IF
! error
                    print '("error-5 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

                 CASE (EMPTY_CORNER_WALL_LEFT)
                    CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
               
                 CASE (EMPTY_CORNER_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)

              END SELECT

           ELSE IF (ion_to_add(s)%part(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top left cell of the area

              SELECT CASE (c_left_top_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((c_X_area_min-ion_to_add(s)%part(k)%X).LT.(ion_to_add(s)%part(k)%Y-c_Y_area_max)) THEN
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)

                 CASE (FLAT_WALL_LEFT)
                    IF (ion_to_add(s)%part(k)%Y.LE.c_Y_area_max) THEN                 
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                       print '("error-6 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                    ELSE
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                    END IF

                 CASE (SURROUNDED_BY_WALL)
!                    IF ((c_X_area_min-ion_to_add(s)%part(k)%X).LT.(ion_to_add(s)%part(k)%Y-c_Y_area_max)) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    ELSE
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    END IF
! error
                    print '("error-7 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

                 CASE (EMPTY_CORNER_WALL_LEFT)
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)

                 CASE (EMPTY_CORNER_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)

              END SELECT
           END IF
           CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE
        END IF   !### IF (ion_to_add(s)%part(k)%X.LT.c_X_area_min) THEN

        IF (ion_to_add(s)%part(k)%X.GT.c_X_area_max) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the ion by the period length
              ion_to_add(s)%part(k)%X = ion_to_add(s)%part(k)%X - L_period_X
              IF (ion_to_add(s)%part(k)%Y.LT.c_Y_area_min) THEN
! particle is somewhere near the left bottom cell of the area
                 IF (Rank_of_master_below.LT.0) THEN
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-8 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)                       
                 END IF
                 CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              ELSE IF (ion_to_add(s)%part(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
                 IF (Rank_of_master_above.LT.0) THEN
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-9 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                 END IF
                 CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              END IF
              CYCLE
           END IF
           
           IF ( (ion_to_add(s)%part(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
              & (ion_to_add(s)%part(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

              IF (Rank_of_master_right.GE.0) THEN
! right neighbor cluster exists
                 CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
              ELSE
! right neighbor cluster does not exist
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                 print '("error-10 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF

           ELSE IF (ion_to_add(s)%part(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the right bottom cell of the area

              SELECT CASE (c_right_bottom_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((ion_to_add(s)%part(k)%X-c_X_area_max).LT.(c_Y_area_min-ion_to_add(s)%part(k)%Y)) THEN
                       CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                
                 CASE (FLAT_WALL_RIGHT)
                    IF (ion_to_add(s)%part(k)%Y.GE.c_Y_area_min) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                       print '("error-11 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

                    ELSE
                       CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)                       
                    END IF

                 CASE (SURROUNDED_BY_WALL)
!                    IF ((ion_to_add(s)%part(k)%X-c_X_area_max).LT.(c_Y_area_min-ion_to_add(s)%part(k)%Y)) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    ELSE
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    END IF
! error
                       print '("error-12 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

                 CASE (EMPTY_CORNER_WALL_RIGHT)
                    CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)

                 CASE (EMPTY_CORNER_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)

              END SELECT

           ELSE IF (ion_to_add(s)%part(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top right cell of the area

              SELECT CASE (c_right_top_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((ion_to_add(s)%part(k)%X-c_X_area_max).LT.(ion_to_add(s)%part(k)%Y-c_Y_area_max)) THEN
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)

                 CASE (FLAT_WALL_RIGHT)
                    IF (ion_to_add(s)%part(k)%Y.LE.c_Y_area_max) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                       print '("error-13 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                    ELSE
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)                    
                    END IF

                 CASE (SURROUNDED_BY_WALL)
!                    IF ((ion_to_add(s)%part(k)%X-c_X_area_max).LT.(ion_to_add(s)%part(k)%Y-c_Y_area_max)) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    ELSE
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    END IF
! error
                       print '("error-14 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

                 CASE (EMPTY_CORNER_WALL_RIGHT)
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)

                 CASE (EMPTY_CORNER_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)

              END SELECT

           END IF
           CALL REMOVE_ION_FROM_ADD_LIST(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE
        END IF        !###  IF (ion_to_add(s)%part(k)%X.GT.c_X_area_max) THEN

! since we are here, c_X_area_min <= ion_to_add(s)&part(k)%Y <= c_X_area_max

        IF (ion_to_add(s)%part(k)%Y.GT.c_Y_area_max) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
              IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
                 CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
              ELSE
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                 print '("error-15 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
              CALL REMOVE_ION_FROM_ADD_LIST(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              CYCLE
           END IF

           IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
              CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
           ELSE
! neighbor cluster above does not exist
              IF ((ion_to_add(s)%part(k)%X.GE.(c_X_area_min+1.0_8)).AND.(ion_to_add(s)%part(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                 print '("error-16 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE IF (ion_to_add(s)%part(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left top corner
                 IF (c_left_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                 ELSE
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-17 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 END IF
              ELSE IF (ion_to_add(s)%part(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right top corner
                 IF (c_right_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX,  ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                 ELSE
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-18 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 END IF
              END IF
           END IF
           CALL REMOVE_ION_FROM_ADD_LIST(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE
        END IF    !### IF (ion_to_add(s)%part(k)%Y.GT.c_Y_area_max) THEN

        IF (ion_to_add(s)%part(k)%Y.LT.c_Y_area_min) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
              IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster above exists
                 CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
              ELSE
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                 print '("error-19 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
              CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              CYCLE
           END IF

           IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster below exists, remove particle and prepare to send it to the neighbor below
              CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
           ELSE
! neighbor cluster below does not exist
              IF ((ion_to_add(s)%part(k)%X.GE.(c_X_area_min+1.0_8)).AND.(ion_to_add(s)%part(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                 print '("error-20 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE IF (ion_to_add(s)%part(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left bottom corner
                 IF (c_left_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                 ELSE
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-21 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 END IF
              ELSE IF (ion_to_add(s)%part(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right bottom corner
                 IF (c_right_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
                 ELSE
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-22 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                END IF
              END IF
           END IF
           CALL REMOVE_ION_FROM_ADD_LIST(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE
        END IF   !### IF (ion_to_add(s)%part(k)%Y.LT.c_Y_area_min) THEN

     END DO   !### DO WHILE (k.LT.N_ions_to_add(s))
  END DO    !### DO s = 1, N_spec

END SUBROUTINE FIND_ALIENS_IN_ION_ADD_LIST

!----------------------------------------
!
SUBROUTINE FIND_INNER_OBJECT_COLL_IN_ION_ADD_LIST

!  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles, ONLY : N_ions_to_add, ion_to_add, N_spec

  IMPLICIT NONE

  INTEGER s, k, n 

  IF (N_of_inner_objects.EQ.0) RETURN

  DO s = 1, N_spec
  
! find, process, and exclude ions which collided with inner objects
     k=0
     DO WHILE (k.LT.N_ions_to_add(s))
        k = k+1
        DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (ion_to_add(s)%part(k)%X.LE.whole_object(n)%Xmin) CYCLE
           IF (ion_to_add(s)%part(k)%X.GE.whole_object(n)%Xmax) CYCLE
           IF (ion_to_add(s)%part(k)%Y.LE.whole_object(n)%Ymin) CYCLE
           IF (ion_to_add(s)%part(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! collision detected
           CALL TRY_ION_COLL_WITH_INNER_OBJECT( s, &
                                              & ion_to_add(s)%part(k)%X, &
                                              & ion_to_add(s)%part(k)%Y, &
                                              & ion_to_add(s)%part(k)%VX, &
                                              & ion_to_add(s)%part(k)%VY, &
                                              & ion_to_add(s)%part(k)%VZ, &
                                              & ion_to_add(s)%part(k)%tag )  !, &
!                                                  & whole_object(n) )
           CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions_to_add(s) = N_ions_to_add(s) - 1 and k = k-1
           EXIT
        END DO
     END DO

  END DO

END SUBROUTINE FIND_INNER_OBJECT_COLL_IN_ION_ADD_LIST

!----------------------------------------
!
SUBROUTINE PROCESS_ADDED_IONS

  USE ParallelOperationValues
  USE CurrentProblemValues

  USE IonParticles, ONLY : N_ions_to_add, N_ions, max_N_ions, ion, ion_to_add, N_spec

  IMPLICIT NONE

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER s, k, current_N

  DO s = 1, N_spec
  
     IF (N_ions_to_add(s).GT.(max_N_ions(s)-N_ions(s))) THEN
! increase the size of the main ion array
        current_N = max_N_ions(s)
        ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
        DO k = 1, current_N
           bufer(k)%X   = ion(s)%part(k)%X
           bufer(k)%Y   = ion(s)%part(k)%Y
           bufer(k)%VX  = ion(s)%part(k)%VX
           bufer(k)%VY  = ion(s)%part(k)%VY
           bufer(k)%VZ  = ion(s)%part(k)%VZ
           bufer(k)%tag = ion(s)%part(k)%tag
        END DO
        DEALLOCATE(ion(s)%part, STAT=DEALLOC_ERR)
        !NULLIFY(ion(s)%part)
        max_N_ions(s) = max_N_ions(s) + MAX(N_ions_to_add(s)-(max_N_ions(s)-N_ions(s)), max_N_ions(s)/10)
        ALLOCATE(ion(s)%part(1:max_N_ions(s)), STAT=ALLOC_ERR)
        DO k = 1, current_N
           ion(s)%part(k)%X   = bufer(k)%X
           ion(s)%part(k)%Y   = bufer(k)%Y
           ion(s)%part(k)%VX  = bufer(k)%VX
           ion(s)%part(k)%VY  = bufer(k)%VY
           ion(s)%part(k)%VZ  = bufer(k)%VZ
           ion(s)%part(k)%tag = bufer(k)%tag
        END DO
        DEALLOCATE(bufer, STAT=DEALLOC_ERR)
     END IF
  
     DO k = 1, N_ions_to_add(s)
        ion(s)%part(k+N_ions(s))%X   = ion_to_add(s)%part(k)%X
        ion(s)%part(k+N_ions(s))%Y   = ion_to_add(s)%part(k)%Y
        ion(s)%part(k+N_ions(s))%VX  = ion_to_add(s)%part(k)%VX
        ion(s)%part(k+N_ions(s))%VY  = ion_to_add(s)%part(k)%VY
        ion(s)%part(k+N_ions(s))%VZ  = ion_to_add(s)%part(k)%VZ
        ion(s)%part(k+N_ions(s))%tag = ion_to_add(s)%part(k)%tag
     END DO

! update electron counter
     N_ions(s) = N_ions(s) + N_ions_to_add(s)

  END DO

! clear counter of particles to be added to this process
  N_ions_to_add = 0  

END SUBROUTINE PROCESS_ADDED_IONS
!------------------------------
!
SUBROUTINE GATHER_ION_CHARGE_DENSITY

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries, ONLY : indx_x_min, indx_x_max, indx_y_min, indx_y_max
!  USE Diagnostics

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

!  REAL(8), ALLOCATABLE :: c_rho_i(:,:)

  INTEGER n1  ! number of nodes in the y-direction
  INTEGER n2  !
  INTEGER n3  ! number of nodes in the x-direction

  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR
  INTEGER bufsize

!  INTEGER npc ! n-umber of p-robe in c-luster
!  INTEGER npa ! n-umber of p-robe a-ll [in the global list of probes]

  INTEGER s, i, j, k
  INTEGER pos_i_j, pos_ip1_j, pos_i_jp1, pos_ip1_jp1
  REAL(8) ax_ip1, ax_i, ay_jp1, ay_j
  REAL(8) vij, vip1j, vijp1

  INTEGER pos

  INTEGER nio, position_flag

  IF ( (cluster_rank_key.NE.0) &                     ! this array is needed only temporarily for non-master processes
     & .OR. &
     & (periodicity_flag.EQ.PERIODICITY_NONE) &      ! if PETSc is used, master processes don't need this array outside this sub
     & .OR. &
     & (periodicity_flag.EQ.PERIODICITY_X_PETSC) &
     & .OR. &
     & (periodicity_flag.EQ.PERIODICITY_X_Y) ) &
     & THEN
     ALLOCATE(c_rho_i(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
  END IF

  n1 = c_indx_y_max - c_indx_y_min + 1
! was  n2 = -c_indx_x_min + 1 - c_indx_y_min * n1
  n3 = c_indx_x_max - c_indx_x_min + 1
  n2 = -c_indx_x_min + 1 - c_indx_y_min * n3

  bufsize = n1 * n3
  ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

  rbufer = 0.0_8

!!  IF (N_of_probes_cluster.GT.0) probe_Ni_cluster = 0.0_8

  DO s = 1, N_spec
     DO k = 1, N_ions(s)
     
        i = INT(ion(s)%part(k)%X)
        j = INT(ion(s)%part(k)%Y)
        IF (ion(s)%part(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
        IF (ion(s)%part(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

        if ((i.lt.c_indx_x_min).or.(i.gt.(c_indx_x_max-1)).or.(j.lt.c_indx_y_min).or.(j.gt.(c_indx_y_max-1))) then
           print '("Process ",i4," : Error-1 in GATHER_ION_CHARGE_DENSITY : index out of bounds")', Rank_of_process
           print '("Process ",i4," : k/s/N_ions(s) : ",i8,2x,i8)', Rank_of_process, k, s, N_ions(s)
           print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag
           print '("Process ",i4," : minx/maxx/miny/maxy : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        end if

!     pos = i - c_indx_x_min + 1 + (j - c_indx_y_min) * (c_indx_x_max - c_indx_x_min + 1)

        pos_i_j     = i + j * n3 + n2
        pos_ip1_j   = pos_i_j + 1
        pos_i_jp1   = pos_i_j + n3
        pos_ip1_jp1 = pos_i_jp1 + 1

        ax_ip1 = ion(s)%part(k)%X - DBLE(i)
        ax_i   = 1.0_8 - ax_ip1

        ay_jp1 = ion(s)%part(k)%Y - DBLE(j)
        ay_j = 1.0_8 - ay_jp1

        vij   = ax_i   * ay_j
        vip1j = ax_ip1 * ay_j
        vijp1 = ax_i   * ay_jp1

!! probe diagnostics
!        DO npc = 1, N_of_probes_cluster
!           npa = List_of_probes_cluster(npc)
!           IF (i.EQ.Probe_position(1,npa)) THEN
!              IF (j.EQ.Probe_position(2,npa)) THEN
!                 probe_Ni_cluster(npc, s) = probe_Ni_cluster(npc, s) + vij
!              ELSE IF ((j+1).EQ.Probe_position(2,npa)) THEN
!                 probe_Ni_cluster(npc, s) = probe_Ni_cluster(npc, s) + vijp1
!              END IF
!           ELSE IF ((i+1).EQ.Probe_position(1,npa)) THEN
!              IF (j.EQ.Probe_position(2,npa)) THEN
!                 probe_Ni_cluster(npc, s) = probe_Ni_cluster(npc, s) + vip1j
!              ELSE IF ((j+1).EQ.Probe_position(2,npa)) THEN
!                 probe_Ni_cluster(npc, s) = probe_Ni_cluster(npc, s) + 1.0_8 - vij - vip1j - vijp1
!              END IF
!           END IF
!        END DO

        SELECT CASE(Qs(s))
! assume that there may be positive and negative ions with charges up to +/-3e
! thus we use addition instead of multiplication which is faster, right?
           CASE(1)
              rbufer(pos_i_j)     = rbufer(pos_i_j)     + vij                           !ax_i   * ay_j
              rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   + vip1j                         !ax_ip1 * ay_j
              rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   + vijp1                         !ax_i   * ay_jp1
              rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) + 1.0_8 - vij - vip1j - vijp1   !ax_ip1 * ay_jp1
           CASE(2)
              rbufer(pos_i_j)     = rbufer(pos_i_j)     + vij + vij                                           !ax_i   * ay_j
              rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   + vip1j + vip1j                                       !ax_ip1 * ay_j
              rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   + vijp1 + vijp1                                       !ax_i   * ay_jp1
              rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) + 2.0_8 - vij - vip1j - vijp1 - vij - vip1j - vijp1   !ax_ip1 * ay_jp1
           CASE(3)
              rbufer(pos_i_j)     = rbufer(pos_i_j)     + vij + vij + vij                                                           !ax_i   * ay_j
              rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   + vip1j + vip1j + vip1j                                                     !ax_ip1 * ay_j
              rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   + vijp1 + vijp1 + vijp1                                                     !ax_i   * ay_jp1
              rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) + 3.0_8 - vij - vip1j - vijp1 - vij - vip1j - vijp1 - vij - vip1j - vijp1   !ax_ip1 * ay_jp1
           CASE(-1)
              rbufer(pos_i_j)     = rbufer(pos_i_j)     - vij                           !ax_i   * ay_j
              rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   - vip1j                         !ax_ip1 * ay_j
              rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   - vijp1                         !ax_i   * ay_jp1
              rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) - 1.0_8 + vij + vip1j + vijp1   !ax_ip1 * ay_jp1
           CASE(-2)
              rbufer(pos_i_j)     = rbufer(pos_i_j)     - vij - vij                                           !ax_i   * ay_j
              rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   - vip1j - vip1j                                       !ax_ip1 * ay_j
              rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   - vijp1 - vijp1                                       !ax_i   * ay_jp1
              rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) - 2.0_8 + vij + vip1j + vijp1 + vij + vip1j + vijp1   !ax_ip1 * ay_jp1
           CASE(-3)
              rbufer(pos_i_j)     = rbufer(pos_i_j)     - vij - vij - vij                                                           !ax_i   * ay_j
              rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   - vip1j - vip1j - vip1j                                                     !ax_ip1 * ay_j
              rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   - vijp1 - vijp1 - vijp1                                                     !ax_i   * ay_jp1
              rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) - 3.0_8 + vij + vip1j + vijp1 + vij + vip1j + vijp1 + vij + vip1j + vijp1   !ax_ip1 * ay_jp1
        END SELECT

     END DO
  END DO

! collect densities from all processes in a cluster
  CALL MPI_REDUCE(rbufer, c_rho_i, bufsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_CLUSTER, ierr)
  
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! now cluster masters exchange information about densities in overlapping nodes
  IF (cluster_rank_key.EQ.0) THEN

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! special case of self-connected X-periodic cluster
        DO j = c_indx_y_min, c_indx_y_max
           c_rho_i(c_indx_x_min+1, j) = c_rho_i(c_indx_x_min+1, j) + c_rho_i(c_indx_x_max, j) 
           c_rho_i(c_indx_x_max-1, j) = c_rho_i(c_indx_x_max-1, j) + c_rho_i(c_indx_x_min, j)
           c_rho_i(c_indx_x_min, j) = c_rho_i(c_indx_x_max-1, j)
           c_rho_i(c_indx_x_max, j) = c_rho_i(c_indx_x_min+1, j)
        END DO
     END IF

     IF (WHITE_CLUSTER) THEN  
! "white processes"

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

        IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right densities in the right edge
           rbufer(1:n1) = c_rho_i(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left densities in the left edge
           rbufer(1:n1) = c_rho_i(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left densities in the vertical line next to the left edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho_i(c_indx_x_min+1, j) = c_rho_i(c_indx_x_min+1, j) + rbufer(j-c_indx_y_min+1)
           END DO
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right densities in the vertical line next to the right edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho_i(c_indx_x_max-1, j) = c_rho_i(c_indx_x_max-1, j) + rbufer(j-c_indx_y_min+1)
           END DO           
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

        IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up densities in the top edge
           rbufer(1:n3) = c_rho_i(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down densities in the bottom edge
           rbufer(1:n3) = c_rho_i(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below densities in the vertical line above the bottom line
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho_i(i, c_indx_y_min+1) = c_rho_i(i, c_indx_y_min+1) + rbufer(i-c_indx_x_min+1)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above densities in the vertical line under the top line
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho_i(i, c_indx_y_max-1) = c_rho_i(i, c_indx_y_max-1) + rbufer(i-c_indx_x_min+1)
           END DO
        END IF

     ELSE
! "black" processes

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

        IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left densities in the vertical line next to the left edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho_i(c_indx_x_min+1, j) = c_rho_i(c_indx_x_min+1, j) + rbufer(j-c_indx_y_min+1)
           END DO
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right densities in the vertical line next to the right edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho_i(c_indx_x_max-1, j) = c_rho_i(c_indx_x_max-1, j) + rbufer(j-c_indx_y_min+1)
           END DO           
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right densities in the right edge
           rbufer(1:n1) = c_rho_i(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left densities in the left edge
           rbufer(1:n1) = c_rho_i(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

        IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below densities in the vertical line above the bottom line
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho_i(i, c_indx_y_min+1) = c_rho_i(i, c_indx_y_min+1) + rbufer(i-c_indx_x_min+1)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above densities in the vertical line under the top line
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho_i(i, c_indx_y_max-1) = c_rho_i(i, c_indx_y_max-1) + rbufer(i-c_indx_x_min+1)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up densities in the top edge
           rbufer(1:n3) = c_rho_i(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down densities in the bottom edge
           rbufer(1:n3) = c_rho_i(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

     END IF

! adjust densities at the boundaries with material walls

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! special case of self-connected X-periodic cluster
        IF (Rank_of_master_above.LT.0) THEN
           DO i = c_indx_x_min, c_indx_x_max
              c_rho_i(i, c_indx_y_max) = 2.0_8 * c_rho_i(i, c_indx_y_max)
           END DO
        END IF
  
        IF (Rank_of_master_below.LT.0) THEN
           DO i = c_indx_x_min, c_indx_x_max
              c_rho_i(i, c_indx_y_min) = 2.0_8 * c_rho_i(i, c_indx_y_min)
           END DO
        END IF

     ELSE
 
        IF (Rank_of_master_left.LT.0) THEN
           DO j = c_indx_y_min+1, c_indx_y_max-1
              c_rho_i(c_indx_x_min, j) = 2.0_8 * c_rho_i(c_indx_x_min, j)
           END DO
        END IF

        IF (Rank_of_master_right.LT.0) THEN
           DO j = c_indx_y_min+1, c_indx_y_max-1
              c_rho_i(c_indx_x_max, j) = 2.0_8 * c_rho_i(c_indx_x_max, j)
           END DO
        END IF
  
        IF (Rank_of_master_above.LT.0) THEN
           DO i = c_indx_x_min+1, c_indx_x_max-1
              c_rho_i(i, c_indx_y_max) = 2.0_8 * c_rho_i(i, c_indx_y_max)
           END DO
        END IF
  
        IF (Rank_of_master_below.LT.0) THEN
           DO i = c_indx_x_min+1, c_indx_x_max-1
              c_rho_i(i, c_indx_y_min) = 2.0_8 * c_rho_i(i, c_indx_y_min)
           END DO
        END IF

        SELECT CASE (c_left_bottom_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho_i(c_indx_x_min, c_indx_y_min) = 4.0_8 * c_rho_i(c_indx_x_min, c_indx_y_min) 
           CASE (EMPTY_CORNER_WALL_LEFT)
              c_rho_i(c_indx_x_min, c_indx_y_min+1) = 0.66666666666666_8 * c_rho_i(c_indx_x_min, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_BELOW)
              c_rho_i(c_indx_x_min+1, c_indx_y_min) = 0.66666666666666_8 * c_rho_i(c_indx_x_min+1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
        END SELECT
     
        SELECT CASE (c_left_top_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho_i(c_indx_x_min, c_indx_y_max) = 4.0_8 * c_rho_i(c_indx_x_min, c_indx_y_max) 
           CASE (EMPTY_CORNER_WALL_LEFT)
              c_rho_i(c_indx_x_min, c_indx_y_max-1) = 0.66666666666666_8 * c_rho_i(c_indx_x_min, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_ABOVE)
              c_rho_i(c_indx_x_min+1, c_indx_y_max) = 0.66666666666666_8 * c_rho_i(c_indx_x_min+1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
        END SELECT

        SELECT CASE (c_right_bottom_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho_i(c_indx_x_max, c_indx_y_min) = 4.0_8 * c_rho_i(c_indx_x_max, c_indx_y_min) 
           CASE (EMPTY_CORNER_WALL_RIGHT)
              c_rho_i(c_indx_x_max, c_indx_y_min+1) = 0.66666666666666_8 * c_rho_i(c_indx_x_max, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_BELOW)
              c_rho_i(c_indx_x_max-1, c_indx_y_min) = 0.66666666666666_8 * c_rho_i(c_indx_x_max-1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
        END SELECT
     
        SELECT CASE (c_right_top_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho_i(c_indx_x_max, c_indx_y_max) = 4.0_8 * c_rho_i(c_indx_x_max, c_indx_y_max) 
           CASE (EMPTY_CORNER_WALL_RIGHT)
              c_rho_i(c_indx_x_max, c_indx_y_max-1) = 0.66666666666666_8 * c_rho_i(c_indx_x_max, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_ABOVE)
              c_rho_i(c_indx_x_max-1, c_indx_y_max) = 0.66666666666666_8 * c_rho_i(c_indx_x_max-1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
        END SELECT

     END IF

  END IF

  IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

     IF (cluster_rank_key.EQ.0) THEN

! prepare and send charge density to field calculators
        DO k = 2, cluster_N_blocks
           bufsize = (field_calculator(k)%indx_x_max - field_calculator(k)%indx_x_min + 1) * &
                   & (field_calculator(k)%indx_y_max - field_calculator(k)%indx_y_min + 1)
           IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
           ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)
           pos=0
           DO j = field_calculator(k)%indx_y_min, field_calculator(k)%indx_y_max
              DO i = field_calculator(k)%indx_x_min, field_calculator(k)%indx_x_max
                 pos = pos+1
                 rbufer(pos) = c_rho_i(i,j)

!                 CALL FIND_INNER_OBJECT_CONTAINING_POINT(i,j,nio,position_flag)
!                 SELECT CASE (position_flag)
!                    CASE (1,3,5,7)
!! point at the corner of a dielectric object
!                       rbufer(pos) = (1.0_8) * c_rho_i(i,j)
!                    CASE (2,4,6,8)
!! point on the surface of a dielectric object
!! the factor used here allows to use common factor -1/4*N_of_particles_cell in SOLVE_POTENTIAL_WITH_PETSC
!! it also allows to use uncorrected values of volume charge density in nodes on the flat material surface
!                       rbufer(pos) = 2.0_8 * c_rho_i(i,j) / (1.0_8 + whole_object(nio)%eps_diel)
!                 END SELECT

              END DO
           END DO
           CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
        END DO

! cluster master is a field calaculator too, prepare its own charge density
        DO j = indx_y_min, indx_y_max
!           rho_i(indx_x_min:indx_x_max, j) = c_rho_i(indx_x_min:indx_x_max, j)
           DO i = indx_x_min, indx_x_max
              rho_i(i, j) = c_rho_i(i, j)

!              CALL FIND_INNER_OBJECT_CONTAINING_POINT(i,j,nio,position_flag)
!              SELECT CASE (position_flag)
!                 CASE (1,3,5,7)
!! point at the corner of a dielectric object
!                    rho_i(i, j) = (1.0_8) * c_rho_i(i,j)
!                 CASE (2,4,6,8)
!! point on the surface of a dielectric object
!! the factor used here allows to use common factor -1/4*N_of_particles_cell in SOLVE_POTENTIAL_WITH_PETSC
!! it also allows to use uncorrected values of volume charge density in nodes on the flat material surface
!                    rho_i(i, j) = 2.0_8 * c_rho_i(i,j) / (1.0_8 + whole_object(nio)%eps_diel)
!              END SELECT

           END DO
        END DO

     ELSE

        bufsize = (indx_x_max - indx_x_min + 1) * (indx_y_max - indx_y_min + 1)
        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)
        
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, field_master, field_master, MPI_COMM_WORLD, stattus, ierr)
        
        pos = 0
        DO j = indx_y_min, indx_y_max
           DO i = indx_x_min, indx_x_max
              pos = pos+1
              rho_i(i,j) = rbufer(pos)
           END DO
        END DO
     
     END IF

  END IF
  
  IF (ALLOCATED(rbufer))  DEALLOCATE(rbufer,  STAT=ALLOC_ERR)

  IF ( (cluster_rank_key.NE.0).OR.(periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN
     IF (ALLOCATED(c_rho_i)) DEALLOCATE(c_rho_i, STAT=ALLOC_ERR)
  END IF

END SUBROUTINE GATHER_ION_CHARGE_DENSITY
