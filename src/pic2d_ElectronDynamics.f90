!---------------------------------------
!
SUBROUTINE ADVANCE_ELECTRONS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER k
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
  N_e_to_send_left = 0
  N_e_to_send_right = 0
  N_e_to_send_above = 0
  N_e_to_send_below = 0

! clear counters of particles that hit boundary objects
  whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count = 0
  e_colls_with_bo(1:N_of_boundary_and_inner_objects)%N_of_saved_parts = 0

! clear counters of particles emitted by boundary objects in order to account for the secondary electron emission
  whole_object(1:N_of_boundary_and_inner_objects)%electron_emit_count = 0  !### ?????

!print *, "enter", Rank_of_process, N_electrons, max_N_electrons

  k=0
  DO WHILE (k.LT.N_electrons)

     k = k + 1

     if ( (electron(k)%X.lt.c_X_area_min).or. &
        & (electron(k)%X.gt.c_X_area_max).or. &
        & (electron(k)%Y.lt.c_Y_area_min).or. &
        & (electron(k)%Y.gt.c_Y_area_max) ) then
        print '("Process ",i4," : Error-1 in ADVANCE_ELECTRONS : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
        print '("Process ",i4," : k/N_electrons : ",i8,2x,i8)', Rank_of_process, k, N_electrons
        print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     end if

! interpolate electric field

     i = INT(electron(k)%X)
     j = INT(electron(k)%Y)

     IF (electron(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (electron(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

if ((i.lt.c_indx_x_min).or.(i.gt.(c_indx_x_max-1)).or.(j.lt.c_indx_y_min).or.(j.gt.(c_indx_y_max-1))) then
   print '("Process ",i4," : Error-2 in ADVANCE_ELECTRONS : index out of bounds")', Rank_of_process
   print '("Process ",i4," : k/N_electrons : ",i8,2x,i8)', Rank_of_process, k, N_electrons
   print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag
   print '("Process ",i4," : minx/maxx/miny/maxy : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
   CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
end if

     ax_ip1 = electron(k)%X - DBLE(i)
     ax_i   = 1.0_8 - ax_ip1

     ay_jp1 = electron(k)%Y - DBLE(j)
     ay_j = 1.0_8 - ay_jp1

     E_X = EX(i,j) * ax_i * ay_j + EX(i+1,j) * ax_ip1 * ay_j + EX(i,j+1) * ax_i * ay_jp1 + EX(i+1,j+1) * ax_ip1 * ay_jp1
     E_Y = EY(i,j) * ax_i * ay_j + EY(i+1,j) * ax_ip1 * ay_j + EY(i,j+1) * ax_i * ay_jp1 + EY(i+1,j+1) * ax_ip1 * ay_jp1
     E_Z = Ez(electron(k)%X, electron(k)%Y)

! calculate magnetic field factors

     alfa_x = -0.5_8 * Bx(electron(k)%X, electron(k)%Y)
     alfa_y = -0.5_8 * By(electron(k)%X, electron(k)%Y)
     alfa_z = -0.5_8 * Bz(electron(k)%X, electron(k)%Y)

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

     VX_minus = electron(k)%VX - 0.5_8 * E_X
     VY_minus = electron(k)%VY - 0.5_8 * E_Y
     VZ_minus = electron(k)%VZ - 0.5_8 * E_Z

! velocity advance: rotation in the magnetic field

     VX_plus = K11 * VX_minus + K12 * VY_minus + K13 * VZ_minus
     VY_plus = K21 * VX_minus + K22 * VY_minus + K23 * VZ_minus
     VZ_plus = K31 * VX_minus + K32 * VY_minus + K33 * VZ_minus

! velocity advance: second half-acceleration due to electric field

     electron(k)%VX = VX_plus - 0.5_8 * E_X
     electron(k)%VY = VY_plus - 0.5_8 * E_Y
     electron(k)%VZ = VZ_plus - 0.5_8 * E_Z

! coordinate advance

!###     IF (ABS(electron(k)%VX).GT.1.0_8) THEN
!        PRINT '("Error in process ",i4," too high VX of electron ",i9," :: X/Y/VX/VY/VZ are ",5(2x,e12.5))', &
!             & Rank_of_process, k, electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ
!        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
!###     END IF

!###     IF(ABS(electron(k)%VY).GT.1.0_8) THEN
!        PRINT '("Error in process ",i4," too high VY of electron ",i9," :: X/Y/VX/VY/VZ are ",5(2x,e12.5))', &
!             & Rank_of_process, k, electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ
!        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
!###     END IF

     electron(k)%X = electron(k)%X + electron(k)%VX
     electron(k)%Y = electron(k)%Y + electron(k)%VY

! a particle crossed symmetry plane, reflect it
     IF (symmetry_plane_X_left) THEN
        IF (electron(k)%X.LT.c_X_area_min) THEN
           electron(k)%X = MAX(c_X_area_min, c_X_area_min + c_X_area_min - electron(k)%X)
           electron(k)%VX = -electron(k)%VX
!###          electron(k)%VZ = -electron(k)%VZ   !###??? do we not have to do this when BY is on ??? 
        END IF
     END IF

! check whether a collision with an inner object occurred
     collision_with_inner_object_occurred = .FALSE.
     DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
        IF (electron(k)%X.LE.whole_object(n)%Xmin) CYCLE
        IF (electron(k)%X.GE.whole_object(n)%Xmax) CYCLE
        IF (electron(k)%Y.LE.whole_object(n)%Ymin) CYCLE
        IF (electron(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! collision detected
        CALL TRY_ELECTRON_COLL_WITH_INNER_OBJECT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag) !, whole_object(n))
        CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
        collision_with_inner_object_occurred = .TRUE.
        EXIT
     END DO

     IF (collision_with_inner_object_occurred) CYCLE

! most probable situation when the particle remains inside the area
     IF ( (electron(k)%X.GE.c_X_area_min) .AND. &
        & (electron(k)%X.LE.c_X_area_max) .AND. &
        & (electron(k)%Y.GE.c_Y_area_min) .AND. &
        & (electron(k)%Y.LE.c_Y_area_max) ) CYCLE

! since we are here, a particle did not collide with an inner object but crossed an area boundary
! note, in a periodic system the particle still may collide with an inner object after transfer to the other domain
! therefore it is still necessary to check add list for collisions after process receives particles from neighbors

     IF (electron(k)%X.LT.c_X_area_min) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the electron by the period length
           electron(k)%X = electron(k)%X + L_period_X
           IF (electron(k)%Y.LT.c_Y_area_min) THEN
! particle is somewhere near the left bottom cell of the area
              IF (Rank_of_master_below.LT.0) THEN
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)                       
              END IF
              CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
           ELSE IF (electron(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
              IF (Rank_of_master_above.LT.0) THEN
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              END IF
              CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
           END IF
           CYCLE
        END IF

        IF ( (electron(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
           & (electron(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

           IF (Rank_of_master_left.GE.0) THEN
! left neighbor cluster exists
              CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           ELSE
! left neighbor cluster does not exist
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX,  electron(k)%VY, electron(k)%VZ, electron(k)%tag)   ! left
           END IF

        ELSE IF (electron(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the left bottom cell of the area

           SELECT CASE (c_left_bottom_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((c_X_area_min-electron(k)%X).LT.(c_Y_area_min-electron(k)%Y)) THEN
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 END IF

              CASE (FLAT_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)

              CASE (FLAT_WALL_LEFT)
                 IF (electron(k)%Y.GE.c_Y_area_min) THEN                 
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)                       
                 END IF

              CASE (SURROUNDED_BY_WALL)
                 IF ((c_X_area_min-electron(k)%X).LT.(c_Y_area_min-electron(k)%Y)) THEN
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 END IF

              CASE (EMPTY_CORNER_WALL_LEFT)
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
               
              CASE (EMPTY_CORNER_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)

           END SELECT

        ELSE IF (electron(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top left cell of the area

           SELECT CASE (c_left_top_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((c_X_area_min-electron(k)%X).LT.(electron(k)%Y-c_Y_area_max)) THEN
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 END IF

              CASE (FLAT_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)

              CASE (FLAT_WALL_LEFT)
                 IF (electron(k)%Y.LE.c_Y_area_max) THEN                 
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 END IF

              CASE (SURROUNDED_BY_WALL)
                 IF ((c_X_area_min-electron(k)%X).LT.(electron(k)%Y-c_Y_area_max)) THEN
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 END IF

              CASE (EMPTY_CORNER_WALL_LEFT)
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)

              CASE (EMPTY_CORNER_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)

           END SELECT

        ELSE
! ERROR, we shouldn't be here
           PRINT '("ERROR-1 in ADVANCE_ELECTRONS: we should not be here")'
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        END IF
        CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
        CYCLE
     END IF

     IF (electron(k)%X.GT.c_X_area_max) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the electron by the period length
           electron(k)%X = electron(k)%X - L_period_X
           IF (electron(k)%Y.LT.c_Y_area_min) THEN
! particle is somewhere near the left bottom cell of the area
              IF (Rank_of_master_below.LT.0) THEN
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)                       
              END IF
              CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
           ELSE IF (electron(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
              IF (Rank_of_master_above.LT.0) THEN
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              END IF
              CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
           END IF
           CYCLE
        END IF

        IF ( (electron(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
           & (electron(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

           IF (Rank_of_master_right.GE.0) THEN
! right neighbor cluster exists
              CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           ELSE
! right neighbor cluster does not exist
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           END IF

        ELSE IF (electron(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the right bottom cell of the area

           SELECT CASE (c_right_bottom_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((electron(k)%X-c_X_area_max).LT.(c_Y_area_min-electron(k)%Y)) THEN
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 END IF

              CASE (FLAT_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                
              CASE (FLAT_WALL_RIGHT)
                 IF (electron(k)%Y.GE.c_Y_area_min) THEN
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)                       
                 END IF

              CASE (SURROUNDED_BY_WALL)
                 IF ((electron(k)%X-c_X_area_max).LT.(c_Y_area_min-electron(k)%Y)) THEN
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 END IF

              CASE (EMPTY_CORNER_WALL_RIGHT)
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)

              CASE (EMPTY_CORNER_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)

           END SELECT

        ELSE IF (electron(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top right cell of the area

           SELECT CASE (c_right_top_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((electron(k)%X-c_X_area_max).LT.(electron(k)%Y-c_Y_area_max)) THEN
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 END IF

              CASE (FLAT_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)

              CASE (FLAT_WALL_RIGHT)
                 IF (electron(k)%Y.LE.c_Y_area_max) THEN
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)                    
                 END IF

              CASE (SURROUNDED_BY_WALL)
                 IF ((electron(k)%X-c_X_area_max).LT.(electron(k)%Y-c_Y_area_max)) THEN
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 ELSE
                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
                 END IF

              CASE (EMPTY_CORNER_WALL_RIGHT)
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)

              CASE (EMPTY_CORNER_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)

           END SELECT

        ELSE
! ERROR, we shouldn't be here
           PRINT '("ERROR-2 in ADVANCE_ELECTRONS: we should not be here")'
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        END IF
        CALL REMOVE_ELECTRON(k)  !       this subroutine does  N_electrons = N_electrons - 1
        CYCLE
     END IF

! since we are here, c_X_area_min <= electron(k)%Y <= c_X_area_max

     IF (electron(k)%Y.GT.c_Y_area_max) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
           IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
              CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           ELSE
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           END IF
           CALL REMOVE_ELECTRON(k)  !       this subroutine does  N_electrons = N_electrons - 1
           CYCLE
        END IF

        IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
           CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
        ELSE
! neighbor cluster above does not exist
           IF ((electron(k)%X.GE.(c_X_area_min+1.0_8)).AND.(electron(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           ELSE IF (electron(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left top corner
              IF (c_left_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              END IF
           ELSE IF (electron(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right top corner
              IF (c_right_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX,  electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              END IF
           END IF
        END IF
        CALL REMOVE_ELECTRON(k)  !       this subroutine does  N_electrons = N_electrons - 1
        CYCLE
     END IF

     IF (electron(k)%Y.LT.c_Y_area_min) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
           IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster above exists
              CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           ELSE
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           END IF
           CALL REMOVE_ELECTRON(k)  !       this subroutine does  N_electrons = N_electrons - 1
           CYCLE
        END IF

        IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster below exists, remove particle and prepare to send it to the neighbor below
           CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
        ELSE
! neighbor cluster below does not exist
           IF ((electron(k)%X.GE.(c_X_area_min+1.0_8)).AND.(electron(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           ELSE IF (electron(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left bottom corner
              IF (c_left_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              END IF
           ELSE IF (electron(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right bottom corner
              IF (c_right_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              END IF
           END IF
        END IF
        CALL REMOVE_ELECTRON(k)  !       this subroutine does  N_electrons = N_electrons - 1
        CYCLE
     END IF
! 

  END DO

!print '("Proc ",i4," will add ",i8," send l/r/a/b ",4(2x,i8))', Rank_of_process, N_e_to_add, N_e_to_send_left, N_e_to_send_right, N_e_to_send_above, N_e_to_send_below

END SUBROUTINE ADVANCE_ELECTRONS

!----------------------------------------
!
SUBROUTINE REMOVE_ELECTRON(k)

  USE ParallelOperationValues
  USE ElectronParticles

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(INOUT) :: k

  IF ((k.LT.1).OR.(k.GT.N_electrons)) THEN
     PRINT '("Process ",i6," : ERROR in REMOVE_ELECTRON : index k invalid")', Rank_of_process
     PRINT '("Process ",i6," : k= ", i7," N_electrons= ",i7)', Rank_of_process, k, N_electrons
     PRINT '("Process ",i6," : PROGRAM TERMINATED")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  IF (k.LT.N_electrons) THEN

     electron(k)%X   = electron(N_electrons)%X
     electron(k)%Y   = electron(N_electrons)%Y
     electron(k)%VX  = electron(N_electrons)%VX
     electron(k)%VY  = electron(N_electrons)%VY
     electron(k)%VZ  = electron(N_electrons)%VZ
     electron(k)%tag = electron(N_electrons)%tag

  END IF

  N_electrons = N_electrons - 1
  k = k-1                  ! to ensure that the new k-particle is processed

!print '("Process ",i4," called REMOVE_ELECTRON(k), T_cntr= ",i7," k= ",i4)', Rank_of_process, T_cntr, k

END SUBROUTINE REMOVE_ELECTRON

!------------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_SEND_LEFT(x, y, vx, vy, vz, tag)

  USE ElectronParticles, ONLY : N_e_to_send_left, max_N_e_to_send_left, electron_to_send_left
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_X_left, L_period_x

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

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
  
  N_e_to_send_left = N_e_to_send_left + 1

  IF (N_e_to_send_left.GT.max_N_e_to_send_left) THEN
! increase the size of the list array
     current_N = max_N_e_to_send_left
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron_to_send_left(k)%X
        bufer(k)%Y   = electron_to_send_left(k)%Y
        bufer(k)%VX  = electron_to_send_left(k)%VX
        bufer(k)%VY  = electron_to_send_left(k)%VY
        bufer(k)%VZ  = electron_to_send_left(k)%VZ
        bufer(k)%tag = electron_to_send_left(k)%tag
     END DO
     IF (ALLOCATED(electron_to_send_left)) DEALLOCATE(electron_to_send_left, STAT=DEALLOC_ERR)
     max_N_e_to_send_left = max_N_e_to_send_left + MAX(50, max_N_e_to_send_left/10)
     ALLOCATE(electron_to_send_left(1:max_N_e_to_send_left), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        electron_to_send_left(k)%X   = bufer(k)%X
        electron_to_send_left(k)%Y   = bufer(k)%Y
        electron_to_send_left(k)%VX  = bufer(k)%VX
        electron_to_send_left(k)%VY  = bufer(k)%VY
        electron_to_send_left(k)%VZ  = bufer(k)%VZ
        electron_to_send_left(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  IF (periodic_boundary_X_left) THEN
     electron_to_send_left(N_e_to_send_left)%X = x + L_period_X
  ELSE
     electron_to_send_left(N_e_to_send_left)%X = x
  END IF
  electron_to_send_left(N_e_to_send_left)%Y   = y
  electron_to_send_left(N_e_to_send_left)%VX  = vx
  electron_to_send_left(N_e_to_send_left)%VY  = vy
  electron_to_send_left(N_e_to_send_left)%VZ  = vz
  electron_to_send_left(N_e_to_send_left)%tag = tag

!print '("Process ",i4," called ADD_ELECTRON_TO_SEND_LEFT, T_cntr= ",i7)', Rank_of_process, T_cntr

END SUBROUTINE ADD_ELECTRON_TO_SEND_LEFT

!------------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_SEND_RIGHT(x, y, vx, vy, vz, tag)

  USE ElectronParticles, ONLY : N_e_to_send_right, max_N_e_to_send_right, electron_to_send_right
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_X_right, L_period_x

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

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
  
  N_e_to_send_right = N_e_to_send_right + 1

  IF (N_e_to_send_right.GT.max_N_e_to_send_right) THEN
! increase the size of the list array
     current_N = max_N_e_to_send_right
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron_to_send_right(k)%X
        bufer(k)%Y   = electron_to_send_right(k)%Y
        bufer(k)%VX  = electron_to_send_right(k)%VX
        bufer(k)%VY  = electron_to_send_right(k)%VY
        bufer(k)%VZ  = electron_to_send_right(k)%VZ
        bufer(k)%tag = electron_to_send_right(k)%tag
     END DO
     IF (ALLOCATED(electron_to_send_right)) DEALLOCATE(electron_to_send_right, STAT=DEALLOC_ERR)
     max_N_e_to_send_right = max_N_e_to_send_right + MAX(50, max_N_e_to_send_right/10)
     ALLOCATE(electron_to_send_right(1:max_N_e_to_send_right), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        electron_to_send_right(k)%X   = bufer(k)%X
        electron_to_send_right(k)%Y   = bufer(k)%Y
        electron_to_send_right(k)%VX  = bufer(k)%VX
        electron_to_send_right(k)%VY  = bufer(k)%VY
        electron_to_send_right(k)%VZ  = bufer(k)%VZ
        electron_to_send_right(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  IF (periodic_boundary_X_right) THEN
     electron_to_send_right(N_e_to_send_right)%X = x - L_period_x
  ELSE
     electron_to_send_right(N_e_to_send_right)%X = x
  END IF
  electron_to_send_right(N_e_to_send_right)%Y   = y
  electron_to_send_right(N_e_to_send_right)%VX  = vx
  electron_to_send_right(N_e_to_send_right)%VY  = vy
  electron_to_send_right(N_e_to_send_right)%VZ  = vz
  electron_to_send_right(N_e_to_send_right)%tag = tag

!print '("Process ",i4," called ADD_ELECTRON_TO_SEND_RIGHT, T_cntr= ",i7)', Rank_of_process, T_cntr
!print '("Process ",i4," : ADD_ELECTRON_TO_SEND_RIGHT : x/y/vx/vy/vz/tag ",5(2x,e14.7),2x,i4)', Rank_of_process, x, y, vx, vy, vz,tag

END SUBROUTINE ADD_ELECTRON_TO_SEND_RIGHT

!------------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_SEND_ABOVE(x, y, vx, vy, vz, tag)

  USE ElectronParticles, ONLY : N_e_to_send_above, max_N_e_to_send_above, electron_to_send_above
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_Y_above, L_period_y

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

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
  
  N_e_to_send_above = N_e_to_send_above + 1

  IF (N_e_to_send_above.GT.max_N_e_to_send_above) THEN
! increase the size of the list array
     current_N = max_N_e_to_send_above
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron_to_send_above(k)%X
        bufer(k)%Y   = electron_to_send_above(k)%Y
        bufer(k)%VX  = electron_to_send_above(k)%VX
        bufer(k)%VY  = electron_to_send_above(k)%VY
        bufer(k)%VZ  = electron_to_send_above(k)%VZ
        bufer(k)%tag = electron_to_send_above(k)%tag
     END DO
     IF (ALLOCATED(electron_to_send_above)) DEALLOCATE(electron_to_send_above, STAT=DEALLOC_ERR)
     max_N_e_to_send_above = max_N_e_to_send_above + MAX(50, max_N_e_to_send_above/10)
     ALLOCATE(electron_to_send_above(1:max_N_e_to_send_above), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        electron_to_send_above(k)%X   = bufer(k)%X
        electron_to_send_above(k)%Y   = bufer(k)%Y
        electron_to_send_above(k)%VX  = bufer(k)%VX
        electron_to_send_above(k)%VY  = bufer(k)%VY
        electron_to_send_above(k)%VZ  = bufer(k)%VZ
        electron_to_send_above(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  electron_to_send_above(N_e_to_send_above)%X   = x
  IF (periodic_boundary_Y_above) THEN
     electron_to_send_above(N_e_to_send_above)%Y = y - L_period_y
  ELSE
     electron_to_send_above(N_e_to_send_above)%Y = y
  END IF
  electron_to_send_above(N_e_to_send_above)%VX  = vx
  electron_to_send_above(N_e_to_send_above)%VY  = vy
  electron_to_send_above(N_e_to_send_above)%VZ  = vz
  electron_to_send_above(N_e_to_send_above)%tag = tag

!print '("Process ",i4," called ADD_ELECTRON_TO_SEND_ABOVE, T_cntr= ",i7)', Rank_of_process, T_cntr

END SUBROUTINE ADD_ELECTRON_TO_SEND_ABOVE

!------------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_SEND_BELOW(x, y, vx, vy, vz, tag)

  USE ElectronParticles, ONLY : N_e_to_send_below, max_N_e_to_send_below, electron_to_send_below
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_Y_below, L_period_y

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

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
  
  N_e_to_send_below = N_e_to_send_below + 1

  IF (N_e_to_send_below.GT.max_N_e_to_send_below) THEN
! increase the size of the list array
     current_N = max_N_e_to_send_below
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron_to_send_below(k)%X
        bufer(k)%Y   = electron_to_send_below(k)%Y
        bufer(k)%VX  = electron_to_send_below(k)%VX
        bufer(k)%VY  = electron_to_send_below(k)%VY
        bufer(k)%VZ  = electron_to_send_below(k)%VZ
        bufer(k)%tag = electron_to_send_below(k)%tag
     END DO
     IF (ALLOCATED(electron_to_send_below)) DEALLOCATE(electron_to_send_below, STAT=DEALLOC_ERR)
     max_N_e_to_send_below = max_N_e_to_send_below + MAX(50, max_N_e_to_send_below/10)
     ALLOCATE(electron_to_send_below(1:max_N_e_to_send_below), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        electron_to_send_below(k)%X   = bufer(k)%X
        electron_to_send_below(k)%Y   = bufer(k)%Y
        electron_to_send_below(k)%VX  = bufer(k)%VX
        electron_to_send_below(k)%VY  = bufer(k)%VY
        electron_to_send_below(k)%VZ  = bufer(k)%VZ
        electron_to_send_below(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  electron_to_send_below(N_e_to_send_below)%X   = x
  IF (periodic_boundary_Y_below) THEN
     electron_to_send_below(N_e_to_send_below)%Y = y + L_period_y
  ELSE
     electron_to_send_below(N_e_to_send_below)%Y = y
  END IF
  electron_to_send_below(N_e_to_send_below)%VX  = vx
  electron_to_send_below(N_e_to_send_below)%VY  = vy
  electron_to_send_below(N_e_to_send_below)%VZ  = vz
  electron_to_send_below(N_e_to_send_below)%tag = tag

!print '("Process ",i4," called ADD_ELECTRON_TO_SEND_BELOW, T_cntr= ",i7)', Rank_of_process, T_cntr
!print '("Process ",i4," : ADD_ELECTRON_TO_SEND_BELOW : x/y/vx/vy/vz/tag ",5(2x,e14.7),2x,i4)', Rank_of_process, x, y, vx, vy, vz,tag

END SUBROUTINE ADD_ELECTRON_TO_SEND_BELOW

!----------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)

  USE ElectronParticles, ONLY : N_e_to_add, max_N_e_to_add, electron_to_add
!  USE ParallelOperationValues

  IMPLICIT NONE

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
  
  N_e_to_add = N_e_to_add + 1

  IF (N_e_to_add.GT.max_N_e_to_add) THEN
! increase the size of the list array
     current_N = max_N_e_to_add
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron_to_add(k)%X
        bufer(k)%Y   = electron_to_add(k)%Y
        bufer(k)%VX  = electron_to_add(k)%VX
        bufer(k)%VY  = electron_to_add(k)%VY
        bufer(k)%VZ  = electron_to_add(k)%VZ
        bufer(k)%tag = electron_to_add(k)%tag
     END DO
     IF (ALLOCATED(electron_to_add)) DEALLOCATE(electron_to_add, STAT=DEALLOC_ERR)
     max_N_e_to_add = max_N_e_to_add + MAX(50, max_N_e_to_add/10)
     ALLOCATE(electron_to_add(1:max_N_e_to_add), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        electron_to_add(k)%X   = bufer(k)%X
        electron_to_add(k)%Y   = bufer(k)%Y
        electron_to_add(k)%VX  = bufer(k)%VX
        electron_to_add(k)%VY  = bufer(k)%VY
        electron_to_add(k)%VZ  = bufer(k)%VZ
        electron_to_add(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  electron_to_add(N_e_to_add)%X   = x
  electron_to_add(N_e_to_add)%Y   = y
  electron_to_add(N_e_to_add)%VX  = vx
  electron_to_add(N_e_to_add)%VY  = vy
  electron_to_add(N_e_to_add)%VZ  = vz
  electron_to_add(N_e_to_add)%tag = tag

!print '("Process ",i4," called ADD_ELECTRON_TO_ADD_LIST, T_cntr= ",i7)', Rank_of_process, T_cntr
!print '("Process ",i4," : ADD_ELECTRON_TO_ADD_LIST : x/y/vx/vy/vz/tag ",5(2x,e14.7),2x,i4)', Rank_of_process, x, y, vx, vy, vz,tag

END SUBROUTINE ADD_ELECTRON_TO_ADD_LIST

!----------------------------------------
!
SUBROUTINE REMOVE_ELECTRON_FROM_ADD_LIST(k)

  USE ParallelOperationValues
  USE ElectronParticles

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(INOUT) :: k

  IF ((k.LT.1).OR.(k.GT.N_e_to_add)) THEN
     PRINT '("Process ",i6," : ERROR in REMOVE_ELECTRON_FROM_ADD_LIST : index k invalid")', Rank_of_process
     PRINT '("Process ",i6," : k= ", i7," N_e_to_add= ",i7)', Rank_of_process, k, N_e_to_add
     PRINT '("Process ",i6," : PROGRAM TERMINATED")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  IF (k.LT.N_e_to_add) THEN

     electron_to_add(k)%X   = electron_to_add(N_e_to_add)%X
     electron_to_add(k)%Y   = electron_to_add(N_e_to_add)%Y
     electron_to_add(k)%VX  = electron_to_add(N_e_to_add)%VX
     electron_to_add(k)%VY  = electron_to_add(N_e_to_add)%VY
     electron_to_add(k)%VZ  = electron_to_add(N_e_to_add)%VZ
     electron_to_add(k)%tag = electron_to_add(N_e_to_add)%tag

  END IF

  N_e_to_add = N_e_to_add - 1
  k = k-1                  ! to ensure that the new k-particle is processed

!print '("Process ",i4," called REMOVE_ELECTRON_FROM_ADD_LIST(k), T_cntr= ",i7," k= ",i4)', Rank_of_process, T_cntr, k

END SUBROUTINE REMOVE_ELECTRON_FROM_ADD_LIST

!----------------------------------------
! This subroutine is called after ADVANCE_ELECTRONS, before exchange of electrons takes place.
! It removes particles which do not belong to the cluster domain.
! Such particles may appear after emission from surface of inner objects.
! It is expected that the emission is never directed into boundary objects aligned along the main simulation domain's boundary.
! The algorithm below, which places alien particles into proper SEND* lists is the same as in ADVANCE_ELECTRONS.
! However, collision with a boundary object at this stage is not expected and is considered an error.
!
SUBROUTINE FIND_ALIENS_IN_ELECTRON_ADD_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE ElectronParticles, ONLY : N_e_to_add, electron_to_add

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER k, n 

  IF (N_of_inner_objects.EQ.0) RETURN

! find, process, and exclude electrons which collided with inner objects
  k=0
  DO WHILE (k.LT.N_e_to_add)
     k = k+1

     IF (symmetry_plane_X_left) THEN
        IF (electron_to_add(k)%X.LT.c_X_area_min) THEN
! since we are here, a particle to be added is beyond the symmetry plane at X=0
! this is a rare but possible situation, for example:
! a particle crosses the symmetry plane and collides with an inner object at the same time
! such a particle is reflected before trying for collision, but the crossing point for the reflected particle may still be at X<0
! [because for collisions with inner objects we do "ray tracing" to find the point of collision]
! if the collision is followed by emission of a secondary electron, then the emitted electron will have X<0
! so we simply move this particle rightward, symmetrically relative to plane X=0
           electron_to_add(k)%X = MAX(c_X_area_min, c_X_area_min + c_X_area_min - electron_to_add(k)%X)
! unlike in ADVANCE_ELECTRONS, here we do not change velocity because it is random
!           electron(k)%VX = -electron(k)%VX
!!###          electron(k)%VZ = -electron(k)%VZ   !###??? do we not have to do this when BY is on ??? 
        END IF
     END IF


! most probable situation when the particle is inside the area
     IF ( (electron_to_add(k)%X.GE.c_X_area_min) .AND. &
        & (electron_to_add(k)%X.LE.c_X_area_max) .AND. &
        & (electron_to_add(k)%Y.GE.c_Y_area_min) .AND. &
        & (electron_to_add(k)%Y.LE.c_Y_area_max) ) CYCLE

! since we are here, particle k is outside the domain of this cluster

     IF (electron_to_add(k)%X.LT.c_X_area_min) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the electron by the period length
           electron_to_add(k)%X = electron_to_add(k)%X + L_period_X
           IF (electron_to_add(k)%Y.LT.c_Y_area_min) THEN
! particle is below the bottom side of the area
              IF (Rank_of_master_below.LT.0) THEN
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-1 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)                       
              END IF
              CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           ELSE IF (electron_to_add(k)%Y.GT.c_Y_area_max) THEN
! particle is above the top side of the area
              IF (Rank_of_master_above.LT.0) THEN
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-2 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
              END IF
              CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           END IF
           CYCLE
        END IF

        IF ( (electron_to_add(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
           & (electron_to_add(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

           IF (Rank_of_master_left.GE.0) THEN
! left neighbor cluster exists
              CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
           ELSE
! left neighbor cluster does not exist
!#              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX,  electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)   ! left
! error
                 print '("error-3 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF

        ELSE IF (electron_to_add(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the left bottom cell of the area

           SELECT CASE (c_left_bottom_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((c_X_area_min-electron_to_add(k)%X).LT.(c_Y_area_min-electron_to_add(k)%Y)) THEN
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
                 END IF

              CASE (FLAT_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

              CASE (FLAT_WALL_LEFT)
                 IF (electron_to_add(k)%Y.GE.c_Y_area_min) THEN                 
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                    print '("error-4 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)                       
                 END IF

              CASE (SURROUNDED_BY_WALL)
!                 IF ((c_X_area_min-electron_to_add(k)%X).LT.(c_Y_area_min-electron_to_add(k)%Y)) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 ELSE
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 END IF
! error
                 print '("error-5 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

              CASE (EMPTY_CORNER_WALL_LEFT)
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

              CASE (EMPTY_CORNER_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

           END SELECT

        ELSE IF (electron_to_add(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top left cell of the area

           SELECT CASE (c_left_top_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((c_X_area_min-electron_to_add(k)%X).LT.(electron_to_add(k)%Y-c_Y_area_max)) THEN
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
                 END IF

              CASE (FLAT_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

              CASE (FLAT_WALL_LEFT)
                 IF (electron_to_add(k)%Y.LE.c_Y_area_max) THEN                 
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                    print '("error-6 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
                 END IF

              CASE (SURROUNDED_BY_WALL)
!                 IF ((c_X_area_min-electron_to_add(k)%X).LT.(electron_to_add(k)%Y-c_Y_area_max)) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 ELSE
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 END IF
! error
                 print '("error-7 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

              CASE (EMPTY_CORNER_WALL_LEFT)
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

              CASE (EMPTY_CORNER_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

           END SELECT
        END IF
        CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
        CYCLE
     END IF         !### IF (electron_to_add(k)%X.LT.c_X_area_min) THEN

     IF (electron_to_add(k)%X.GT.c_X_area_max) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the electron by the period length
           electron_to_add(k)%X = electron_to_add(k)%X - L_period_X
           IF (electron_to_add(k)%Y.LT.c_Y_area_min) THEN
! particle is below the bottom side of the area
              IF (Rank_of_master_below.LT.0) THEN
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-8 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)                       
              END IF
              CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           ELSE IF (electron_to_add(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
              IF (Rank_of_master_above.LT.0) THEN
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-9 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
              END IF
              CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           END IF
           CYCLE
        END IF

        IF ( (electron_to_add(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
           & (electron_to_add(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

           IF (Rank_of_master_right.GE.0) THEN
! right neighbor cluster exists
              CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
           ELSE
! right neighbor cluster does not exist
!              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
              print '("error-10 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF

        ELSE IF (electron_to_add(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the right bottom cell of the area

           SELECT CASE (c_right_bottom_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((electron_to_add(k)%X-c_X_area_max).LT.(c_Y_area_min-electron_to_add(k)%Y)) THEN
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
                 END IF

              CASE (FLAT_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
                
              CASE (FLAT_WALL_RIGHT)
                 IF (electron_to_add(k)%Y.GE.c_Y_area_min) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                    print '("error-11 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)                       
                 END IF

              CASE (SURROUNDED_BY_WALL)
!                 IF ((electron_to_add(k)%X-c_X_area_max).LT.(c_Y_area_min-electron_to_add(k)%Y)) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 ELSE
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 END IF
! error
                    print '("error-12 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

              CASE (EMPTY_CORNER_WALL_RIGHT)
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

              CASE (EMPTY_CORNER_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

           END SELECT

        ELSE IF (electron_to_add(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top right cell of the area

           SELECT CASE (c_right_top_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((electron_to_add(k)%X-c_X_area_max).LT.(electron_to_add(k)%Y-c_Y_area_max)) THEN
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
                 END IF

              CASE (FLAT_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

              CASE (FLAT_WALL_RIGHT)
                 IF (electron_to_add(k)%Y.LE.c_Y_area_max) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                    print '("error-13 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)                    
                 END IF

              CASE (SURROUNDED_BY_WALL)
!                 IF ((electron_to_add(k)%X-c_X_area_max).LT.(electron_to_add(k)%Y-c_Y_area_max)) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 ELSE
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 END IF
! error
                    print '("error-14 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

              CASE (EMPTY_CORNER_WALL_RIGHT)
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

              CASE (EMPTY_CORNER_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)

           END SELECT

        END IF
        CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
        CYCLE
     END IF         !### IF (electron_to_add(k)%X.GT.c_X_area_max) THEN

! since we are here, c_X_area_min <= electron_to_add(k)%Y <= c_X_area_max

     IF (electron_to_add(k)%Y.GT.c_Y_area_max) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
           IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
              CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
           ELSE
!              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
              print '("error-15 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF
           CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           CYCLE
        END IF

        IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
           CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
        ELSE
! neighbor cluster above does not exist
           IF ((electron_to_add(k)%X.GE.(c_X_area_min+1.0_8)).AND.(electron_to_add(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
!              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
              print '("error-16 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           ELSE IF (electron_to_add(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left top corner
              IF (c_left_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
              ELSE
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-17 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
           ELSE IF (electron_to_add(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right top corner
              IF (c_right_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX,  electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
              ELSE
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-18 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
           END IF
        END IF
        CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
        CYCLE
     END IF     !### IF (electron_to_add(k)%Y.GT.c_Y_area_max) THEN

     IF (electron_to_add(k)%Y.LT.c_Y_area_min) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
           IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster above exists
              CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
           ELSE
!              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
              print '("error-19 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF
           CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           CYCLE
        END IF

        IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster below exists, remove particle and prepare to send it to the neighbor below
           CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
        ELSE
! neighbor cluster below does not exist
           IF ((electron_to_add(k)%X.GE.(c_X_area_min+1.0_8)).AND.(electron_to_add(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
!              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
              print '("error-20 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           ELSE IF (electron_to_add(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left bottom corner
              IF (c_left_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
              ELSE
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-21 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
           ELSE IF (electron_to_add(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right bottom corner
              IF (c_right_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
              ELSE
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-22 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
           END IF
        END IF
        CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
        CYCLE
     END IF     !### IF (electron_to_add(k)%Y.LT.c_Y_area_min) THEN

  END DO   !### DO WHILE (k.LT.N_e_to_add)

END SUBROUTINE FIND_ALIENS_IN_ELECTRON_ADD_LIST

!----------------------------------------
!
SUBROUTINE FIND_INNER_OBJECT_COLL_IN_ELECTRON_ADD_LIST

!  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles, ONLY : N_e_to_add, electron_to_add

  IMPLICIT NONE

  INTEGER k, n 

  IF (N_of_inner_objects.EQ.0) RETURN

! find, process, and exclude electrons which collided with inner objects
  k=0
  DO WHILE (k.LT.N_e_to_add)
     k = k+1
     DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
        IF (electron_to_add(k)%X.LE.whole_object(n)%Xmin) CYCLE
        IF (electron_to_add(k)%X.GE.whole_object(n)%Xmax) CYCLE
        IF (electron_to_add(k)%Y.LE.whole_object(n)%Ymin) CYCLE
        IF (electron_to_add(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! collision detected
        CALL TRY_ELECTRON_COLL_WITH_INNER_OBJECT( electron_to_add(k)%X, &
                                                & electron_to_add(k)%Y, &
                                                & electron_to_add(k)%VX, &
                                                & electron_to_add(k)%VY, &
                                                & electron_to_add(k)%VZ, &
                                                & electron_to_add(k)%tag) !, &
!                                                & whole_object(n) )
        CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
        EXIT
     END DO
  END DO

END SUBROUTINE FIND_INNER_OBJECT_COLL_IN_ELECTRON_ADD_LIST

!----------------------------------------
!
SUBROUTINE PROCESS_ADDED_ELECTRONS

  USE ParallelOperationValues
  USE CurrentProblemValues

  USE ElectronParticles, ONLY : N_e_to_add, N_electrons, max_N_electrons, electron, electron_to_add

  USE rng_wrapper  !###???

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

  INTEGER ALLOC_ERR
  INTEGER k, current_N

  INTEGER random_j
  INTEGER temptag
  REAL(8) tempX, tempY, tempVX, tempVY, tempVZ
  
  IF (N_e_to_add.GT.(max_N_electrons-N_electrons)) THEN
! increase the size of the main electron array
     current_N = max_N_electrons
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron(k)%X
        bufer(k)%Y   = electron(k)%Y
        bufer(k)%VX  = electron(k)%VX
        bufer(k)%VY  = electron(k)%VY
        bufer(k)%VZ  = electron(k)%VZ
        bufer(k)%tag = electron(k)%tag
     END DO
     DEALLOCATE(electron, STAT=ALLOC_ERR)
     max_N_electrons = max_N_electrons + MAX(N_e_to_add-(max_N_electrons-N_electrons), max_N_electrons/10)
     ALLOCATE(electron(1:max_N_electrons), STAT=ALLOC_ERR)
     DO k = 1, current_N
        electron(k)%X   = bufer(k)%X
        electron(k)%Y   = bufer(k)%Y
        electron(k)%VX  = bufer(k)%VX
        electron(k)%VY  = bufer(k)%VY
        electron(k)%VZ  = bufer(k)%VZ
        electron(k)%tag = bufer(k)%tag
     END DO
     DEALLOCATE(bufer, STAT=ALLOC_ERR)
  END IF
  
  DO k = 1, N_e_to_add
     electron(k+N_electrons)%X   = electron_to_add(k)%X
     electron(k+N_electrons)%Y   = electron_to_add(k)%Y
     electron(k+N_electrons)%VX  = electron_to_add(k)%VX
     electron(k+N_electrons)%VY  = electron_to_add(k)%VY
     electron(k+N_electrons)%VZ  = electron_to_add(k)%VZ
     electron(k+N_electrons)%tag = electron_to_add(k)%tag
  END DO

! update electron counter
  N_electrons = N_electrons + N_e_to_add

!???????????? shuffle the added and the available particles ?????????????
  DO k = N_electrons - N_e_to_add + 1, N_electrons
!     random_j = MAX(1, MIN( N_electrons - N_e_to_add, INT(well_random_number() * (N_electrons-N_e_to_add))) )
     random_j = MAX(1, MIN( N_electrons, INT(well_random_number() * N_electrons)))

     IF (random_j.EQ.k) CYCLE

     tempX   = electron(random_j)%X
     tempY   = electron(random_j)%Y
     tempVX  = electron(random_j)%VX
     tempVY  = electron(random_j)%VY
     tempVZ  = electron(random_j)%VZ
     temptag = electron(random_j)%tag

     electron(random_j)%X   = electron(k)%X
     electron(random_j)%Y   = electron(k)%Y
     electron(random_j)%VX  = electron(k)%VX
     electron(random_j)%VY  = electron(k)%VY
     electron(random_j)%VZ  = electron(k)%VZ
     electron(random_j)%tag = electron(k)%tag

     electron(k)%X   = tempX
     electron(k)%Y   = tempY
     electron(k)%VX  = tempVX
     electron(k)%VY  = tempVY
     electron(k)%VZ  = tempVZ
     electron(k)%tag = temptag
  END DO

!print '("Process ",i4," called PROCESS_ADDED_ELECTRONS, T_cntr= ",i7," N_e_to_add= ",i8," N_electrons= ",i8," max_N_electrons= ",i8)', &
!     & Rank_of_process, T_cntr, N_e_to_add, N_electrons, max_N_electrons

! clear counter of particles to be added to this process
  N_e_to_add = 0  

END SUBROUTINE PROCESS_ADDED_ELECTRONS

!------------------------------
!
SUBROUTINE GATHER_ELECTRON_CHARGE_DENSITY

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries, ONLY : indx_x_min, indx_x_max, indx_y_min, indx_y_max
!  USE Diagnostics

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

!  REAL(8), ALLOCATABLE :: c_rho(:,:)

  INTEGER n1  ! number of nodes in the y-direction
  INTEGER n2  !
  INTEGER n3  ! number of nodes in the x-direction

  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR
  INTEGER bufsize

!  INTEGER npc   ! number of probe in a cluster list of probes
!  INTEGER npa   ! number of probe in the global list of probes

  INTEGER i, j, k
  INTEGER pos_i_j, pos_ip1_j, pos_i_jp1, pos_ip1_jp1
  REAL(8) ax_ip1, ax_i, ay_jp1, ay_j
  REAL(8) vij, vip1j, vijp1

  INTEGER pos

  INTEGER nio, position_flag
! function
  REAL(8) Get_Surface_Charge_Inner_Object

  IF ( (cluster_rank_key.NE.0) &                     ! this array is needed only temporarily for non-master processes
     & .OR. &
     & (periodicity_flag.EQ.PERIODICITY_NONE) &      ! if PETSc is used, master processes don't need this array outside this sub
     & .OR. &
     & (periodicity_flag.EQ.PERIODICITY_X_PETSC) &
     & .OR. &
     & (periodicity_flag.EQ.PERIODICITY_X_Y) ) &
     & THEN
     ALLOCATE(c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
  END IF

  n1 = c_indx_y_max - c_indx_y_min + 1
! was  n2 = -c_indx_x_min + 1 - c_indx_y_min * n1
  n3 = c_indx_x_max - c_indx_x_min + 1
  n2 = -c_indx_x_min + 1 - c_indx_y_min * n3

  bufsize = n1 * n3
  ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

  rbufer = 0.0_8

  DO k = 1, N_electrons
     
     i = INT(electron(k)%X)
     j = INT(electron(k)%Y)
     IF (electron(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (electron(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

if ((i.lt.c_indx_x_min).or.(i.gt.(c_indx_x_max-1)).or.(j.lt.c_indx_y_min).or.(j.gt.(c_indx_y_max-1))) then
   print '("Process ",i4," : Error-1 in GATHER_ELECTRON_CHARGE_DENSITY : index out of bounds")', Rank_of_process
   print '("Process ",i4," : k/N_electrons : ",i8,2x,i8)', Rank_of_process, k, N_electrons
   print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag
   print '("Process ",i4," : minx/maxx/miny/maxy : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
   CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
end if

!     pos = i - c_indx_x_min + 1 + (j - c_indx_y_min) * (c_indx_x_max - c_indx_x_min + 1)

     pos_i_j     = i + j * n3 + n2
     pos_ip1_j   = pos_i_j + 1
     pos_i_jp1   = pos_i_j + n3
     pos_ip1_jp1 = pos_i_jp1 + 1

     ax_ip1 = electron(k)%X - DBLE(i)
     ax_i   = 1.0_8 - ax_ip1

     ay_jp1 = electron(k)%Y - DBLE(j)
     ay_j = 1.0_8 - ay_jp1

     vij   = ax_i   * ay_j
     vip1j = ax_ip1 * ay_j
     vijp1 = ax_i   * ay_jp1

if ((pos_i_j.gt.bufsize)) then
   print '(2x,8(2x,i10))', Rank_of_process, bufsize, pos_i_j, i, j, k, n1, n2
end if

     rbufer(pos_i_j)     = rbufer(pos_i_j)     + vij                         !ax_i   * ay_j
     rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   + vip1j                       !ax_ip1 * ay_j
     rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   + vijp1                       !ax_i   * ay_jp1
     rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) + 1.0_8 - vij - vip1j - vijp1 !ax_ip1 * ay_jp1

  END DO

! collect densities from all processes in a cluster
  CALL MPI_REDUCE(rbufer, c_rho, bufsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_CLUSTER, ierr)
  
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! now cluster masters exchange information about densities in overlapping nodes
  IF (cluster_rank_key.EQ.0) THEN

!! ################ diagnostics, electron density ##################
!! since GATHER_ELECTRON_CHARGE_DENSITY is called at every time step,
!! one can collect diagnostics data only at required time steps
!     IF (T_cntr.EQ.Save_probes_data_T_cntr) THEN
!        DO npc = 1, N_of_probes_cluster
!           npa = List_of_probes_cluster(npc)
!           i = Probe_position(1,npa)
!           j = Probe_position(2,npa)
!           probe_Ne_cluster(npc) = c_rho(i,j)   ! add scale factor ###################
!        END DO
!     END IF

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! special case of self-connected X-periodic cluster
        DO j = c_indx_y_min, c_indx_y_max
           c_rho(c_indx_x_min+1, j) = c_rho(c_indx_x_min+1, j) + c_rho(c_indx_x_max, j) 
           c_rho(c_indx_x_max-1, j) = c_rho(c_indx_x_max-1, j) + c_rho(c_indx_x_min, j)
           c_rho(c_indx_x_min, j) = c_rho(c_indx_x_max-1, j)
           c_rho(c_indx_x_max, j) = c_rho(c_indx_x_min+1, j)
        END DO
     END IF

     IF (WHITE_CLUSTER) THEN  
! "white processes"

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

        IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right densities in the right edge
           rbufer(1:n1) = c_rho(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left densities in the left edge
           rbufer(1:n1) = c_rho(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left densities in the vertical line next to the left edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_min+1, j) = c_rho(c_indx_x_min+1, j) + rbufer(j-c_indx_y_min+1)
           END DO
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right densities in the vertical line next to the right edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_max-1, j) = c_rho(c_indx_x_max-1, j) + rbufer(j-c_indx_y_min+1)
           END DO           
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

        IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up densities in the top edge
           rbufer(1:n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down densities in the bottom edge
           rbufer(1:n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below densities in the vertical line above the bottom line
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_min+1) = c_rho(i, c_indx_y_min+1) + rbufer(i-c_indx_x_min+1)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above densities in the vertical line under the top line
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_max-1) = c_rho(i, c_indx_y_max-1) + rbufer(i-c_indx_x_min+1)
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
              c_rho(c_indx_x_min+1, j) = c_rho(c_indx_x_min+1, j) + rbufer(j-c_indx_y_min+1)
           END DO
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right densities in the vertical line next to the right edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_max-1, j) = c_rho(c_indx_x_max-1, j) + rbufer(j-c_indx_y_min+1)
           END DO           
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right densities in the right edge
           rbufer(1:n1) = c_rho(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left densities in the left edge
           rbufer(1:n1) = c_rho(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

        IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below densities in the vertical line above the bottom line
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_min+1) = c_rho(i, c_indx_y_min+1) + rbufer(i-c_indx_x_min+1)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above densities in the vertical line under the top line
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_max-1) = c_rho(i, c_indx_y_max-1) + rbufer(i-c_indx_x_min+1)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up densities in the top edge
           rbufer(1:n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down densities in the bottom edge
           rbufer(1:n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

     END IF

! adjust densities at the boundaries with material walls

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! special case of self-connected X-periodic cluster
        IF (Rank_of_master_above.LT.0) THEN
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_max) = 2.0_8 * c_rho(i, c_indx_y_max)
           END DO
        END IF
  
        IF (Rank_of_master_below.LT.0) THEN
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_min) = 2.0_8 * c_rho(i, c_indx_y_min)
           END DO
        END IF

     ELSE
 
        IF (Rank_of_master_left.LT.0) THEN
           DO j = c_indx_y_min+1, c_indx_y_max-1
              c_rho(c_indx_x_min, j) = 2.0_8 * c_rho(c_indx_x_min, j)
           END DO
        END IF

        IF (Rank_of_master_right.LT.0) THEN
           DO j = c_indx_y_min+1, c_indx_y_max-1
              c_rho(c_indx_x_max, j) = 2.0_8 * c_rho(c_indx_x_max, j)
           END DO
        END IF
  
        IF (Rank_of_master_above.LT.0) THEN
           DO i = c_indx_x_min+1, c_indx_x_max-1
              c_rho(i, c_indx_y_max) = 2.0_8 * c_rho(i, c_indx_y_max)
           END DO
        END IF
  
        IF (Rank_of_master_below.LT.0) THEN
           DO i = c_indx_x_min+1, c_indx_x_max-1
              c_rho(i, c_indx_y_min) = 2.0_8 * c_rho(i, c_indx_y_min)
           END DO
        END IF

        SELECT CASE (c_left_bottom_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho(c_indx_x_min, c_indx_y_min) = 4.0_8 * c_rho(c_indx_x_min, c_indx_y_min) 
           CASE (EMPTY_CORNER_WALL_LEFT)
              c_rho(c_indx_x_min, c_indx_y_min+1) = 0.66666666666666_8 * c_rho(c_indx_x_min, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_BELOW)
              c_rho(c_indx_x_min+1, c_indx_y_min) = 0.66666666666666_8 * c_rho(c_indx_x_min+1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
        END SELECT
     
        SELECT CASE (c_left_top_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho(c_indx_x_min, c_indx_y_max) = 4.0_8 * c_rho(c_indx_x_min, c_indx_y_max) 
           CASE (EMPTY_CORNER_WALL_LEFT)
              c_rho(c_indx_x_min, c_indx_y_max-1) = 0.66666666666666_8 * c_rho(c_indx_x_min, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_ABOVE)
              c_rho(c_indx_x_min+1, c_indx_y_max) = 0.66666666666666_8 * c_rho(c_indx_x_min+1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
        END SELECT

        SELECT CASE (c_right_bottom_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho(c_indx_x_max, c_indx_y_min) = 4.0_8 * c_rho(c_indx_x_max, c_indx_y_min) 
           CASE (EMPTY_CORNER_WALL_RIGHT)
              c_rho(c_indx_x_max, c_indx_y_min+1) = 0.66666666666666_8 * c_rho(c_indx_x_max, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_BELOW)
              c_rho(c_indx_x_max-1, c_indx_y_min) = 0.66666666666666_8 * c_rho(c_indx_x_max-1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
        END SELECT
     
        SELECT CASE (c_right_top_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho(c_indx_x_max, c_indx_y_max) = 4.0_8 * c_rho(c_indx_x_max, c_indx_y_max) 
           CASE (EMPTY_CORNER_WALL_RIGHT)
              c_rho(c_indx_x_max, c_indx_y_max-1) = 0.66666666666666_8 * c_rho(c_indx_x_max, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_ABOVE)
              c_rho(c_indx_x_max-1, c_indx_y_max) = 0.66666666666666_8 * c_rho(c_indx_x_max-1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
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
                 rbufer(pos) = c_rho(i,j)

                 DO nio = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
                    CALL CHECK_IF_INNER_OBJECT_CONTAINS_POINT(whole_object(nio), i, j, position_flag)
                    rbufer(pos) = rbufer(pos) - Get_Surface_Charge_Inner_Object(i, j, position_flag, whole_object(nio))
                 END DO

!                 CALL FIND_INNER_OBJECT_CONTAINING_POINT(i,j,nio,position_flag)
!                 SELECT CASE (position_flag)
!! for points on surface of dielectric objects, combine surface charge density and volume charge density
!                    CASE (1,3,5,7)
!! point at the corner of a dielectric object
!                       rbufer(pos) = (1.0_8) * c_rho(i,j) - Get_Surface_Charge_Inner_Object(i,j,position_flag,whole_object(nio))
!                    CASE (2,4,6,8)
!! point on the surface of a dielectric object
!! the factor used here allows to use common factor -1/4*N_of_particles_cell in SOLVE_POTENTIAL_WITH_PETSC
!! it also allows to use uncorrected values of volume charge density in nodes on the flat material surface
!                       rbufer(pos) = 2.0_8 * (c_rho(i,j) - Get_Surface_Charge_Inner_Object(i,j,position_flag,whole_object(nio))) / (1.0_8 + whole_object(nio)%eps_diel)
!                 END SELECT

              END DO
           END DO
           CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
        END DO

! cluster master is a field calaculator too, prepare its own charge density
        DO j = indx_y_min, indx_y_max
           DO i = indx_x_min, indx_x_max
              rho_e(i, j) = c_rho(i, j)

              DO nio = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
                 CALL CHECK_IF_INNER_OBJECT_CONTAINS_POINT(whole_object(nio), i, j, position_flag)
                 rho_e(i, j) = rho_e(i, j) - Get_Surface_Charge_Inner_Object(i, j, position_flag, whole_object(nio))
              END DO

!              CALL FIND_INNER_OBJECT_CONTAINING_POINT(i,j,nio,position_flag)
!              SELECT CASE (position_flag)
!! for points on surface of dielectric objects, combine surface charge density and volume charge density
!                 CASE (1,3,5,7)
!! point at the corner of a dielectric object
!                    rho_e(i, j) = (1.0_8) * c_rho(i,j) - Get_Surface_Charge_Inner_Object(i,j,position_flag,whole_object(nio))
!                 CASE (2,4,6,8)
!! point on the surface of a dielectric object
!! the factor used here allows to use common factor -1/4*N_of_particles_cell in SOLVE_POTENTIAL_WITH_PETSC
!! it also allows to use uncorrected values of volume charge density in nodes on the flat material surface
!                    rho_e(i, j) = 2.0_8 * (c_rho(i,j) - Get_Surface_Charge_Inner_Object(i,j,position_flag,whole_object(nio))) / (1.0_8 + whole_object(nio)%eps_diel)
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
              rho_e(i,j) = rbufer(pos)
           END DO
        END DO

     END IF
     
  END IF
  
  IF (ALLOCATED(rbufer))  DEALLOCATE(rbufer,  STAT=ALLOC_ERR)

  IF ((cluster_rank_key.NE.0).OR.(periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN
     IF (ALLOCATED(c_rho)) DEALLOCATE(c_rho, STAT=ALLOC_ERR)
  END IF

END SUBROUTINE GATHER_ELECTRON_CHARGE_DENSITY
