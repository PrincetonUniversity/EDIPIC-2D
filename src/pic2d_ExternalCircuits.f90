!------------------------------------------------
!
SUBROUTINE CALCULATE_OBJECT_POTENTIALS_2D

  USE ParallelOperationValues
  USE ExternalCircuit
  USE CurrentProblemValues
  USE BlockAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER nn, n, noi, i, j

  IF (N_of_object_potentials_to_solve.EQ.0) RETURN

! set charge density equal to zero everywhere (this includes the surface charge density which normally is incorporated into rho_e)
  rho_i = 0.0_8
  rho_e = 0.0_8

  DO nn = 1, N_of_object_potentials_to_solve

! set zero potential of all metal electrodes
     DO n = 1, N_of_boundary_and_inner_objects
        IF (whole_object(n)%object_type.NE.METAL_WALL) CYCLE
        whole_object(n)%phi = 0.0_8
     END DO

     noi = object_charge_calculation(nn)%noi   ! number of object of interest

! set potential of the electrode of interest to 1   !### it is dimensionless 1 for now
     whole_object(noi)%phi = 1.0_8

! set potentials across vacuum gaps (same as in UPDATE_WALL_POTENTIALS ### make it a separate routine???)
     DO n = 1, N_of_boundary_objects
        IF (whole_object(n)%object_type.NE.VACUUM_GAP) CYCLE
     
        IF (whole_object(n)%segment(1)%jstart.NE.whole_object(n)%segment(1)%jend) THEN
! the gap is vertical (stretched along y)
           DO j = whole_object(n)%segment(1)%jstart+1, whole_object(n)%segment(1)%jend-1
              whole_object(n)%phi_profile(j) = whole_object(whole_object(n)%n_connected_to_start)%phi +  &
                                             & (whole_object(whole_object(n)%n_connected_to_end)%phi - whole_object(whole_object(n)%n_connected_to_start)%phi) * &
                                             & DBLE(j - whole_object(n)%segment(1)%jstart) / DBLE(whole_object(n)%segment(1)%jend - whole_object(n)%segment(1)%jstart)
           END DO
           whole_object(n)%phi_profile(whole_object(n)%segment(1)%jstart) = whole_object(whole_object(n)%n_connected_to_start)%phi
           whole_object(n)%phi_profile(whole_object(n)%segment(1)%jend)   = whole_object(whole_object(n)%n_connected_to_end)%phi
        ELSE IF (whole_object(n)%segment(1)%istart.NE.whole_object(n)%segment(1)%iend) THEN
! the gap is horizontal (stretched along x)
           DO i = whole_object(n)%segment(1)%istart+1, whole_object(n)%segment(1)%iend-1
              whole_object(n)%phi_profile(i) = whole_object(whole_object(n)%n_connected_to_start)%phi +  &
                                             & (whole_object(whole_object(n)%n_connected_to_end)%phi - whole_object(whole_object(n)%n_connected_to_start)%phi) * &
                                             & DBLE(i - whole_object(n)%segment(1)%istart) / DBLE(whole_object(n)%segment(1)%iend - whole_object(n)%segment(1)%istart)
           END DO
           whole_object(n)%phi_profile(whole_object(n)%segment(1)%istart) = whole_object(whole_object(n)%n_connected_to_start)%phi
           whole_object(n)%phi_profile(whole_object(n)%segment(1)%iend)   = whole_object(whole_object(n)%n_connected_to_end)%phi
        END IF
     END DO

     CALL SOLVE_POTENTIAL_WITH_PETSC

     DO j = indx_y_min, indx_y_max
        phi_due_object(indx_x_min:indx_x_max, j, nn) = phi(indx_x_min:indx_x_max, j)
     END DO

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  END DO

END SUBROUTINE CALCULATE_OBJECT_POTENTIALS_2D

!------------------------------------------------
!
SUBROUTINE CALCULATE_OBJECT_POTENTIAL_CHARGE_COEFFS

  USE ParallelOperationValues
  USE ExternalCircuit
  USE CurrentProblemValues
  USE BlockAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER nn      ! index of the object in object_charge_calculation array
  INTEGER noi     ! index of the object in whole_object array
  INTEGER pos     ! index of an object surface point in object_charge_calculation(nn)%control array
  INTEGER nseg    ! index of a segment of the object

  INTEGER segment_jstart, segment_jend
  INTEGER segment_start_type, segment_end_type
  INTEGER jstart, jend
  INTEGER segment_istart, segment_iend
  INTEGER istart, iend
  INTEGER i, j, m, nnode
  REAL(8) eps_ishifted, eps_jshifted, eps_shifted_4

  INTEGER bufsize
  REAL(8), ALLOCATABLE :: rbufer(:), rbufer2(:)
  INTEGER ALLOC_ERR

  IF (N_of_object_potentials_to_solve.EQ.0) RETURN

  object_charge_coeff = 0.0_8

  DO nn = 1, N_of_object_potentials_to_solve

     noi = object_charge_calculation(nn)%noi   ! number of object of interest

print *, Rank_of_process, noi

     object_charge_calculation(nn)%N_of_points_to_process = 0
     pos = 0

     IF (noi.LE.N_of_boundary_objects) THEN
! objects along the boundary

        DO nseg = 1, whole_object(noi)%number_of_segments

           IF (whole_object(noi)%segment(nseg)%istart.EQ.whole_object(noi)%segment(nseg)%iend) THEN
! vertical segment
              IF (whole_object(noi)%segment(nseg)%istart.LT.indx_x_min) CYCLE  ! segment not inside the block
              IF (whole_object(noi)%segment(nseg)%istart.GT.indx_x_max) CYCLE  ! segment not inside the block
              
! note, it is expected that jstart<=jend and istart<=iend

              segment_jstart = whole_object(noi)%segment(nseg)%jstart  !MIN(whole_object(noi)%segment(nseg)%jstart, whole_object(noi)%segment(nseg)%jend)
              segment_jend   = whole_object(noi)%segment(nseg)%jend    !MAX(whole_object(noi)%segment(nseg)%jstart, whole_object(noi)%segment(nseg)%jend)

              segment_start_type = whole_object(noi)%segment(nseg)%start_type
              segment_end_type   = whole_object(noi)%segment(nseg)%end_type

              IF (segment_jstart.GE.indx_y_max) CYCLE  ! segment not inside the block
              IF (segment_jend.LE.indx_y_min) CYCLE    ! segment not inside the block

              IF (Rank_of_process_below.LT.0) THEN
! there is no other block below this block
                 jstart = MAX(segment_jstart, indx_y_min)
              ELSE
! there is another block below this block, avoid overlapping point
                 jstart = MAX(segment_jstart, indx_y_min+1)
              END IF

              IF (Rank_of_process_above.LT.0) THEN
! there is no other block above this block
                 jend = MIN(segment_jend, indx_y_max)
              ELSE
! there is another block above this block, avoid overlapping point
                 jend = MIN(segment_jend, indx_y_max-1)
              END IF

              IF (jend.LT.jstart) CYCLE

              CALL ADD_MORE_CONTROL_POINTS(nn, MAX(0, jend-jstart+1))

              IF (whole_object(noi)%segment(nseg)%istart.EQ.indx_x_min) THEN  !-------------------------------------------
! wall on the left

                 i = indx_x_min

                 IF (jstart.EQ.segment_jstart) THEN
! block includes start end of the segment
                    pos = pos+1
                    j = jstart
                    SELECT CASE (segment_start_type)
                       CASE (END_FLAT)

                          CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted) ! i + 1/2                          !       eps_C = eps_B
                          CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) + 0.25_8, eps_shifted_4)          ! i+1/2, j+1/4, eps_D

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 5
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (8.0_8 * eps_shifted_4 + 3.0_8 * eps_ishifted) / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted * 3.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) *  eps_ishifted / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(4) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(4) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(4) = -DBLE(N_of_particles_cell) *  eps_shifted_4 * 6.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(5) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(5) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(5) =  DBLE(N_of_particles_cell) *  (eps_ishifted - 2.0_8 * eps_shifted_4) / 16.0_8

                       CASE (END_CORNER_CONCAVE)

                          CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) + 0.25_8, eps_shifted_4)          ! i+1/2, j+1/4, eps_D

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 3
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_shifted_4 * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.375_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.125_8

                       CASE (CONCAVE_CORNER)

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 0    !  node coordinates (indices) stored in the first element
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.125_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(1:1), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(1:1), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j

                       CASE (CONVEX_CORNER)

                          CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted)  ! i + 1/2
                          CALL SET_EPS_JSHIFTED(i,   j, eps_jshifted)  ! j - 1/2

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.375_8   ! the density is corrected

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                    END SELECT

                    jstart = jstart+1
                 END IF   !### IF (jstart.EQ.segment_jstart) THEN

                 IF (jend.EQ.segment_jend) THEN
! block includes last end of the segment
                    pos=pos+1
                    j = jend
                    SELECT CASE (segment_end_type)
                       CASE (END_FLAT)

                          CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted) ! i + 1/2                          !       eps_G = eps_H
                          CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) - 0.25_8, eps_shifted_4)          ! i+1/2, j-1/4, eps_F

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 5
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (8.0_8 * eps_shifted_4 + 3.0_8 * eps_ishifted) / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted * 3.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) *  eps_ishifted / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(4) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(4) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(4) = -DBLE(N_of_particles_cell) *  eps_shifted_4 * 6.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(5) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(5) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(5) =  DBLE(N_of_particles_cell) *  (eps_ishifted - 2.0_8 * eps_shifted_4) / 16.0_8

                       CASE (END_CORNER_CONCAVE)

                          CALL GET_EPS_IN_POINT(DBLE(i)+0.5_8, DBLE(j) - 0.25_8, eps_shifted_4)          ! i+1/2, j-1/4, eps_F

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 3
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_shifted_4 * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.375_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.125_8

                       CASE (CONCAVE_CORNER)

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 0    !  node coordinates (indices) stored in the first element
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.125_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(1:1), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(1:1), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j

                       CASE (CONVEX_CORNER)

                          CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted) ! i + 1/2
                          CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted) ! j + 1/2

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.375_8  ! the density is corrected

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                    END SELECT

                    jend = jend-1
                 END IF   !### IF (jend.EQ.segment_jend) THEN

                 DO j = jstart, jend
                  
                    CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted) ! i + 1/2

                    pos = pos+1

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_ishifted

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i+1
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                 END DO   !### DO j = jstart, jend

              ELSE IF (whole_object(noi)%segment(nseg)%istart.EQ.indx_x_max) THEN  !-------------------------------------------
! wall on the right

                 i = indx_x_max

                 IF (jstart.EQ.segment_jstart) THEN
! block includes start end of the segment
                    pos = pos+1
                    j = jstart
                    SELECT CASE (segment_start_type)
                       CASE (END_FLAT)

                          CALL SET_EPS_ISHIFTED(i, j, eps_ishifted)   ! i - 1/2                          !       eps_G = eps_H
                          CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) + 0.25_8, eps_shifted_4)          ! i-1/2, j+1/4, eps_F

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 5
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (8.0_8 * eps_shifted_4 + 3.0_8 * eps_ishifted) / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted * 3.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) *  eps_ishifted / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(4) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(4) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(4) = -DBLE(N_of_particles_cell) *  eps_shifted_4 * 6.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(5) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(5) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(5) =  DBLE(N_of_particles_cell) *  (eps_ishifted - 2.0_8 * eps_shifted_4) / 16.0_8

                       CASE (END_CORNER_CONCAVE)

                          CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) + 0.25_8, eps_shifted_4)          ! i-1/2, j+1/4, eps_F

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 3
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_shifted_4 * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.375_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.125_8

                       CASE (CONCAVE_CORNER)

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 0    !  node coordinates (indices) stored in the first element
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.125_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(1:1), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(1:1), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j

                       CASE (CONVEX_CORNER)

                          CALL SET_EPS_ISHIFTED(i, j, eps_ishifted)  ! i - 1/2
                          CALL SET_EPS_JSHIFTED(i, j, eps_jshifted)  ! j - 1/2

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.375_8  ! the density is corrected

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                    END SELECT

                    jstart = jstart+1
                 END IF   !### IF (jstart.EQ.segment_jstart) THEN

                 IF (jend.EQ.segment_jend) THEN
! block includes last end of the segment
                    pos=pos+1
                    j = jend
                    SELECT CASE (segment_end_type)
                       CASE (END_FLAT)

                          CALL SET_EPS_ISHIFTED(i, j, eps_ishifted)   ! i - 1/2                          !       eps_C = eps_B
                          CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) - 0.25_8, eps_shifted_4)          ! i-1/2, j-1/4, eps_D

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 5
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (8.0_8 * eps_shifted_4 + 3.0_8 * eps_ishifted) / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted * 3.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) *  eps_ishifted / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(4) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(4) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(4) = -DBLE(N_of_particles_cell) *  eps_shifted_4 * 6.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(5) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(5) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(5) =  DBLE(N_of_particles_cell) *  (eps_ishifted - 2.0_8 * eps_shifted_4) / 16.0_8

                       CASE (END_CORNER_CONCAVE)
                          CALL GET_EPS_IN_POINT(DBLE(i)-0.5_8, DBLE(j) - 0.25_8, eps_shifted_4)          ! i-1/2, j-1/4, eps_D

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 3
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_shifted_4 * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.375_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.125_8

                       CASE (CONCAVE_CORNER)

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 0    !  node coordinates (indices) stored in the first element
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.125_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(1:1), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(1:1), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j

                       CASE (CONVEX_CORNER)

                          CALL SET_EPS_ISHIFTED(i,   j, eps_ishifted) ! i - 1/2
                          CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted) ! j + 1/2

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.375_8   ! the density is corrected

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                    END SELECT

                    jend = jend-1
                 END IF   !### IF (jend.EQ.segment_jend) THEN

                 DO j = jstart, jend
                  
                    CALL SET_EPS_ISHIFTED(i, j, eps_ishifted) ! i - 1/2

                    pos = pos+1

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_ishifted

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i-1
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                 END DO   !### DO j = jstart, jend

              ELSE    ! IF (whole_object(noi)%segment(nseg)%istart.EQ.indx_x_min) THEN
! error
                 PRINT '("CALCULATE_OBJECT_POTENTIAL_CHARGE_COEFFS :: improbable ERROR-1")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF  ! IF (whole_object(noi)%segment(nseg)%istart.EQ.indx_x_min) THEN

           ELSE IF (whole_object(noi)%segment(nseg)%jstart.EQ.whole_object(noi)%segment(nseg)%jend) THEN   !-------------------------------------------------------------------
! horizontal segment

              IF (whole_object(noi)%segment(nseg)%jstart.LT.indx_y_min) CYCLE  ! segment not inside the block
              IF (whole_object(noi)%segment(nseg)%jstart.GT.indx_y_max) CYCLE  ! segment not inside the block
              
! note, it is expected that jstart<=jend and istart<=iend

              segment_istart = whole_object(noi)%segment(nseg)%istart  !MIN(whole_object(noi)%segment(nseg)%jstart, whole_object(noi)%segment(nseg)%jend)
              segment_iend   = whole_object(noi)%segment(nseg)%iend    !MAX(whole_object(noi)%segment(nseg)%jstart, whole_object(noi)%segment(nseg)%jend)

              segment_start_type = whole_object(noi)%segment(nseg)%start_type
              segment_end_type   = whole_object(noi)%segment(nseg)%end_type

              IF (segment_istart.GE.indx_x_max) CYCLE  ! segment not inside the block
              IF (segment_iend.LE.indx_x_min) CYCLE    ! segment not inside the block

              IF (Rank_of_process_left.LT.0) THEN
! there is no other block left from this block
                 istart = MAX(segment_istart, indx_x_min)
              ELSE
! there is another block left from this block, avoid overlapping point
                 istart = MAX(segment_istart, indx_x_min+1)
              END IF

              IF (Rank_of_process_right.LT.0) THEN
! there is no other block right from this block
                 iend = MIN(segment_iend, indx_x_max)
              ELSE
! there is another block above this block, avoid overlapping point
                 iend = MIN(segment_iend, indx_x_max-1)
              END IF

              IF (iend.LT.istart) CYCLE

              CALL ADD_MORE_CONTROL_POINTS(nn, MAX(0, iend-istart+1))

              IF (whole_object(noi)%segment(nseg)%jstart.EQ.indx_y_min) THEN  !-------------------------------------------
! wall on the bottom

                 j = indx_y_min

                 IF (istart.EQ.segment_istart) THEN
! block includes start end of the segment
                    pos = pos+1
                    i = istart
                    SELECT CASE (segment_start_type)
                       CASE (END_FLAT)

                          CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted) ! j + 1/2                          !       eps_G = eps_H
                          CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j) + 0.5_8, eps_shifted_4)          ! i+1/4, j+1/2, eps_F

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 5
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (8.0_8 * eps_shifted_4 + 3.0_8 * eps_jshifted) / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted * 3.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) *  eps_jshifted / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(4) = i
                          object_charge_calculation(nn)%control(pos)%use_j(4) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(4) = -DBLE(N_of_particles_cell) *  eps_shifted_4 * 6.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(5) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(5) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(5) =  DBLE(N_of_particles_cell) *  (eps_jshifted - 2.0_8 * eps_shifted_4) / 16.0_8

                       CASE (END_CORNER_CONCAVE)

                          CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j) + 0.5_8, eps_shifted_4)          ! i+1/4, j+1/2, eps_F

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 3
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_shifted_4 * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.375_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.125_8

                       CASE (CONCAVE_CORNER)

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 0    !  node coordinates (indices) stored in the first element
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.125_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(1:1), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(1:1), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j

                       CASE (CONVEX_CORNER)

                          CALL SET_EPS_ISHIFTED(i,   j, eps_ishifted)  ! i - 1/2
                          CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted)  ! j + 1/2

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.375_8   ! the density is corrected

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                    END SELECT

                    istart = istart+1
                 END IF   !### IF (istart.EQ.segment_istart) THEN

                 IF (iend.EQ.segment_iend) THEN
! block includes last end of the segment
                    pos=pos+1
                    i = iend
                    SELECT CASE (segment_end_type)
                       CASE (END_FLAT)
 
                          CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted) ! j + 1/2                          !       eps_C = eps_B
                          CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) + 0.5_8, eps_shifted_4)          ! i-1/4, j+1/2, eps_D

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 5
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (8.0_8 * eps_shifted_4 + 3.0_8 * eps_jshifted) / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted * 3.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) *  eps_jshifted / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(4) = i
                          object_charge_calculation(nn)%control(pos)%use_j(4) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(4) = -DBLE(N_of_particles_cell) *  eps_shifted_4 * 6.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(5) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(5) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(5) =  DBLE(N_of_particles_cell) *  (eps_jshifted - 2.0_8 * eps_shifted_4) / 16.0_8

                       CASE (END_CORNER_CONCAVE)

                          CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) + 0.5_8, eps_shifted_4)          ! i-1/4, j+1/2, eps_D

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 3
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_shifted_4 * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.375_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.125_8

                       CASE (CONCAVE_CORNER)

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 0    !  node coordinates (indices) stored in the first element
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.125_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(1:1), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(1:1), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j

                       CASE (CONVEX_CORNER)

                          CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted) ! i + 1/2
                          CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted) ! j + 1/2

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.375_8   ! the density is corrected

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j+1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                    END SELECT

                    iend = iend-1
                 END IF   !### IF (iend.EQ.segment_iend) THEN

                 DO i = istart, iend
                  
                    CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted) ! j + 1/2

                    pos = pos+1

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_jshifted

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j+1
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                 END DO   !### DO i = istart, iend

              ELSE IF (whole_object(noi)%segment(nseg)%istart.EQ.indx_x_max) THEN  !-------------------------------------------
! wall above

                 j = indx_y_max

                 IF (istart.EQ.segment_istart) THEN
! block includes start end of the segment
                    pos = pos+1
                    i = istart
                    SELECT CASE (segment_start_type)
                       CASE (END_FLAT)

                          CALL SET_EPS_JSHIFTED(i, j, eps_jshifted)   ! j - 1/2                          !       eps_C = eps_B
                          CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j) - 0.5_8, eps_shifted_4)          ! i+1/4, j-1/2, eps_D

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 5
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (8.0_8 * eps_shifted_4 + 3.0_8 * eps_jshifted) / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted * 3.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) *  eps_jshifted / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(4) = i
                          object_charge_calculation(nn)%control(pos)%use_j(4) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(4) = -DBLE(N_of_particles_cell) *  eps_shifted_4 * 6.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(5) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(5) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(5) =  DBLE(N_of_particles_cell) *  (eps_jshifted - 2.0_8 * eps_shifted_4) / 16.0_8

                       CASE (END_CORNER_CONCAVE)

                          CALL GET_EPS_IN_POINT(DBLE(i)+0.25_8, DBLE(j) - 0.5_8, eps_shifted_4)          ! i+1/4, j-1/2, eps_D

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 3
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_shifted_4 * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.375_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.125_8

                       CASE (CONCAVE_CORNER)

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 0    !  node coordinates (indices) stored in the first element
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.125_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(1:1), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(1:1), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j

                       CASE (CONVEX_CORNER)

                          CALL SET_EPS_ISHIFTED(i, j, eps_ishifted)  ! i - 1/2
                          CALL SET_EPS_JSHIFTED(i, j, eps_jshifted)  ! j - 1/2

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.375_8   ! the density is corrected

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                    END SELECT

                    istart = istart+1
                 END IF   !### IF (istart.EQ.segment_istart) THEN

                 IF (iend.EQ.segment_iend) THEN
! block includes last end of the segment
                    pos=pos+1
                    i = iend
                    SELECT CASE (segment_end_type)
                       CASE (END_FLAT)
                    
                          CALL SET_EPS_JSHIFTED(i, j, eps_jshifted)   ! j - 1/2                          !       eps_G = eps_H
                          CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) - 0.5_8, eps_shifted_4)          ! i-1/4, j-1/2, eps_F

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 5
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (8.0_8 * eps_shifted_4 + 3.0_8 * eps_jshifted) / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted * 3.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i+1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) *  eps_jshifted / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(4) = i
                          object_charge_calculation(nn)%control(pos)%use_j(4) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(4) = -DBLE(N_of_particles_cell) *  eps_shifted_4 * 6.0_8 / 16.0_8

                          object_charge_calculation(nn)%control(pos)%use_i(5) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(5) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(5) =  DBLE(N_of_particles_cell) *  (eps_jshifted - 2.0_8 * eps_shifted_4) / 16.0_8

                       CASE (END_CORNER_CONCAVE)

                          CALL GET_EPS_IN_POINT(DBLE(i)-0.25_8, DBLE(j) - 0.5_8, eps_shifted_4)          ! i-1/4, j-1/2, eps_F

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 3
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.25_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_shifted_4 * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.375_8

                          object_charge_calculation(nn)%control(pos)%use_i(3) = i-1
                          object_charge_calculation(nn)%control(pos)%use_j(3) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(3) = -DBLE(N_of_particles_cell) * eps_shifted_4 * 0.125_8

                       CASE (CONCAVE_CORNER)

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 0    !  node coordinates (indices) stored in the first element
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.125_8

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(1:1), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(1:1), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j

                       CASE (CONVEX_CORNER)

                          CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted) ! i + 1/2
                          CALL SET_EPS_JSHIFTED(i,   j, eps_jshifted) ! j - 1/2

                          object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                          object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.375_8   ! the density is corrected

                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                          ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                          object_charge_calculation(nn)%control(pos)%use_i(1) = i
                          object_charge_calculation(nn)%control(pos)%use_j(1) = j
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                          object_charge_calculation(nn)%control(pos)%use_i(2) = i
                          object_charge_calculation(nn)%control(pos)%use_j(2) = j-1
                          object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                    END SELECT

                    iend = iend-1
                 END IF   !### IF (iend.EQ.segment_iend) THEN

                 DO i = istart, iend
                  
                    CALL SET_EPS_JSHIFTED(i, j, eps_jshifted) ! j - 1/2

                    pos = pos+1

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_jshifted

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j-1
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                 END DO   !### DO j = jstart, jend

              ELSE    ! IF (whole_object(noi)%segment(nseg)%jstart.EQ.indx_y_min) THEN
! error
                 PRINT '("CALCULATE_OBJECT_POTENTIAL_CHARGE_COEFFS :: improbable ERROR-2")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF  ! IF (whole_object(noi)%segment(nseg)%jstart.EQ.indx_y_min) THEN

           ELSE     !### IF (whole_object(noi)%segment(nseg)%istart.EQ.whole_object(noi)%segment(nseg)%iend) THEN
! error
              PRINT '("CALCULATE_OBJECT_POTENTIAL_CHARGE_COEFFS :: improbable ERROR-3")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF   !### IF (whole_object(noi)%segment(nseg)%istart.EQ.whole_object(noi)%segment(nseg)%iend) THEN
        
        END DO  !### DO nseg = 1, whole_object(noi)%number_of_segments

     ELSE    !### IF (noi.LE.N_of_boundary_objects) THEN
! inner object
! for inner object, ends of each segment are CONVEX_CORNER

        DO nseg = 1, whole_object(noi)%number_of_segments

           SELECT CASE (nseg)

              CASE (1)
! segment 1 is VERTICAL and is the left side of the inner object [wall on the right] -------------------------------------------------------------------------------------------

! here use GE and LE to make sure that a point on object surface is conisidered by one block only
                 IF (whole_object(noi)%segment(nseg)%istart.LE.indx_x_min) CYCLE  ! segment not inside the block
                 IF (whole_object(noi)%segment(nseg)%istart.GE.indx_x_max) CYCLE  ! segment not inside the block
              
                 segment_jstart = whole_object(noi)%segment(nseg)%jstart
                 segment_jend   = whole_object(noi)%segment(nseg)%jend

                 jstart = MAX(segment_jstart, indx_y_min+1)
                 jend   = MIN(segment_jend,   indx_y_max-1)

                 IF (jend.LT.jstart) CYCLE

                 CALL ADD_MORE_CONTROL_POINTS(nn, MAX(0, jend-jstart+1))

                 i = whole_object(noi)%segment(nseg)%istart

                 IF (jstart.EQ.segment_jstart) THEN
! block includes start end of the segment ::::::: LEFT VERTICAL EDGE, BOTTOM CORNER
                    pos = pos + 1
                    j = jstart
           
                    CALL SET_EPS_ISHIFTED(i, j, eps_ishifted)  ! i - 1/2
                    CALL SET_EPS_JSHIFTED(i, j, eps_jshifted)  ! j - 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2 ! 3
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8 ! -0.375_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i-1
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                    jstart = jstart + 1
                 END IF

                 IF (jend.EQ.segment_jend) THEN
! block includes last end of the segment ::::::: LEFT VERTICAL EDGE, TOP CORNER
                    pos = pos + 1
                    j = jend
           
                    CALL SET_EPS_ISHIFTED(i, j,   eps_ishifted)  ! i - 1/2
                    CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted)  ! j + 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8 ! -0.375_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i-1
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                    jend = jend - 1
                 END IF

                 DO j = jstart, jend
! inner points :::::::::::::::::: LEFT VERTICAL EDGE

                    pos = pos+1
                    CALL SET_EPS_ISHIFTED(i, j, eps_ishifted) ! i - 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -1.0_8  ! -0.5_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_ishifted

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i-1
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                 END DO   !### DO j = jstart, jend

              CASE (2)
! segment 2 is HORIZONTAL and is the top side of the inner object [wall on the bottom] -------------------------------------------------------------------------------------------

! here use GE and LE to make sure that a point on object surface is conisidered by one block only
                 IF (whole_object(noi)%segment(nseg)%jstart.LE.indx_y_min) CYCLE  ! segment not inside the block
                 IF (whole_object(noi)%segment(nseg)%jstart.GE.indx_y_max) CYCLE  ! segment not inside the block
              
                 segment_istart = whole_object(noi)%segment(nseg)%istart
                 segment_iend   = whole_object(noi)%segment(nseg)%iend

                 istart = MAX(segment_istart, indx_x_min+1)
                 iend   = MIN(segment_iend,   indx_x_max-1)

                 IF (iend.LT.istart) CYCLE

                 CALL ADD_MORE_CONTROL_POINTS(nn, MAX(0, iend-istart+1))

                 j = whole_object(noi)%segment(nseg)%jstart

                 IF (istart.EQ.segment_istart) THEN
! block includes start end of the segment ::::::::::::: TOP HORIZONTAL EDGE, LEFT CORNER
                    pos = pos + 1
                    i = istart
  
                    CALL SET_EPS_ISHIFTED(i,   j, eps_ishifted)  ! i - 1/2
                    CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted)  ! j + 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8 ! -0.375_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j+1
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                    istart = istart + 1
                 END IF

                 IF (iend.EQ.segment_iend) THEN
! block includes last end of the segment ::::::::::::: TOP HORIZONTAL EDGE, RIGHT CORNER
                    pos = pos + 1
                    i = iend
  
                    CALL SET_EPS_ISHIFTED(i+1,   j, eps_ishifted)  ! i + 1/2
                    CALL SET_EPS_JSHIFTED(i,   j+1, eps_jshifted)  ! j + 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8 ! -0.375_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j+1
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                    iend = iend - 1
                 END IF

                 DO i = istart, iend
! inner points :::::::::::::::::: TOP HORIZONTAL EDGE

                    pos = pos+1
                    CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted) ! j + 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -1.0_8  ! -0.5_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_jshifted

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j+1
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                 END DO   !### DO i = istart, iend

              CASE (3)
! segment 3 is VERTICAL and is the right side of the inner object [wall on the left] -------------------------------------------------------------------------------------------

! here use GE and LE to make sure that a point on object surface is conisidered by one block only
                 IF (whole_object(noi)%segment(nseg)%istart.LE.indx_x_min) CYCLE  ! segment not inside the block
                 IF (whole_object(noi)%segment(nseg)%istart.GE.indx_x_max) CYCLE  ! segment not inside the block
              
                 segment_jstart = whole_object(noi)%segment(nseg)%jstart
                 segment_jend   = whole_object(noi)%segment(nseg)%jend

                 jstart = MAX(segment_jstart, indx_y_min+1)
                 jend   = MIN(segment_jend,   indx_y_max-1)

                 IF (jend.LT.jstart) CYCLE

                 CALL ADD_MORE_CONTROL_POINTS(nn, MAX(0, jend-jstart+1))

                 i = whole_object(noi)%segment(nseg)%istart

                 IF (jstart.EQ.segment_jstart) THEN
! block includes start end of the segment ::::::: RIGHT VERTICAL EDGE, BOTTOM CORNER
                    pos = pos + 1
                    j = jstart
           
                    CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted)  ! i + 1/2
                    CALL SET_EPS_JSHIFTED(i,   j, eps_jshifted)  ! j - 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8 ! -0.375_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i+1
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                    jstart = jstart + 1
                 END IF

                 IF (jend.EQ.segment_jend) THEN
! block includes last end of the segment ::::::: RIGHT VERTICAL EDGE, TOP CORNER
                    pos = pos + 1
                    j = jend
           
                    CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted)  ! i + 1/2
                    CALL SET_EPS_JSHIFTED(i, j+1, eps_jshifted)  ! j + 1/2
                    
                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8 ! -0.375_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i+1
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                    jend = jend - 1
                 END IF

                 DO j = jstart, jend
! inner points :::::::::::::::::: RIGHT VERTICAL EDGE

                    pos = pos + 1
                    CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted) ! i + 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -1.0_8  ! -0.5_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_ishifted

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i+1
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_ishifted

                 END DO   !### DO j = jstart, jend

              CASE (4)
! segment 4 is HORIZONTAL and is the bottom side of the inner object [wall on the top] -------------------------------------------------------------------------------------------

! here use GE and LE to make sure that a point on object surface is conisidered by one block only
                 IF (whole_object(noi)%segment(nseg)%jstart.LE.indx_y_min) CYCLE  ! segment not inside the block
                 IF (whole_object(noi)%segment(nseg)%jstart.GE.indx_y_max) CYCLE  ! segment not inside the block
              
                 segment_istart = whole_object(noi)%segment(nseg)%istart
                 segment_iend   = whole_object(noi)%segment(nseg)%iend

                 istart = MAX(segment_istart, indx_x_min+1)
                 iend   = MIN(segment_iend,   indx_x_max-1)

                 IF (iend.LT.istart) CYCLE

                 CALL ADD_MORE_CONTROL_POINTS(nn, MAX(0, iend-istart+1))

                 j = whole_object(noi)%segment(nseg)%jstart

                 IF (istart.EQ.segment_istart) THEN
! block includes start end of the segment ::::::::::::: BOTTOM HORIZONTAL EDGE, LEFT CORNER
                    pos = pos + 1
                    i = istart
  
                    CALL SET_EPS_ISHIFTED(i, j, eps_ishifted)  ! i - 1/2
                    CALL SET_EPS_JSHIFTED(i, j, eps_jshifted)  ! j - 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8 ! -0.375_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j-1
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                    istart = istart + 1
                 END IF

                 IF (iend.EQ.segment_iend) THEN
! block includes last end of the segment ::::::::::::: BOTTOM HORIZONTAL EDGE, RIGHT CORNER
                    pos = pos + 1
                    i = iend
  
                    CALL SET_EPS_ISHIFTED(i+1, j, eps_ishifted)  ! i + 1/2
                    CALL SET_EPS_JSHIFTED(i,   j, eps_jshifted)  ! j - 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -0.5_8 ! -0.375_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * (eps_ishifted + eps_jshifted) * 0.5_8

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j-1
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                    iend = iend - 1
                 END IF

                 DO i = istart, iend
! inner points :::::::::::::::::: BOTTOM HORIZONTAL EDGE

                    pos = pos + 1
                    CALL SET_EPS_JSHIFTED(i, j, eps_jshifted) ! j - 1/2

                    object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use = 2
                    object_charge_calculation(nn)%control(pos)%use_alpha_rho = -1.0_8  ! -0.5_8 if the density is corrected

                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_i(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_j(        1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)
                    ALLOCATE(object_charge_calculation(nn)%control(pos)%use_alpha_phi(1:object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use), STAT = ALLOC_ERR)

                    object_charge_calculation(nn)%control(pos)%use_i(1) = i
                    object_charge_calculation(nn)%control(pos)%use_j(1) = j
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(1) =  DBLE(N_of_particles_cell) * eps_jshifted

                    object_charge_calculation(nn)%control(pos)%use_i(2) = i
                    object_charge_calculation(nn)%control(pos)%use_j(2) = j-1
                    object_charge_calculation(nn)%control(pos)%use_alpha_phi(2) = -DBLE(N_of_particles_cell) * eps_jshifted

                 END DO   !### DO j = jstart, jend

           END SELECT   !### SELECT CASE (nseg)

        END DO   !###   DO nseg = 1, whole_object(noi)%number_of_segments

     END IF  !### IF (noi.LE.N_of_boundary_objects) THEN

     DO m = 1, N_of_object_potentials_to_solve
        DO pos = 1, object_charge_calculation(nn)%N_of_points_to_process
           DO nnode = 1, object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use
              object_charge_coeff(m,nn) = object_charge_coeff(m,nn) + &
                         &                object_charge_calculation(nn)%control(pos)%use_alpha_phi(nnode) * &
                         & phi_due_object(object_charge_calculation(nn)%control(pos)%use_i(nnode), &
                         &                object_charge_calculation(nn)%control(pos)%use_j(nnode), m)
           END DO
        END DO
     END DO

  END DO   !### DO nn = 1, N_of_object_potentials_to_solve

! assemble contributions from all processes in the zero-rank process

  bufsize = N_of_object_potentials_to_solve * N_of_object_potentials_to_solve
  ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer2(1:bufsize), STAT=ALLOC_ERR)
  pos=0
  DO nn = 1, N_of_object_potentials_to_solve
     DO m = 1, N_of_object_potentials_to_solve
        pos = pos+1
        rbufer(pos) = object_charge_coeff(m,nn)
     END DO
  END DO

  CALL MPI_REDUCE(rbufer, rbufer2, bufsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  IF (Rank_of_process.EQ.0) THEN
     pos=0
     DO nn = 1, N_of_object_potentials_to_solve
        DO m = 1, N_of_object_potentials_to_solve
           pos = pos+1
           object_charge_coeff(m,nn) = rbufer2(pos)
        END DO
     END DO
  END IF

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)

END SUBROUTINE CALCULATE_OBJECT_POTENTIAL_CHARGE_COEFFS

!----------------------------------------
! this procedure is invoked only a few times in the beginning of the simulation, so it does not have to be pretty
SUBROUTINE ADD_MORE_CONTROL_POINTS(nn, N_to_add)

  USE ExternalCircuit

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(IN) :: nn, N_to_add

  INTEGER ibuflength, dbuflength

  INTEGER ALLOC_ERR
  INTEGER N_points, i
  INTEGER, ALLOCATABLE :: ibufer(:)
  REAL(8), ALLOCATABLE :: dbufer(:)

  INTEGER ipos, dpos
  INTEGER N_to_use
  
  IF (N_to_add.LE.0) RETURN

  IF (.NOT.ALLOCATED(object_charge_calculation(nn)%control)) THEN
     ALLOCATE(object_charge_calculation(nn)%control(1:N_to_add), STAT = ALLOC_ERR)
     object_charge_calculation(nn)%N_of_points_to_process = N_to_add
     RETURN
  END IF

  N_points = object_charge_calculation(nn)%N_of_points_to_process

! save present array into the buffers

! define size of buffers
  ibuflength = 0
  dbuflength = 0
  DO i = 1, N_points
     ibuflength = ibuflength + 1 + 2 * MAX(1, object_charge_calculation(nn)%control(i)%N_of_nodes_to_use)   ! if N_of_nodes_to_use=0 we still have to save 3 integers
     dbuflength = dbuflength + 1 + object_charge_calculation(nn)%control(i)%N_of_nodes_to_use
  END DO

  ALLOCATE(ibufer(1:ibuflength), STAT=ALLOC_ERR)
  ALLOCATE(dbufer(1:dbuflength), STAT=ALLOC_ERR)

  ipos = 1
  dpos = 1
  DO i = 1, N_points

     N_to_use = object_charge_calculation(nn)%control(i)%N_of_nodes_to_use

! fool proof
     IF ((ipos+2*MAX(N_to_use,1)).GT.ibuflength) THEN
! error
        PRINT '("ADD_MORE_CONTROL_POINTS :: ERROR ipos")'
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF

     IF ((dpos+N_to_use).GT.ibuflength) THEN
! error
        PRINT '("ADD_MORE_CONTROL_POINTS :: ERROR dpos")'
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF

     ibufer(ipos) = N_to_use
     ibufer(ipos+1                :ipos+MAX(N_to_use,1)                ) = object_charge_calculation(nn)%control(i)%use_i(1:MAX(N_to_use,1))
     ibufer(ipos+MAX(N_to_use,1)+1:ipos+MAX(N_to_use,1)+MAX(N_to_use,1)) = object_charge_calculation(nn)%control(i)%use_j(1:MAX(N_to_use,1))
     ipos = ipos + 2 * MAX(N_to_use,1) + 1

     dbufer(dpos) = object_charge_calculation(nn)%control(i)%use_alpha_rho
     IF (N_to_use.GT.0) dbufer(dpos+1:dpos+N_to_use) = object_charge_calculation(nn)%control(i)%use_alpha_phi(1:N_to_use)
     dpos = dpos + 1 + N_to_use

  END DO

! resize the array

! delete the array
  DO i = 1, N_points
     IF (ALLOCATED(object_charge_calculation(nn)%control(i)%use_i))         DEALLOCATE(object_charge_calculation(nn)%control(i)%use_i, STAT=ALLOC_ERR)
     IF (ALLOCATED(object_charge_calculation(nn)%control(i)%use_j))         DEALLOCATE(object_charge_calculation(nn)%control(i)%use_j, STAT=ALLOC_ERR)
     IF (ALLOCATED(object_charge_calculation(nn)%control(i)%use_alpha_phi)) DEALLOCATE(object_charge_calculation(nn)%control(i)%use_alpha_phi, STAT=ALLOC_ERR)
  END DO
  DEALLOCATE(object_charge_calculation(nn)%control, STAT = ALLOC_ERR)

  object_charge_calculation(nn)%N_of_points_to_process = object_charge_calculation(nn)%N_of_points_to_process + N_to_add
  ALLOCATE(object_charge_calculation(nn)%control(1:object_charge_calculation(nn)%N_of_points_to_process), STAT = ALLOC_ERR)

! restore data from the buffers

  ipos = 1
  dpos = 1
  DO i = 1, N_points
     N_to_use = ibufer(ipos)
     object_charge_calculation(nn)%control(i)%N_of_nodes_to_use  = N_to_use

     ALLOCATE(object_charge_calculation(nn)%control(i)%use_i(1:MAX(1,N_to_use)), STAT=ALLOC_ERR)
     ALLOCATE(object_charge_calculation(nn)%control(i)%use_j(1:MAX(1,N_to_use)), STAT=ALLOC_ERR)
     IF (N_to_use.GT.0) ALLOCATE(object_charge_calculation(nn)%control(i)%use_alpha_phi(1:N_to_use), STAT=ALLOC_ERR)

     object_charge_calculation(nn)%control(i)%use_i(1:MAX(N_to_use,1)) = ibufer(ipos+1                :ipos+MAX(N_to_use,1)                )
     object_charge_calculation(nn)%control(i)%use_j(1:MAX(N_to_use,1)) = ibufer(ipos+MAX(N_to_use,1)+1:ipos+MAX(N_to_use,1)+MAX(N_to_use,1))
     ipos = ipos + 2 * MAX(N_to_use,1) + 1

     object_charge_calculation(nn)%control(i)%use_alpha_rho = dbufer(dpos)
     IF (N_to_use.GT.0) object_charge_calculation(nn)%control(i)%use_alpha_phi(1:N_to_use) = dbufer(dpos+1:dpos+N_to_use)
     dpos = dpos + 1 + N_to_use

  END DO

  DEALLOCATE(ibufer, STAT = ALLOC_ERR)
  DEALLOCATE(dbufer, STAT = ALLOC_ERR)

  RETURN

END SUBROUTINE ADD_MORE_CONTROL_POINTS

!------------------------------------------------
!
SUBROUTINE SOLVE_EXTERNAL_CONTOUR

  USE ParallelOperationValues
  USE ExternalCircuit
  USE CurrentProblemValues
  USE BlockAndItsBoundaries
  USE IonParticles, ONLY : N_spec, Qs
  USE Diagnostics, ONLY : Save_probes_e_data_T_cntr, N_of_probes_block, Probe_params_block_list, probe_F_block

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER nn  ! index of object
  INTEGER pos ! index of point on the surface of the object
  INTEGER i, j !, s
  INTEGER nnode

  INTEGER ALLOC_ERR
  REAL(8), ALLOCATABLE :: rbufer(:), rbufer2(:)

  REAL(8), ALLOCATABLE :: a(:,:)
  REAL(8), ALLOCATABLE :: rhs(:)
  REAL(8), ALLOCATABLE :: dQ_full(:)

  INTEGER circuit_type

  REAL(8) factor_C
  REAL(8) dU_source

  REAL(8) d0, d1, d2

!  INTEGER noi
  INTEGER m
  INTEGER npb

! function
  REAL(8) ECPS_Voltage

  IF (N_of_object_potentials_to_solve.EQ.0) RETURN

! create coefficients of the linear system

  DO nn = 1, N_of_object_potentials_to_solve

     object_charge_coeff(0,nn) = 0.0_8

     DO pos = 1, object_charge_calculation(nn)%N_of_points_to_process

        i = object_charge_calculation(nn)%control(pos)%use_i(1)   ! coordinate indices of point pos (point on surface  of object)
        j = object_charge_calculation(nn)%control(pos)%use_j(1)   ! 

        object_charge_coeff(0,nn) = object_charge_coeff(0,nn) + &
                                  & object_charge_calculation(nn)%control(pos)%use_alpha_rho * &
                                  & (rho_i(i,j) - rho_e(i,j))

        DO nnode = 1, object_charge_calculation(nn)%control(pos)%N_of_nodes_to_use
           
           i = object_charge_calculation(nn)%control(pos)%use_i(nnode)
           j = object_charge_calculation(nn)%control(pos)%use_j(nnode)

           object_charge_coeff(0,nn) = object_charge_coeff(0,nn) + &
                                     & object_charge_calculation(nn)%control(pos)%use_alpha_phi(nnode) * phi(i, j)  ! here phi is 2d potential obtained with 
                                                                                                                    ! zero potentials of objects 1:N_of_object_potentials_to_solve
        END DO
     END DO
  END DO

! assemble contributions from all processes in the zero-rank process

  ALLOCATE(rbufer(1:N_of_object_potentials_to_solve), STAT=ALLOC_ERR)
  ALLOCATE(rbufer2(1:N_of_object_potentials_to_solve), STAT=ALLOC_ERR)
  DO nn = 1, N_of_object_potentials_to_solve
     rbufer(nn) = object_charge_coeff(0,nn)
  END DO

  CALL MPI_REDUCE(rbufer, rbufer2, N_of_object_potentials_to_solve, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  IF (Rank_of_process.EQ.0) THEN

! presently the circuit equation is solved by process with rank zero

     DO nn = 1, N_of_object_potentials_to_solve
        object_charge_coeff(0,nn) = rbufer2(nn)
     END DO

! now we have charge_of_object_nn = object_charge_coeff(0,nn) + 
!                                   object_charge_coeff(1,nn) * potential_of_object(1) + 
!                                   object_charge_coeff(2,nn) * potential_of_object(2) + 
!                                   object_charge_coeff(3,nn) * potential_of_object(3) + ...

! create coefficients of system of equations to find potential_of_object(1:N_of_object_potentials_to_solve)
!
! 
! !! a(1,1) a(1,2) a(1,3) ...  !! f(1)   rhs(1)
! !! a(2,1) a(2,2) a(2,3) ...  !! f(2) = rhs(2)
! !! a(3,1) a(3,2) a(3,3) ...  !! f(3)   rhs(3)
! !! ...    ...    ...    ...  !! ...    ...

! here it is done for one circuit equation
! d phi_1 / dt = d U_source / dt - (1/C) * (d Q_full_1 / dt - d Q_plasma_1 / dt)
! where phi_1 is potential of electrode 1
! U_source is the voltage of the power source (given function of time)
! C is the capacity of the capacitor
! Q_full_1 is the full charge of electrode 1
! Q_plasma_1 is th echarge of electrode 1 due to plasma electrons and ions collided with this electrode and emission of electrons
! the system is a rf discharge with two electrodes, 
! one electrode (#1) is connected to the rf voltage source, capacitor, and then the ground
! the other electrode is simply grounded
! N_of_object_potentials_to_solve = 1


     ALLOCATE(      a(1:N_of_object_potentials_to_solve,1:N_of_object_potentials_to_solve), STAT=ALLOC_ERR)
     ALLOCATE(    rhs(1:N_of_object_potentials_to_solve), STAT=ALLOC_ERR)
     ALLOCATE(dQ_full(1:N_of_object_potentials_to_solve), STAT=ALLOC_ERR)
  
     circuit_type = 1   ! or 2 or 3 etc later make this a part of an input file

     SELECT CASE (circuit_type)
        CASE (1)
! one rf electrode connected in series to rf voltage source, capacitor, ground
! one electrode grounded

           factor_C = eps_0_Fm / (capacitor_C_F(1) * DBLE(N_of_particles_cell))

           a(1,1) = 1.0_8 + factor_C * object_charge_coeff(1,1)

!           dU_source = source_U * ( SIN(source_omega * T_cntr + source_phase) - &     ! source voltage at t^n
!                                  & SIN(source_omega * (T_cntr-1) + source_phase) )   ! source voltage at t^{n-1}
           dU_source = ECPS_Voltage(1, T_cntr) - ECPS_Voltage(1, T_cntr-1)

!     noi = object_charge_calculation(1)%noi

           rhs(1) = potential_of_object(1) + &                                        ! here potential_of_object(1) is at t^{n-1}
                  & dU_source - &
                  & factor_C * (object_charge_coeff(0,1) - charge_of_object(1)) + &   ! here charge_of_object(1) is at t^{n-1}
                  & factor_C * dQ_plasma_of_object(1)                                 ! charge of plasma particles deposited at the object during one timestep

! solve the linear system

           IF (a(1,1).NE.0.0_8) THEN
              potential_of_object(1) = rhs(1) / a(1,1)
           ELSE
! error
              PRINT '("error zero a(1,1)")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF

        CASE (2)
! one electrode with floating potential
! potentials of all other electrodes given

!           charge_of_object(1) + dQ_plasma_of_object(1) = object_charge_coeff(0,1) + object_charge_coeff(1,1) * potential_of_object(1)

           a(1,1) = object_charge_coeff(1,1)
           rhs(1) = charge_of_object(1) + dQ_plasma_of_object(1) - object_charge_coeff(0,1)

           IF (a(1,1).NE.0.0_8) THEN
              potential_of_object(1) = rhs(1) / a(1,1)
           ELSE
! error
              PRINT '("error zero a(1,1)")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF

        CASE (3)
! two electrodes with floating potential
! potentials of all other electrodes given

!           charge_of_object(1) + dQ_plasma_of_object(1) = object_charge_coeff(0,1) + object_charge_coeff(1,1) * potential_of_object(1) + object_charge_coeff(2,1) * potential_of_object(2)
!           charge_of_object(2) + dQ_plasma_of_object(2) = object_charge_coeff(0,2) + object_charge_coeff(1,2) * potential_of_object(1) + object_charge_coeff(2,2) * potential_of_object(2)

           a(1,1) = object_charge_coeff(1,1)
           a(1,2) = object_charge_coeff(2,1)
           rhs(1) = charge_of_object(1) + dQ_plasma_of_object(1) - object_charge_coeff(0,1)

           a(2,1) = object_charge_coeff(1,2)
           a(2,2) = object_charge_coeff(2,2)
           rhs(2) = charge_of_object(2) + dQ_plasma_of_object(2) - object_charge_coeff(0,2)

! solve the linear system

           d0 = a(1,1) * a(2,2) - a(1,2) * a(2,1)
           d1 = rhs(1) * a(2,2) - a(1,2) * rhs(2)
           d2 = a(1,1) * rhs(2) - rhs(1) * a(2,1)
           IF (d0.NE.0.0_8) THEN
              potential_of_object(1) = d1 / d0
              potential_of_object(2) = d2 / d0
           ELSE
! error
              PRINT '("error zero d0")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF

     END SELECT

!     dQ_full(1) = -charge_of_object(1)
     dQ_full = -charge_of_object

! update the full object charge
     DO nn = 1, N_of_object_potentials_to_solve
        charge_of_object(nn) = object_charge_coeff(0,nn)
        DO m = 1, N_of_object_potentials_to_solve
           charge_of_object(nn) = charge_of_object(nn) + &
                             & object_charge_coeff(m,nn) * potential_of_object(m)
        END DO
     END DO

!     dQ_full(1) = charge_of_object(1) + dQ_full(1)
     dQ_full = charge_of_object + dQ_full

     OPEN (21, FILE = 'history_ext_circuit.dat', POSITION = 'APPEND')
     WRITE (21, '(2x,i9,8(2x,e14.7))') &
          & T_cntr, &                                                             ! 1
          & T_cntr * delta_t_s * 1.0d9, &                                         ! 2
!          & source_U * SIN(source_omega * T_cntr + source_phase) * F_scale_V, &        ! 3
          & ECPS_Voltage(1, T_cntr) * F_scale_V, &        ! 3
          & potential_of_object(1) * F_scale_V, &                                      ! 4
          & charge_of_object(1), &                      ! 5
          & dQ_full(1), &                               ! 6
          & dQ_plasma_of_object(1), &                             ! 7
          & (dQ_full(1) - dQ_plasma_of_object(1)) * (e_Cl * N_plasma_m3 * delta_x_m**2 / (DBLE(N_of_particles_cell) * delta_t_s)), &    ! 8 ! current in external circuit
          &             - dQ_plasma_of_object(1)  * (e_Cl * N_plasma_m3 * delta_x_m**2 / (DBLE(N_of_particles_cell) * delta_t_s))       ! 9 ! current in plasma
     CLOSE (21, STATUS = 'KEEP')

  END IF   !### IF (Rank_of_process.EQ.0) THEN

  CALL MPI_BCAST(potential_of_object, N_of_object_potentials_to_solve, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

! at this stage we know the potentials of the electrodes
! include their contribution to the full 2d potential profile

  DO nn = 1, N_of_object_potentials_to_solve
     DO j = indx_y_min, indx_y_max
        DO i = indx_x_min, indx_x_max
           phi(i,j) = phi(i,j) + potential_of_object(nn) * phi_due_object(i,j,nn)
        END DO
     END DO
  END DO

! ################ diagnostics, electrostatic potential #################
  IF (T_cntr.EQ.Save_probes_e_data_T_cntr) THEN
     DO npb = 1, N_of_probes_block
        i = Probe_params_block_list(1,npb)
        j = Probe_params_block_list(2,npb)
        probe_F_block(npb) = phi(i,j)
     END DO
  END IF

! cleanup
  IF (ALLOCATED(a))   DEALLOCATE(a, STAT = ALLOC_ERR)
  IF (ALLOCATED(rhs)) DEALLOCATE(rhs, STAT = ALLOC_ERR)
  IF (ALLOCATED(dQ_full))   DEALLOCATE(dQ_full, STAT = ALLOC_ERR)

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)

END SUBROUTINE SOLVE_EXTERNAL_CONTOUR
