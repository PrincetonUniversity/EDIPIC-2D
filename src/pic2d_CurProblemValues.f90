!------------------------------------------
!
SUBROUTINE SET_PHYSICAL_CONSTANTS

  USE CurrentProblemValues, ONLY : true_eps_0_Fm, eps_0_Fm
  USE ParallelOperationValues, ONLY : Rank_of_process

  IMPLICIT NONE

  LOGICAL exists
  INTEGER IOS
  CHARACTER(1) buf
  REAL(8) temp

  eps_0_Fm = true_eps_0_Fm   ! default value

  INQUIRE (FILE = 'init_physconstants.dat', EXIST = exists)

  IF (.NOT.exists) THEN
     IF (Rank_of_process.EQ.0) PRINT '("###### File init_physconstants.dat not found, use default/true eps_0 = ",e15.8," F/m ######")', eps_0_Fm
     RETURN
  END IF
 
  OPEN (9, FILE = 'init_physconstants.dat')

  READ (9, '(A1)', IOSTAT=IOS) buf ! the vacuum permittivity constant is in the line below (the true value is 8.854188e-12 F/m)

  IF (IOS.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("###### Error while reading from file init_physconstants.dat, use default/true eps_0 = ",e15.8," F/m ######")', eps_0_Fm
     CLOSE (9, STATUS = 'KEEP')
     RETURN
  END IF

  READ (9, *, IOSTAT=IOS) temp

  IF (IOS.NE.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("###### Error while reading from file init_physconstants.dat, use default/true eps_0 = ",e15.8," F/m ######")', eps_0_Fm
     CLOSE (9, STATUS = 'KEEP')
     RETURN
  END IF

  CLOSE (9, STATUS = 'KEEP')

  IF (temp.LE.0.0_8) THEN
! omit negative or zero values
     IF (Rank_of_process.EQ.0) PRINT '("###### The eps_0 value acquired from file init_physconstants.dat is not positive, use default/true eps_0 = ",e15.8," F/m ######")', eps_0_Fm
     RETURN
  END IF

! since we are here use what is in the file
  eps_0_Fm = temp
  IF (Rank_of_process.EQ.0) PRINT '("###### WARNING ####### The value of vacuum permittivity is obtained from init_physconstants.dat and is eps_0 = ",e15.8," F/m ###### WARNING ######")', eps_0_Fm

END SUBROUTINE SET_PHYSICAL_CONSTANTS

!------------------------------------------
!
!
SUBROUTINE INITIATE_PARAMETERS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE IonParticles
  USE LoadBalancing
  USE BlockAndItsBoundaries, ONLY : block_row, block_column
  USE ClusterAndItsBoundaries, ONLY : c_row, c_column, c_N_of_local_object_parts, c_indx_x_min, c_indx_x_max, c_indx_y_min, c_indx_y_max
  USE Diagnostics, ONLY : N_of_probes, N_of_probes_cluster
  USE BlockAndItsBoundaries, ONLY : indx_x_min, indx_x_max, indx_y_min, indx_y_max
  USE ExternalFields
  USE Checkpoints
  USE SetupValues, ONLY : ht_grid_requested, F_grid, grid_j

  USE rng_wrapper

  USE PETSc_Solver
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL exists

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  CHARACTER(1) buf
  INTEGER ALLOC_ERR
  INTEGER i, n, m, itmp, s, j

  LOGICAL config_inconsistent

!  INTEGER n_connected_to_start
!  INTEGER n_connected_to_end
  INTEGER nn

  INTEGER metal_object_counter

  INTEGER, ALLOCATABLE :: processed_flag(:)
  INTEGER iperiodx, jperiody
  INTEGER nio, nio2
  INTEGER n_of_left_boundary, n_of_right_boundary
  INTEGER n_of_bottom_boundary, n_of_top_boundary

  INTEGER mm

  REAL(8) t_sim_ns

  LOGICAL period_x_found, period_y_found

  INTEGER boundary_master

  INTEGER init_random_seed
  REAL(8) myran

  CHARACTER(32) bo_filename     ! boundary_object_NNN_segments.dat
                                ! ----*----I----*----I----*----I--

  CHARACTER(24) iobox_filename  ! inner_object_NNN_box.dat
  CHARACTER(24) BxBy_filename   ! proc_NNNN_BxBy_vs_xy.dat
                                ! ----x----I----x----I----

  INTEGER bufsize
  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER pos1, pos2

! functions
  REAL(8) Bx, By
  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface
  INTEGER convert_logical_to_int

! default values
  given_F_double_period_sys = 1000000.0_8        !
  i_given_F_double_period_sys = -7777777   ! should be sufficient to not to trigger accidentally the node with given potential
  j_given_F_double_period_sys = -7777777   ! for a double-periodic system without given potential metal boundaries

  INQUIRE (FILE = 'init_configuration.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (.NOT.exists) THEN    
     PRINT '(2x,"Process ",i5," : ERROR : init_configuration.dat not found. Program terminated")', Rank_of_process
     CALL MPI_FINALIZE(ierr)
     STOP
  END IF

  IF (Rank_of_process.EQ.0) THEN
     PRINT '(2x,"Process ",i5," : init_configuration.dat is found. Reading the data file...")', Rank_of_process
  END IF

  OPEN (9, FILE = 'init_configuration.dat')

  READ (9, '(A1)') buf !"---dddd.ddddddd----- scale electron temperature [eV]")')
  READ (9, '(3x,f12.7)') T_e_eV
  READ (9, '(A1)') buf !"---+d.dddddddE+dd--- scale electron density [m^-3]")')
  READ (9, '(3x,e14.7)') N_plasma_m3
  READ (9, '(A1)') buf !"---ddd---------- number of cells per scale electron Debye length")')
  READ (9, '(3x,i3)') N_of_cells_debye
  READ (9, '(A1)') buf !"---ddd---------- maximal expected velocity [units of scale thermal electron velocity]")')
  READ (9, '(3x,i3)') N_max_vel
  READ (9, '(A1)') buf !"---ddd---------- number of blocks (processes) along the X (horizontal) direction")')
  READ (9, '(3x,i3)') N_blocks_x
  READ (9, '(A1)') buf !"---ddd---------- number of blocks (processes) along the Y (vertical) direction")')
  READ (9, '(3x,i3)') N_blocks_y
  READ (9, '(A1)') buf !"---ddd---------- number of cells along the X-direction in a block")')
  READ (9, '(3x,i3)') N_grid_block_x
  READ (9, '(A1)') buf !"---ddd---------- number of cells along the Y-direction in a block")')
  READ (9, '(3x,i3)') N_grid_block_y
  READ (9, '(A1)') buf !"--dddd---------- number of macroparticles per cell for the scale density")')
  READ (9, '(2x,i4)') N_of_particles_cell
  READ (9, '(A1)') buf !"-----d---------- number of blocks in a cluster along the X-direction")')
  READ (9, '(5x,i1)') cluster_N_blocks_x
  READ (9, '(A1)') buf !"-----d---------- number of blocks in a cluster along the Y-direction")')
  READ (9, '(5x,i1)') cluster_N_blocks_y
  READ (9, '(A1)') buf !"---ddd---ddd---- number of objects along domain boundary // number of material inner objects (>=0), each inner objects is a rectangle")')
  READ (9, '(3x,i3,3x,i3)') N_of_boundary_objects, N_of_inner_objects

! configuration consistency check (rectangular domain)
! N_blocks_x * N_blocks_y = N_of_processes
! (N_blocks_x / cluster_N_blocks_x) * cluster_N_blocks_x = N_blocks_x 
! (N_blocks_y / cluster_N_blocks_y) * cluster_N_blocks_y = N_blocks_y

  config_inconsistent = .FALSE.

  IF (N_of_processes.NE.(N_blocks_x*N_blocks_y)) THEN
     config_inconsistent = .TRUE.
     IF (Rank_of_process.EQ.0) PRINT '("@@@ INCONSISTENT CONFIGURATION ERROR-1, N_blocks_x * N_blocks_y .NE. N_of_processes :: ",i4," * ",i4," = ",i4," instead of ",i4," @@@")', &
          & N_blocks_x, N_blocks_y, N_blocks_x*N_blocks_y, N_of_processes
  END IF

  IF (N_blocks_x.NE.((N_blocks_x / cluster_N_blocks_x) * cluster_N_blocks_x)) THEN
     config_inconsistent = .TRUE.
     IF (Rank_of_process.EQ.0) PRINT '("@@@ INCONSISTENT CONFIGURATION ERROR-2, N_blocks_x .NE. (N_blocks_x / cluster_N_blocks_x) * cluster_N_blocks_x :: (",i4," / ",i4,") * ",i4," = ",i4," instead of ",i4," @@@")', &
          & N_blocks_x, cluster_N_blocks_x, cluster_N_blocks_x, (N_blocks_x / cluster_N_blocks_x) * cluster_N_blocks_x, N_blocks_x
  END IF

  IF (N_blocks_y.NE.((N_blocks_y / cluster_N_blocks_y) * cluster_N_blocks_y)) THEN
     config_inconsistent = .TRUE.
     IF (Rank_of_process.EQ.0) PRINT '("@@@ INCONSISTENT CONFIGURATION ERROR-3, N_blocks_y .NE. (N_blocks_y / cluster_N_blocks_y) * cluster_N_blocks_y :: (",i4," / ",i4,") * ",i4," = ",i4," instead of ",i4," @@@")', &
          & N_blocks_y, cluster_N_blocks_y, cluster_N_blocks_y, (N_blocks_y / cluster_N_blocks_y) * cluster_N_blocks_y, N_blocks_y
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (config_inconsistent) THEN
     CALL MPI_FINALIZE(ierr)
     STOP
  END IF
     
  IF (Rank_of_process.EQ.0) THEN
     IF (N_of_inner_objects.EQ.0) THEN
        PRINT '(2x,"### No inner material objects requested ###")'
     ELSE
        PRINT '(2x,"### ",i3," inner material objects requested ###")', N_of_inner_objects
     END IF
  END IF

  N_of_boundary_and_inner_objects = N_of_boundary_objects + N_of_inner_objects

  ALLOCATE(whole_object(1:N_of_boundary_and_inner_objects), STAT=ALLOC_ERR)

  ALLOCATE(ion_colls_with_bo(1:N_of_boundary_and_inner_objects), STAT=ALLOC_ERR)
  ALLOCATE(  e_colls_with_bo(1:N_of_boundary_and_inner_objects), STAT=ALLOC_ERR)

  DO n = 1, N_of_boundary_objects
     READ (9, '(A1)') buf !"===dd===dd=== object type, number of segments")')
     READ (9, '(3x,i2,3x,i2)') whole_object(n)%object_type, &
                             & whole_object(n)%number_of_segments
     READ (9, '(A1)') buf !"---dddddd---dddddd---dddddd---dddddd---segment start X/Y end X/Y [global node index]")')
     ALLOCATE(whole_object(n)%segment(1:whole_object(n)%number_of_segments), STAT=ALLOC_ERR)
     DO m = 1, whole_object(n)%number_of_segments
        READ (9, '(3x,i6,3x,i6,3x,i6,3x,i6)') whole_object(n)%segment(m)%istart, &
                                            & whole_object(n)%segment(m)%jstart, &
                                            & whole_object(n)%segment(m)%iend, &
                                            & whole_object(n)%segment(m)%jend
! sort endpoints so that start-end always goes either left-righ or bottom-top
! also check that the segment is either vertical or horizontal
        IF (whole_object(n)%segment(m)%istart.EQ.whole_object(n)%segment(m)%iend) THEN
! vertical segment
           itmp = MAX(whole_object(n)%segment(m)%jstart, whole_object(n)%segment(m)%jend)
           whole_object(n)%segment(m)%jstart = MIN(whole_object(n)%segment(m)%jstart, whole_object(n)%segment(m)%jend)
           whole_object(n)%segment(m)%jend   = itmp
        ELSE IF (whole_object(n)%segment(m)%jstart.EQ.whole_object(n)%segment(m)%jend) THEN
! horizontal segment
           itmp = MAX(whole_object(n)%segment(m)%istart, whole_object(n)%segment(m)%iend)
           whole_object(n)%segment(m)%istart = MIN(whole_object(n)%segment(m)%istart, whole_object(n)%segment(m)%iend)
           whole_object(n)%segment(m)%iend   = itmp
        ELSE
! error
           PRINT '("Process ",i4," ERROR in INITIATE_PARAMETERS : boundary object ",i4," has segment ",i4," with inconsistent start i/j ",2(2x,i6)," and end i/j ",2(2x,i6))', &
                & Rank_of_process, n, m, whole_object(n)%segment(m)%istart, whole_object(n)%segment(m)%jstart, whole_object(n)%segment(m)%iend, whole_object(n)%segment(m)%jend
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        END IF
     END DO
  END DO

  DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects  !N_of_inner_objects
     READ (9, '(A1)') buf !"===dd=== object type")')
     READ (9, '(3x,i2)') whole_object(n)%object_type
     IF ((whole_object(n)%object_type.NE.METAL_WALL).AND.(whole_object(n)%object_type.NE.DIELECTRIC)) THEN
        IF (Rank_of_process.EQ.0) PRINT '("Error, inner material object ",i2," has type ",i3," which is not permitted")', n, whole_object(n)%object_type
        CALL MPI_FINALIZE(ierr)
        STOP
     END IF
     READ (9, '(A1)') buf !"---dddddd---dddddd---dddddd---dddddd--- coordinates of left bottom X/Y corner and right top X/Y corners [global node index]")')
     READ (9, '(3x,i6,3x,i6,3x,i6,3x,i6)') whole_object(n)%ileft, whole_object(n)%jbottom, whole_object(n)%iright, whole_object(n)%jtop

     whole_object(n)%number_of_segments = 4
     ALLOCATE(whole_object(n)%segment(1:whole_object(n)%number_of_segments), STAT=ALLOC_ERR)
! left side
     whole_object(n)%segment(1)%istart = whole_object(n)%ileft
     whole_object(n)%segment(1)%jstart = whole_object(n)%jbottom
     whole_object(n)%segment(1)%iend = whole_object(n)%ileft
     whole_object(n)%segment(1)%jend = whole_object(n)%jtop
     ALLOCATE(whole_object(n)%segment(1)%cell_is_covered(whole_object(n)%jbottom:(whole_object(n)%jtop-1)), STAT=ALLOC_ERR)
     whole_object(n)%segment(1)%cell_is_covered = .FALSE.
! top side
     whole_object(n)%segment(2)%istart = whole_object(n)%ileft
     whole_object(n)%segment(2)%jstart = whole_object(n)%jtop
     whole_object(n)%segment(2)%iend = whole_object(n)%iright
     whole_object(n)%segment(2)%jend = whole_object(n)%jtop
     ALLOCATE(whole_object(n)%segment(2)%cell_is_covered(whole_object(n)%ileft:(whole_object(n)%iright-1)), STAT=ALLOC_ERR)
     whole_object(n)%segment(2)%cell_is_covered = .FALSE.
! right side
     whole_object(n)%segment(3)%istart = whole_object(n)%iright
     whole_object(n)%segment(3)%jstart = whole_object(n)%jbottom
     whole_object(n)%segment(3)%iend = whole_object(n)%iright
     whole_object(n)%segment(3)%jend = whole_object(n)%jtop
     ALLOCATE(whole_object(n)%segment(3)%cell_is_covered(whole_object(n)%jbottom:(whole_object(n)%jtop-1)), STAT=ALLOC_ERR)
     whole_object(n)%segment(3)%cell_is_covered = .FALSE.
! bottom side
     whole_object(n)%segment(4)%istart = whole_object(n)%ileft
     whole_object(n)%segment(4)%jstart = whole_object(n)%jbottom
     whole_object(n)%segment(4)%iend = whole_object(n)%iright
     whole_object(n)%segment(4)%jend = whole_object(n)%jbottom
     ALLOCATE(whole_object(n)%segment(4)%cell_is_covered(whole_object(n)%ileft:(whole_object(n)%iright-1)), STAT=ALLOC_ERR)
     whole_object(n)%segment(4)%cell_is_covered = .FALSE.

     whole_object(n)%Xmin = DBLE(whole_object(n)%ileft)
     whole_object(n)%Xmax = DBLE(whole_object(n)%iright)
     whole_object(n)%Ymin = DBLE(whole_object(n)%jbottom)
     whole_object(n)%Ymax = DBLE(whole_object(n)%jtop)

     whole_object(n)%N_boundary_nodes = ((whole_object(n)%jtop - whole_object(n)%jbottom + 1) + (whole_object(n)%iright - whole_object(n)%ileft - 1))*2
!        whole_object(n)%L = whole_object(n)%N_boundary_nodes

     IF (whole_object(n)%object_type.EQ.DIELECTRIC) THEN
        ALLOCATE(whole_object(n)%surface_charge_variation(1:whole_object(n)%N_boundary_nodes), STAT = ALLOC_ERR)
        whole_object(n)%surface_charge_variation = 0.0_8
     END IF

     whole_object(n)%n_connected_to_start = -1
     whole_object(n)%n_connected_to_end = -1

 !! ### todo :: add consistency check    
  END DO   !###  DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! we close this file after we figure out what is the periodicity of the system
!     CLOSE (9, STATUS = 'KEEP')

  DO n = 1, N_of_boundary_and_inner_objects
     whole_object(n)%object_id_number = n
     whole_object(n)%object_copy_periodic_X_right = -1  ! default values
     whole_object(n)%object_copy_periodic_Y_above = -1
     whole_object(n)%object_does_NOT_cross_symmetry_plane_X = .TRUE.
  END DO
        
! configure vacuum gaps if they are included
  DO n = 1, N_of_boundary_objects

     whole_object(n)%n_connected_to_start = -1
     whole_object(n)%n_connected_to_end = -1
 
     IF (whole_object(n)%object_type.NE.VACUUM_GAP) CYCLE

! note that a vacuum gap has only one segment
     IF (whole_object(n)%segment(1)%jstart.NE.whole_object(n)%segment(1)%jend) THEN
! vertical gap (along y)
        ALLOCATE(whole_object(n)%phi_profile(whole_object(n)%segment(1)%jstart:whole_object(n)%segment(1)%jend), STAT=ALLOC_ERR)
     ELSE IF (whole_object(n)%segment(1)%istart.NE.whole_object(n)%segment(1)%iend) THEN
! horizontal gap (along x)
        ALLOCATE(whole_object(n)%phi_profile(whole_object(n)%segment(1)%istart:whole_object(n)%segment(1)%iend), STAT=ALLOC_ERR)
     END IF

! check whether the start point of the gap coinsides with a start or an end of any segment of any other object
     DO nn = 1, N_of_boundary_objects
        IF (nn.EQ.n) CYCLE
        DO m = 1, whole_object(nn)%number_of_segments
           IF ( ( (whole_object(nn)%segment(m)%istart.EQ.whole_object(n)%segment(1)%istart).AND. &
              &   (whole_object(nn)%segment(m)%jstart.EQ.whole_object(n)%segment(1)%jstart) ) &
              & .OR. &
              & ( (whole_object(nn)%segment(m)%iend.EQ.whole_object(n)%segment(1)%istart).AND. &
              &   (whole_object(nn)%segment(m)%jend.EQ.whole_object(n)%segment(1)%jstart) ) ) THEN 
              whole_object(n)%n_connected_to_start = nn
              EXIT
           END IF
        END DO
        IF (whole_object(n)%n_connected_to_start.GT.0) EXIT
     END DO

! check whether the end point of the gap coinsides with a start or an end of any segment of any other object
     DO nn = 1, N_of_boundary_objects
        IF (nn.EQ.n) CYCLE
        DO m = 1, whole_object(nn)%number_of_segments
           IF ( ( (whole_object(nn)%segment(m)%istart.EQ.whole_object(n)%segment(1)%iend).AND. &
              &   (whole_object(nn)%segment(m)%jstart.EQ.whole_object(n)%segment(1)%jend) ) &
              & .OR. &
              & ( (whole_object(nn)%segment(m)%iend.EQ.whole_object(n)%segment(1)%iend).AND. &
              &   (whole_object(nn)%segment(m)%jend.EQ.whole_object(n)%segment(1)%jend) ) ) THEN 
              whole_object(n)%n_connected_to_end = nn
              EXIT
           END IF
        END DO
        IF (whole_object(n)%n_connected_to_end.GT.0) EXIT
     END DO

     IF (whole_object(n)%n_connected_to_start.EQ.-1) THEN
        IF (Rank_of_process.EQ.0) PRINT '("Error, start point of vacuum gap boundary object ",i2," is not connected to any other boundary object")', n
        CALL MPI_FINALIZE(ierr)
        STOP
     END IF

     IF (whole_object(n)%n_connected_to_end.EQ.-1) THEN
        IF (Rank_of_process.EQ.0) PRINT '("Error, end point of vacuum gap boundary object ",i2," is not connected to any other boundary object")', n
        CALL MPI_FINALIZE(ierr)
        STOP
     END IF

     IF (Rank_of_process.EQ.0) PRINT '("Vacuum gap object ",i2," is connected to objects ",i2," (start) and ",i2," (end)")', n, whole_object(n)%n_connected_to_start, whole_object(n)%n_connected_to_end

  END DO   !###   DO n = 1, N_of_boundary_objects

! identify whether periodicity was requested
  periodicity_flag = PERIODICITY_NONE
  period_x_found = .FALSE.
  period_y_found = .FALSE.
  DO n = 1, N_of_boundary_objects
     IF (whole_object(n)%object_type.EQ.PERIODIC_PIPELINE_X) period_x_found = .TRUE.
     IF (whole_object(n)%object_type.EQ.-PERIODIC_PIPELINE_X) THEN   !### if the flag is negative, use PETSc based solver
        period_x_found = .TRUE.
        periodicity_flag = PERIODICITY_X_PETSC
        whole_object(n)%object_type = PERIODIC_PIPELINE_X            !### don't forget to fix the value
     END IF
     IF (whole_object(n)%object_type.EQ.PERIODIC_PIPELINE_Y) period_y_found = .TRUE.
  END DO
  IF (period_x_found) THEN
     IF (periodicity_flag.EQ.PERIODICITY_NONE) periodicity_flag = PERIODICITY_X   ! if flag was set to PERIODICITY_X_PETSC above, this is skipped
  END IF
  IF (period_x_found.AND.period_y_found) periodicity_flag = PERIODICITY_X_Y       ! if flag was set to PERIODICITY_X_PETSC above, it is redefined here
  IF (period_y_found.AND.(.NOT.period_x_found)) THEN
! error, case with periodicity in y, no periodicity in x is not included
     PRINT '(2x,"Process ",i3," : ERROR : periodicity in Y requested without periodicity in X")', Rank_of_process
     CALL MPI_FINALIZE(ierr)
     STOP
  END IF

!! identify whether periodicity was requested
!  periodicity_flag = PERIODICITY_NONE
!  period_x_found = .FALSE.
!  period_y_found = .FALSE.
!  DO n = 1, N_of_boundary_objects
!     IF (whole_object(n)%object_type.EQ.PERIODIC_PIPELINE_X) period_x_found = .TRUE.
!     IF (whole_object(n)%object_type.EQ.PERIODIC_PIPELINE_Y) period_y_found = .TRUE.
!  END DO
!  IF (period_x_found)                    periodicity_flag = PERIODICITY_X
!  IF (period_x_found.AND.period_y_found) periodicity_flag = PERIODICITY_X_Y
!  IF (period_y_found.AND.(.NOT.period_x_found)) THEN
!! error, case with periodicity in y, no periodicity in x is not included
!     PRINT '(2x,"Process ",i3," : ERROR : periodicity in Y requested without periodicity in X")', Rank_of_process
!     CALL MPI_FINALIZE(ierr)
!     STOP
!  END IF

  IF (periodicity_flag.EQ.PERIODICITY_X_Y) THEN
! in a system periodic in both X and Y directions
! if there is no objects with given pontential inside
! expect that the user specifies a node (i,j) and a potential in that node 
     metal_object_counter = 0
     DO n = 1, N_of_boundary_and_inner_objects
        IF (whole_object(n)%object_type.EQ.METAL_WALL) metal_object_counter = metal_object_counter + 1
     END DO
     IF (metal_object_counter.EQ.0) THEN
        READ (9, '(A1)') buf !"---dddd---dddd--- for a system periodic along both X and Y directions, without objects with given potential, specify X/Y (global node index i/j) of a point with given potential ")')
        READ (9, '(3x,i4,3x,i4)') i_given_F_double_period_sys, j_given_F_double_period_sys
        READ (9, '(A1)') buf !"---ddddd.d------- value of the potential at this point [V]
        READ (9, '(3x,f7.1)')  given_F_double_period_sys
        IF (Rank_of_process.EQ.0) PRINT '("### system periodic alog X and Y directions, no objects with given potential found, potential of point i/j ",i4,"/",i4," will be set to ",f7.1," V")', &
             & i_given_F_double_period_sys, j_given_F_double_period_sys, given_F_double_period_sys
     END IF
  END IF

  CLOSE (9, STATUS = 'KEEP')   !##### closed file init_configuration.dat

  IF (N_of_inner_objects.GT.0) THEN
     IF (periodicity_flag.EQ.PERIODICITY_X) THEN 
        periodicity_flag = PERIODICITY_X_PETSC
        IF (Rank_of_process.EQ.0) PRINT '("### System periodic along the X-direction with ",i3," inner objects, FFT-based field solver is off, PETSc-based solver is used instead")', N_of_inner_objects
     END IF
  END IF

!if ((periodicity_flag.EQ.PERIODICITY_X).and.(N_of_boundary_and_inner_objects.GT.4)) then
!   periodicity_flag = PERIODICITY_X_PETSC
!   IF (Rank_of_process.EQ.0) PRINT '("### System periodic along the X-direction with ",i3," boundary objects, FFT-based field solver is off, PETSc-based solver is used instead")', N_of_boundary_and_inner_objects
!end if

! for periodic systems, check whether there are inner objects crossing periodic boundaries =============================================================================================================================
! if there are, identify matching objects (copies) at opposite periodic boundaries
! the code below should work if there are several segments of periodicity along the X-border
! this may be used when the whole simulation domain is not a rectangle

  ALLOCATE(processed_flag(1:N_of_boundary_objects), STAT=ALLOC_ERR)

! X-boundary
  processed_flag = 0
  DO n = 1, N_of_boundary_objects
     IF (whole_object(n)%object_type.NE.PERIODIC_PIPELINE_X) CYCLE
! boundary object n is a X-periodic boundary
! skip already processed boundary object
     IF (processed_flag(n).EQ.1) CYCLE 
     processed_flag(n) = 1
! find the opposite X-perioidic boundary
     iperiodx = -1
     DO nn = 1, N_of_boundary_objects
        IF (nn.EQ.n) CYCLE
        IF (whole_object(nn)%object_type.NE.PERIODIC_PIPELINE_X) CYCLE
        IF (processed_flag(nn).EQ.1) CYCLE 
        IF (MIN(whole_object(n)%segment(1)%jstart, whole_object(n)%segment(1)%jend).NE.MIN(whole_object(nn)%segment(1)%jstart, whole_object(nn)%segment(1)%jend)) CYCLE
        IF (MAX(whole_object(n)%segment(1)%jstart, whole_object(n)%segment(1)%jend).NE.MAX(whole_object(nn)%segment(1)%jstart, whole_object(nn)%segment(1)%jend)) CYCLE
! ends of boundary objects have same Y-coordinates, therefore these two are the matching X-periodic boundaries
        processed_flag(nn) = 1
! find which one is left which one is right
        IF (whole_object(n)%segment(1)%istart.LT.whole_object(nn)%segment(1)%istart) THEN
           n_of_left_boundary = n
           n_of_right_boundary = nn
        ELSE ! IF (whole_object(nn)%segment(1)%istart.LT.whole_object(n)%segment(1)%istart) THEN
           n_of_left_boundary = nn
           n_of_right_boundary = n
        END IF
        iperiodx = whole_object(n_of_right_boundary)%segment(1)%istart - whole_object(n_of_left_boundary)%segment(1)%istart - 1  ! exclude overlapping
        EXIT
     END DO   !### DO nn = 1, N_of_boundary_objects
     IF (iperiodx.LT.0) THEN
! error
        PRINT '("error, X-periodic boundary ",i3," has no match")', n
        CALL MPI_FINALIZE(ierr)
        STOP
     END IF
! scan through inner objects, find one which crosses the left boundary
     DO nio = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
        IF (whole_object(nio)%ileft.GT.whole_object(n_of_left_boundary)%segment(1)%istart) CYCLE
        IF (whole_object(nio)%iright.LT.whole_object(n_of_left_boundary)%segment(1)%istart) CYCLE
        IF (whole_object(nio)%jbottom.GE.MAX(whole_object(n_of_left_boundary)%segment(1)%jstart, whole_object(n_of_left_boundary)%segment(1)%jend)) CYCLE
        IF (whole_object(nio)%jtop.LE.MIN(whole_object(n_of_left_boundary)%segment(1)%jstart, whole_object(n_of_left_boundary)%segment(1)%jend)) CYCLE
! the object crosses left boundary, let's find its matching copy
        DO nio2 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (whole_object(nio2)%ileft.NE.(whole_object(nio)%ileft+iperiodx)) CYCLE
           IF (whole_object(nio2)%iright.NE.(whole_object(nio)%iright+iperiodx)) CYCLE
           IF (whole_object(nio2)%jbottom.NE.whole_object(nio)%jbottom) CYCLE
           IF (whole_object(nio2)%jtop.NE.whole_object(nio)%jtop) CYCLE
! this is the copy
           whole_object(nio)%object_copy_periodic_X_right = nio2
           EXIT
        END DO
! fool proof
        IF (whole_object(nio)%object_copy_periodic_X_right.LT.0) THEN
! error
           PRINT '("error, matching object for inner object ",i3," which crosses X-periodic boundary ",i3," not found")', nio, n_of_left_boundary
           CALL MPI_FINALIZE(ierr)
           STOP
        END IF
     END DO   !### DO nio = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
  END DO   !### DO n = 1, N_of_boundary_objects

! Y-boundary
  processed_flag = 0
  DO n = 1, N_of_boundary_objects
     IF (whole_object(n)%object_type.NE.PERIODIC_PIPELINE_Y) CYCLE
! boundary object n is a Y-periodic boundary
! skip already processed boundary object
     IF (processed_flag(n).EQ.1) CYCLE 
     processed_flag(n) = 1
! find the opposite Y-perioidic boundary
     jperiody = -1
     DO nn = 1, N_of_boundary_objects
        IF (nn.EQ.n) CYCLE
        IF (whole_object(nn)%object_type.NE.PERIODIC_PIPELINE_Y) CYCLE
        IF (processed_flag(nn).EQ.1) CYCLE 
        IF (MIN(whole_object(n)%segment(1)%istart, whole_object(n)%segment(1)%iend).NE.MIN(whole_object(nn)%segment(1)%istart, whole_object(nn)%segment(1)%iend)) CYCLE
        IF (MAX(whole_object(n)%segment(1)%istart, whole_object(n)%segment(1)%iend).NE.MAX(whole_object(nn)%segment(1)%istart, whole_object(nn)%segment(1)%iend)) CYCLE
! ends of boundary objects have same X-coordinates, therefore these two are the matching Y-periodic boundaries
        processed_flag(nn) = 1
! find which one is bottom which one is top
        IF (whole_object(n)%segment(1)%jstart.LT.whole_object(nn)%segment(1)%jstart) THEN
           n_of_bottom_boundary = n
           n_of_top_boundary = nn
        ELSE ! IF (whole_object(nn)%segment(1)%jstart.LT.whole_object(n)%segment(1)%jstart) THEN
           n_of_bottom_boundary = nn
           n_of_top_boundary = n
        END IF
        jperiody = whole_object(n_of_top_boundary)%segment(1)%jstart - whole_object(n_of_bottom_boundary)%segment(1)%jstart - 1  ! exclude overlapping
        EXIT
     END DO   !### DO nn = 1, N_of_boundary_objects
     IF (jperiody.LT.0) THEN
! error
        PRINT '("error, Y-periodic boundary ",i3," has no match")', n
        CALL MPI_FINALIZE(ierr)
        STOP
     END IF
! scan through inner objects, find one which crosses the bottom boundary
     DO nio = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
        IF (whole_object(nio)%jbottom.GT.whole_object(n_of_bottom_boundary)%segment(1)%jstart) CYCLE
        IF (whole_object(nio)%jtop.LT.whole_object(n_of_bottom_boundary)%segment(1)%jstart) CYCLE
        IF (whole_object(nio)%ileft.GE.MAX(whole_object(n_of_bottom_boundary)%segment(1)%istart, whole_object(n_of_bottom_boundary)%segment(1)%iend)) CYCLE
        IF (whole_object(nio)%iright.LE.MIN(whole_object(n_of_bottom_boundary)%segment(1)%istart, whole_object(n_of_bottom_boundary)%segment(1)%iend)) CYCLE
! the object crosses bottom boundary, let's find its matching copy
        DO nio2 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (whole_object(nio2)%ileft.NE.whole_object(nio)%ileft) CYCLE
           IF (whole_object(nio2)%iright.NE.whole_object(nio)%iright) CYCLE
           IF (whole_object(nio2)%jbottom.NE.(whole_object(nio)%jbottom+jperiody)) CYCLE
           IF (whole_object(nio2)%jtop.NE.(whole_object(nio)%jtop+jperiody)) CYCLE
! this is the copy
           whole_object(nio)%object_copy_periodic_Y_above = nio2
           EXIT
        END DO
! fool proof
        IF (whole_object(nio)%object_copy_periodic_Y_above.LT.0) THEN
! error
           PRINT '("error, matching object for inner object ",i3," which crosses Y-periodic boundary ",i3," not found")', nio, n_of_bottom_boundary
           CALL MPI_FINALIZE(ierr)
           STOP
        END IF
     END DO   !### DO nio = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
  END DO   !### DO n = 1, N_of_boundary_objects

  IF (ALLOCATED(processed_flag)) DEALLOCATE( processed_flag, STAT=ALLOC_ERR)

! identify covered parts of segments of inner boundary objects ========================================================================================================================================================

  DO nio = N_of_boundary_objects + 1, N_of_boundary_and_inner_objects

! left and right vertical segments (#1 and #3)
     DO m = 1, 3, 2
        DO j = whole_object(nio)%segment(m)%jstart, whole_object(nio)%segment(m)%jend-1

! first, check if there are any ordinary boundary objects attached
           DO n = 1, N_of_boundary_objects
              IF (whole_object(n)%object_type.EQ.PERIODIC_PIPELINE_X) CYCLE
              IF (whole_object(n)%object_type.EQ.PERIODIC_PIPELINE_Y) CYCLE
              IF (whole_object(n)%object_type.EQ.SYMMETRY_PLANE) CYCLE
              DO mm = 1, whole_object(n)%number_of_segments
                 IF (whole_object(n)%segment(mm)%istart.NE.whole_object(n)%segment(mm)%iend) CYCLE    ! skip horizontal segments
                 IF (whole_object(n)%segment(mm)%istart.NE.whole_object(nio)%segment(m)%istart) CYCLE ! skip vertical segment with different i
                 IF ((j.GE.whole_object(n)%segment(mm)%jstart).AND.(j.LT.whole_object(n)%segment(mm)%jend)) THEN
                    whole_object(nio)%segment(m)%cell_is_covered(j) = .TRUE.
                    EXIT
                 END IF
              END DO   !### DO mm = 1, whole_object(n)%number_of_segments
              IF (whole_object(nio)%segment(m)%cell_is_covered(j)) EXIT
           END DO   !###  DO n = 1, N_of_boundary_objects

           IF (whole_object(nio)%segment(m)%cell_is_covered(j)) CYCLE

! second, check if there is attachment or overlapping with another inner object
           DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
              IF (n.EQ.nio) CYCLE
              IF (whole_object(nio)%segment(m)%istart.LT.whole_object(n)%ileft) CYCLE
              IF (whole_object(nio)%segment(m)%istart.GT.whole_object(n)%iright) CYCLE
              IF ((j.GE.whole_object(n)%jbottom).AND.(j.LT.whole_object(n)%jtop)) THEN
                 whole_object(nio)%segment(m)%cell_is_covered(j) = .TRUE.
                 EXIT
              END IF
           END DO   !### DO n = 1, N_of_boundary_objects+1, N_of_boundary_and_inner_objects

        END DO   !### DO j = whole_object(nio)%segment(m)%jstart, whole_object(nio)%segment(m)%jend-1
     END DO   !### DO m = 1, 3, 2

! top and bottom horizontal segments (#2 and #4)
     DO m = 2, 4, 2
        DO i = whole_object(nio)%segment(m)%istart, whole_object(nio)%segment(m)%iend-1

! first, check if there are any ordinary boundary objects attached
           DO n = 1, N_of_boundary_objects
              IF (whole_object(n)%object_type.EQ.PERIODIC_PIPELINE_X) CYCLE
              IF (whole_object(n)%object_type.EQ.PERIODIC_PIPELINE_Y) CYCLE
              IF (whole_object(n)%object_type.EQ.SYMMETRY_PLANE) CYCLE
              DO mm = 1, whole_object(n)%number_of_segments
                 IF (whole_object(n)%segment(mm)%jstart.NE.whole_object(n)%segment(mm)%jend) CYCLE    ! skip vertical segments
                 IF (whole_object(n)%segment(mm)%jstart.NE.whole_object(nio)%segment(m)%jstart) CYCLE ! skip horizontal segment with different j
                 IF ((i.GE.whole_object(n)%segment(mm)%istart).AND.(i.LT.whole_object(n)%segment(mm)%iend)) THEN
                    whole_object(nio)%segment(m)%cell_is_covered(i) = .TRUE.
                    EXIT
                 END IF
              END DO   !### DO mm = 1, whole_object(n)%number_of_segments
              IF (whole_object(nio)%segment(m)%cell_is_covered(i)) EXIT
           END DO   !###  DO n = 1, N_of_boundary_objects

           IF (whole_object(nio)%segment(m)%cell_is_covered(i)) CYCLE

! second, check if there is attachment or overlapping with another inner object
           DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
              IF (n.EQ.nio) CYCLE
              IF (whole_object(nio)%segment(m)%jstart.LT.whole_object(n)%jbottom) CYCLE
              IF (whole_object(nio)%segment(m)%jstart.GT.whole_object(n)%jtop) CYCLE
              IF ((i.GE.whole_object(n)%ileft).AND.(i.LT.whole_object(n)%iright)) THEN
                 whole_object(nio)%segment(m)%cell_is_covered(i) = .TRUE.
                 EXIT
              END IF
           END DO   !### DO n = 1, N_of_boundary_objects+1, N_of_boundary_and_inner_objects

        END DO   !### DO j = whole_object(nio)%segment(m)%istart, whole_object(nio)%segment(m)%iend-1
     END DO   !### DO m = 2, 4, 2

  END DO   !### DO nio = N_of_boundary_objects + 1, N_of_boundary_and_inner_objects

! find whether any inner object is positioned across symmetry plane at x=0 (if there is one, of course) ================================================================================================================
  DO n = 1, N_of_boundary_objects
     IF (whole_object(n)%object_type.EQ.SYMMETRY_PLANE) THEN
! a symmetry plane x=0 is found
        DO nio = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF ((whole_object(nio)%ileft.LT.0).AND.(whole_object(nio)%iright.GT.0)) THEN
              IF (whole_object(nio)%ileft.NE.-whole_object(nio)%iright) THEN
                 PRINT '("Error, inner object ",i3," is placed across the symmetry plane but is not symmetric")', nio
                 CALL MPI_FINALIZE(ierr)
                 STOP
              END IF
              whole_object(nio)%object_does_NOT_cross_symmetry_plane_X = .FALSE.
           END IF
        END DO
        EXIT
     END IF
  END DO   !### DO n = 1, N_of_boundary_objects

! calculate length of not covered part of the inner object, to be used later when calculating number of particles to inject  ===========================================================================================

  DO nio = N_of_boundary_objects + 1, N_of_boundary_and_inner_objects
     whole_object(nio)%L = 0
! vertical segments # 1 and 3
     DO m = 1, 3, 2
        DO j = whole_object(nio)%segment(m)%jstart, whole_object(nio)%segment(m)%jend-1
           IF (whole_object(nio)%segment(m)%cell_is_covered(j)) CYCLE
           whole_object(nio)%L = whole_object(nio)%L + 1
        END DO
     END DO
! horizontal segments # 2 and 4
     DO m = 2, 4, 2
        DO i = whole_object(nio)%segment(m)%istart, whole_object(nio)%segment(m)%iend-1
           IF (whole_object(nio)%segment(m)%cell_is_covered(i)) CYCLE
           whole_object(nio)%L = whole_object(nio)%L + 1
        END DO
     END DO
  END DO

  v_Te_ms     = SQRT(2.0_8 * T_e_eV * e_Cl / m_e_kg)
  W_plasma_s1 = SQRT(N_plasma_m3 * e_Cl**2 / (eps_0_Fm * m_e_kg))
  L_debye_m   = v_Te_ms / W_plasma_s1
  delta_x_m   = L_debye_m / N_of_cells_debye    
  delta_t_s   = delta_x_m / (N_max_vel * v_Te_ms)

! save geometry of inner objects (note that we know delta_x_m now :)

  IF (Rank_of_process.EQ.0) THEN
     DO nn = 1, N_of_boundary_objects
        bo_filename = 'boundary_object_NNN_segments.dat'
!                      ----*----I----*----I----*----I--
        bo_filename(17:19) = convert_int_to_txt_string(nn, 3)

        open (19, file = bo_filename)
        write (19, '("# segments of boundary object ",i3)') nn
        write (19, '("# column 1 is x-index of the grid node in the segment end [dim-less]")')
        write (19, '("# column 2 is y-index of the grid node in the segment end [dim-less]")')
        write (19, '("# column 3 is x-coordinate of the segment end [m]")')
        write (19, '("# column 4 is y-coordinate of the segment end [m]")')

        DO m = 1, whole_object(nn)%number_of_segments
           write (19, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') whole_object(nn)%segment(m)%istart,             whole_object(nn)%segment(m)%jstart, &
                                                       & whole_object(nn)%segment(m)%istart * delta_x_m, whole_object(nn)%segment(m)%jstart * delta_x_m
           write (19, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') whole_object(nn)%segment(m)%iend,             whole_object(nn)%segment(m)%jend, &
                                                       & whole_object(nn)%segment(m)%iend * delta_x_m, whole_object(nn)%segment(m)%jend * delta_x_m
           write (19, '(" ")')
        END DO
        close (19, status = 'keep')
        print '("### created boundary object segments file ",A32)', bo_filename
     END DO

     DO nio = N_of_boundary_objects + 1, N_of_boundary_and_inner_objects
        iobox_filename = 'inner_object_NNN_box.dat'
!                         ----*----I----*----I----
        iobox_filename(14:16) = convert_int_to_txt_string(nio, 3)

        open (19, file = iobox_filename)
        write (19, '("# spatial box of inner boundary object ",i3)') nio
        write (19, '("# column 1 is x-index of the grid node in the box corner [dim-less]")')
        write (19, '("# column 2 is y-index of the grid node in the box corner [dim-less]")')
        write (19, '("# column 3 is x-coordinate of the box corner [m]")')
        write (19, '("# column 4 is y-coordinate of the box corner [m]")')
! left-right along X, up-down along Y
        write (19, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') whole_object(nio)%ileft,  whole_object(nio)%jbottom,  whole_object(nio)%ileft * delta_x_m, whole_object(nio)%jbottom * delta_x_m  ! left bottom
        write (19, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') whole_object(nio)%ileft,  whole_object(nio)%jtop,     whole_object(nio)%ileft * delta_x_m, whole_object(nio)%jtop    * delta_x_m  ! left top
        write (19, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') whole_object(nio)%iright, whole_object(nio)%jtop,    whole_object(nio)%iright * delta_x_m, whole_object(nio)%jtop    * delta_x_m  ! right top
        write (19, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') whole_object(nio)%iright, whole_object(nio)%jbottom, whole_object(nio)%iright * delta_x_m, whole_object(nio)%jbottom * delta_x_m  ! right bottom
        write (19, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') whole_object(nio)%ileft,  whole_object(nio)%jbottom, whole_object(nio)%ileft  * delta_x_m, whole_object(nio)%jbottom * delta_x_m  ! left bottom
        close (19, status = 'keep')
        print '("### created inner object box file ",A24)', iobox_filename
        
        iobox_filename = 'inner_object_NNN_map.dat'
!                         ----*----I----*----I----
        iobox_filename(14:16) = convert_int_to_txt_string(nio, 3)
        open (19, file = iobox_filename)
        write (19, '("# map of covered/open surface cells of inner boundary object ",i3)') nio

        write (19, '("# column 1 is x-index of the middle of the cell, dimensionless, may be half-integer")')
        write (19, '("# column 2 is y-index of the middle of the cell, dimensionless, may be half-integer")')
        write (19, '("# column 3 is x-coordinate of the middle of the cell [m]")')
        write (19, '("# column 4 is y-coordinate of the middle of the cell [m]")')
        write (19, '("# column 5 is integer flag showing whether the cell is open (0) or covered (1)")')
        write (19, '("# column 6 is the coorresponding value of cell_is_covered array (logical)")')
! left side, segment 1
        i = whole_object(nio)%ileft
        DO j = whole_object(nio)%jbottom, whole_object(nio)%jtop-1
           WRITE (19, '(2x,f6.1,2x,f6.1,2x,f12.9,2x,f12.9,2x,i1,2x,L1)') REAL(i), REAL(j)+0.5, i * delta_x_m, (REAL(j)+0.5) * delta_x_m, &
                & convert_logical_to_int(whole_object(nio)%segment(1)%cell_is_covered(j)), &
                & whole_object(nio)%segment(1)%cell_is_covered(j)
        END DO
! top side, segment 2
        j = whole_object(nio)%jtop
        DO i = whole_object(nio)%ileft, whole_object(nio)%iright-1
           WRITE (19, '(2x,f6.1,2x,f6.1,2x,f12.9,2x,f12.9,2x,i1,2x,L1)') REAL(i)+0.5, REAL(j), (REAL(i)+0.5) * delta_x_m, j * delta_x_m, &
                & convert_logical_to_int(whole_object(nio)%segment(2)%cell_is_covered(i)), &
                & whole_object(nio)%segment(2)%cell_is_covered(i)
        END DO
! right side, segment 3
        i = whole_object(nio)%iright
        DO j = whole_object(nio)%jtop-1, whole_object(nio)%jbottom, -1
           WRITE (19, '(2x,f6.1,2x,f6.1,2x,f12.9,2x,f12.9,2x,i1,2x,L1)') REAL(i), REAL(j)+0.5, i * delta_x_m, (REAL(j)+0.5) * delta_x_m, &
                & convert_logical_to_int(whole_object(nio)%segment(3)%cell_is_covered(j)), &
                & whole_object(nio)%segment(3)%cell_is_covered(j)
        END DO
! bottom side, segment 4
        j = whole_object(nio)%jbottom
        DO i = whole_object(nio)%iright-1, whole_object(nio)%ileft, -1
           WRITE (19, '(2x,f6.1,2x,f6.1,2x,f12.9,2x,f12.9,2x,i1,2x,L1)') REAL(i)+0.5, REAL(j), (REAL(i)+0.5) * delta_x_m, j * delta_x_m, &
                & convert_logical_to_int(whole_object(nio)%segment(4)%cell_is_covered(i)), &
                & whole_object(nio)%segment(4)%cell_is_covered(i)
        END DO
        close (19, status = 'keep')
        print '("### created inner object surface map file ",A24)', iobox_filename
        
     END DO   !###     DO nio = N_of_boundary_objects + 1, N_of_boundary_and_inner_objects
  END IF   !###  IF (Rank_of_process.EQ.0) THEN

  INQUIRE (FILE = 'init_simcontrol.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (.NOT.exists) THEN
     PRINT '(2x,"Process ",i5," : ERROR : init_simcontrol.dat not found. Program terminated")', Rank_of_process
     CALL MPI_FINALIZE(ierr)
     STOP
  END IF

  IF (Rank_of_process.EQ.0) THEN
     PRINT '(2x,"Process ",i5," : init_simcontrol.dat is found. Reading the data file...")', Rank_of_process
  END IF

  OPEN (9, FILE = 'init_simcontrol.dat')

  READ (9, '(A1)') buf !"---ddddddd.ddd----- simulation time [ns]")')
  READ (9, '(3x,f11.3)') t_sim_ns
  READ (9, '(A1)') buf !"--------dd--------- number of electron sub-cycles per ion cycle (odd)")')
  READ (9, '(8x,i2)') N_subcycles
  READ (9, '(A1)') buf !"------dddd--------- number of ion cycles between internal cluster load balancing events")')
  READ (9, '(6x,i4)') dT_cluster_load_balance
  READ (9, '(A1)') buf !"------dddd--------- number of internal cluster load balancing events between global load balancing events")')
  READ (9, '(6x,i4)') dT_global_load_balance
  READ (9, '(A1)') buf !"---ddddddd--------- number of ion cycles between checkpoints (no checkpoints if 0, MPIIO if >0, POSIX if <0)")')
  READ (9, '(3x,i7)') dT_save_checkpoint
  READ (9, '(A1)') buf !"---------d--------- use checkpoint (2/1/0 = Yes, to start/Yes, to continue/No)")')
  READ (9, '(9x,i1)') use_checkpoint
  READ (9, '(A1)') buf !"--dddddddd--------- time step when the checkpoint to be used was saved (dim-less)")')
  READ (9, '(2x,i8)') T_cntr_to_continue

  CLOSE (9, STATUS = 'KEEP')

  Max_T_cntr = t_sim_ns / (1.0d9 * delta_t_s)

! fool proof
  IF (MOD(N_subcycles,2).EQ.0) THEN
     PRINT '("ERROR :: requested even number of electron subcycles per ion cycle")', N_subcycles
     CALL MPI_FINALIZE(ierr)
     STOP
  END IF

  dT_cluster_load_balance = dT_cluster_load_balance * N_subcycles             ! must be an integer number of N_subcycles
  dT_global_load_balance  = dT_global_load_balance * dT_cluster_load_balance  ! must be an integer number of dT_cluster_load_balance

  use_mpiio_checkpoint = .TRUE.
  IF (dT_save_checkpoint.LT.0) THEN
! in this case a checkpoint will consist of many files (one per process)
     use_mpiio_checkpoint = .FALSE.
     dT_save_checkpoint = ABS(dT_save_checkpoint)
  END IF

  dT_save_checkpoint = dT_save_checkpoint * N_subcycles                       ! must be an integer number of N_subcycles

  IF (dT_save_checkpoint.EQ.0) dT_save_checkpoint = -N_subcycles              ! zero dT_save_checkpoint produces unnecessary checkpoint at 
                                                                              ! T_cntr=0 or T_cntr=T_cntr_save_checkpoint

  IF ((use_checkpoint.EQ.0).OR.(use_checkpoint.EQ.2)) THEN
     T_cntr_save_checkpoint = dT_save_checkpoint                              ! default value
  ELSE IF (use_checkpoint.EQ.1) THEN
     T_cntr_save_checkpoint = T_cntr_to_continue + dT_save_checkpoint         ! next save if restart from checkpoint
  END IF

  IF (Rank_of_process.EQ.0) THEN
     PRINT '("### Maximal number of electron time steps ",i8," ###")', Max_T_cntr
     PRINT '("### Internal balancing of cluster load occurs every ",i7," electron time steps ###")', dT_cluster_load_balance
     PRINT '("### Global load balancing occurs every ",i7," electron time steps ###")', dT_global_load_balance
     PRINT '("### Checkpoints will be saves every ",i7," electron time steps ###")', dT_save_checkpoint
     IF (use_checkpoint.EQ.2) THEN
        PRINT '("### Simulation STARTS at saved checkpoint")'
     ELSE IF (use_checkpoint.EQ.1) THEN
        PRINT '("### Simulation CONTINUES at saved checkpoint")'
     ELSE
        PRINT '("### Fresh start, checkpoint NOT USED")'
     END IF
  END IF

  V_scale_ms = N_max_vel * v_Te_ms
  E_scale_Vm = m_e_kg * V_scale_ms / (e_Cl * delta_t_s)
  B_scale_T  = m_e_kg / (e_Cl * delta_t_s)
  F_scale_V  = e_Cl * N_plasma_m3 * delta_x_m**2 / eps_0_Fm

  N_scale_part_m3 = N_plasma_m3 / N_of_particles_cell
  current_factor_Am2 = e_Cl * V_scale_ms * N_scale_part_m3
  energy_factor_eV = 0.5_8 * m_e_kg * V_scale_ms**2 / e_Cl

  temperature_factor_eV = m_e_kg * V_scale_ms**2 / e_Cl
  heat_flow_factor_Wm2 = 0.5_8 * m_e_kg * V_scale_ms**2 * N_scale_part_m3

  given_F_double_period_sys = given_F_double_period_sys / F_scale_V

!  DO n = 1, N_of_boundary_objects
!     whole_object(n)%phi = whole_object(n)%phi / F_scale_V
!     IF (whole_object(n)%object_type.EQ.VACUUM_GAP) whole_object(n)%phi_profile = whole_object(n)%phi_profile / F_scale_V
!  END DO

  cluster_N_blocks = cluster_N_blocks_x * cluster_N_blocks_y
  global_maximal_i = N_blocks_x * N_grid_block_x + 1
  global_maximal_j = N_blocks_y * N_grid_block_y + 1 

  CALL PREPARE_EXTERNAL_FIELDS

! calculate total length of the object (will be used to calculate number of injected particles in each cluster)
  IF ((periodicity_flag.EQ.PERIODICITY_X).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC)) THEN
! account that for periodic systems the boundary in the direction of periodicity has one extra cell for periodic overlapping
     DO n = 1, N_of_boundary_objects
        whole_object(n)%L = 0
        DO m = 1, whole_object(n)%number_of_segments
           IF (whole_object(n)%segment(m)%jend.EQ.whole_object(n)%segment(m)%jstart) THEN
! horizontal segment
              whole_object(n)%L = whole_object(n)%L + ABS( MIN(global_maximal_i-1, whole_object(n)%segment(m)%iend) - &  ! exclude the periodic overlapping cell
                                                         & MIN(global_maximal_i-1, whole_object(n)%segment(m)%istart) )  ! should work when istart>iend
           ELSE IF (whole_object(n)%segment(m)%iend.EQ.whole_object(n)%segment(m)%istart) THEN
! vertical segment
              whole_object(n)%L = whole_object(n)%L + ABS(whole_object(n)%segment(m)%jend - whole_object(n)%segment(m)%jstart)
           END IF
        END DO
     END DO
! the IF branch above should work for periodic systems with segmented electrodes
  ELSE
     DO n = 1, N_of_boundary_objects
        whole_object(n)%L = 0
        DO m = 1, whole_object(n)%number_of_segments
           whole_object(n)%L = whole_object(n)%L + ABS(whole_object(n)%segment(m)%jend - whole_object(n)%segment(m)%jstart) + &
                                                 & ABS(whole_object(n)%segment(m)%iend - whole_object(n)%segment(m)%istart)
        END DO
     END DO
  END IF 

  N_clusters_x = N_blocks_x/cluster_N_blocks_x 
  N_clusters_y = N_blocks_y/cluster_N_blocks_y

  IF (Rank_of_process.EQ.0) THEN
     PRINT '("electron thermal velocity for the scale temperature v_Te_ms = ",e14.7," m/s")', v_Te_ms
     PRINT '("circular electron plasma frequency for the scale density W_plasma_s1 = ",e14.7," s^-1")', W_plasma_s1
     PRINT '("electron Debye length for the scale density/temperature L_debye_m = ",e14.7," m")', L_debye_m
     PRINT '("mesh size delta_x_m = ",e14.7," m")', delta_x_m
     PRINT '("electron time step delta_t_s = ",e14.7," s")', delta_t_s
     PRINT '("scale value of electric field E_scale_Vm = ",e14.7," V/m")', E_scale_Vm
     PRINT '("scale value of magnetic field B_scale_T = ",e14.7," T")', B_scale_T
     PRINT '("scale value of velocity V_scale_ms = ",e14.7," m/s")', V_scale_ms
  END IF

  block_row    = 1 + Rank_of_process / N_blocks_x
  block_column = 1 + Rank_of_process - N_blocks_x * (block_row - 1)

  c_row    = 1 + (block_row  - 1) / cluster_N_blocks_y
  c_column = 1 + (block_column-1) / cluster_N_blocks_x

  CALL IDENTIFY_BLOCK_NEIGHBOURS ! index limits for field-related arrays are set here

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
if (Rank_of_process.eq.0) print *, "IDENTIFY_BLOCK_NEIGHBOURS done"

  CALL IDENTIFY_BLOCK_BOUNDARIES ! 

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
if (Rank_of_process.eq.0) print *, "IDENTIFY_BLOCK_BOUNDARIES done"

  CALL INCLUDE_BLOCK_PERIODICITY

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
if (Rank_of_process.eq.0) print *, "INCLUDE_BLOCK_PERIODICITY done"

  CALL CALCULATE_BLOCK_OFFSET
  
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
if (Rank_of_process.eq.0) print *, "CALCULATE_BLOCK_OFFSET done"

  Rank_of_bottom_left_cluster_master = 0  ! in MPI_COMM_WORLD, ### in future, this may be some other process, fixmeplease!!!

  CALL SET_CLUSTER_STRUCTURE

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
if (Rank_of_process.eq.0) print *, "SET_CLUSTER_STRUCTURE done"

  INQUIRE (FILE = 'init_particles.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (.NOT.exists) THEN   
     PRINT '(2x,"Process ",i3," : ERROR : init_particles.dat not found. Program terminated")', Rank_of_process
     CALL MPI_FINALIZE(ierr)
     STOP
  END IF

  IF (Rank_of_process.EQ.0) THEN
     PRINT '(2x,"Process ",i3," : init_particles.dat is found. Reading the data file...")', Rank_of_process
  END IF

  OPEN (9, FILE = 'init_particles.dat')

  READ (9, '(A1)') buf !"---dddd.ddd----- initial electron temperature [eV]")')
  READ (9, '(3x,f8.3)') init_Te_eV
  READ (9, '(A1)') buf !"---+d.dddE+dd--- initial electron density [m^-3]")')
  READ (9, '(3x,e10.3)') init_Ne_m3
  READ (9, '(A1)') buf !"------d--------- number of ion species")')
  READ (9, '(6x,i1)') N_spec

  ALLOCATE(Qs(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(M_i_amu(1:N_spec), STAT = ALLOC_ERR)

  ALLOCATE(init_Ti_eV(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(init_NiNe(1:N_spec), STAT = ALLOC_ERR)

  READ (9, '(A1)') buf !"---ddd.d---+d---dddd.ddd---ddd.dd-- ion mass [amu] / charge [e] / initial temperature [eV] / initial relative concentration [%]")')
  DO s = 1, N_spec
     READ (9, '(3x,f5.1,3x,i2,3x,f8.3,3x,f6.2)') M_i_amu(s), Qs(s), init_Ti_eV(s), init_NiNe(s)
     init_NiNe(s) = 0.01_8 * init_NiNe(s)   ! convert % to fraction
  END DO

  IF (N_spec.EQ.1) THEN
     init_NiNe(1) = 1.0_8
  END IF
     
  READ (9, '(A1)') buf !"---ddddddddd---- seed for random number generator for process with rank zero")')
  READ (9, '(3x,i9)') init_random_seed

  CLOSE (9, STATUS = 'KEEP')

! now we know N_spec, so we can allocate ion hit counters in wall objects
  DO n = 1, N_of_boundary_and_inner_objects
     ALLOCATE(whole_object(n)%ion_hit_count(1:N_spec), STAT = ALLOC_ERR)
     whole_object(n)%ion_hit_count(1:N_spec) = 0
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL SET_ION_SPECIES

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL PREPARE_SETUP_VALUES                   ! <<<<<<<< SETUP <<<<<<<<< additional processes are here <<<<<<<<<<<<
                                              ! also calls PREPARE_WAVEFORMS
                                              ! also calls PREPARE_EXTERNAL_CIRCUIT

  CALL PREPARE_WALL_MATERIALS                 ! the secondary electron emission is initialized here
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL PREPARE_HT_SETUP_VALUES                ! <<<<<<<< SETUP <<<<<<<<< additional processes are here <<<<<<<<<<<<

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL IDENTIFY_CLUSTER_BOUNDARIES    ! masters only, sets connections between parts of dielectric boundary objects
                                      ! allocates surface_charge for inner objects

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! create communicator for collecting boundary-particle interaction data
! this communicator will never be changed
! only boundary cluster-master processes will be used

  boundary_master = 0
  IF ((cluster_rank_key.EQ.0).AND.(c_N_of_local_object_parts.GT.0)) boundary_master=1

  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, boundary_master, Rank_of_process, COMM_BOUNDARY, ierr)
  CALL MPI_COMM_RANK(COMM_BOUNDARY, Rank_boundary, ierr)
  CALL MPI_COMM_SIZE(COMM_BOUNDARY, N_processes_boundary, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL INCLUDE_CLUSTER_PERIODICITY

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL ANALYZE_CORNERS                ! masters only, for self-connected X-periodic clusters, left/right neighbor master ranks are set to -1 here

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (periodicity_flag.EQ.PERIODICITY_X) THEN

     CALL PREPARE_FFT_X

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL PREPARE_SYS_Y

  ELSE

     CALL SolverInitialization(MPI_COMM_WORLD)

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (use_checkpoint.GT.0) THEN
     CALL READ_A_CHECKPOINT      ! may change particle_master, reads random number generator state, sets all particles
                                 ! calls either READ_CHECKPOINT_MPIIO_2 or READ_CHECKPOINT_POSIX
  END IF

  IF (use_checkpoint.EQ.2) THEN
     T_cntr_global_load_balance = T_cntr_global_load_balance - Start_T_cntr
     T_cntr_cluster_load_balance = T_cntr_cluster_load_balance - Start_T_cntr
     Start_T_cntr = 0
  END IF

  CALL SET_COMMUNICATIONS

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF ((use_checkpoint.NE.0).AND.(cluster_rank_key.EQ.0)) THEN
! checkpoint has permanent surface charge arrays for inner objects, but only the master of the left bottom cluster reads that
! now this master must share the inner object surface charge data with all other masters
     bufsize = 0
     DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
        IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE 
        bufsize = bufsize + whole_object(n)%N_boundary_nodes
     END DO
     IF (bufsize.GT.0) THEN
        ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)
        IF (Rank_horizontal.EQ.0) THEN
! this process must have Rank_of_process==Rank_of_bottom_left_cluster_master
           pos2=0
           DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
              IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE
              pos1 = pos2+1
              pos2 = pos2+whole_object(n)%N_boundary_nodes
              rbufer(pos1:pos2) = whole_object(n)%surface_charge(1:whole_object(n)%N_boundary_nodes)
           END DO
        END IF
        CALL MPI_BCAST(rbufer, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_HORIZONTAL, ierr)
        IF (Rank_horizontal.GT.0) THEN
           pos2=0
           DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
              IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE
              pos1 = pos2+1
              pos2 = pos2+whole_object(n)%N_boundary_nodes
              whole_object(n)%surface_charge(1:whole_object(n)%N_boundary_nodes) = rbufer(pos1:pos2)
           END DO
        END IF
        DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     END IF   !### IF (bufsize.GT.0) THEN
  END IF  !### IF ((use_checkpoint.NE.0).AND.(cluster_rank_key.EQ.0)) THEN

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL COLLECT_MASTERS_INFO           ! masters only

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  N_of_probes_cluster = 0                    ! ensures correct call of DISTRIBUTE_CLUSTER_PARAMETERS below

  CALL DISTRIBUTE_CLUSTER_PARAMETERS         ! particle calculators set accumulated fields to zero here

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (use_checkpoint.EQ.0) THEN

! initialize random numbers generators
     IF (Rank_of_process.EQ.0) THEN 
        PRINT '(/2x,"Process ",i4," : Seed for well_random_seed: ",i12)', Rank_of_process, init_random_seed

        CALL well_random_seed(init_random_seed)

        DO j = 1, 1000000
           myran = well_random_number()
        END DO

        DO i = 1, N_of_processes - 1
           itmp = INT(2000000000.0_8 * well_random_number())                       
           CALL MPI_SEND(itmp, 1, MPI_INTEGER, i, 101, MPI_COMM_WORLD, ierr)
        END DO
     ELSE
        CALL MPI_RECV(init_random_seed, 1, MPI_INTEGER, 0, 101, MPI_COMM_WORLD, stattus, ierr)

        PRINT '(/2x,"Process ",i4," : Seed for well_random_seed: ",i12)', Rank_of_process, init_random_seed

        CALL well_random_seed(init_random_seed)

        DO j = 1, 1000000
           myran = well_random_number()
        END DO

     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL DISTRIBUTE_PARTICLES
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  END IF

  N_e_to_add = 0

  max_N_e_to_add        = 0       !       
  max_N_e_to_send_left  = 0       ! these values will be adjusted
  max_N_e_to_send_right = 0       ! and the proper arrays will be allocated once
  max_N_e_to_send_above = 0       ! a particle will be added to the list 
  max_N_e_to_send_below = 0       ! 

  N_ions_to_add = 0

  max_N_ions_to_add        = 0  ! arrays of values for different ion species
  max_N_ions_to_send_left  = 0  ! same as above,
  max_N_ions_to_send_right = 0  ! these values will be adjusted
  max_N_ions_to_send_above = 0  ! and the proper arrays (ion_todosomething(s)%part) will be allocated once
  max_N_ions_to_send_below = 0  ! a particle will be added to the list
  
  ALLOCATE(ion_to_add(1:N_spec), STAT=ALLOC_ERR)
  ALLOCATE(ion_to_send_left(1:N_spec), STAT=ALLOC_ERR)
  ALLOCATE(ion_to_send_right(1:N_spec), STAT=ALLOC_ERR)
  ALLOCATE(ion_to_send_above(1:N_spec), STAT=ALLOC_ERR)
  ALLOCATE(ion_to_send_below(1:N_spec), STAT=ALLOC_ERR)

  whole_object%phi_const  = whole_object%phi_const / F_scale_V
  whole_object%phi_var    = whole_object%phi_var   / F_scale_V
  whole_object%omega      = whole_object%omega * 2.0d6 * pi * delta_t_s    ! convert from MHz
  whole_object%phase      = whole_object%phase * pi /180.0_8

! do this here to set whole_object%phi which may be used below to calculate external electric field
  CALL UPDATE_WALL_POTENTIALS(0)

  IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN
     phi=0.0_8
  END IF

! save vectors of the external magnetic field
  IF (Rank_cluster.EQ.0) THEN
     BxBy_filename = 'proc_NNNN_BxBy_vs_xy.dat'
     BxBy_filename(6:9) = convert_int_to_txt_string(Rank_of_process,4)
     OPEN (10, FILE = BxBy_filename)
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           WRITE (10, '(2x,i5,2x,i5,2x,f12.9,2x,f12.9,2x,e14.7,2x,e14.7)') &
                & i, &
                & j, &
                & i * delta_x_m, &
                & j * delta_x_m, &
                & B_scale_T * Bx(DBLE(i), DBLE(j)), &
                & B_scale_T * By(DBLE(i), DBLE(j))
        END DO
        WRITE (10, '(" ")')
     END DO
     CLOSE (10, STATUS = 'KEEP')
     PRINT '("file ",A24," is ready")', BxBy_filename
  END IF

END SUBROUTINE INITIATE_PARAMETERS

!--------------------------------------------------------
INTEGER FUNCTION convert_logical_to_int(logvar)
  IMPLICIT NONE
  LOGICAL logvar
  IF (logvar) THEN
     convert_logical_to_int = 1
  ELSE
     convert_logical_to_int = 0
  END IF
  RETURN
END FUNCTION convert_logical_to_int

!--------------------------------------------------------
! must be preceded by a call to IDENTIFY_BLOCK_NEIGHBOURS
! which defines block_row and block_column
!
SUBROUTINE SET_CLUSTER_STRUCTURE

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries, ONLY : block_row, block_column, indx_x_min, indx_x_max, indx_y_min, indx_y_max

  IMPLICIT NONE

  INTEGER n, m
  INTEGER kx, ky, k
  INTEGER ALLOC_ERR
  
  character(20) boxproc_filename ! master_proc_NNNN.dat
  integer devid  ! id of device for writing the file
  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface

! moved to INITIATE_PARAMETERS
!  c_row    = 1 + (block_row  - 1) / cluster_N_blocks_y
!  c_column = 1 + (block_column-1) / cluster_N_blocks_x

  n = 1 + (c_row    - 1) * cluster_N_blocks_y
  m = 1 + (c_column - 1) * cluster_N_blocks_x

  IF ((m.EQ.block_column).AND.(n.EQ.block_row)) THEN
! process is in the left bottom corner of a cluster and thus is declared a master     
     particle_master = Rank_of_process  ! identifies processes belonging to COMM_CLUSTER

     cluster_rank_key = 0          ! this ensures that the master gets zero rank in COMM_CLUSTER
                                   ! the zero value of this key is the only official sign of a particle master 
                                   ! to be used outside of COMM_CLUSTER

     c_indx_x_min = (block_column - 1) * N_grid_block_x
     c_indx_x_max = (block_column - 1 + cluster_N_blocks_x) * N_grid_block_x + 1 
     c_indx_y_min = (block_row    - 1) * N_grid_block_y
     c_indx_y_max = (block_row    - 1 + cluster_N_blocks_y) * N_grid_block_y + 1 

! note that for non-master processes these arrays will be allocated/reallocated in subroutine DISTRIBUTE_CLUSTER_PARAMETERS
! which will be called every time the load balancing takes place and the set of processes in the cluster changes 
     ALLOCATE(EX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(EY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(acc_EX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(acc_EY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     acc_EX=0.0_8
     acc_EY=0.0_8

     IF (periodicity_flag.EQ.PERIODICITY_X) THEN
        ALLOCATE(c_rho(  c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(c_rho_i(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(c_phi(  c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     END IF

     c_X_area_min = DBLE(c_indx_x_min)
     c_X_area_max = DBLE(c_indx_x_max)
     c_Y_area_min = DBLE(c_indx_y_min)
     c_Y_area_max = DBLE(c_indx_y_max)        

! find ranks of master processes in neighbor clusters
     IF (block_column.EQ.1) THEN
        Rank_of_master_left  = -1
     ELSE
        Rank_of_master_left  = Rank_of_process - cluster_N_blocks_x
     END IF

     IF (block_column.EQ.(N_blocks_x - cluster_N_blocks_x + 1)) THEN
        Rank_of_master_right = -1
     ELSE
        Rank_of_master_right = Rank_of_process + cluster_N_blocks_x
     END IF

     IF (block_row.EQ.1) THEN
        Rank_of_master_below = -1
     ELSE
        Rank_of_master_below = Rank_of_process - N_blocks_x * cluster_N_blocks_y
     END IF

     IF (block_row.EQ.(N_blocks_y - cluster_N_blocks_y + 1)) THEN
        Rank_of_master_above = -1
     ELSE
        Rank_of_master_above = Rank_of_process + N_blocks_x * cluster_N_blocks_y
     END IF

! identify color of the cluster
     IF (MOD(c_row,2).EQ.1) THEN            ! rows    1, 3, 5, ...
        IF (MOD(c_column,2).EQ.1) THEN      ! columns 1, 3, 5, ...
           WHITE_CLUSTER = .TRUE.
        ELSE
           WHITE_CLUSTER = .FALSE.
        END IF
     ELSE                                         ! rows    2, 4, 6, ...
        IF (MOD(c_column,2).EQ.0) THEN      ! columns 2, 4, 6, ...
           WHITE_CLUSTER = .TRUE.
        ELSE
           WHITE_CLUSTER = .FALSE.
        END IF
     END IF
        
     IF (WHITE_CLUSTER) THEN
        PRINT '("Process ",i4," : WHITE master ; cluster grid : row ",i4," column ",i4," ; neighbour masters : left ",i4," right ",i4," above ",i4," below ",i4)', &
             & Rank_of_process, c_row, c_column, Rank_of_master_left, Rank_of_master_right, Rank_of_master_above, Rank_of_master_below
     ELSE
        PRINT '("Process ",i4," : BLACK master ; cluster grid : row ",i4," column ",i4," ; neighbour masters : left ",i4," right ",i4," above ",i4," below ",i4)', &
             & Rank_of_process, c_row, c_column, Rank_of_master_left, Rank_of_master_right, Rank_of_master_above, Rank_of_master_below
     END IF
     
! create list of field calculators associated with the cluster, this list will never be changed
! for non-periodic setup the initial blocks do the field calculation
     ALLOCATE(field_calculator(2:cluster_N_blocks), STAT=ALLOC_ERR)
     DO ky = 1, cluster_N_blocks_y
        DO kx = 1, cluster_N_blocks_x
           IF ((kx.EQ.1).AND.(ky.EQ.1)) CYCLE
           k = (ky-1) * cluster_N_blocks_x + kx   ! k starts with 2
           
           field_calculator(k)%rank = Rank_of_process + (ky - 1) * N_blocks_x + kx - 1
           
           field_calculator(k)%indx_x_min = c_indx_x_min + (kx-1) * N_grid_block_x
           field_calculator(k)%indx_x_max = c_indx_x_min +  kx    * N_grid_block_x + 1
           
           field_calculator(k)%indx_y_min = c_indx_y_min + (ky-1) * N_grid_block_y
           field_calculator(k)%indx_y_max = c_indx_y_min +  ky    * N_grid_block_y + 1
        END DO
     END DO
     
     boxproc_filename = 'mfield_calc_NNNN.dat'
     boxproc_filename(13:16) = convert_int_to_txt_string(Rank_of_process, 4)
     devid = 9+Rank_of_process
     open (devid, file=boxproc_filename)
     write (devid, '(3(2x,i5))') indx_x_min, indx_y_min, Rank_of_process
     write (devid, '(3(2x,i5))') indx_x_min, indx_y_max, Rank_of_process
     write (devid, '(3(2x,i5))') indx_x_max, indx_y_max, Rank_of_process
     write (devid, '(3(2x,i5))') indx_x_max, indx_y_min, Rank_of_process
     write (devid, '(3(2x,i5))') indx_x_min, indx_y_min, Rank_of_process
     write (devid, '(" ")')
     write (devid, '(" ")')
     do k = 2, cluster_N_blocks
        write (devid, '(3(2x,i5))') field_calculator(k)%indx_x_min, field_calculator(k)%indx_y_min, field_calculator(k)%rank
        write (devid, '(3(2x,i5))') field_calculator(k)%indx_x_min, field_calculator(k)%indx_y_max, field_calculator(k)%rank
        write (devid, '(3(2x,i5))') field_calculator(k)%indx_x_max, field_calculator(k)%indx_y_max, field_calculator(k)%rank
        write (devid, '(3(2x,i5))') field_calculator(k)%indx_x_max, field_calculator(k)%indx_y_min, field_calculator(k)%rank
        write (devid, '(3(2x,i5))') field_calculator(k)%indx_x_min, field_calculator(k)%indx_y_min, field_calculator(k)%rank
        write (devid, '(" ")')
        write (devid, '(" ")')
     end do
     close (devid, status = 'keep')
     print '("Master process ",i4," : SET_CLUSTER_STRUCTURE : file ",A20," is ready")', Rank_of_process, boxproc_filename
     
!### for setup with periodicity, the initial blocks are still used to parallelize calculation of charge density and electric field
!### in addition, they parallelize the field solver, which is why the field_calculator structure is expanded at PREPARE_FFT_X ###

     boxproc_filename = 'master_proc_NNNN.dat'
     boxproc_filename(13:16) = convert_int_to_txt_string(Rank_of_process, 4)
     devid = 9+Rank_of_process
     open (devid, file=boxproc_filename)
     write (devid, '(3(2x,i5),2(2x,e14.7))') c_indx_x_min, c_indx_y_min, Rank_of_process, c_X_area_min * delta_x_m, c_Y_area_min * delta_x_m
     write (devid, '(3(2x,i5),2(2x,e14.7))') c_indx_x_min, c_indx_y_max, Rank_of_process, c_X_area_min * delta_x_m, c_Y_area_max * delta_x_m
     write (devid, '(3(2x,i5),2(2x,e14.7))') c_indx_x_max, c_indx_y_max, Rank_of_process, c_X_area_max * delta_x_m, c_Y_area_max * delta_x_m
     write (devid, '(3(2x,i5),2(2x,e14.7))') c_indx_x_max, c_indx_y_min, Rank_of_process, c_X_area_max * delta_x_m, c_Y_area_min * delta_x_m
     write (devid, '(3(2x,i5),2(2x,e14.7))') c_indx_x_min, c_indx_y_min, Rank_of_process, c_X_area_min * delta_x_m, c_Y_area_min * delta_x_m
     close (devid, status = 'keep')
     print '("Master process ",i4," : SET_CLUSTER_STRUCTURE : file ",A20," is ready")', Rank_of_process, boxproc_filename

  ELSE
! process is inside a cluster with c_row and c_column, but is not a master of this cluster

! note that m and n calculated in the beginning of the subroutine are the block_column (m) and the block_row (n) of the master process of this cluster

     field_master = N_blocks_x * (n-1) + m - 1    ! this value will never be changed

     particle_master = field_master  ! identifies processes belonging to COMM_CLUSTER
                                     ! for a process which is not a master this value may change
                                     ! initially assign blocks which are inside a cluster to this cluster

     cluster_rank_key = Rank_of_process      ! must be >0 but !=0, 
                                             ! this ensures that the master gets zero rank in COMM_CLUSTER
                                             ! while COMM_CLUSTER ranks of particle calculators will be sorted 
                                             ! in ascending order relative to their MPI_COMM_WORLD ranks
                                             ! the non-zero value of this key is the only official sign of a process
                                             ! which is not a particle master to be used outside of COMM_CLUSTER
                                             ! this value never changes

     PRINT '("Process ",i4," : NOT a master ; my field master : ",i4," with grid column ",i4," row ",i4," MY grid column ",i4," row ",i4)', &
          & Rank_of_process, field_master, m, n, block_column, block_row

  END IF

END SUBROUTINE SET_CLUSTER_STRUCTURE

!------------------------------------------
!
SUBROUTINE SET_ION_SPECIES

  USE CurrentProblemValues
  USE IonParticles

  IMPLICIT NONE

  INTEGER s
  INTEGER ALLOC_ERR

! for simplicity, start with one H+ ion species

  ALLOCATE(      Ms(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(    QM2s(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(QM2sNsub(1:N_spec), STAT = ALLOC_ERR)

  ALLOCATE(    N_ions(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(max_N_ions(1:N_spec), STAT = ALLOC_ERR)

  ALLOCATE(    N_ions_to_add(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(max_N_ions_to_add(1:N_spec), STAT = ALLOC_ERR)

  ALLOCATE(N_ions_to_send_left(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(N_ions_to_send_right(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(N_ions_to_send_above(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(N_ions_to_send_below(1:N_spec), STAT = ALLOC_ERR)

  ALLOCATE(max_N_ions_to_send_left(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(max_N_ions_to_send_right(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(max_N_ions_to_send_above(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(max_N_ions_to_send_below(1:N_spec), STAT = ALLOC_ERR)

  ALLOCATE(ion_to_send_left(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(ion_to_send_right(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(ion_to_send_above(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(ion_to_send_below(1:N_spec), STAT = ALLOC_ERR)

  DO s = 1, N_spec
     Ms(s) = (amu_kg * M_i_amu(s)) / m_e_kg
     QM2s(s) = 0.5_8 * Qs(s) / Ms(s)
     QM2sNsub(s) = N_subcycles * QM2s(s)
  END DO

END SUBROUTINE SET_ION_SPECIES

!-------------------------------------------------
! this subroutine works only with master processes
!
SUBROUTINE IDENTIFY_CLUSTER_BOUNDARIES
  
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER n, m
  INTEGER jstart, jend, istart, iend
  INTEGER minimal_index, maximal_index

  INTEGER ALLOC_ERR
  INTEGER nn, nwo, nsg

  INTEGER ibufer(4)

  IF (cluster_rank_key.NE.0) RETURN

! if there are inner objects, allocate permanent surface charge storage for dielectric inner objects
  DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
     IF (whole_object(n)%object_type.EQ.DIELECTRIC) THEN
        ALLOCATE(whole_object(n)%surface_charge(1:whole_object(n)%N_boundary_nodes), STAT = ALLOC_ERR)
        whole_object(n)%surface_charge = 0.0_8
     END IF
  END DO

! most probable default assumption - no boundary objects

  connect_left  = .FALSE.
  connect_right = .FALSE.
  connect_above = .FALSE.
  connect_below = .FALSE.

  symmetry_plane_X_left = .FALSE.

  n_left = 0
  n_right = 0
  n_below = 0
  n_above = 0

  c_N_of_local_object_parts = 0

  c_N_of_local_object_parts_left  = 0
  c_N_of_local_object_parts_right = 0
  c_N_of_local_object_parts_above = 0
  c_N_of_local_object_parts_below = 0

  IF ( (Rank_of_master_left.GE.0).AND. &
       (Rank_of_master_right.GE.0).AND. &
       (Rank_of_master_below.GE.0).AND. &
       (Rank_of_master_above.GE.0) ) THEN
     PRINT '(" Master process ",i4," is not connected to any boundary object")', Rank_of_process
     RETURN
  END IF

! since we are here, there is at least one boundary object

! check the left edge of the cluster
  IF (Rank_of_master_left.EQ.-1) THEN
     DO n = 1, N_of_boundary_objects
        DO m = 1, whole_object(n)%number_of_segments
           IF ( (whole_object(n)%segment(m)%istart.EQ.c_indx_x_min) .AND. &
              & (whole_object(n)%segment(m)%iend.EQ.c_indx_x_min) ) THEN
              
              jstart = MAX(whole_object(n)%segment(m)%jstart, c_indx_y_min)
              jend   = MIN(whole_object(n)%segment(m)%jend, c_indx_y_max)

              IF (jstart.LT.jend) THEN
                 c_N_of_local_object_parts = c_N_of_local_object_parts + 1
                 IF (c_N_of_local_object_parts.GT.c_max_N_of_local_object_parts) THEN
                    PRINT '("Process ",i4," : ERROR-1 in IDENTIFY_CLUSTER_BOUNDARIES : maximal number of boundary parts exceeded : ",i4)', Rank_of_process, c_N_of_local_object_parts
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 END IF
                 c_local_object_part(c_N_of_local_object_parts)%object_number = n
                 c_local_object_part(c_N_of_local_object_parts)%segment_number = m
                 c_local_object_part(c_N_of_local_object_parts)%istart = c_indx_x_min
                 c_local_object_part(c_N_of_local_object_parts)%jstart = jstart 
                 c_local_object_part(c_N_of_local_object_parts)%iend   = c_indx_x_min
                 c_local_object_part(c_N_of_local_object_parts)%jend   = jend

                 IF (whole_object(n)%object_type.EQ.DIELECTRIC) THEN
                    ALLOCATE(c_local_object_part(c_N_of_local_object_parts)%surface_charge(jstart:jend), STAT = ALLOC_ERR)
                    c_local_object_part(c_N_of_local_object_parts)%surface_charge(jstart:jend) = 0.0_8
                 END IF

                 IF (whole_object(n)%object_type.EQ.SYMMETRY_PLANE) symmetry_plane_X_left = .TRUE.

                 c_N_of_local_object_parts_left = c_N_of_local_object_parts_left+1
                 c_index_of_local_object_part_left(c_N_of_local_object_parts_left) = c_N_of_local_object_parts
              END IF

           END IF
        END DO
     END DO
  END IF

  IF (symmetry_plane_X_left) THEN
     PRINT '("Proc ",i4," (cluster master) has a symmetry plane on its left edge" )', Rank_of_process
  END IF

! check the top edge of the cluster
  IF (Rank_of_master_above.EQ.-1) THEN
     DO n = 1, N_of_boundary_objects
        DO m = 1, whole_object(n)%number_of_segments
           IF ( (whole_object(n)%segment(m)%jstart.EQ.c_indx_y_max) .AND. &
              & (whole_object(n)%segment(m)%jend.EQ.c_indx_y_max) ) THEN
              
              istart = MAX(whole_object(n)%segment(m)%istart, c_indx_x_min)
              iend   = MIN(whole_object(n)%segment(m)%iend, c_indx_x_max)

              IF (istart.LT.iend) THEN
                 c_N_of_local_object_parts = c_N_of_local_object_parts + 1
                 IF (c_N_of_local_object_parts.GT.c_max_N_of_local_object_parts) THEN
                    PRINT '("Process ",i4," : ERROR-2 in IDENTIFY_CLUSTER_BOUNDARIES : maximal number of boundary parts exceeded : ",i4)', Rank_of_process, c_N_of_local_object_parts
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 END IF
                 c_local_object_part(c_N_of_local_object_parts)%object_number = n
                 c_local_object_part(c_N_of_local_object_parts)%segment_number = m
                 c_local_object_part(c_N_of_local_object_parts)%istart = istart
                 c_local_object_part(c_N_of_local_object_parts)%jstart = c_indx_y_max
                 c_local_object_part(c_N_of_local_object_parts)%iend   = iend
                 c_local_object_part(c_N_of_local_object_parts)%jend   = c_indx_y_max

                 IF (whole_object(n)%object_type.EQ.DIELECTRIC) THEN
                    ALLOCATE(c_local_object_part(c_N_of_local_object_parts)%surface_charge(istart:iend), STAT = ALLOC_ERR)
                    c_local_object_part(c_N_of_local_object_parts)%surface_charge(istart:iend) = 0.0_8
                 END IF

                 c_N_of_local_object_parts_above = c_N_of_local_object_parts_above+1
                 c_index_of_local_object_part_above(c_N_of_local_object_parts_above) = c_N_of_local_object_parts
              END IF

           END IF
        END DO
     END DO
  END IF

! check the right edge of the cluster
  IF (Rank_of_master_right.EQ.-1) THEN
     DO n = 1, N_of_boundary_objects
        DO m = 1, whole_object(n)%number_of_segments
           IF ( (whole_object(n)%segment(m)%istart.EQ.c_indx_x_max) .AND. &
              & (whole_object(n)%segment(m)%iend.EQ.c_indx_x_max) ) THEN
              
              jstart = MAX(whole_object(n)%segment(m)%jstart, c_indx_y_min)
              jend   = MIN(whole_object(n)%segment(m)%jend, c_indx_y_max)

              IF (jstart.LT.jend) THEN
                 c_N_of_local_object_parts = c_N_of_local_object_parts + 1
                 IF (c_N_of_local_object_parts.GT.c_max_N_of_local_object_parts) THEN
                    PRINT '("Process ",i4," : ERROR-3 in IDENTIFY_CLUSTER_BOUNDARIES : maximal number of boundary parts exceeded : ",i4)', Rank_of_process, c_N_of_local_object_parts
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 END IF
                 c_local_object_part(c_N_of_local_object_parts)%object_number = n
                 c_local_object_part(c_N_of_local_object_parts)%segment_number = m
                 c_local_object_part(c_N_of_local_object_parts)%istart = c_indx_x_max
                 c_local_object_part(c_N_of_local_object_parts)%jstart = jstart 
                 c_local_object_part(c_N_of_local_object_parts)%iend   = c_indx_x_max
                 c_local_object_part(c_N_of_local_object_parts)%jend   = jend

                 IF (whole_object(n)%object_type.EQ.DIELECTRIC) THEN
                    ALLOCATE(c_local_object_part(c_N_of_local_object_parts)%surface_charge(jstart:jend), STAT = ALLOC_ERR)
                    c_local_object_part(c_N_of_local_object_parts)%surface_charge(jstart:jend) = 0.0_8
                 END IF

                 c_N_of_local_object_parts_right = c_N_of_local_object_parts_right+1
                 c_index_of_local_object_part_right(c_N_of_local_object_parts_right) = c_N_of_local_object_parts
              END IF

           END IF
        END DO
     END DO
  END IF

! check the bottom edge of the cluster
  IF (Rank_of_master_below.EQ.-1) THEN
     DO n = 1, N_of_boundary_objects
        DO m = 1, whole_object(n)%number_of_segments
           IF ( (whole_object(n)%segment(m)%jstart.EQ.c_indx_y_min) .AND. &
              & (whole_object(n)%segment(m)%jend.EQ.c_indx_y_min) ) THEN
              
              istart = MAX(whole_object(n)%segment(m)%istart, c_indx_x_min)
              iend   = MIN(whole_object(n)%segment(m)%iend, c_indx_x_max)

              IF (istart.LT.iend) THEN
                 c_N_of_local_object_parts = c_N_of_local_object_parts + 1
                 IF (c_N_of_local_object_parts.GT.c_max_N_of_local_object_parts) THEN
                    PRINT '("Process ",i4," : ERROR-4 in IDENTIFY_CLUSTER_BOUNDARIES : maximal number of boundary parts exceeded : ",i4)', Rank_of_process, c_N_of_local_object_parts
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 END IF
                 c_local_object_part(c_N_of_local_object_parts)%object_number = n
                 c_local_object_part(c_N_of_local_object_parts)%segment_number = m
                 c_local_object_part(c_N_of_local_object_parts)%istart = istart
                 c_local_object_part(c_N_of_local_object_parts)%jstart = c_indx_y_min
                 c_local_object_part(c_N_of_local_object_parts)%iend   = iend
                 c_local_object_part(c_N_of_local_object_parts)%jend   = c_indx_y_min

                 IF (whole_object(n)%object_type.EQ.DIELECTRIC) THEN
                    ALLOCATE(c_local_object_part(c_N_of_local_object_parts)%surface_charge(istart:iend), STAT = ALLOC_ERR)
                    c_local_object_part(c_N_of_local_object_parts)%surface_charge(istart:iend) = 0.0_8
                 END IF

                 c_N_of_local_object_parts_below = c_N_of_local_object_parts_below+1
                 c_index_of_local_object_part_below(c_N_of_local_object_parts_below) = c_N_of_local_object_parts
              END IF

           END IF
        END DO
     END DO
  END IF

! setup connections between clusters that share dielectric walls

! left edge of the cluster
  DO nn = 1, c_N_of_local_object_parts_left
     n = c_index_of_local_object_part_left(nn)
     nwo = c_local_object_part(n)%object_number
     IF (whole_object(nwo)%object_type.NE.DIELECTRIC) CYCLE
     nsg = c_local_object_part(n)%segment_number
     IF ((c_local_object_part(n)%jstart.EQ.c_indx_y_min).AND.(c_local_object_part(n)%jstart.GT.whole_object(nwo)%segment(nsg)%jstart)) THEN
        connect_below = .TRUE.
        n_below(1) = n
     END IF
     IF ((c_local_object_part(n)%jend.EQ.c_indx_y_max).AND.(c_local_object_part(n)%jend.LT.whole_object(nwo)%segment(nsg)%jend)) THEN
        connect_above = .TRUE.
        n_above(1) = n
     END IF
  END DO

! right edge of the cluster
  DO nn = 1, c_N_of_local_object_parts_right
     n = c_index_of_local_object_part_right(nn)
     nwo = c_local_object_part(n)%object_number
     IF (whole_object(nwo)%object_type.NE.DIELECTRIC) CYCLE
     nsg = c_local_object_part(n)%segment_number
     IF ((c_local_object_part(n)%jstart.EQ.c_indx_y_min).AND.(c_local_object_part(n)%jstart.GT.whole_object(nwo)%segment(nsg)%jstart)) THEN
        connect_below = .TRUE.
        n_below(2) = n
     END IF
     IF ((c_local_object_part(n)%jend.EQ.c_indx_y_max).AND.(c_local_object_part(n)%jend.LT.whole_object(nwo)%segment(nsg)%jend)) THEN
        connect_above = .TRUE.
        n_above(2) = n
     END IF
  END DO

! bottom edge of the cluster
  DO nn = 1, c_N_of_local_object_parts_below
     n = c_index_of_local_object_part_below(nn)
     nwo = c_local_object_part(n)%object_number
     IF (whole_object(nwo)%object_type.NE.DIELECTRIC) CYCLE
     nsg = c_local_object_part(n)%segment_number
     IF ( ((c_local_object_part(n)%istart.EQ.c_indx_x_min).AND.(c_local_object_part(n)%istart.GT.whole_object(nwo)%segment(nsg)%istart)) &
        & .OR. &
        & ((periodicity_flag.EQ.PERIODICITY_X).AND.(c_local_object_part(n)%istart.EQ.0)) ) THEN
        connect_left = .TRUE.
        n_left(1) = n
     END IF
     IF ( ((c_local_object_part(n)%iend.EQ.c_indx_x_max).AND.(c_local_object_part(n)%iend.LT.whole_object(nwo)%segment(nsg)%iend)) &
        & .OR. &
        & ((periodicity_flag.EQ.PERIODICITY_X).AND.(c_local_object_part(n)%iend.EQ.global_maximal_i)) ) THEN
        connect_right = .TRUE.
        n_right(1) = n
     END IF
  END DO

! top edge of the cluster
  DO nn = 1, c_N_of_local_object_parts_above
     n = c_index_of_local_object_part_above(nn)
     nwo = c_local_object_part(n)%object_number
     IF (whole_object(nwo)%object_type.NE.DIELECTRIC) CYCLE
     nsg = c_local_object_part(n)%segment_number
     IF ( ((c_local_object_part(n)%istart.EQ.c_indx_x_min).AND.(c_local_object_part(n)%istart.GT.whole_object(nwo)%segment(nsg)%istart)) &
        & .OR. &
        & ((periodicity_flag.EQ.PERIODICITY_X).AND.(c_local_object_part(n)%istart.EQ.0)) ) THEN
        connect_left = .TRUE.
        n_left(2) = n
     END IF
     IF ( ((c_local_object_part(n)%iend.EQ.c_indx_x_max).AND.(c_local_object_part(n)%iend.LT.whole_object(nwo)%segment(nsg)%iend)) &
        & .OR. &
        & ((periodicity_flag.EQ.PERIODICITY_X).AND.(c_local_object_part(n)%iend.EQ.global_maximal_i)) ) THEN
        connect_right = .TRUE.
        n_right(2) = n
     END IF
  END DO

END SUBROUTINE IDENTIFY_CLUSTER_BOUNDARIES

!-------------------------------------------------
! this subroutine works mostly with master processes
! but other processes are involved as well
! in order to create new temporary communicators
!
SUBROUTINE INCLUDE_CLUSTER_PERIODICITY
  
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER periodic_message(1:13)
  INTEGER n, m, nwo

  INTEGER, ALLOCATABLE :: all_periodic_messages(:,:)
  INTEGER ALLOC_ERR
  LOGICAL pair_found
  INTEGER ibufer(1:4)
  
! this subroutine will be called BEFORE first call of SET_COMMUNICATIONS
! therefore here we create temporary horizontal communicators which will be destroyed at the end

  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, cluster_rank_key, particle_master, COMM_HORIZONTAL, ierr)
  CALL MPI_COMM_RANK(COMM_HORIZONTAL, Rank_horizontal, ierr)
  CALL MPI_COMM_SIZE(COMM_HORIZONTAL, N_processes_horizontal, ierr)

  IF (cluster_rank_key.EQ.0) THEN
! check the presence of periodic boundaries

     periodic_message(1:12)=-1
     periodic_message(13) = Rank_of_process       ! the global rank

     periodic_boundary_X_left  = .FALSE.
     periodic_boundary_X_right = .FALSE.
     periodic_boundary_Y_below = .FALSE.
     periodic_boundary_Y_above = .FALSE.

     i_period_x = 0
     i_period_y = 0
     L_period_X = DBLE(i_period_x)
     L_period_Y = DBLE(i_period_y)

! left boundary (periodicity in X)
     DO n = 1, c_N_of_local_object_parts_left
        m = c_index_of_local_object_part_left(n)
        nwo = c_local_object_part(m)%object_number
        IF (whole_object(nwo)%object_type.EQ.PERIODIC_PIPELINE_X) THEN
           periodic_message(1) = nwo   ! 1<=nwo<=N_of_boundary_objects
           periodic_message(2) = c_local_object_part(m)%jstart
           periodic_message(3) = c_local_object_part(m)%jend
        END IF
     END DO

! right boundary (periodicity in X)
     DO n = 1, c_N_of_local_object_parts_right
        m = c_index_of_local_object_part_right(n)
        nwo = c_local_object_part(m)%object_number
        IF (whole_object(nwo)%object_type.EQ.PERIODIC_PIPELINE_X) THEN
           periodic_message(4) = nwo   ! 1<=nwo<=N_of_boundary_objects
           periodic_message(5) = c_local_object_part(m)%jstart
           periodic_message(6) = c_local_object_part(m)%jend
        END IF
     END DO

! bottom boundary (periodicity in Y)
     DO n = 1, c_N_of_local_object_parts_below
        m = c_index_of_local_object_part_below(n)
        nwo = c_local_object_part(m)%object_number
        IF (whole_object(nwo)%object_type.EQ.PERIODIC_PIPELINE_Y) THEN
           periodic_message(7) = nwo   ! 1<=nwo<=N_of_boundary_objects
           periodic_message(8) = c_local_object_part(m)%istart
           periodic_message(9) = c_local_object_part(m)%iend
        END IF
     END DO

! top boundary (periodicity in Y)
     DO n = 1, c_N_of_local_object_parts_above
        m = c_index_of_local_object_part_above(n)
        nwo = c_local_object_part(m)%object_number
        IF (whole_object(nwo)%object_type.EQ.PERIODIC_PIPELINE_Y) THEN
           periodic_message(10) = nwo   ! 1<=nwo<=N_of_boundary_objects
           periodic_message(11) = c_local_object_part(m)%istart
           periodic_message(12) = c_local_object_part(m)%iend
        END IF
     END DO

! send/receive-process periodic_message to the cluster with rank zero, use COMM_HORIZONTAL
     IF (Rank_horizontal.EQ.0) THEN
        ALLOCATE(all_periodic_messages(1:15,0:N_processes_horizontal-1), STAT=ALLOC_ERR)
        all_periodic_messages = -1
        all_periodic_messages(1:12,0) = periodic_message(1:12)
        all_periodic_messages(13,0) = Rank_of_process
! receive periodic_message from all other clusters
        DO n = 1, N_processes_horizontal-1
           CALL MPI_RECV(all_periodic_messages(1:13,n), 13, MPI_INTEGER, n, n, COMM_HORIZONTAL, stattus, ierr)        
        END DO

! find pairs of clusters with periodic boundary in X
        DO n = 0, N_processes_horizontal-1
! find one such a cluster
           IF (all_periodic_messages(1,n).GT.0) THEN
! find a matching pair
              pair_found=.FALSE.
              DO m = 0, N_processes_horizontal-1
                 IF ((all_periodic_messages(4,m).GT.0).AND.(all_periodic_messages(5,m).EQ.all_periodic_messages(2,n)).AND.(all_periodic_messages(6,m).EQ.all_periodic_messages(3,n))) THEN
                    pair_found=.TRUE.
                    all_periodic_messages(14,n) = all_periodic_messages(13,m)  ! MPI_COMM_WORLD rank of process m from COMM_HORIZONTAL
                    all_periodic_messages(14,m) = all_periodic_messages(13,n)  ! MPI_COMM_WORLD rank of process n from COMM_HORIZONTAL
                    EXIT
                 END IF
              END DO
              IF (.NOT.pair_found) THEN
! error message if the pair is not found
                 PRINT '("Process ",i4," :: Error in INCLUDE_PERIODICITY :: cannot find a matching pair for X-periodic boudnary of process ",i4," (horizontal rank ",i4,")")', Rank_of_process, all_periodic_messages(13,n), n
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
           END IF
        END DO

! find pairs of clusters with periodic boundary in Y
        DO n = 0, N_processes_horizontal-1
! find one such a cluster
           IF (all_periodic_messages(7,n).GT.0) THEN
! find a matching pair
              pair_found=.FALSE.
              DO m = 0, N_processes_horizontal-1
                 IF ((all_periodic_messages(10,m).GT.0).AND.(all_periodic_messages(11,m).EQ.all_periodic_messages(8,n)).AND.(all_periodic_messages(12,m).EQ.all_periodic_messages(9,n))) THEN
                    pair_found=.TRUE.
                    all_periodic_messages(15,n) = all_periodic_messages(13,m)  ! MPI_COMM_WORLD rank of process m from COMM_HORIZONTAL
                    all_periodic_messages(15,m) = all_periodic_messages(13,n)  ! MPI_COMM_WORLD rank of process n from COMM_HORIZONTAL
                    EXIT
                 END IF
              END DO
              IF (.NOT.pair_found) THEN
! error message if the pair is not found
                 PRINT '("Process ",i4," :: Error in INCLUDE_PERIODICITY :: cannot find a matching pair for Y-periodic boundary of process ",i4," (horizontal rank ",i4,")")', Rank_of_process, all_periodic_messages(13,n), n
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
           END IF
        END DO


print '(4(2x,i4))', 0, all_periodic_messages(13:15,0)
! send  information about periodicity pairs to each process
        DO n = 1, N_processes_horizontal-1
print '(4(2x,i4))', n, all_periodic_messages(13:15,n)
           CALL MPI_SEND(all_periodic_messages(14:15,n), 2, MPI_INTEGER, n, Rank_horizontal, COMM_HORIZONTAL, request, ierr)        
        END DO
     
! process own periodicity pair info
        IF (all_periodic_messages(14,0).GT.0) THEN
           IF (Rank_of_master_left.LT.0) THEN
              Rank_of_master_left  = all_periodic_messages(14,0)
              periodic_boundary_X_left = .TRUE.
           END IF
           IF (Rank_of_master_right.LT.0) THEN
              Rank_of_master_right = all_periodic_messages(14,0)
              periodic_boundary_X_right = .TRUE.
           END IF
        END IF
        IF (all_periodic_messages(15,0).GT.0) THEN
           IF (Rank_of_master_below.LT.0) THEN
              Rank_of_master_below = all_periodic_messages(15,0)
              periodic_boundary_Y_below = .TRUE.
           END IF
           IF (Rank_of_master_above.LT.0) THEN
              Rank_of_master_above = all_periodic_messages(15,0)
              periodic_boundary_Y_above = .TRUE.
           END IF
        END IF

        DEALLOCATE(all_periodic_messages, STAT=ALLOC_ERR)
 
     ELSE
! send periodic_message to zero-rank cluster
        CALL MPI_SEND(periodic_message, 13, MPI_INTEGER, 0, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
! receive information about periodicity pairs
        CALL MPI_RECV(ibufer(1:2), 2, MPI_INTEGER, 0, 0, COMM_HORIZONTAL, stattus, ierr) 
        IF (ibufer(1).GE.0) THEN
           IF (Rank_of_master_left.LT.0) THEN
              Rank_of_master_left  = ibufer(1)
              periodic_boundary_X_left = .TRUE.
           END IF
           IF (Rank_of_master_right.LT.0) THEN
              Rank_of_master_right = ibufer(1)
              periodic_boundary_X_right = .TRUE.
           END IF
        END IF
        IF (ibufer(2).GE.0) THEN
           IF (Rank_of_master_below.LT.0) THEN
              Rank_of_master_below = ibufer(2)
              periodic_boundary_Y_below = .TRUE.
           END IF
           IF (Rank_of_master_above.LT.0) THEN
              Rank_of_master_above = ibufer(2)
              periodic_boundary_Y_above = .TRUE.
           END IF
        END IF
     END IF

  END IF

print '("Process ",i4," :: P.B. L/R/B/A :: ",4(1x,L1))', Rank_of_process, periodic_boundary_X_left, periodic_boundary_X_right, periodic_boundary_Y_below, periodic_boundary_Y_above

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! all processes must participate in destroying the communicator, I guess...
  CALL MPI_COMM_FREE(COMM_HORIZONTAL, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! let master processes with periodic boundaries communicate with each other to get the distances between the boundaries

  IF (cluster_rank_key.EQ.0) THEN

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
        i_period_x = c_indx_x_max - c_indx_x_min - 1
        L_period_X = DBLE(i_period_x)
!     Rank_of_master_left = -1      ! this is peformed in the end of ANALYZE_CORNERS to allow reasonable values of corner types
!     Rank_of_master_right = -1     ! although it is overkill since the corner types are not used in this case
        PRINT '("INCLUDE_PERIODICITY :: master process ",i4," is periodically [X] connected to itself, L_period_X = ",f10.4)', Rank_of_process, L_period_X
     ELSE IF (periodic_boundary_X_left) THEN
! send 
        CALL MPI_SEND(c_indx_x_min, 1, MPI_INTEGER, Rank_of_master_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! receive 
        CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, Rank_of_master_left, Rank_of_master_left, MPI_COMM_WORLD, stattus, ierr) 
        i_period_x = ibufer(1) - c_indx_x_min - 1
        L_period_X = DBLE(i_period_x)
! report
        PRINT '("INCLUDE_PERIODICITY :: master process ",i4," communicated with process ",i4," and got L_period_X = ",f10.4)', Rank_of_process, Rank_of_master_left, L_period_X
     ELSE IF (periodic_boundary_X_right) THEN
! receive
        CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, Rank_of_master_right, Rank_of_master_right, MPI_COMM_WORLD, stattus, ierr) 
        i_period_x = c_indx_x_max - ibufer(1) - 1   ! last cell is overlapping
        L_period_X = DBLE(i_period_x)
! send
        CALL MPI_SEND(c_indx_x_max, 1, MPI_INTEGER, Rank_of_master_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! report
        PRINT '("INCLUDE_PERIODICITY :: master process ",i4," communicated with process ",i4," and got L_period_X = ",f10.4)', Rank_of_process, Rank_of_master_right, L_period_X
     END IF

!  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! ######## future work - needs case when a single cluster is periodically connected to itself in both X- and Y-directions #######

     IF (periodic_boundary_Y_below) THEN
! send 
        CALL MPI_SEND(c_indx_y_min, 1, MPI_INTEGER, Rank_of_master_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! receive 
        CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, Rank_of_master_below, Rank_of_master_below, MPI_COMM_WORLD, stattus, ierr) 
        i_period_y = ibufer(1) - c_indx_y_min - 1
        L_period_Y = DBLE(i_period_y)
! report
        PRINT '("INCLUDE_PERIODICITY :: master process ",i4," communicated with process ",i4," and got L_period_Y = ",f10.4)', Rank_of_process, Rank_of_master_below, L_period_Y
     ELSE IF (periodic_boundary_Y_above) THEN
! receive
        CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, Rank_of_master_above, Rank_of_master_above, MPI_COMM_WORLD, stattus, ierr) 
        i_period_y = c_indx_y_max - ibufer(1) - 1   ! last cell is overlapping
        L_period_Y = DBLE(i_period_y)   ! last cell is overlapping
! send
        CALL MPI_SEND(c_indx_y_max, 1, MPI_INTEGER, Rank_of_master_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! report
        PRINT '("INCLUDE_PERIODICITY :: master process ",i4," communicated with process ",i4," and got L_period_Y = ",f10.4)', Rank_of_process, Rank_of_master_above, L_period_Y
     END IF

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  IF (cluster_rank_key.NE.0) RETURN

! check that periodicity along X on the left boundary is not mixed with the symmetry plane
  IF (periodic_boundary_X_left.AND.symmetry_plane_X_left) THEN
     PRINT '("Proc ",i4," :: Error in INCLUDE_CLUSTER_PERIODICITY, both periodicity and symmetry plane at the left edge of a cluster detected, program terminated")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

! check connections 

  IF (WHITE_CLUSTER) THEN  
! "white processes"

     IF (connect_right) THEN
! ## 1 ## send right data on the boundary piece(s) that will participate in surface charge exchange 
        ibufer(1:4) = 0
        DO n = 1, 2
           IF (n_right(n).GT.0) THEN
              ibufer(2*n-1) = c_local_object_part(n_right(n))%object_number
              ibufer(2*n)   = c_local_object_part(n_right(n))%segment_number
           END IF
        END DO
        CALL MPI_SEND(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_left) THEN
! ## 2 ## send left data on the boundary piece(s) that will participate in surface charge exchange 
        ibufer(1:4) = 0
        DO n = 1, 2
           IF (n_left(n).GT.0) THEN
              ibufer(2*n-1) = c_local_object_part(n_left(n))%object_number
              ibufer(2*n)   = c_local_object_part(n_left(n))%segment_number
           END IF
        END DO
        CALL MPI_SEND(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_left) THEN
! ## 3 ## receive from left data on the boundary piece that will participate in surface charge exchange 
        CALL MPI_RECV(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_left, Rank_of_master_left, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_left(n).GT.0) THEN
              IF (ibufer(2*n-1).NE.c_local_object_part(n_left(n))%object_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its left neighbor ",i4," status = ERROR, object number does not match")', Rank_of_process, Rank_of_master_left
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE IF (ibufer(2*n).NE.c_local_object_part(n_left(n))%segment_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its left neighbor ",i4," status = ERROR, segment number does not match")', Rank_of_process, Rank_of_master_left
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE
                 PRINT '("Cluster ",i4," received connection info from its left neighbor ",i4," status = SUCCESS")', Rank_of_process, Rank_of_master_left
              END IF
           END IF
        END DO
     END IF

     IF (connect_right) THEN
! ## 4 ## receive from right data on the boundary piece that will participate in surface charge exchange 
        CALL MPI_RECV(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_right, Rank_of_master_right, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_right(n).GT.0) THEN
              IF (ibufer(2*n-1).NE.c_local_object_part(n_right(n))%object_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its right neighbor ",i4," status = ERROR, object number does not match")', Rank_of_process, Rank_of_master_right
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE IF (ibufer(2*n).NE.c_local_object_part(n_right(n))%segment_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its right neighbor ",i4," status = ERROR, segment number does not match")', Rank_of_process, Rank_of_master_right
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE
                 PRINT '("Cluster ",i4," received connection info from its right neighbor ",i4," status = SUCCESS")', Rank_of_process, Rank_of_master_right
              END IF
           END IF
        END DO
     END IF

     IF (connect_above) THEN
! ## 5 ## send above data on the boundary piece that will participate in surface charge exchange 
        ibufer(1:4) = 0
        DO n = 1, 2
           IF (n_above(n).GT.0) THEN
              ibufer(2*n-1) = c_local_object_part(n_above(n))%object_number
              ibufer(2*n)   = c_local_object_part(n_above(n))%segment_number
           END IF
        END DO
        CALL MPI_SEND(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_below) THEN
! ## 6 ## send below data on the boundary piece that will participate in surface charge exchange 
        ibufer(1:4) = 0
        DO n = 1, 2
           IF (n_below(n).GT.0) THEN
              ibufer(2*n-1) = c_local_object_part(n_below(n))%object_number
              ibufer(2*n)   = c_local_object_part(n_below(n))%segment_number
           END IF
        END DO
        CALL MPI_SEND(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_below) THEN
! ## 7 ## receive from below data on the boundary piece that will participate in surface charge exchange 
        CALL MPI_RECV(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_below, Rank_of_master_below, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_below(n).GT.0) THEN
              IF (ibufer(2*n-1).NE.c_local_object_part(n_below(n))%object_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its neighbor below ",i4," status = ERROR, object number does not match")', Rank_of_process, Rank_of_master_below
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE IF (ibufer(2*n).NE.c_local_object_part(n_below(n))%segment_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its neighbor below ",i4," status = ERROR, segment number does not match")', Rank_of_process, Rank_of_master_below
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE
                 PRINT '("Cluster ",i4," received connection info from its neighbor below ",i4," status = SUCCESS")', Rank_of_process, Rank_of_master_below
              END IF
           END IF
        END DO
     END IF

     IF (connect_above) THEN
! ## 8 ## receive from above data on the boundary piece that will participate in surface charge exchange 
        CALL MPI_RECV(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_above, Rank_of_master_above, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_above(n).GT.0) THEN
              IF (ibufer(2*n-1).NE.c_local_object_part(n_above(n))%object_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its neighbor above ",i4," status = ERROR, object number does not match")', Rank_of_process, Rank_of_master_above
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE IF (ibufer(2*n).NE.c_local_object_part(n_above(n))%segment_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its neighbor above ",i4," status = ERROR, segment number does not match")', Rank_of_process, Rank_of_master_above
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE
                 PRINT '("Cluster ",i4," received connection info from its neighbor above ",i4," status = SUCCESS")', Rank_of_process, Rank_of_master_above
              END IF
           END IF
        END DO
     END IF

  ELSE
! "black" processes

     IF (connect_left) THEN
! ## 1 ## receive from left 
        CALL MPI_RECV(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_left, Rank_of_master_left, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_left(n).GT.0) THEN
              IF (ibufer(2*n-1).NE.c_local_object_part(n_left(n))%object_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its left neighbor ",i4," status = ERROR, object number does not match")', Rank_of_process, Rank_of_master_left
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE IF (ibufer(2*n).NE.c_local_object_part(n_left(n))%segment_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its left neighbor ",i4," status = ERROR, segment number does not match")', Rank_of_process, Rank_of_master_left
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE
                 PRINT '("Cluster ",i4," received connection info from its left neighbor ",i4," status = SUCCESS")', Rank_of_process, Rank_of_master_left
              END IF
           END IF
        END DO
     END IF

     IF (connect_right) THEN
! ## 2 ## receive from right 
        CALL MPI_RECV(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_right, Rank_of_master_right, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_right(n).GT.0) THEN
              IF (ibufer(2*n-1).NE.c_local_object_part(n_right(n))%object_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its right neighbor ",i4," status = ERROR, object number does not match")', Rank_of_process, Rank_of_master_right
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE IF (ibufer(2*n).NE.c_local_object_part(n_right(n))%segment_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its right neighbor ",i4," status = ERROR, segment number does not match")', Rank_of_process, Rank_of_master_right
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE
                 PRINT '("Cluster ",i4," received connection info from its right neighbor ",i4," status = SUCCESS")', Rank_of_process, Rank_of_master_right
              END IF
           END IF
        END DO
     END IF

     IF (connect_right) THEN
! ## 3 ## send right 
        ibufer(1:4) = 0
        DO n = 1, 2
           IF (n_right(n).GT.0) THEN
              ibufer(2*n-1) = c_local_object_part(n_right(n))%object_number
              ibufer(2*n)   = c_local_object_part(n_right(n))%segment_number
           END IF
        END DO
        CALL MPI_SEND(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_left) THEN
! ## 4 ## send left 
        ibufer(1:4) = 0
        DO n = 1, 2
           IF (n_left(n).GT.0) THEN
              ibufer(2*n-1) = c_local_object_part(n_left(n))%object_number
              ibufer(2*n)   = c_local_object_part(n_left(n))%segment_number
           END IF
        END DO
        CALL MPI_SEND(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_below) THEN
! ## 5 ## receive from below 
        CALL MPI_RECV(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_below, Rank_of_master_below, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_below(n).GT.0) THEN
              IF (ibufer(2*n-1).NE.c_local_object_part(n_below(n))%object_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its neighbor below ",i4," status = ERROR, object number does not match")', Rank_of_process, Rank_of_master_below
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE IF (ibufer(2*n).NE.c_local_object_part(n_below(n))%segment_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its neighbor below ",i4," status = ERROR, segment number does not match")', Rank_of_process, Rank_of_master_below
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE
                 PRINT '("Cluster ",i4," received connection info from its neighbor below ",i4," status = SUCCESS")', Rank_of_process, Rank_of_master_below
              END IF
           END IF
        END DO
     END IF

     IF (connect_above) THEN
! ## 6 ## receive from above 
        CALL MPI_RECV(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_above, Rank_of_master_above, MPI_COMM_WORLD, stattus, ierr)
        DO n = 1, 2
           IF (n_above(n).GT.0) THEN
              IF (ibufer(2*n-1).NE.c_local_object_part(n_above(n))%object_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its neighbor above ",i4," status = ERROR, object number does not match")', Rank_of_process, Rank_of_master_above
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE IF (ibufer(2*n).NE.c_local_object_part(n_above(n))%segment_number) THEN
                 PRINT '("Cluster ",i4," received connection info from its neighbor above ",i4," status = ERROR, segment number does not match")', Rank_of_process, Rank_of_master_above
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
               ELSE
                 PRINT '("Cluster ",i4," received connection info from its neighbor above ",i4," status = SUCCESS")', Rank_of_process, Rank_of_master_above
              END IF
           END IF
        END DO
     END IF

     IF (connect_above) THEN
! ## 7 ## send above 
        ibufer(1:4) = 0
        DO n = 1, 2
           IF (n_above(n).GT.0) THEN
              ibufer(2*n-1) = c_local_object_part(n_above(n))%object_number
              ibufer(2*n)   = c_local_object_part(n_above(n))%segment_number
           END IF
        END DO
        CALL MPI_SEND(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (connect_below) THEN
! ## 8 ## send below 
        ibufer(1:4) = 0
        DO n = 1, 2
           IF (n_below(n).GT.0) THEN
              ibufer(2*n-1) = c_local_object_part(n_below(n))%object_number
              ibufer(2*n)   = c_local_object_part(n_below(n))%segment_number
           END IF
        END DO
        CALL MPI_SEND(ibufer(1:4), 4, MPI_INTEGER, Rank_of_master_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

  END IF

END SUBROUTINE INCLUDE_CLUSTER_PERIODICITY

!-------------------------------------------------
! this subroutine works only with master processes
!
SUBROUTINE ANALYZE_CORNERS
  
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INTEGER n, m
!  INTEGER jstart, jend, istart, iend
  INTEGER minimal_index, maximal_index

  IF (cluster_rank_key.NE.0) RETURN

! most probable default assumption - no boundary objects

  c_left_bottom_corner_type  = HAS_TWO_NEIGHBORS
  c_left_top_corner_type     = HAS_TWO_NEIGHBORS
  c_right_bottom_corner_type = HAS_TWO_NEIGHBORS
  c_right_top_corner_type    = HAS_TWO_NEIGHBORS

  IF ( (Rank_of_master_left.GE.0).AND. &
       (Rank_of_master_right.GE.0).AND. &
       (Rank_of_master_below.GE.0).AND. &
       (Rank_of_master_above.GE.0) ) THEN
     PRINT '(" Master process ",i4," is not connected to any boundary object")', Rank_of_process
     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
        Rank_of_master_left = -1
        Rank_of_master_right = -1
        PRINT '("ANALYZE_CORNERS :: master process ",i4," is periodically [X] connected to itself, ranks of left/right neighbors set to -1")', Rank_of_process
     END IF
     RETURN
  END IF

! since we are here, there is at least one boundary object

! analyze configuration in the corners

! left bottom corner
  IF ((Rank_of_master_left.LT.0).AND.(Rank_of_master_below.LT.0)) THEN
     c_left_bottom_corner_type = SURROUNDED_BY_WALL           !###
  ELSE IF ((Rank_of_master_left.LT.0).AND.(Rank_of_master_below.GE.0)) THEN
     c_left_bottom_corner_type = FLAT_WALL_LEFT               !###
! check for the empty corner
     minimal_index=c_indx_y_max
     DO n = 1, c_N_of_local_object_parts_left
        IF ( c_local_object_part(c_index_of_local_object_part_left(n))%jstart.LT.minimal_index ) THEN
           minimal_index = c_local_object_part(c_index_of_local_object_part_left(n))%jstart
        END IF
     END DO
     IF (minimal_index.GT.c_indx_y_min) THEN
        c_left_bottom_corner_type = EMPTY_CORNER_WALL_LEFT    !###
     END IF
  ELSE IF ((Rank_of_master_left.GE.0).AND.(Rank_of_master_below.LT.0)) THEN
     c_left_bottom_corner_type = FLAT_WALL_BELOW              !###
! check for the empty corner
     minimal_index=c_indx_x_max
     DO n = 1, c_N_of_local_object_parts_below
        IF ( c_local_object_part(c_index_of_local_object_part_below(n))%istart.LT.minimal_index ) THEN
           minimal_index = c_local_object_part(c_index_of_local_object_part_below(n))%istart
        END IF
     END DO
     IF (minimal_index.GT.c_indx_x_min) THEN
        c_left_bottom_corner_type = EMPTY_CORNER_WALL_BELOW   !###
     END IF
  END IF

! left top corner
  IF ((Rank_of_master_left.LT.0).AND.(Rank_of_master_above.LT.0)) THEN
     c_left_top_corner_type = SURROUNDED_BY_WALL           !###
  ELSE IF ((Rank_of_master_left.LT.0).AND.(Rank_of_master_above.GE.0)) THEN
     c_left_top_corner_type = FLAT_WALL_LEFT               !###
! check for the empty corner
     maximal_index=c_indx_y_min
     DO n = 1, c_N_of_local_object_parts_left
        IF ( c_local_object_part(c_index_of_local_object_part_left(n))%jend.GT.maximal_index ) THEN
           maximal_index = c_local_object_part(c_index_of_local_object_part_left(n))%jend
        END IF
     END DO
     IF (maximal_index.LT.c_indx_y_max) THEN
        c_left_top_corner_type = EMPTY_CORNER_WALL_LEFT    !###
     END IF
  ELSE IF ((Rank_of_master_left.GE.0).AND.(Rank_of_master_above.LT.0)) THEN
     c_left_top_corner_type = FLAT_WALL_ABOVE              !###
! check for the empty corner
     minimal_index=c_indx_x_max
     DO n = 1, c_N_of_local_object_parts_above
        IF ( c_local_object_part(c_index_of_local_object_part_above(n))%istart.LT.minimal_index ) THEN
           minimal_index = c_local_object_part(c_index_of_local_object_part_above(n))%istart
        END IF
     END DO
     IF (minimal_index.GT.c_indx_x_min) THEN
        c_left_top_corner_type = EMPTY_CORNER_WALL_ABOVE   !###
     END IF
  END IF

! right bottom corner
  IF ((Rank_of_master_right.LT.0).AND.(Rank_of_master_below.LT.0)) THEN
     c_right_bottom_corner_type = SURROUNDED_BY_WALL           !###
  ELSE IF ((Rank_of_master_right.LT.0).AND.(Rank_of_master_below.GE.0)) THEN
     c_right_bottom_corner_type = FLAT_WALL_RIGHT              !###
! check for the empty corner
     minimal_index=c_indx_y_max
     DO n = 1, c_N_of_local_object_parts_right
        IF ( c_local_object_part(c_index_of_local_object_part_right(n))%jstart.LT.minimal_index ) THEN
           minimal_index = c_local_object_part(c_index_of_local_object_part_right(n))%jstart
        END IF
     END DO
     IF (minimal_index.GT.c_indx_y_min) THEN
        c_right_bottom_corner_type = EMPTY_CORNER_WALL_RIGHT   !###
     END IF
  ELSE IF ((Rank_of_master_right.GE.0).AND.(Rank_of_master_below.LT.0)) THEN
     c_right_bottom_corner_type = FLAT_WALL_BELOW              !###
! check for the empty corner
     maximal_index=c_indx_x_min
     DO n = 1, c_N_of_local_object_parts_below
        IF ( c_local_object_part(c_index_of_local_object_part_below(n))%iend.GT.maximal_index ) THEN
           maximal_index = c_local_object_part(c_index_of_local_object_part_below(n))%iend
        END IF
     END DO
     IF (maximal_index.LT.c_indx_x_max) THEN
        c_right_bottom_corner_type = EMPTY_CORNER_WALL_BELOW   !###
     END IF
  END IF

! right top corner
  IF ((Rank_of_master_right.LT.0).AND.(Rank_of_master_above.LT.0)) THEN
     c_right_top_corner_type = SURROUNDED_BY_WALL           !###
  ELSE IF ((Rank_of_master_right.LT.0).AND.(Rank_of_master_above.GE.0)) THEN
     c_right_top_corner_type = FLAT_WALL_RIGHT              !###
! check for the empty corner
     maximal_index=c_indx_y_min
     DO n = 1, c_N_of_local_object_parts_right
        IF ( c_local_object_part(c_index_of_local_object_part_right(n))%jend.GT.maximal_index ) THEN
           maximal_index = c_local_object_part(c_index_of_local_object_part_right(n))%jend
        END IF
     END DO
     IF (maximal_index.LT.c_indx_y_max) THEN
        c_right_top_corner_type = EMPTY_CORNER_WALL_RIGHT   !###
     END IF
  ELSE IF ((Rank_of_master_right.GE.0).AND.(Rank_of_master_above.LT.0)) THEN
     c_right_top_corner_type = FLAT_WALL_ABOVE              !###
! check for the empty corner
     maximal_index=c_indx_x_min
     DO n = 1, c_N_of_local_object_parts_above
        IF ( c_local_object_part(c_index_of_local_object_part_above(n))%iend.GT.maximal_index ) THEN
           maximal_index = c_local_object_part(c_index_of_local_object_part_above(n))%iend
        END IF
     END DO
     IF (maximal_index.LT.c_indx_x_max) THEN
        c_right_top_corner_type = EMPTY_CORNER_WALL_ABOVE   !###
     END IF
  END IF

  IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! we are here if the cluster is periodic in x and connected to itself
! note, the electron dynamics procedure ADVANCE_ELECTRONS does not use corner types in this case
! so we kept the non-negative left/right neighbor master ranks till now to get
! reasonable values of the corner types just in case 
     Rank_of_master_left = -1
     Rank_of_master_right = -1
     PRINT '("ANALYZE_CORNERS :: master process ",i4," is periodically [X] connected to itself, ranks of left/right neighbors set to -1")', Rank_of_process
  END IF

print '("ANALYZE_CORNERS :: proc ",i4," corners LT/RT/LB/RB ",4(1x,i3))', Rank_of_process, c_left_top_corner_type, c_right_top_corner_type, c_left_bottom_corner_type, c_right_bottom_corner_type

END SUBROUTINE ANALYZE_CORNERS

!-------------------------------------------
!
SUBROUTINE SET_COMMUNICATIONS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER ALLOC_ERR
  INTEGER, ALLOCATABLE :: horizontal_rank_in_cluster_level(:)
  INTEGER, ALLOCATABLE :: horizontal_rank_in_cluster_level_left(:)
  INTEGER, ALLOCATABLE :: horizontal_rank_in_cluster_level_right(:)
  INTEGER, ALLOCATABLE :: horizontal_rank_in_cluster_level_above(:)
  INTEGER, ALLOCATABLE :: horizontal_rank_in_cluster_level_below(:)

  INTEGER ibufer(8)

  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, particle_master, cluster_rank_key, COMM_CLUSTER, ierr)
  CALL MPI_COMM_RANK(COMM_CLUSTER, Rank_cluster, ierr)
  CALL MPI_COMM_SIZE(COMM_CLUSTER, N_processes_cluster, ierr)
! COMM_CLUSTER will be used for communications between levels (processes) within a single cluster
! note that N_processes_cluster >= 1 always (1 is the cluster has a master only, no additional processes assigned)

  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, Rank_cluster, particle_master, COMM_HORIZONTAL, ierr)
  CALL MPI_COMM_RANK(COMM_HORIZONTAL, Rank_horizontal, ierr)
  CALL MPI_COMM_SIZE(COMM_HORIZONTAL, N_processes_horizontal, ierr)
! COMM_HORIZONTAL will be used for communications between processes from different clusters which are at the same level

! fool-proof check
  IF (Rank_of_process.EQ.Rank_of_bottom_left_cluster_master) THEN
     IF (Rank_horizontal.NE.0) THEN
        PRINT '("Error-000 in SET_COMMUNICATIONS ",3(2x,i4))', Rank_of_process, Rank_of_bottom_left_cluster_master, Rank_horizontal
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF
  END IF

! collect values of Rank_horizontal from all levels (processes) in the master process
  ALLOCATE(horizontal_rank_in_cluster_level(1:N_processes_cluster), STAT=ALLOC_ERR)
  CALL MPI_GATHER(Rank_horizontal, 1, MPI_INTEGER, horizontal_rank_in_cluster_level, 1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

! default assumption is that a cluster has no neighbors
  N_processes_cluster_left = -1
  N_processes_cluster_right = -1
  N_processes_cluster_above = -1
  N_processes_cluster_below = -1

  Rank_horizontal_left=-1
  Rank_horizontal_right=-1
  Rank_horizontal_above=-1
  Rank_horizontal_below=-1

! send the collected values to neighbor clusters, use the main communicator
  
  IF (cluster_rank_key.EQ.0) THEN
     IF (WHITE_CLUSTER) THEN  
! "white processes"

        IF (Rank_of_master_right.GE.0) THEN
! ##  1 ## send right the number of levels ----------------------------------
           CALL MPI_SEND(N_processes_cluster, 1, MPI_INTEGER, Rank_of_master_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! ##  2 ## send right the ranks
           CALL MPI_SEND(horizontal_rank_in_cluster_level, N_processes_cluster, MPI_INTEGER, Rank_of_master_right, Rank_of_process+SHIFT1, MPI_COMM_WORLD, request, ierr)     
        END IF

        IF (Rank_of_master_left.GE.0) THEN
! ##  3 ## send left the number of levels
           CALL MPI_SEND(N_processes_cluster, 1, MPI_INTEGER, Rank_of_master_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! ##  4 ## send left the ranks
           CALL MPI_SEND(horizontal_rank_in_cluster_level, N_processes_cluster, MPI_INTEGER, Rank_of_master_left, Rank_of_process+SHIFT1, MPI_COMM_WORLD, request, ierr)     
        END IF

        IF (Rank_of_master_left.GE.0) THEN
! ##  5 ## receive from left the number of levels ===========================
           CALL MPI_RECV(N_processes_cluster_left, 1, MPI_INTEGER, Rank_of_master_left, Rank_of_master_left, MPI_COMM_WORLD, stattus, ierr)
! ##  6 ## receive from left the ranks
           ALLOCATE(horizontal_rank_in_cluster_level_left(1:N_processes_cluster_left), STAT=ALLOC_ERR)
           CALL MPI_RECV(horizontal_rank_in_cluster_level_left, N_processes_cluster_left, MPI_INTEGER, Rank_of_master_left, Rank_of_master_left+SHIFT1, MPI_COMM_WORLD, stattus, ierr)     
        END IF

        IF (Rank_of_master_right.GE.0) THEN
! ##  7 ## receive from right the number of levels
           CALL MPI_RECV(N_processes_cluster_right, 1, MPI_INTEGER, Rank_of_master_right, Rank_of_master_right, MPI_COMM_WORLD, stattus, ierr)
! ##  8 ## receive from right the ranks
           ALLOCATE(horizontal_rank_in_cluster_level_right(1:N_processes_cluster_right), STAT=ALLOC_ERR)
           CALL MPI_RECV(horizontal_rank_in_cluster_level_right, N_processes_cluster_right, MPI_INTEGER, Rank_of_master_right, Rank_of_master_right+SHIFT1, MPI_COMM_WORLD, stattus, ierr)     
        END IF

        IF (Rank_of_master_above.GE.0) THEN
! ##  9 ## send up the number of levels -------------------------------------
           CALL MPI_SEND(N_processes_cluster, 1, MPI_INTEGER, Rank_of_master_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! ## 10 ## send up the ranks
           CALL MPI_SEND(horizontal_rank_in_cluster_level, N_processes_cluster, MPI_INTEGER, Rank_of_master_above, Rank_of_process+SHIFT1, MPI_COMM_WORLD, request, ierr)     
        END IF

        IF (Rank_of_master_below.GE.0) THEN
! ## 11 ## send down the number of levels
           CALL MPI_SEND(N_processes_cluster, 1, MPI_INTEGER, Rank_of_master_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! ## 12 ## send down the ranks
           CALL MPI_SEND(horizontal_rank_in_cluster_level, N_processes_cluster, MPI_INTEGER, Rank_of_master_below, Rank_of_process+SHIFT1, MPI_COMM_WORLD, request, ierr)     
        END IF

        IF (Rank_of_master_below.GE.0) THEN
! ## 13 ## receive from below the number of levels ==========================
           CALL MPI_RECV(N_processes_cluster_below, 1, MPI_INTEGER, Rank_of_master_below, Rank_of_master_below, MPI_COMM_WORLD, stattus, ierr)
! ## 14 ## receive from below the ranks
           ALLOCATE(horizontal_rank_in_cluster_level_below(1:N_processes_cluster_below), STAT=ALLOC_ERR)
           CALL MPI_RECV(horizontal_rank_in_cluster_level_below, N_processes_cluster_below, MPI_INTEGER, Rank_of_master_below, Rank_of_master_below+SHIFT1, MPI_COMM_WORLD, stattus, ierr)     
        END IF

        IF (Rank_of_master_above.GE.0) THEN
! ## 15 ## receive from above the number of levels ==========================
           CALL MPI_RECV(N_processes_cluster_above, 1, MPI_INTEGER, Rank_of_master_above, Rank_of_master_above, MPI_COMM_WORLD, stattus, ierr)
! ## 16 ## receive from above the ranks
           ALLOCATE(horizontal_rank_in_cluster_level_above(1:N_processes_cluster_above), STAT=ALLOC_ERR)
           CALL MPI_RECV(horizontal_rank_in_cluster_level_above, N_processes_cluster_above, MPI_INTEGER, Rank_of_master_above, Rank_of_master_above+SHIFT1, MPI_COMM_WORLD, stattus, ierr)     
        END IF

     ELSE
! "black" processes

        IF (Rank_of_master_left.GE.0) THEN
! ##  1 ## receive from left the number of levels ===========================
           CALL MPI_RECV(N_processes_cluster_left, 1, MPI_INTEGER, Rank_of_master_left, Rank_of_master_left, MPI_COMM_WORLD, stattus, ierr)
! ##  2 ## receive from left the ranks
           ALLOCATE(horizontal_rank_in_cluster_level_left(1:N_processes_cluster_left), STAT=ALLOC_ERR)
           CALL MPI_RECV(horizontal_rank_in_cluster_level_left, N_processes_cluster_left, MPI_INTEGER, Rank_of_master_left, Rank_of_master_left+SHIFT1, MPI_COMM_WORLD, stattus, ierr)     
        END IF

        IF (Rank_of_master_right.GE.0) THEN
! ##  3 ## receive from right the number of levels
           CALL MPI_RECV(N_processes_cluster_right, 1, MPI_INTEGER, Rank_of_master_right, Rank_of_master_right, MPI_COMM_WORLD, stattus, ierr)
! ##  4 ## receive from right the ranks
           ALLOCATE(horizontal_rank_in_cluster_level_right(1:N_processes_cluster_right), STAT=ALLOC_ERR)
           CALL MPI_RECV(horizontal_rank_in_cluster_level_right, N_processes_cluster_right, MPI_INTEGER, Rank_of_master_right, Rank_of_master_right+SHIFT1, MPI_COMM_WORLD, stattus, ierr)     
        END IF

        IF (Rank_of_master_right.GE.0) THEN
! ##  5 ## send right the number of levels ----------------------------------
           CALL MPI_SEND(N_processes_cluster, 1, MPI_INTEGER, Rank_of_master_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! ##  6 ## send right the ranks
           CALL MPI_SEND(horizontal_rank_in_cluster_level, N_processes_cluster, MPI_INTEGER, Rank_of_master_right, Rank_of_process+SHIFT1, MPI_COMM_WORLD, request, ierr)     
        END IF

        IF (Rank_of_master_left.GE.0) THEN
! ##  7 ## send left the number of levels
           CALL MPI_SEND(N_processes_cluster, 1, MPI_INTEGER, Rank_of_master_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! ##  8 ## send left the ranks
           CALL MPI_SEND(horizontal_rank_in_cluster_level, N_processes_cluster, MPI_INTEGER, Rank_of_master_left, Rank_of_process+SHIFT1, MPI_COMM_WORLD, request, ierr)     
        END IF

        IF (Rank_of_master_below.GE.0) THEN
! ##  9 ## receive from below the number of levels ==========================
           CALL MPI_RECV(N_processes_cluster_below, 1, MPI_INTEGER, Rank_of_master_below, Rank_of_master_below, MPI_COMM_WORLD, stattus, ierr)
! ## 10 ## receive from below the ranks
           ALLOCATE(horizontal_rank_in_cluster_level_below(1:N_processes_cluster_below), STAT=ALLOC_ERR)
           CALL MPI_RECV(horizontal_rank_in_cluster_level_below, N_processes_cluster_below, MPI_INTEGER, Rank_of_master_below, Rank_of_master_below+SHIFT1, MPI_COMM_WORLD, stattus, ierr)     
        END IF

        IF (Rank_of_master_above.GE.0) THEN
! ## 11 ## receive from above the number of levels ==========================
           CALL MPI_RECV(N_processes_cluster_above, 1, MPI_INTEGER, Rank_of_master_above, Rank_of_master_above, MPI_COMM_WORLD, stattus, ierr)
! ## 12 ## receive from above the ranks
           ALLOCATE(horizontal_rank_in_cluster_level_above(1:N_processes_cluster_above), STAT=ALLOC_ERR)
           CALL MPI_RECV(horizontal_rank_in_cluster_level_above, N_processes_cluster_above, MPI_INTEGER, Rank_of_master_above, Rank_of_master_above+SHIFT1, MPI_COMM_WORLD, stattus, ierr)     
        END IF

        IF (Rank_of_master_above.GE.0) THEN
! ## 13 ## send up the number of levels -------------------------------------
           CALL MPI_SEND(N_processes_cluster, 1, MPI_INTEGER, Rank_of_master_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! ## 14 ## send up the ranks
           CALL MPI_SEND(horizontal_rank_in_cluster_level, N_processes_cluster, MPI_INTEGER, Rank_of_master_above, Rank_of_process+SHIFT1, MPI_COMM_WORLD, request, ierr)     
        END IF

        IF (Rank_of_master_below.GE.0) THEN
! ## 15 ## send down the number of levels
           CALL MPI_SEND(N_processes_cluster, 1, MPI_INTEGER, Rank_of_master_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! ## 16 ## send down the ranks
           CALL MPI_SEND(horizontal_rank_in_cluster_level, N_processes_cluster, MPI_INTEGER, Rank_of_master_below, Rank_of_process+SHIFT1, MPI_COMM_WORLD, request, ierr)     
        END IF

     END IF
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! distribute numbers of levels (processes) in neighbor clusters 
! and ranks of neighbor masters among the levels of this cluster
  IF (Rank_cluster.EQ.0) THEN
     ibufer(1) = N_processes_cluster_left
     ibufer(2) = N_processes_cluster_right
     ibufer(3) = N_processes_cluster_above
     ibufer(4) = N_processes_cluster_below
     ibufer(5) = Rank_of_master_left
     ibufer(6) = Rank_of_master_right
     ibufer(7) = Rank_of_master_above
     ibufer(8) = Rank_of_master_below
     CALL MPI_BCAST(ibufer, 8, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
  ELSE
     CALL MPI_BCAST(ibufer, 8, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
     N_processes_cluster_left  = ibufer(1)
     N_processes_cluster_right = ibufer(2)
     N_processes_cluster_above = ibufer(3)
     N_processes_cluster_below = ibufer(4)
     Rank_of_master_left  = ibufer(5)
     Rank_of_master_right = ibufer(6)
     Rank_of_master_above = ibufer(7)
     Rank_of_master_below = ibufer(8)

     ALLOCATE(horizontal_rank_in_cluster_level_left( 1:N_processes_cluster_left), STAT=ALLOC_ERR)    
     ALLOCATE(horizontal_rank_in_cluster_level_right(1:N_processes_cluster_right), STAT=ALLOC_ERR)
     ALLOCATE(horizontal_rank_in_cluster_level_above(1:N_processes_cluster_above), STAT=ALLOC_ERR)
     ALLOCATE(horizontal_rank_in_cluster_level_below(1:N_processes_cluster_below), STAT=ALLOC_ERR)
  END IF

! distribute vectors with intercommunication ranks of levels (processes) from neighbor clusters among the levels of this cluster
! since all processes know the values of N_processes_cluster_* and have arrays horizontal_rank_in_cluster_level_* ready, use MPI_BCAST
! the IFs below help exclude boundaries
  IF (N_processes_cluster_left.GT.0)  CALL MPI_BCAST(horizontal_rank_in_cluster_level_left,  N_processes_cluster_left,  MPI_INTEGER, 0, COMM_CLUSTER, ierr)
  IF (N_processes_cluster_right.GT.0) CALL MPI_BCAST(horizontal_rank_in_cluster_level_right, N_processes_cluster_right, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
  IF (N_processes_cluster_above.GT.0) CALL MPI_BCAST(horizontal_rank_in_cluster_level_above, N_processes_cluster_above, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
  IF (N_processes_cluster_below.GT.0) CALL MPI_BCAST(horizontal_rank_in_cluster_level_below, N_processes_cluster_below, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

! each level must identify whether each neighbor cluster has the same levels as itself
! note that a level, which is also an index in arrays horizontal_rank_in_cluster_level_*, is defined as Rank_cluster+1, i.e. level >=1
  IF (Rank_cluster.LT.N_processes_cluster_left)  Rank_horizontal_left  = horizontal_rank_in_cluster_level_left( Rank_cluster+1)
  IF (Rank_cluster.LT.N_processes_cluster_right) Rank_horizontal_right = horizontal_rank_in_cluster_level_right(Rank_cluster+1)
  IF (Rank_cluster.LT.N_processes_cluster_above) Rank_horizontal_above = horizontal_rank_in_cluster_level_above(Rank_cluster+1)
  IF (Rank_cluster.LT.N_processes_cluster_below) Rank_horizontal_below = horizontal_rank_in_cluster_level_below(Rank_cluster+1)

! arrays with horizontal levels are not necessary any more, deallocate them
  DEALLOCATE(horizontal_rank_in_cluster_level, STAT=ALLOC_ERR)
  IF (ALLOCATED(horizontal_rank_in_cluster_level_left))  DEALLOCATE(horizontal_rank_in_cluster_level_left,  STAT=ALLOC_ERR)
  IF (ALLOCATED(horizontal_rank_in_cluster_level_right)) DEALLOCATE(horizontal_rank_in_cluster_level_right, STAT=ALLOC_ERR)
  IF (ALLOCATED(horizontal_rank_in_cluster_level_above)) DEALLOCATE(horizontal_rank_in_cluster_level_above, STAT=ALLOC_ERR)
  IF (ALLOCATED(horizontal_rank_in_cluster_level_below)) DEALLOCATE(horizontal_rank_in_cluster_level_below, STAT=ALLOC_ERR)

! at this stage each level, including the master process, knows the following:
!
!  COMM_HORIZONTAL
!
!  N_processes_cluster_left
!  N_processes_cluster_right
!  N_processes_cluster_above
!  N_processes_cluster_below
!
!  Rank_horizontal_left
!  Rank_horizontal_right
!  Rank_horizontal_above
!  Rank_horizontal_below
!
!  Rank_of_master_left    ! keeping in memory both Rank_horizontal_* and Rank_of_master_* is not a redundancy
!  Rank_of_master_right   ! while Rank_horizontal_* are used for particle exchange
!  Rank_of_master_above   ! Rank_of_master_* are used in particle advance procedure to distinguish between
!  Rank_of_master_below   ! wall/no wall situations and properly place escaping particles 
!                         ! even when no horizontal neighbors at the level of th eprocess are available
  
print '("Process ",i4," particle_master ",i4," Rank_cluster ",i4," Rank_horizontal ",i4," neighbor horizontal ranks [L/R/A/B] ",4(2x,i4))', &
     & Rank_of_process, particle_master, Rank_cluster, Rank_horizontal, Rank_horizontal_left, Rank_horizontal_right, Rank_horizontal_above, Rank_horizontal_below

END SUBROUTINE SET_COMMUNICATIONS

!--------------------------------------------------
! this subroutine works only with master processes
!
SUBROUTINE COLLECT_MASTERS_INFO

  USE ParallelOperationValues
  USE LoadBalancing

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER myibufer(2)
  INTEGER ALLOC_ERR
  INTEGER i

  IF (cluster_rank_key.NE.0) RETURN

  ALLOCATE(ibufer(0:2*N_processes_horizontal-1), STAT=ALLOC_ERR)

  myibufer(1) = particle_master
  myibufer(2) = N_processes_cluster

  CALL MPI_GATHER(myibufer, 2, MPI_INTEGER, ibufer, 2, MPI_INTEGER, 0, COMM_HORIZONTAL, ierr)
  
  IF (Rank_horizontal.EQ.0) THEN
     ALLOCATE(cluster(0:N_processes_horizontal-1), STAT=ALLOC_ERR)
     DO i = 0, N_processes_horizontal-1
        cluster(i)%particle_master = ibufer(2*i)
        cluster(i)%N_processes     = ibufer(2*i+1)
     END DO
  END IF

  DEALLOCATE(ibufer, STAT=ALLOC_ERR)

END SUBROUTINE COLLECT_MASTERS_INFO

!--------------------------------------------
!
SUBROUTINE DISTRIBUTE_CLUSTER_PARAMETERS

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues
  USE Diagnostics, ONLY : N_of_probes_cluster, List_of_probes_cluster !, probe_Ni_cluster
  USE IonParticles, ONLY : N_spec
  USE SetupValues

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER bufer_length
  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  INTEGER pos, n

  IF (Rank_cluster.EQ.0) THEN

     bufer_length = 14 + &
                  & 4 + &
                  & 2 + &
                  & 2 + &
                  & 6 * c_N_of_local_object_parts + &
                  & c_N_of_local_object_parts_left + &
                  & c_N_of_local_object_parts_right + &
                  & c_N_of_local_object_parts_above + &
                  & c_N_of_local_object_parts_below + &
                  & 1

     CALL MPI_BCAST(bufer_length, 1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

     ALLOCATE(ibufer(1:bufer_length), STAT=ALLOC_ERR)

     ibufer(1) = c_indx_x_min
     ibufer(2) = c_indx_x_max
     ibufer(3) = c_indx_y_min
     ibufer(4) = c_indx_y_max

     ibufer(5) = c_left_bottom_corner_type
     ibufer(6) = c_left_top_corner_type
     ibufer(7) = c_right_bottom_corner_type
     ibufer(8) = c_right_top_corner_type

     ibufer(9) = c_N_of_local_object_parts

     ibufer(10) = c_N_of_local_object_parts_left
     ibufer(11) = c_N_of_local_object_parts_right
     ibufer(12) = c_N_of_local_object_parts_above
     ibufer(13) = c_N_of_local_object_parts_below

     IF (WHITE_CLUSTER) THEN
        ibufer(14) = 1
     ELSE
        ibufer(14) = 0
     END IF

     ibufer(15) = 0
     ibufer(16) = 0
     ibufer(17) = 0
     ibufer(18) = 0

     IF (periodic_boundary_X_left) THEN  !### periodic_boundary_X_left and symmetry_plane_X_left cannot be .TRUE. at the same time, but both can be .FALSE.
        ibufer(15) = 1
     ELSE IF (symmetry_plane_X_left) THEN
        ibufer(15) = -1
     END IF

     IF (periodic_boundary_X_right) ibufer(16) = 1
     IF (periodic_boundary_Y_below) ibufer(17) = 1
     IF (periodic_boundary_Y_above) ibufer(18) = 1

     ibufer(19) = i_period_x
     ibufer(20) = i_period_y

     IF (ht_use_ionization_source) THEN
        ibufer(21) = c_j_ion_source_1
        ibufer(22) = c_j_ion_source_2
     ELSE
        ibufer(21) = 0
        ibufer(22) = 0        
     END IF

     pos=22
     DO n = 1, c_N_of_local_object_parts 
        ibufer(pos+1) = c_local_object_part(n)%object_number
        ibufer(pos+2) = c_local_object_part(n)%segment_number
        ibufer(pos+3) = c_local_object_part(n)%istart
        ibufer(pos+4) = c_local_object_part(n)%jstart
        ibufer(pos+5) = c_local_object_part(n)%iend
        ibufer(pos+6) = c_local_object_part(n)%jend
        pos=pos+6
     END DO

     DO n = 1, c_N_of_local_object_parts_left
        ibufer(pos+n) = c_index_of_local_object_part_left(n)
     END DO

     pos = pos + c_N_of_local_object_parts_left

     DO n = 1, c_N_of_local_object_parts_right
        ibufer(pos+n) = c_index_of_local_object_part_right(n)
     END DO

     pos = pos + c_N_of_local_object_parts_right

     DO n = 1, c_N_of_local_object_parts_above
        ibufer(pos+n) = c_index_of_local_object_part_above(n)
     END DO

     pos = pos + c_N_of_local_object_parts_above
     
     DO n = 1, c_N_of_local_object_parts_below
        ibufer(pos+n) = c_index_of_local_object_part_below(n)
     END DO

     pos = pos + c_N_of_local_object_parts_below
     ibufer(pos+1) = N_of_probes_cluster
     
!if ((pos+1).ne.bufer_length) 
!print '(">>>",4(2x,i4))', Rank_of_process, bufer_length, pos+1, N_of_probes_cluster

     CALL MPI_BCAST(ibufer, bufer_length, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

     IF (N_of_probes_cluster.GT.0) THEN

        CALL MPI_BCAST(List_of_probes_cluster, N_of_probes_cluster, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

     END IF
    
     IF (ht_use_ionization_source) THEN
        CALL MPI_BCAST(yi(0:c_R_max), c_R_max+1, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)        
     END IF

  ELSE

     CALL MPI_BCAST(bufer_length, 1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
     
     ALLOCATE(ibufer(1:bufer_length), STAT=ALLOC_ERR)

     CALL MPI_BCAST(ibufer, bufer_length, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

     c_indx_x_min = ibufer(1)
     c_indx_x_max = ibufer(2)
     c_indx_y_min = ibufer(3)
     c_indx_y_max = ibufer(4)

     IF (ALLOCATED(EX))     DEALLOCATE(EX,     STAT=ALLOC_ERR)
     IF (ALLOCATED(EY))     DEALLOCATE(EY,     STAT=ALLOC_ERR)
     IF (ALLOCATED(acc_EX)) DEALLOCATE(acc_EX, STAT=ALLOC_ERR)
     IF (ALLOCATED(acc_EY)) DEALLOCATE(acc_EY, STAT=ALLOC_ERR)

     ALLOCATE(EX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(EY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(acc_EX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE(acc_EY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     
!#########
!     EX=0.0_8 !-36.491_8 / E_scale_Vm !0.0_8
!     EY=0.0_8 !-1.0_8 / E_scale_Vm
     acc_EX=0.0_8 !-36.491_8 / E_scale_Vm !0.0_8
     acc_EY=0.0_8 !-1.0_8 / E_scale_Vm
!#########

     c_X_area_min = DBLE(c_indx_x_min)
     c_X_area_max = DBLE(c_indx_x_max)
     c_Y_area_min = DBLE(c_indx_y_min)
     c_Y_area_max = DBLE(c_indx_y_max)        

     c_left_bottom_corner_type  = ibufer(5)
     c_left_top_corner_type     = ibufer(6)
     c_right_bottom_corner_type = ibufer(7)
     c_right_top_corner_type    = ibufer(8)

     c_N_of_local_object_parts = ibufer(9)

     c_N_of_local_object_parts_left  = ibufer(10)
     c_N_of_local_object_parts_right = ibufer(11)
     c_N_of_local_object_parts_above = ibufer(12)
     c_N_of_local_object_parts_below = ibufer(13)
     
     IF (ibufer(14).EQ.1) THEN
        WHITE_CLUSTER = .TRUE.
     ELSE
        WHITE_CLUSTER = .FALSE.
     END IF

     periodic_boundary_X_left  = .FALSE.
     periodic_boundary_X_right = .FALSE.
     periodic_boundary_Y_below = .FALSE.
     periodic_boundary_Y_above = .FALSE.

     symmetry_plane_X_left = .FALSE.

     IF (ibufer(15).EQ.1) THEN
        periodic_boundary_X_left  = .TRUE.
     ELSE IF (ibufer(15).EQ.-1) THEN 
        symmetry_plane_X_left = .TRUE.
     END IF

     IF (ibufer(16).EQ.1) periodic_boundary_X_right = .TRUE.
     IF (ibufer(17).EQ.1) periodic_boundary_Y_below = .TRUE.
     IF (ibufer(18).EQ.1) periodic_boundary_Y_above = .TRUE.

     i_period_x = ibufer(19)
     i_period_y = ibufer(20)

     L_period_X = DBLE(i_period_x)
     L_period_Y = DBLE(i_period_y)
    
     IF (ht_use_ionization_source) THEN
        c_j_ion_source_1 = ibufer(21)
        c_j_ion_source_2 = ibufer(22)
     ELSE
     END IF

     pos=22
     DO n = 1, c_N_of_local_object_parts 
        c_local_object_part(n)%object_number  = ibufer(pos+1)
        c_local_object_part(n)%segment_number = ibufer(pos+2)
        c_local_object_part(n)%istart         = ibufer(pos+3)
        c_local_object_part(n)%jstart         = ibufer(pos+4)
        c_local_object_part(n)%iend           = ibufer(pos+5)
        c_local_object_part(n)%jend           = ibufer(pos+6)

! for dielectric walls
        IF (whole_object(c_local_object_part(n)%object_number)%object_type.EQ.DIELECTRIC) THEN
! remove previously allocated array of surface charge if necessary
           IF (ALLOCATED(c_local_object_part(n)%surface_charge)) DEALLOCATE(c_local_object_part(n)%surface_charge, STAT = ALLOC_ERR)
           IF ((c_local_object_part(n)%istart.EQ.c_local_object_part(n)%iend).AND.(c_local_object_part(n)%jstart.LT.c_local_object_part(n)%jend)) THEN
! this part is vertical (along the Y-direction)
              ALLOCATE(c_local_object_part(n)%surface_charge(c_local_object_part(n)%jstart:c_local_object_part(n)%jend), STAT = ALLOC_ERR)
              c_local_object_part(n)%surface_charge(c_local_object_part(n)%jstart:c_local_object_part(n)%jend) = 0.0_8
           ELSE IF ((c_local_object_part(n)%jstart.EQ.c_local_object_part(n)%jend).AND.(c_local_object_part(n)%istart.LT.c_local_object_part(n)%iend)) THEN
! this part is horizontal (along the X-direction)
              ALLOCATE(c_local_object_part(n)%surface_charge(c_local_object_part(n)%istart:c_local_object_part(n)%iend), STAT = ALLOC_ERR)
              c_local_object_part(n)%surface_charge(c_local_object_part(n)%istart:c_local_object_part(n)%iend) = 0.0_8
           ELSE
! error
              PRINT '("Proc ",i4," Error in DISTRIBUTE_CLUSTER_PARAMETERS")', Rank_of_process
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF
        END IF

        pos=pos+6
     END DO

     DO n = 1, c_N_of_local_object_parts_left
        c_index_of_local_object_part_left(n) = ibufer(pos+n)
     END DO

     pos = pos + c_N_of_local_object_parts_left

     DO n = 1, c_N_of_local_object_parts_right
        c_index_of_local_object_part_right(n) = ibufer(pos+n)
     END DO

     pos = pos + c_N_of_local_object_parts_right

     DO n = 1, c_N_of_local_object_parts_above
        c_index_of_local_object_part_above(n) = ibufer(pos+n)
     END DO

     pos = pos + c_N_of_local_object_parts_above

     DO n = 1, c_N_of_local_object_parts_below
        c_index_of_local_object_part_below(n) = ibufer(pos+n)
     END DO

     pos = pos + c_N_of_local_object_parts_below

     N_of_probes_cluster = ibufer(pos+1)

     IF (ALLOCATED(List_of_probes_cluster)) DEALLOCATE(List_of_probes_cluster, STAT=ALLOC_ERR)

     IF (N_of_probes_cluster.GT.0) THEN
        ALLOCATE(List_of_probes_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        CALL MPI_BCAST(List_of_probes_cluster, N_of_probes_cluster, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
     END IF

     IF (ht_use_ionization_source) THEN
        CALL MPI_BCAST(yi(0:c_R_max), c_R_max+1, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)        
     END IF

  END IF

  DEALLOCATE(ibufer, STAT=ALLOC_ERR)

END SUBROUTINE DISTRIBUTE_CLUSTER_PARAMETERS

!--------------------------------------------
! initially, particles will be given to the master processes only
! later (at the very first timestep) a standard procedure of global and inner load balancing will be applied
! to ensure that the system is balanced from the beginning
!
SUBROUTINE DISTRIBUTE_PARTICLES

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE IonParticles
  USE ClusterAndItsBoundaries

  USE rng_wrapper

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(8) x_limit_left, x_limit_right
  REAL(8) y_limit_bot, y_limit_top
  INTEGER n_limit_x, n_limit_y

  REAL(8) myEy_Vm, myEz_Vm
  REAL(8) myBx_T, myBy_T, myBz_T, myB2

  INTEGER ALLOC_ERR
  REAL(8) factor_convert
  REAL(8) vx_drift, vy_drift, vz_drift, v

  INTEGER k, n, s, pos
  INTEGER sum_Ni

! functions
  REAL(8) Bx, By, Bz, Ez

  IF (Rank_cluster.EQ.0) THEN

! initialize particles

! fixmeplease
     myEy_Vm = 0.0_8 !!!### turned it off, was like this ((-whole_object(2)%phi * F_scale_V) / (global_maximal_j * delta_x_m))

!! electrons
!     N_electrons = INT(DBLE(N_of_particles_cell) * (init_Ne_m3 / N_plasma_m3) * N_grid_block_x * N_grid_block_y * cluster_N_blocks)
!!NEW ######################################################################
!! in clusters adjacent to the upper wall along y (maximal y) one needs to add one row of cells
!     IF (Rank_of_master_above.LT.0) THEN
!        N_electrons = INT(DBLE(N_of_particles_cell) * (init_Ne_m3 / N_plasma_m3) * (N_grid_block_x * N_grid_block_y * cluster_N_blocks + N_grid_block_x * cluster_N_blocks_x))
!     END IF

     x_limit_left = c_X_area_min
     x_limit_right = c_X_area_max - 1.0_8
     n_limit_x = c_indx_x_max - c_indx_x_min - 1
     IF (Rank_of_master_right.LT.0) THEN
        x_limit_right = c_X_area_max
        n_limit_x = c_indx_x_max - c_indx_x_min
     END IF

     y_limit_bot = c_Y_area_min
     y_limit_top = c_Y_area_max - 1.0_8
     n_limit_y = c_indx_y_max - c_indx_y_min - 1
     IF (Rank_of_master_above.LT.0) THEN
        y_limit_top = c_Y_area_max
        n_limit_y = c_indx_y_max - c_indx_y_min
     END IF

     N_electrons =  INT( DBLE(N_of_particles_cell) * (init_Ne_m3 / N_plasma_m3) * n_limit_x * n_limit_y )

     PRINT '(2x,"Master process ",i3," : ",i8," electron macroparticles")', Rank_of_process, N_electrons

     max_N_electrons = N_electrons+50

     ALLOCATE(electron(1:max_N_electrons), STAT=ALLOC_ERR)

     factor_convert = SQRT(init_Te_eV / T_e_eV) / N_max_vel

     DO k = 1, N_electrons

        electron(k)%X = MAX(c_X_area_min, MIN(c_X_area_max, x_limit_left + (x_limit_right - x_limit_left) * well_random_number()))

        electron(k)%Y = MAX(c_Y_area_min, MIN(c_Y_area_max, y_limit_bot + (y_limit_top - y_limit_bot) * well_random_number()))

!        electron(k)%X = MIN(c_X_area_max-1.0_8, c_X_area_min + (c_X_area_max-1.0_8-c_X_area_min) * well_random_number())
!!NEW #######################################################################
!        IF (Rank_of_master_above.GE.0) THEN
!! old version valid for clusters not adjacent to the upper wall
!           electron(k)%Y = MIN(c_Y_area_max-1.0_8, c_Y_area_min + (c_Y_area_max-1.0_8-c_Y_area_min) * well_random_number())
!        ELSE
!! in clusters adjacent to the upper wall expand range of y by one cell
!           electron(k)%Y = MIN(c_Y_area_max, c_Y_area_min + (c_Y_area_max-c_Y_area_min) * well_random_number())
!        END IF

        myBx_T = B_scale_T * Bx(electron(k)%X, electron(k)%Y)
        myBy_T = B_scale_T * By(electron(k)%X, electron(k)%Y)
        myBz_T = B_scale_T * Bz(electron(k)%X, electron(k)%Y)

        myB2 = (myBx_T**2 + myBy_T**2 + myBz_T**2) * V_scale_ms

        IF (myB2.GT.0.0_8) THEN
! account possible ExB drifts for electrons
           myEz_Vm = E_scale_Vm * Ez(electron(k)%X, electron(k)%Y)
           vx_drift = ( myEy_Vm * myBz_T - myEz_Vm * myBy_T) / myB2
           vy_drift = ( myEz_Vm * myBx_T) / myB2
           vz_drift = (-myEy_Vm * myBx_T) / myB2
        ELSE
           vx_drift = 0.0_8
           vy_drift = 0.0_8
           vz_drift = 0.0_8
        END IF

        CALL GetMaxwellVelocity(v)
        electron(k)%VX = v * factor_convert + vx_drift

        CALL GetMaxwellVelocity(v)
        electron(k)%VY = v * factor_convert + vy_drift

        CALL GetMaxwellVelocity(v)
        electron(k)%VZ = v * factor_convert + vz_drift
        electron(k)%tag = 0

     END DO

! remove all electrons which were placed inside inner objects
     k=0
     DO WHILE (k.LT.N_electrons)
        IF (N_electrons.EQ.0) EXIT
        k = k+1
        DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (electron(k)%X.LE.whole_object(n)%Xmin) CYCLE
           IF (electron(k)%X.GE.whole_object(n)%Xmax) CYCLE
           IF (electron(k)%Y.LE.whole_object(n)%Ymin) CYCLE
           IF (electron(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! particle inside an inner object
           CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
           EXIT
        END DO
     END DO

! ions
     sum_Ni = 0
     DO s = 1, N_spec-1
        N_ions(s) = INT(DBLE(N_electrons) * init_NiNe(s))
        sum_Ni = sum_Ni + N_ions(s)
        PRINT '(2x,"Master process ",i3," : ",i8," macroparticles of ion species ",i3)', Rank_of_process, N_ions(s), s
        max_N_ions(s) = N_ions(s)+50
     END DO

! process the last ion species separately to make sure that the total number of electrons and ions is the same
     s = N_spec
     N_ions(s) = N_electrons - sum_Ni
     PRINT '(2x,"Master process ",i3," : ",i8," macroparticles of ion species ",i3)', Rank_of_process, N_ions(s), s
     max_N_ions(s) = N_ions(s)+50

     ALLOCATE(ion(1:N_spec), STAT = ALLOC_ERR)
     DO s = 1, N_spec
        ALLOCATE(ion(s)%part(1:max_N_ions(s)), STAT = ALLOC_ERR)
     END DO

     pos = 0
     DO s = 1, N_spec
! since electrons are distributed randomly and are not ordered, we simply use their coordinates

        factor_convert = SQRT(init_Ti_eV(s) / T_e_eV) / (N_max_vel * SQRT(Ms(s)))

        DO k = 1, N_ions(s)
           pos = pos+1
           ion(s)%part(k)%X = electron(pos)%X
           ion(s)%part(k)%Y = electron(pos)%Y

           CALL GetMaxwellVelocity(v)
           ion(s)%part(k)%VX = v * factor_convert

           CALL GetMaxwellVelocity(v)
           ion(s)%part(k)%VY = v * factor_convert

           ion(s)%part(k)%VZ = 0.0_8
           ion(s)%part(k)%tag = 0
        END DO
     END DO

  ELSE

! electrons
     N_electrons = 0
     max_N_electrons = 50
     ALLOCATE(electron(1:max_N_electrons), STAT=ALLOC_ERR)

! ions
     DO s = 1, N_spec
! ###### so far this works for a single ion species only        
        N_ions(s) = 0
        max_N_ions(s) = N_ions(s)+50
     END DO

     ALLOCATE(ion(1:N_spec), STAT = ALLOC_ERR)
     DO s = 1, N_spec
        ALLOCATE(ion(s)%part(1:max_N_ions(s)), STAT = ALLOC_ERR)
     END DO

  END IF

END SUBROUTINE DISTRIBUTE_PARTICLES

