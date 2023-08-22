!---------------------------------------
!
SUBROUTINE ADVANCE_IONS_PLUS

  USE ParallelOperationValues, ONLY : cluster_rank_key
  USE CurrentProblemValues, ONLY : T_cntr, N_of_boundary_and_inner_objects, whole_object, ion_colls_with_bo
  USE Snapshots, ONLY : cs_N, cs_VX, cs_VY, cs_VZ, cs_WX, cs_WY, cs_WZ, cs_VXVY, cs_VXVZ, cs_VYVZ, cs_QX, cs_QY, cs_QZ
  USE AvgSnapshots, ONLY : current_avgsnap, N_of_all_avgsnaps, avgsnapshot, save_avg_data
  USE ClusterAndItsBoundaries, ONLY : c_indx_x_min, c_indx_x_max, c_indx_y_min, c_indx_y_max
  USE IonParticles, ONLY : N_spec, N_ions_to_send_left, N_ions_to_send_right, N_ions_to_send_above, N_ions_to_send_below
  USE Diagnostics, ONLY : Save_probes_i_data_T_cntr

  use mpi

  IMPLICIT NONE


  INTEGER ierr

  LOGICAL collect_ion_moments_now
  INTEGER ALLOC_ERR

  INTEGER k, s

  collect_ion_moments_now = .FALSE.
  
  IF ((current_avgsnap.GE.1).AND.(current_avgsnap.LE.N_of_all_avgsnaps)) THEN
     IF ( (T_cntr.GE.avgsnapshot(current_avgsnap)%T_cntr_begin).AND. &
        & (T_cntr.LE.avgsnapshot(current_avgsnap)%T_cntr_end) ) THEN

        IF ( save_avg_data(4).OR. &
           & save_avg_data(5).OR. &
           & save_avg_data(6).OR. &
           & save_avg_data(23).OR. &
           & save_avg_data(24).OR. &
           & save_avg_data(25).OR. &
           & save_avg_data(26).OR. &
           & save_avg_data(27).OR. &
           & save_avg_data(28).OR. &
           & save_avg_data(29).OR. &
           & save_avg_data(30).OR. &
           & save_avg_data(31).OR. &
           & save_avg_data(32).OR. &
           & save_avg_data(33).OR. &
           & save_avg_data(34).OR. &
           & save_avg_data(35).OR. &
           & save_avg_data(36).OR. &
           & save_avg_data(37).OR. &
           & save_avg_data(38) ) collect_ion_moments_now = .TRUE.
     END IF
  END IF

  IF (collect_ion_moments_now) THEN

     IF (cluster_rank_key.EQ.0) THEN

! arrays required by procedures collecting moments of electron and ion vdfs

        ALLOCATE(cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

        ALLOCATE(cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

        ALLOCATE(cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

        ALLOCATE(cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

        ALLOCATE(cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ELSE

! the moments of the distribution functions are calculated using MPI_REDUCE which takes sum of values from all particle calculators
! the results are stored in master processes in arrays cs_N, cs_VX, etc
! the particle calculators, in general, don't need these arrays at all, 
! so allocating these arrays in non-master processes is just a waste of memory
! however, the compiler reports an error if the code is compiled with -C flag (check everything)
! this can be avoided if at least some minimal size arrays are allocated in the non-master processes
        ALLOCATE(cs_N(1,1), STAT=ALLOC_ERR)

        ALLOCATE(cs_VX(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_VY(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_VZ(1,1), STAT=ALLOC_ERR)

        ALLOCATE(cs_WX(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_WY(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_WZ(1,1), STAT=ALLOC_ERR)

        ALLOCATE(cs_VXVY(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_VXVZ(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_VYVZ(1,1), STAT=ALLOC_ERR)

        ALLOCATE(cs_QX(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_QY(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_QZ(1,1), STAT=ALLOC_ERR)

     END IF

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

     DO s = 1, N_spec

! clear arrays
        cs_N = 0.0

        cs_VX = 0.0
        cs_VY = 0.0
        cs_VZ = 0.0

        cs_WX = 0.0
        cs_WY = 0.0
        cs_WZ = 0.0

        cs_VXVY = 0.0
        cs_VXVZ = 0.0
        cs_VYVZ = 0.0

        cs_QX = 0.0
        cs_QY = 0.0
        cs_QZ = 0.0

        CALL ADVANCE_IONS_AND_CALCULATE_MOMENTS_2D(s)

        IF (cluster_rank_key.EQ.0) CALL COLLECT_ION_DATA_FOR_AVERAGED_SNAPSHOT(s)

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     END DO

! cleanup

     DEALLOCATE(cs_N, STAT = ALLOC_ERR)

     DEALLOCATE(cs_VX, STAT = ALLOC_ERR)
     DEALLOCATE(cs_VY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_VZ, STAT = ALLOC_ERR)

     DEALLOCATE(cs_WX, STAT = ALLOC_ERR)
     DEALLOCATE(cs_WY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_WZ, STAT = ALLOC_ERR)

     DEALLOCATE(cs_VXVY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_VXVZ, STAT = ALLOC_ERR)
     DEALLOCATE(cs_VYVZ, STAT = ALLOC_ERR)

     DEALLOCATE(cs_QX, STAT = ALLOC_ERR)
     DEALLOCATE(cs_QY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_QZ, STAT = ALLOC_ERR)

  ELSE

     IF (T_cntr.EQ.Save_probes_i_data_T_cntr) THEN
        CALL ADVANCE_IONS_AND_CALCULATE_MOMENTS_PROBES
     ELSE
        CALL ADVANCE_IONS
     END IF

  END IF   !###   IF (collect_ion_moments_now) THEN

END SUBROUTINE ADVANCE_IONS_PLUS


!---------------------------------------
!
SUBROUTINE ADVANCE_IONS_AND_CALCULATE_MOMENTS_2D(s)

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles
  USE ClusterAndItsBoundaries
  
!------------------------------------------>>>
  USE Snapshots, ONLY : cs_N, cs_VX, cs_VY, cs_VZ, cs_WX, cs_WY, cs_WZ, cs_VXVY, cs_VXVZ, cs_VYVZ, cs_QX, cs_QY, cs_QZ
  USE Diagnostics
!------------------------------------------<<<

  use mpi

  IMPLICIT NONE


  INTEGER, INTENT(IN) :: s

  INTEGER errcode,ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER i, j, k
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

!------------------------------------>>>
  INTEGER n1  ! number of nodes in the y-direction
  INTEGER n2  !
  INTEGER n3  ! number of nodes in the x-direction

  REAL, ALLOCATABLE :: rbufer_n(:)

  REAL, ALLOCATABLE :: rbufer_vx(:)
  REAL, ALLOCATABLE :: rbufer_vy(:)
  REAL, ALLOCATABLE :: rbufer_vz(:)

  REAL, ALLOCATABLE :: rbufer_vx2(:)
  REAL, ALLOCATABLE :: rbufer_vy2(:)
  REAL, ALLOCATABLE :: rbufer_vz2(:)

  REAL, ALLOCATABLE :: rbufer_vxvy(:)
  REAL, ALLOCATABLE :: rbufer_vxvz(:)
  REAL, ALLOCATABLE :: rbufer_vyvz(:)

  REAL, ALLOCATABLE :: rbufer_vx3(:)
  REAL, ALLOCATABLE :: rbufer_vy3(:)
  REAL, ALLOCATABLE :: rbufer_vz3(:)

  INTEGER ALLOC_ERR
  INTEGER bufsize
  INTEGER pos_i_j, pos_ip1_j, pos_i_jp1, pos_ip1_jp1
  REAL vij, vip1j, vijp1, vip1jp1

  REAL(8) updVX, updVY, updVZ

  INTEGER npc, npa

  REAL rvx, rvy, rvz
  REAL rvx2, rvy2, rvz2
  REAL rvxvy, rvxvz, rvyvz
  REAL rvx3, rvy3, rvz3
  REAL inv_N, rtemp
!------------------------------------<<<

! functions
  REAL(8) Bx, By, Bz, Ez

! counters of particles to be sent to neighbor processes are cleared before this procedure

! counters of particles that hit the boundary objects are cleared before this procedure

! cycle over ion species is outside this procedure

!------------------------------------>>>
  n1 = c_indx_y_max - c_indx_y_min + 1
  n3 = c_indx_x_max - c_indx_x_min + 1
  n2 = -c_indx_x_min + 1 - c_indx_y_min * n3

  bufsize = n1 * n3
  ALLOCATE(rbufer_n(1:bufsize), STAT=ALLOC_ERR)

  ALLOCATE(rbufer_vx(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vy(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vz(1:bufsize), STAT=ALLOC_ERR)

  ALLOCATE(rbufer_vx2(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vy2(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vz2(1:bufsize), STAT=ALLOC_ERR)

  ALLOCATE(rbufer_vxvy(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vxvz(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vyvz(1:bufsize), STAT=ALLOC_ERR)

  ALLOCATE(rbufer_vx3(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vy3(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vz3(1:bufsize), STAT=ALLOC_ERR)

  rbufer_n = 0.0

  rbufer_vx = 0.0
  rbufer_vy = 0.0
  rbufer_vz = 0.0

  rbufer_vx2 = 0.0
  rbufer_vy2 = 0.0
  rbufer_vz2 = 0.0

  rbufer_vxvy = 0.0
  rbufer_vxvz = 0.0
  rbufer_vyvz = 0.0

  rbufer_vx3 = 0.0
  rbufer_vy3 = 0.0
  rbufer_vz3 = 0.0
!------------------------------------<<<

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
        errcode=300
        CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
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
        errcode=301
        CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
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

!------------------------------------>>>
     pos_i_j     = i + j * n3 + n2
     pos_ip1_j   = pos_i_j + 1
     pos_i_jp1   = pos_i_j + n3
     pos_ip1_jp1 = pos_i_jp1 + 1

     vij   = REAL(ax_i   * ay_j)
     vip1j = REAL(ax_ip1 * ay_j)
     vijp1 = REAL(ax_i   * ay_jp1)
     vip1jp1 = 1.0 - vij - vip1j - vijp1

     rbufer_n(pos_i_j)     = rbufer_n(pos_i_j)     + vij     !ax_i   * ay_j
     rbufer_n(pos_ip1_j)   = rbufer_n(pos_ip1_j)   + vip1j   !ax_ip1 * ay_j
     rbufer_n(pos_i_jp1)   = rbufer_n(pos_i_jp1)   + vijp1   !ax_i   * ay_jp1
     rbufer_n(pos_ip1_jp1) = rbufer_n(pos_ip1_jp1) + vip1jp1 !ax_ip1 * ay_jp1

     updVX = ion(s)%part(k)%VX
     updVY = ion(s)%part(k)%VY
     updVZ = ion(s)%part(k)%VZ
!------------------------------------<<<
     
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

!------------------------------------>>>
     updVX = 0.5_8 * (ion(s)%part(k)%VX + updVX)
     updVY = 0.5_8 * (ion(s)%part(k)%VY + updVY)
     updVZ = 0.5_8 * (ion(s)%part(k)%VZ + updVZ)

!### ion velocities updV* correspond to the same time as the coordinates before update ###
     rvx = REAL(updVX)
     rvy = REAL(updVY)
     rvz = REAL(updVZ)

     rvx2 = rvx * rvx
     rvy2 = rvy * rvy
     rvz2 = rvz * rvz

     rvxvy = rvx * rvy
     rvxvz = rvx * rvz
     rvyvz = rvy * rvz

     rvx3 = (rvx2 + rvy2 + rvz2) * rvx
     rvy3 = (rvx2 + rvy2 + rvz2) * rvy
     rvz3 = (rvx2 + rvy2 + rvz2) * rvz

     rbufer_vx(pos_i_j)     = rbufer_vx(pos_i_j)     + rvx * vij
     rbufer_vx(pos_ip1_j)   = rbufer_vx(pos_ip1_j)   + rvx * vip1j
     rbufer_vx(pos_i_jp1)   = rbufer_vx(pos_i_jp1)   + rvx * vijp1
     rbufer_vx(pos_ip1_jp1) = rbufer_vx(pos_ip1_jp1) + rvx * vip1jp1

     rbufer_vy(pos_i_j)     = rbufer_vy(pos_i_j)     + rvy * vij
     rbufer_vy(pos_ip1_j)   = rbufer_vy(pos_ip1_j)   + rvy * vip1j
     rbufer_vy(pos_i_jp1)   = rbufer_vy(pos_i_jp1)   + rvy * vijp1
     rbufer_vy(pos_ip1_jp1) = rbufer_vy(pos_ip1_jp1) + rvy * vip1jp1

     rbufer_vz(pos_i_j)     = rbufer_vz(pos_i_j)     + rvz * vij
     rbufer_vz(pos_ip1_j)   = rbufer_vz(pos_ip1_j)   + rvz * vip1j
     rbufer_vz(pos_i_jp1)   = rbufer_vz(pos_i_jp1)   + rvz * vijp1
     rbufer_vz(pos_ip1_jp1) = rbufer_vz(pos_ip1_jp1) + rvz * vip1jp1

     rbufer_vx2(pos_i_j)     = rbufer_vx2(pos_i_j)     + rvx2 * vij
     rbufer_vx2(pos_ip1_j)   = rbufer_vx2(pos_ip1_j)   + rvx2 * vip1j
     rbufer_vx2(pos_i_jp1)   = rbufer_vx2(pos_i_jp1)   + rvx2 * vijp1
     rbufer_vx2(pos_ip1_jp1) = rbufer_vx2(pos_ip1_jp1) + rvx2 * vip1jp1

     rbufer_vy2(pos_i_j)     = rbufer_vy2(pos_i_j)     + rvy2 * vij
     rbufer_vy2(pos_ip1_j)   = rbufer_vy2(pos_ip1_j)   + rvy2 * vip1j
     rbufer_vy2(pos_i_jp1)   = rbufer_vy2(pos_i_jp1)   + rvy2 * vijp1
     rbufer_vy2(pos_ip1_jp1) = rbufer_vy2(pos_ip1_jp1) + rvy2 * vip1jp1

     rbufer_vz2(pos_i_j)     = rbufer_vz2(pos_i_j)     + rvz2 * vij
     rbufer_vz2(pos_ip1_j)   = rbufer_vz2(pos_ip1_j)   + rvz2 * vip1j
     rbufer_vz2(pos_i_jp1)   = rbufer_vz2(pos_i_jp1)   + rvz2 * vijp1
     rbufer_vz2(pos_ip1_jp1) = rbufer_vz2(pos_ip1_jp1) + rvz2 * vip1jp1

     rbufer_vxvy(pos_i_j)     = rbufer_vxvy(pos_i_j)     + rvxvy * vij
     rbufer_vxvy(pos_ip1_j)   = rbufer_vxvy(pos_ip1_j)   + rvxvy * vip1j
     rbufer_vxvy(pos_i_jp1)   = rbufer_vxvy(pos_i_jp1)   + rvxvy * vijp1
     rbufer_vxvy(pos_ip1_jp1) = rbufer_vxvy(pos_ip1_jp1) + rvxvy * vip1jp1

     rbufer_vxvz(pos_i_j)     = rbufer_vxvz(pos_i_j)     + rvxvz * vij
     rbufer_vxvz(pos_ip1_j)   = rbufer_vxvz(pos_ip1_j)   + rvxvz * vip1j
     rbufer_vxvz(pos_i_jp1)   = rbufer_vxvz(pos_i_jp1)   + rvxvz * vijp1
     rbufer_vxvz(pos_ip1_jp1) = rbufer_vxvz(pos_ip1_jp1) + rvxvz * vip1jp1

     rbufer_vyvz(pos_i_j)     = rbufer_vyvz(pos_i_j)     + rvyvz * vij
     rbufer_vyvz(pos_ip1_j)   = rbufer_vyvz(pos_ip1_j)   + rvyvz * vip1j
     rbufer_vyvz(pos_i_jp1)   = rbufer_vyvz(pos_i_jp1)   + rvyvz * vijp1
     rbufer_vyvz(pos_ip1_jp1) = rbufer_vyvz(pos_ip1_jp1) + rvyvz * vip1jp1

     rbufer_vx3(pos_i_j)     = rbufer_vx3(pos_i_j)     + rvx3 * vij
     rbufer_vx3(pos_ip1_j)   = rbufer_vx3(pos_ip1_j)   + rvx3 * vip1j
     rbufer_vx3(pos_i_jp1)   = rbufer_vx3(pos_i_jp1)   + rvx3 * vijp1
     rbufer_vx3(pos_ip1_jp1) = rbufer_vx3(pos_ip1_jp1) + rvx3 * vip1jp1

     rbufer_vy3(pos_i_j)     = rbufer_vy3(pos_i_j)     + rvy3 * vij
     rbufer_vy3(pos_ip1_j)   = rbufer_vy3(pos_ip1_j)   + rvy3 * vip1j
     rbufer_vy3(pos_i_jp1)   = rbufer_vy3(pos_i_jp1)   + rvy3 * vijp1
     rbufer_vy3(pos_ip1_jp1) = rbufer_vy3(pos_ip1_jp1) + rvy3 * vip1jp1

     rbufer_vz3(pos_i_j)     = rbufer_vz3(pos_i_j)     + rvz3 * vij
     rbufer_vz3(pos_ip1_j)   = rbufer_vz3(pos_ip1_j)   + rvz3 * vip1j
     rbufer_vz3(pos_i_jp1)   = rbufer_vz3(pos_i_jp1)   + rvz3 * vijp1
     rbufer_vz3(pos_ip1_jp1) = rbufer_vz3(pos_ip1_jp1) + rvz3 * vip1jp1
!------------------------------------<<<

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
           errcode=302
           CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
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
           errcode=303
           CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
        END IF   !###   IF ( (ion(s)%part(k)%Y.GE.(c_Y_area_min+1.0_8)).AND.(ion(s)%part(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN

        CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
        CYCLE

     END IF   !###   IF (ion(s)%part(k)%X.GT.c_X_area_max) THEN

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
     END IF   !###   IF (ion(s)%part(k)%Y.GT.c_Y_area_max) THEN

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
     END IF   !###   IF (ion(s)%part(k)%Y.LT.c_Y_area_min) THEN

  END DO  !   ###   DO WHILE (k.LT.N_ions(s))

!------------------------------------------->>>

! collect moments from all processes in a cluster
  CALL MPI_REDUCE(rbufer_n, cs_N, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_REDUCE(rbufer_vx, cs_VX, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vy, cs_VY, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vz, cs_VZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  
  CALL MPI_REDUCE(rbufer_vx2, cs_WX, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vy2, cs_WY, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vz2, cs_WZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  
  CALL MPI_REDUCE(rbufer_vxvy, cs_VXVY, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vxvz, cs_VXVZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vyvz, cs_VYVZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
    
  CALL MPI_REDUCE(rbufer_vx3, cs_QX, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vy3, cs_QY, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vz3, cs_QZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  DEALLOCATE(rbufer_n, STAT=ALLOC_ERR)

  DEALLOCATE(rbufer_vx, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vy, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vz, STAT=ALLOC_ERR)

  DEALLOCATE(rbufer_vx2, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vy2, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vz2, STAT=ALLOC_ERR)

  DEALLOCATE(rbufer_vxvy, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vxvz, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vyvz, STAT=ALLOC_ERR)

  DEALLOCATE(rbufer_vx3, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vy3, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vz3, STAT=ALLOC_ERR)

! particle calculators which are not masters are done
  IF (cluster_rank_key.NE.0) RETURN

  IF (T_cntr.EQ.Save_probes_i_data_T_cntr) THEN
! save ion VDF moments in probes
! do this before synchronization because otherwise nodes in overlapped areas will not be processed correctly
! contributions from multiple clusters will be assembled in DO_PROBE_DIAGNOSTICS* procedures
     DO npc = 1, N_of_probes_cluster
        npa = List_of_probes_cluster(npc)    ! npa = n-umber of p-robe, a-ll   npc = n-umber of p-robe in c-luster
        i = Probe_position(1, npa)
        j = Probe_position(2, npa)

        probe_Ni_cluster(npc, s) = cs_N(i,j)

        probe_VXi_cluster(npc, s) = cs_VX(i,j)
        probe_VYi_cluster(npc, s) = cs_VY(i,j)
        probe_VZi_cluster(npc, s) = cs_VZ(i,j)

        probe_WXi_cluster(npc, s) = cs_WX(i,j)
        probe_WYi_cluster(npc, s) = cs_WY(i,j)
        probe_WZi_cluster(npc, s) = cs_WZ(i,j)

        probe_VXVYi_cluster(npc, s) = cs_VXVY(i,j)
        probe_VXVZi_cluster(npc, s) = cs_VXVZ(i,j)
        probe_VYVZi_cluster(npc, s) = cs_VYVZ(i,j)

        probe_QXi_cluster(npc, s) = cs_QX(i,j)
        probe_QYi_cluster(npc, s) = cs_QY(i,j)
        probe_QZi_cluster(npc, s) = cs_QZ(i,j)
     END DO
  END IF

! now cluster masters exchange information about densities in overlapping nodes

  CALL SYNCHRONIZE_MOMENTS_IN_OVERLAP_NODES

! calculate average dimensional velocities, energies, and heat flows
  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min, c_indx_x_max
        IF (cs_N(i,j).GT.1.0e-9) THEN    ! note this is small but not zero

           inv_N = 1.0 / cs_N(i,j)

           cs_VX(i, j) = cs_VX(i, j) * inv_N
           cs_VY(i, j) = cs_VY(i, j) * inv_N
           cs_VZ(i, j) = cs_VZ(i, j) * inv_N

           cs_WX(i, j) = cs_WX(i, j) * inv_N
           cs_WY(i, j) = cs_WY(i, j) * inv_N
           cs_WZ(i, j) = cs_WZ(i, j) * inv_N

           cs_VXVY(i, j) = 2.0 * cs_VXVY(i, j) * inv_N
           cs_VXVZ(i, j) = 2.0 * cs_VXVZ(i, j) * inv_N
           cs_VYVZ(i, j) = 2.0 * cs_VYVZ(i, j) * inv_N

           rvx = cs_VX(i, j)
           rvy = cs_VY(i, j)
           rvz = cs_VZ(i, j)

           rvx2 = cs_WX(i, j)
           rvy2 = cs_WY(i, j)
           rvz2 = cs_WZ(i, j)

!           rtemp = 2.0 * (rvx * rvx + rvy * rvy + rvz * rvz) - rvx2 - rvy2 - rvz2
           rtemp = rvx * rvx + rvy * rvy + rvz * rvz
           rtemp = rtemp + rtemp - rvx2 - rvy2 - rvz2

           cs_QX(i,j) = cs_QX(i,j) * inv_N + &
                      & (rtemp - rvx2 - rvx2) * rvx - &
                      & cs_VXVY(i,j) * rvy - & 
                      & cs_VXVZ(i,j) * rvz

           cs_QY(i,j) = cs_QY(i,j) * inv_N - &
                      & cs_VXVY(i,j) * rvx + &
                      & (rtemp - rvy2 - rvy2) * rvy - &
                      & cs_VYVZ(i,j) * rvz

           cs_QZ(i,j) = cs_QZ(i,j) * inv_N - &
                      & cs_VXVZ(i,j) * rvx - &
                      & cs_VYVZ(i,j) * rvy + &
                      & (rtemp - rvz2 - rvz2) * rvz

        ELSE
           cs_N(i,j) = 0.0
           cs_VX(i, j) = 0.0
           cs_VY(i, j) = 0.0
           cs_VZ(i, j) = 0.0
           cs_WX(i, j) = 0.0
           cs_WY(i, j) = 0.0
           cs_WZ(i, j) = 0.0
           cs_VXVY(i, j) = 0.0
           cs_VXVZ(i, j) = 0.0
           cs_VYVZ(i, j) = 0.0
           cs_QX(i, j) = 0.0
           cs_QY(i, j) = 0.0
           cs_QZ(i, j) = 0.0
        END IF

     END DO
  END DO
  
! adjust densities at the boundaries with material walls

  CALL ADJUST_DENSITY_AT_WALL_BOUNDARIES
!-------------------------------------------<<<

END SUBROUTINE ADVANCE_IONS_AND_CALCULATE_MOMENTS_2D


!---------------------------------------
!
SUBROUTINE ADVANCE_IONS_AND_CALCULATE_MOMENTS_PROBES

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles
  USE ClusterAndItsBoundaries
  
  USE Diagnostics
  
  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr

!------------------------------------------>>>
  INTEGER bufsize
  REAL, ALLOCATABLE :: rbufer_local(:)
  REAL, ALLOCATABLE :: rbufer_local_2(:)
  INTEGER ALLOC_ERR
!------------------------------------------<<<

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

!------------------------------------------>>>
  REAL(8) updVX, updVY, updVZ
  LOGICAL no_probe_in_cell_corner
  INTEGER npc, npa, pos
  REAL weight
  REAL rvx, rvy, rvz
  REAL rvx2, rvy2, rvz2
  REAL rvxvy, rvxvz, rvyvz
  REAL rvx3, rvy3, rvz3
  REAL inv_N, rtemp
!------------------------------------------<<<

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

!---------------------------------->>>
  IF (N_of_probes_cluster.GT.0) THEN
     bufsize = 13 * N_of_probes_cluster * N_spec  ! {N,VX,VY,VZ,WX,WY,WZ,VXVY,VXVZ,VYVZ,QX,QY,QZ}
     ALLOCATE(rbufer_local(bufsize), STAT=ALLOC_ERR)
     ALLOCATE(rbufer_local_2(bufsize), STAT=ALLOC_ERR)
     rbufer_local   = 0.0
     rbufer_local_2 = 0.0
  END IF
!----------------------------------<<<

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
           errcode=304
           CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
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
           errcode=305
           CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
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

!--------------------------------->>>
        updVX = ion(s)%part(k)%VX
        updVY = ion(s)%part(k)%VY
        updVZ = ion(s)%part(k)%VZ
!---------------------------------<<<

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

!--------------------------------->>>
! collect particle contribution to the electron moments in all close probes
        DO npc = 1, N_of_probes_cluster

           npa = List_of_probes_cluster(npc)
           no_probe_in_cell_corner = .TRUE.

           IF (i.EQ.Probe_position(1,npa)) THEN
              IF (j.EQ.Probe_position(2,npa)) THEN
! the probe is in the left bottom corner of the cell containing the particle
                 no_probe_in_cell_corner = .FALSE.
                 weight = REAL(ax_i * ay_j)   ! wij
              ELSE IF ((j+1).EQ.Probe_position(2,npa)) THEN
! the probe is in the left top corner of the cell containing the particle
                 no_probe_in_cell_corner = .FALSE.
                 weight = REAL(ax_i * ay_jp1)   ! vijp1
              END IF
           ELSE IF ((i+1).EQ.Probe_position(1,npa)) THEN
              IF (j.EQ.Probe_position(2,npa)) THEN
! the probe is in the right bottom corner of the cell containing the particle
                 no_probe_in_cell_corner = .FALSE.
                 weight = REAL(ax_ip1 * ay_j)   ! vip1j
              ELSE IF ((j+1).EQ.Probe_position(2,npa)) THEN
! the probe is in the right top corner of the cell containing the particle
                 no_probe_in_cell_corner = .FALSE.
                 weight = REAL(ax_ip1 * ay_jp1)   ! vip1jp1
              END IF
           END IF

           IF (no_probe_in_cell_corner) CYCLE

           pos = (npc-1) * 13 + (s-1) * N_of_probes_cluster * 13

           rvx = REAL(0.5_8 * (ion(s)%part(k)%VX + updVX))
           rvy = REAL(0.5_8 * (ion(s)%part(k)%VY + updVY))
           rvz = REAL(0.5_8 * (ion(s)%part(k)%VZ + updVZ))

           rvx2 = rvx * rvx
           rvy2 = rvy * rvy
           rvz2 = rvz * rvz

           rvxvy = rvx * rvy
           rvxvz = rvx * rvz
           rvyvz = rvy * rvz

           rvx3 = (rvx2 + rvy2 + rvz2) * rvx
           rvy3 = (rvx2 + rvy2 + rvz2) * rvy
           rvz3 = (rvx2 + rvy2 + rvz2) * rvz

           rbufer_local(pos+1) = rbufer_local(pos+1) + weight

           rbufer_local(pos+2) = rbufer_local(pos+2) + weight * rvx
           rbufer_local(pos+3) = rbufer_local(pos+3) + weight * rvy
           rbufer_local(pos+4) = rbufer_local(pos+4) + weight * rvz

           rbufer_local(pos+5) = rbufer_local(pos+5) + weight * rvx2
           rbufer_local(pos+6) = rbufer_local(pos+6) + weight * rvy2
           rbufer_local(pos+7) = rbufer_local(pos+7) + weight * rvz2

           rbufer_local(pos+8)  = rbufer_local(pos+8)  + weight * rvxvy
           rbufer_local(pos+9)  = rbufer_local(pos+9)  + weight * rvxvz
           rbufer_local(pos+10) = rbufer_local(pos+10) + weight * rvyvz

           rbufer_local(pos+11) = rbufer_local(pos+11) + weight * rvx3
           rbufer_local(pos+12) = rbufer_local(pos+12) + weight * rvy3
           rbufer_local(pos+13) = rbufer_local(pos+13) + weight * rvz3

        END DO   !###   DO npc = 1, N_of_probes_cluster
!--------------------------------->>>

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
              errcode=306
              CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
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
              errcode=307
              CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
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

  IF (N_of_probes_cluster.LE.0) RETURN

! collect moments from all processes in a cluster
  CALL MPI_REDUCE(rbufer_local, rbufer_local_2, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  IF (cluster_rank_key.EQ.0) THEN
! cluster master translates the message and stores data in permanent arrays
     pos=0
     DO s = 1, N_spec
        DO npc = 1, N_of_probes_cluster
           probe_Ni_cluster(npc, s) = rbufer_local_2(pos+1)

           probe_VXi_cluster(npc, s) = rbufer_local_2(pos+2)
           probe_VYi_cluster(npc, s) = rbufer_local_2(pos+3)
           probe_VZi_cluster(npc, s) = rbufer_local_2(pos+4)

           probe_WXi_cluster(npc, s) = rbufer_local_2(pos+5)
           probe_WYi_cluster(npc, s) = rbufer_local_2(pos+6)
           probe_WZi_cluster(npc, s) = rbufer_local_2(pos+7)

! the next three arrays are necessary because there may be probes in the overlapping areas
! and the VDF moments in these probes have not been finalized yet

           probe_VXVYi_cluster(npc, s) = rbufer_local_2(pos+8)
           probe_VXVZi_cluster(npc, s) = rbufer_local_2(pos+9)
           probe_VYVZi_cluster(npc, s) = rbufer_local_2(pos+10)

           probe_QXi_cluster(npc, s) = rbufer_local_2(pos+11)
           probe_QYi_cluster(npc, s) = rbufer_local_2(pos+12)
           probe_QZi_cluster(npc, s) = rbufer_local_2(pos+13)
        
           pos= pos+13
        END DO
     END DO
  END IF

! cleanup
  DEALLOCATE(rbufer_local, STAT = ALLOC_ERR)
  DEALLOCATE(rbufer_local_2, STAT = ALLOC_ERR)

!print '("Proc ",i4," will add ",i8," send l/r/a/b ",4(2x,i8))', Rank_of_process, N_e_to_add, N_ions_to_send_left, N_ions_to_send_right, N_ions_to_send_above, N_ions_to_send_below

END SUBROUTINE ADVANCE_IONS_AND_CALCULATE_MOMENTS_PROBES
