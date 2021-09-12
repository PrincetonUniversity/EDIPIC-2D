!------------------------------
!
SUBROUTINE COLLECT_ION_MOMENTS(s)

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : V_scale_ms, energy_factor_eV, m_e_kg, e_Cl, N_scale_part_m3
  USE IonParticles
  USE ClusterAndItsBoundaries
!  USE BlockAndItsBoundaries, ONLY : indx_x_min, indx_x_max, indx_y_min, indx_y_max
  USE Snapshots

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, INTENT(IN) :: s

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

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
  INTEGER ALLOC_ERR
  INTEGER bufsize

  INTEGER i, j, k
  INTEGER pos_i_j, pos_ip1_j, pos_i_jp1, pos_ip1_jp1
  REAL ax_ip1, ax_i, ay_jp1, ay_j
  REAL vij, vip1j, vijp1, vip1jp1
  REAL v_ij, v_ip1j, v_ijp1, v_ip1jp1

! clear arrays set in MPI_REDUCE (just in case)
  cs_N = 0.0
  cs_VX = 0.0
  cs_VY = 0.0
  cs_VZ = 0.0
  cs_WX = 0.0
  cs_WY = 0.0
  cs_WZ = 0.0

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

  rbufer_n = 0.0
  rbufer_vx = 0.0
  rbufer_vy = 0.0
  rbufer_vz = 0.0
  rbufer_vx2 = 0.0
  rbufer_vy2 = 0.0
  rbufer_vz2 = 0.0

  DO k = 1, N_ions(s)
     
     i = INT(ion(s)%part(k)%X)
     j = INT(ion(s)%part(k)%Y)
     IF (ion(s)%part(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (ion(s)%part(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

if ((i.lt.c_indx_x_min).or.(i.gt.(c_indx_x_max-1)).or.(j.lt.c_indx_y_min).or.(j.gt.(c_indx_y_max-1))) then
   print '("Process ",i4," : Error-1 in COLLECT_ION_MOMENTS : index out of bounds")', Rank_of_process
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

     ax_ip1 = REAL(ion(s)%part(k)%X - DBLE(i))
     ax_i   = 1.0 - ax_ip1

     ay_jp1 = REAL(ion(s)%part(k)%Y - DBLE(j))
     ay_j = 1.0 - ay_jp1

     vij   = ax_i   * ay_j
     vip1j = ax_ip1 * ay_j
     vijp1 = ax_i   * ay_jp1
     vip1jp1 = 1.0 - vij - vip1j - vijp1

     rbufer_n(pos_i_j)     = rbufer_n(pos_i_j)     + vij     !ax_i   * ay_j
     rbufer_n(pos_ip1_j)   = rbufer_n(pos_ip1_j)   + vip1j   !ax_ip1 * ay_j
     rbufer_n(pos_i_jp1)   = rbufer_n(pos_i_jp1)   + vijp1   !ax_i   * ay_jp1
     rbufer_n(pos_ip1_jp1) = rbufer_n(pos_ip1_jp1) + vip1jp1 !ax_ip1 * ay_jp1

     v_ij     = REAL(ion(s)%part(k)%VX) * vij
     v_ip1j   = REAL(ion(s)%part(k)%VX) * vip1j
     v_ijp1   = REAL(ion(s)%part(k)%VX) * vijp1
     v_ip1jp1 = REAL(ion(s)%part(k)%VX) * vip1jp1

     rbufer_vx(pos_i_j)     = rbufer_vx(pos_i_j)     + v_ij
     rbufer_vx(pos_ip1_j)   = rbufer_vx(pos_ip1_j)   + v_ip1j
     rbufer_vx(pos_i_jp1)   = rbufer_vx(pos_i_jp1)   + v_ijp1
     rbufer_vx(pos_ip1_jp1) = rbufer_vx(pos_ip1_jp1) + v_ip1jp1

     rbufer_vx2(pos_i_j)     = rbufer_vx2(pos_i_j)     + REAL(ion(s)%part(k)%VX) * v_ij
     rbufer_vx2(pos_ip1_j)   = rbufer_vx2(pos_ip1_j)   + REAL(ion(s)%part(k)%VX) * v_ip1j
     rbufer_vx2(pos_i_jp1)   = rbufer_vx2(pos_i_jp1)   + REAL(ion(s)%part(k)%VX) * v_ijp1
     rbufer_vx2(pos_ip1_jp1) = rbufer_vx2(pos_ip1_jp1) + REAL(ion(s)%part(k)%VX) * v_ip1jp1

     v_ij     = REAL(ion(s)%part(k)%VY) * vij
     v_ip1j   = REAL(ion(s)%part(k)%VY) * vip1j
     v_ijp1   = REAL(ion(s)%part(k)%VY) * vijp1
     v_ip1jp1 = REAL(ion(s)%part(k)%VY) * vip1jp1

     rbufer_vy(pos_i_j)     = rbufer_vy(pos_i_j)     + v_ij
     rbufer_vy(pos_ip1_j)   = rbufer_vy(pos_ip1_j)   + v_ip1j
     rbufer_vy(pos_i_jp1)   = rbufer_vy(pos_i_jp1)   + v_ijp1
     rbufer_vy(pos_ip1_jp1) = rbufer_vy(pos_ip1_jp1) + v_ip1jp1

     rbufer_vy2(pos_i_j)     = rbufer_vy2(pos_i_j)     + REAL(ion(s)%part(k)%VY) * v_ij
     rbufer_vy2(pos_ip1_j)   = rbufer_vy2(pos_ip1_j)   + REAL(ion(s)%part(k)%VY) * v_ip1j
     rbufer_vy2(pos_i_jp1)   = rbufer_vy2(pos_i_jp1)   + REAL(ion(s)%part(k)%VY) * v_ijp1
     rbufer_vy2(pos_ip1_jp1) = rbufer_vy2(pos_ip1_jp1) + REAL(ion(s)%part(k)%VY) * v_ip1jp1

     v_ij     = REAL(ion(s)%part(k)%VZ) * vij
     v_ip1j   = REAL(ion(s)%part(k)%VZ) * vip1j
     v_ijp1   = REAL(ion(s)%part(k)%VZ) * vijp1
     v_ip1jp1 = REAL(ion(s)%part(k)%VZ) * vip1jp1

     rbufer_vz(pos_i_j)     = rbufer_vz(pos_i_j)     + v_ij
     rbufer_vz(pos_ip1_j)   = rbufer_vz(pos_ip1_j)   + v_ip1j
     rbufer_vz(pos_i_jp1)   = rbufer_vz(pos_i_jp1)   + v_ijp1
     rbufer_vz(pos_ip1_jp1) = rbufer_vz(pos_ip1_jp1) + v_ip1jp1

     rbufer_vz2(pos_i_j)     = rbufer_vz2(pos_i_j)     + REAL(ion(s)%part(k)%VZ) * v_ij
     rbufer_vz2(pos_ip1_j)   = rbufer_vz2(pos_ip1_j)   + REAL(ion(s)%part(k)%VZ) * v_ip1j
     rbufer_vz2(pos_i_jp1)   = rbufer_vz2(pos_i_jp1)   + REAL(ion(s)%part(k)%VZ) * v_ijp1
     rbufer_vz2(pos_ip1_jp1) = rbufer_vz2(pos_ip1_jp1) + REAL(ion(s)%part(k)%VZ) * v_ip1jp1

  END DO

! collect moments from all processes in a cluster
  CALL MPI_REDUCE(rbufer_n, cs_N, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_REDUCE(rbufer_vx, cs_VX, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_REDUCE(rbufer_vy, cs_VY, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_REDUCE(rbufer_vz, cs_VZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  
  CALL MPI_REDUCE(rbufer_vx2, cs_WX, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_REDUCE(rbufer_vy2, cs_WY, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_REDUCE(rbufer_vz2, cs_WZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  DEALLOCATE(rbufer_n, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vx, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vy, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vz, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vx2, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vy2, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vz2, STAT=ALLOC_ERR)

! particle calculators which are not masters are done
  IF (cluster_rank_key.NE.0) RETURN

! now cluster masters exchange information about densities in overlapping nodes

  CALL SYNCHRONIZE_MOMENTS_IN_OVERLAP_NODES

! calculate average dimensional velocities, energies, temperatures, and density
  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min, c_indx_x_max
        IF (cs_N(i,j).GT.0.0) THEN

           cs_VX(i, j) = REAL(V_scale_ms) * (cs_VX(i, j) / cs_N(i,j))
           cs_VY(i, j) = REAL(V_scale_ms) * (cs_VY(i, j) / cs_N(i,j))
           cs_VZ(i, j) = REAL(V_scale_ms) * (cs_VZ(i, j) / cs_N(i,j))

           cs_WX(i, j) = REAL(Ms(s) * energy_factor_eV) * (cs_WX(i, j) / cs_N(i,j))         ! energy_factor_eV =  0.5_8 * m_e_kg * V_scale_ms**2 / e_Cl
           cs_WY(i, j) = REAL(Ms(s) * energy_factor_eV) * (cs_WY(i, j) / cs_N(i,j))
           cs_WZ(i, j) = REAL(Ms(s) * energy_factor_eV) * (cs_WZ(i, j) / cs_N(i,j))

           cs_TX(i, j) = MAX(0.0, 2.0 * cs_WX(i, j) - REAL(Ms(s) * m_e_kg / e_Cl) * cs_VX(i, j)**2)
           cs_TY(i, j) = MAX(0.0, 2.0 * cs_WY(i, j) - REAL(Ms(s) * m_e_kg / e_Cl) * cs_VY(i, j)**2)
           cs_TZ(i, j) = MAX(0.0, 2.0 * cs_WZ(i, j) - REAL(Ms(s) * m_e_kg / e_Cl) * cs_VZ(i, j)**2)

        ELSE
           cs_VX(i, j) = 0.0
           cs_VY(i, j) = 0.0
           cs_VZ(i, j) = 0.0
           cs_WX(i, j) = 0.0
           cs_WY(i, j) = 0.0
           cs_WZ(i, j) = 0.0
           cs_TX(i, j) = 0.0
           cs_TY(i, j) = 0.0
           cs_TZ(i, j) = 0.0
        END IF

        cs_N(i,j) = cs_N(i,j) * REAL(N_scale_part_m3)

     END DO
  END DO
  
! adjust densities at the boundaries with material walls

  CALL ADJUST_DENSITY_AT_WALL_BOUNDARIES

! calculate electric currents
  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min, c_indx_x_max

        cs_JX(i,j) = e_Cl * Qs(s) * cs_N(i,j) * cs_VX(i,j)
        cs_JY(i,j) = e_Cl * Qs(s) * cs_N(i,j) * cs_VY(i,j)
        cs_JZ(i,j) = e_Cl * Qs(s) * cs_N(i,j) * cs_VZ(i,j)

        cs_JXsum(i,j) = cs_JXsum(i,j) + cs_JX(i,j)
        cs_JYsum(i,j) = cs_JYsum(i,j) + cs_JY(i,j)
        cs_JZsum(i,j) = cs_JZsum(i,j) + cs_JZ(i,j)

     END DO
  END DO
  
END SUBROUTINE COLLECT_ION_MOMENTS

!------------------------------
!
SUBROUTINE COLLECT_ION_MOMENTS_IN_CLUSTER_PROBES(ions_moved)

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : T_cntr !V_scale_ms, energy_factor_eV, m_e_kg, e_Cl, N_scale_part_m3
  USE IonParticles
  USE ClusterAndItsBoundaries
  USE Diagnostics

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL, INTENT(IN) :: ions_moved

  INTEGER ierr
!  INTEGER stattus(MPI_STATUS_SIZE)
!  INTEGER request

  INTEGER bufsize
  REAL, ALLOCATABLE :: rbufer_local(:)
  REAL, ALLOCATABLE :: rbufer_local_2(:)
  INTEGER ALLOC_ERR

  INTEGER s, i, j, k
  LOGICAL skip_this_particle
  INTEGER npc   ! probe number in the list of probes belonging to this cluster
  INTEGER npa   ! probe number in the global list of probes

  REAL ax_ip1, ax_i, ay_jp1, ay_j
  REAL vij, vip1j, vijp1, vip1jp1

  REAL(8) updVX, updVY, updVZ

  INTEGER pos

! exit if the current time layer is not assigned for the diagnostic output 
  IF (T_cntr.NE.Save_probes_data_T_cntr) RETURN

  IF (.NOT.ions_moved) RETURN

  IF (N_of_probes_cluster.LE.0) RETURN

  bufsize = 7 * N_of_probes_cluster * N_spec  ! {N,JX,JY,JZ,WX,WY,WZ}

  ALLOCATE(rbufer_local(bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_local_2(bufsize), STAT=ALLOC_ERR)

  rbufer_local   = 0.0
  rbufer_local_2 = 0.0

  DO s = 1, N_spec
     DO k = 1, N_ions(s)
     
        i = INT(ion(s)%part(k)%X)
        j = INT(ion(s)%part(k)%Y)
        IF (ion(s)%part(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
        IF (ion(s)%part(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

        skip_this_particle = .TRUE.

! find whether in any corner of the cell containing the particle there is a probe

        DO npc = 1, N_of_probes_cluster
           npa = List_of_probes_cluster(npc)
           IF (i.EQ.Probe_position(1,npa)) THEN
              IF (j.EQ.Probe_position(2,npa)) THEN
! the probe is in the left top corner of the cell containing the particle
                 skip_this_particle = .FALSE.
                 EXIT
              ELSE IF ((j+1).EQ.Probe_position(2,npa)) THEN
! the probe is in the left top corner of the cell containing the particle
                 skip_this_particle = .FALSE.
                 EXIT
              END IF
           ELSE IF ((i+1).EQ.Probe_position(1,npa)) THEN
              IF (j.EQ.Probe_position(2,npa)) THEN
! the probe is in the right bottom corner of the cell containing the particle
                 skip_this_particle = .FALSE.
                 EXIT
              ELSE IF ((j+1).EQ.Probe_position(2,npa)) THEN
! the probe is in the right top corner of the cell containing the particle
                 skip_this_particle = .FALSE.
                 EXIT
              END IF
           END IF
        END DO

        IF (skip_this_particle) CYCLE

        ax_ip1 = REAL(ion(s)%part(k)%X - DBLE(i))
        ax_i   = 1.0 - ax_ip1

        ay_jp1 = REAL(ion(s)%part(k)%Y - DBLE(j))
        ay_j = 1.0 - ay_jp1

        vij   = ax_i   * ay_j
        vip1j = ax_ip1 * ay_j
        vijp1 = ax_i   * ay_jp1
        vip1jp1 = 1.0 - vij - vip1j - vijp1

        updVX = ion(s)%part(k)%VX
        updVY = ion(s)%part(k)%VY
        updVZ = ion(s)%part(k)%VZ

! collect particle contribution to the electron moments in all close probes
        DO npc = 1, N_of_probes_cluster
           npa = List_of_probes_cluster(npc)
           IF (i.EQ.Probe_position(1,npa)) THEN
              IF (j.EQ.Probe_position(2,npa)) THEN
! the probe is in the left bottom corner of the cell containing the particle

                 pos = (npc-1) * 7 + (s-1) * N_of_probes_cluster * 7

                 rbufer_local(pos+1) = rbufer_local(pos+1) + vij

                 rbufer_local(pos+2) = rbufer_local(pos+2) + vij * updVX
                 rbufer_local(pos+3) = rbufer_local(pos+3) + vij * updVY
                 rbufer_local(pos+4) = rbufer_local(pos+4) + vij * updVZ

                 rbufer_local(pos+5) = rbufer_local(pos+5) + vij * updVX * updVX
                 rbufer_local(pos+6) = rbufer_local(pos+6) + vij * updVY * updVY
                 rbufer_local(pos+7) = rbufer_local(pos+7) + vij * updVZ * updVZ

              ELSE IF ((j+1).EQ.Probe_position(2,npa)) THEN
! the probe is in the left top corner of the cell containing the particle

                 pos = (npc-1) * 7 + (s-1) * N_of_probes_cluster * 7

                 rbufer_local(pos+1) = rbufer_local(pos+1) + vijp1

                 rbufer_local(pos+2) = rbufer_local(pos+2) + vijp1 * updVX
                 rbufer_local(pos+3) = rbufer_local(pos+3) + vijp1 * updVY
                 rbufer_local(pos+4) = rbufer_local(pos+4) + vijp1 * updVZ

                 rbufer_local(pos+5) = rbufer_local(pos+5) + vijp1 * updVX * updVX
                 rbufer_local(pos+6) = rbufer_local(pos+6) + vijp1 * updVY * updVY
                 rbufer_local(pos+7) = rbufer_local(pos+7) + vijp1 * updVZ * updVZ

              END IF

           ELSE IF ((i+1).EQ.Probe_position(1,npa)) THEN
              IF (j.EQ.Probe_position(2,npa)) THEN
! the probe is in the right bottom corner of the cell containing the particle

                 pos = (npc-1) * 7 + (s-1) * N_of_probes_cluster * 7

                 rbufer_local(pos+1) = rbufer_local(pos+1) + vip1j

                 rbufer_local(pos+2) = rbufer_local(pos+2) + vip1j * updVX
                 rbufer_local(pos+3) = rbufer_local(pos+3) + vip1j * updVY
                 rbufer_local(pos+4) = rbufer_local(pos+4) + vip1j * updVZ

                 rbufer_local(pos+5) = rbufer_local(pos+5) + vip1j * updVX * updVX
                 rbufer_local(pos+6) = rbufer_local(pos+6) + vip1j * updVY * updVY
                 rbufer_local(pos+7) = rbufer_local(pos+7) + vip1j * updVZ * updVZ

              ELSE IF ((j+1).EQ.Probe_position(2,npa)) THEN
! the probe is in the right top corner of the cell containing the particle

                 pos = (npc-1) * 7 + (s-1) * N_of_probes_cluster * 7

                 rbufer_local(pos+1) = rbufer_local(pos+1) + vip1jp1

                 rbufer_local(pos+2) = rbufer_local(pos+2) + vip1jp1 * updVX
                 rbufer_local(pos+3) = rbufer_local(pos+3) + vip1jp1 * updVY
                 rbufer_local(pos+4) = rbufer_local(pos+4) + vip1jp1 * updVZ

                 rbufer_local(pos+5) = rbufer_local(pos+5) + vip1jp1 * updVX * updVX
                 rbufer_local(pos+6) = rbufer_local(pos+6) + vip1jp1 * updVY * updVY
                 rbufer_local(pos+7) = rbufer_local(pos+7) + vip1jp1 * updVZ * updVZ

              END IF
           END IF
        END DO   !###   DO npc = 1, N_of_probes_cluster

     END DO   !###   DO k = 1, N_ions(s)
  END DO   !###   DO s = 1, N_spec

! collect moments from all processes in a cluster
  CALL MPI_REDUCE(rbufer_local, rbufer_local_2, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  IF (cluster_rank_key.EQ.0) THEN
! cluster master translates the message and stores data in permanent arrays
     pos=1
     DO s = 1, N_spec
        DO npc = 1, N_of_probes_cluster
           probe_Ni_cluster(npc, s) = rbufer_local_2(pos)

           probe_JXi_cluster(npc, s) = rbufer_local_2(pos+1)
           probe_JYi_cluster(npc, s) = rbufer_local_2(pos+2)
           probe_JZi_cluster(npc, s) = rbufer_local_2(pos+3)

           probe_WXi_cluster(npc, s) = rbufer_local_2(pos+4)
           probe_WYi_cluster(npc, s) = rbufer_local_2(pos+5)
           probe_WZi_cluster(npc, s) = rbufer_local_2(pos+6)
        
           pos= pos+7
        END DO
     END DO
  END IF

  IF (ALLOCATED(rbufer_local))   DEALLOCATE(rbufer_local, STAT=ALLOC_ERR)
  IF (ALLOCATED(rbufer_local_2)) DEALLOCATE(rbufer_local_2, STAT=ALLOC_ERR)
  
  RETURN

END SUBROUTINE COLLECT_ION_MOMENTS_IN_CLUSTER_PROBES
