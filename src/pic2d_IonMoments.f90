!------------------------------
!
SUBROUTINE COLLECT_ION_MOMENTS(s)

  USE ParallelOperationValues
  USE IonParticles
  USE ClusterAndItsBoundaries
  USE Snapshots

  use mpi

  IMPLICIT NONE


  INTEGER, INTENT(IN) :: s

  INTEGER errcode,ierr
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

  REAL, ALLOCATABLE :: rbufer_vxvy(:)
  REAL, ALLOCATABLE :: rbufer_vxvz(:)
  REAL, ALLOCATABLE :: rbufer_vyvz(:)

  REAL, ALLOCATABLE :: rbufer_vx3(:)
  REAL, ALLOCATABLE :: rbufer_vy3(:)
  REAL, ALLOCATABLE :: rbufer_vz3(:)

  INTEGER ALLOC_ERR
  INTEGER bufsize

  INTEGER i, j, k
  INTEGER pos_i_j, pos_ip1_j, pos_i_jp1, pos_ip1_jp1
  REAL ax_ip1, ax_i, ay_jp1, ay_j
  REAL vij, vip1j, vijp1, vip1jp1

  REAL rvx, rvy, rvz
  REAL rvx2, rvy2, rvz2
  REAL rvxvy, rvxvz, rvyvz
  REAL rvx3, rvy3, rvz3
  REAL inv_N, rtemp

! clear arrays set in MPI_REDUCE (just in case)
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
        errcode=310
        CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
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

     rvx = REAL(ion(s)%part(k)%VX)
     rvy = REAL(ion(s)%part(k)%VY)
     rvz = REAL(ion(s)%part(k)%VZ)

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

  END DO

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

END SUBROUTINE COLLECT_ION_MOMENTS

