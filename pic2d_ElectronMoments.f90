!------------------------------
!
SUBROUTINE COLLECT_ELECTRON_MOMENTS

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : EX, EY   ! 
  USE ElectronParticles
  USE ClusterAndItsBoundaries
  USE Snapshots

  IMPLICIT NONE

  INCLUDE 'mpif.h'

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

  REAL(8) dax_ip1, dax_i, day_jp1, day_j
  REAL(8) E_X, E_Y, E_Z
  REAL(8) Bx, By, Bz, Ez       ! functions
  REAL(8) alfa_x, alfa_y, alfa_z
  REAL(8) alfa_x2, alfa_y2, alfa_z2
  REAL(8) theta2, invtheta
  REAL(8) K11, K12, K13, K21, K22, K23, K31, K32, K33
  REAL(8) VX_minus, VY_minus, VZ_minus
  REAL(8) VX_plus, VY_plus, VZ_plus
  REAL(8) updVX, updVY, updVZ

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

  DO k = 1, N_electrons

     i = INT(electron(k)%X)
     j = INT(electron(k)%Y)
     IF (electron(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (electron(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

     if ((i.lt.c_indx_x_min).or.(i.gt.(c_indx_x_max-1)).or.(j.lt.c_indx_y_min).or.(j.gt.(c_indx_y_max-1))) then
        print '("Process ",i4," : Error-1 in COLLECT_ELECTRON_MOMENTS : index out of bounds")', Rank_of_process
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

     ax_ip1 = REAL(electron(k)%X - DBLE(i))
     ax_i   = 1.0 - ax_ip1

     ay_jp1 = REAL(electron(k)%Y - DBLE(j))
     ay_j = 1.0 - ay_jp1

     vij   = ax_i   * ay_j
     vip1j = ax_ip1 * ay_j
     vijp1 = ax_i   * ay_jp1
     vip1jp1 = 1.0 - vij - vip1j - vijp1

     rbufer_n(pos_i_j)     = rbufer_n(pos_i_j)     + vij     !ax_i   * ay_j
     rbufer_n(pos_ip1_j)   = rbufer_n(pos_ip1_j)   + vip1j   !ax_ip1 * ay_j
     rbufer_n(pos_i_jp1)   = rbufer_n(pos_i_jp1)   + vijp1   !ax_i   * ay_jp1
     rbufer_n(pos_ip1_jp1) = rbufer_n(pos_ip1_jp1) + vip1jp1 !ax_ip1 * ay_jp1

! the snapshot is taken after the electric field was calculated but before the electrons were advanced
! therefore, while the electron coordinates correspond to time level n
! the electron velocities correspond to time level n-1/2
! here we advance the electron velocity to time level n (by half-a-time-step)
! this advanced velocity is used to calculate the electron velocity distribution function moments 
! note that the advanced velocity is a temporary value which is not used beyound this subroutine
! the actual velocity vector components of a particle DO NOT CHANGE 

     dax_ip1 = electron(k)%X - DBLE(i)
     dax_i   = 1.0_8 - dax_ip1

     day_jp1 = electron(k)%Y - DBLE(j)
     day_j = 1.0_8 - day_jp1

     E_X = EX(i,j) * dax_i * day_j + EX(i+1,j) * dax_ip1 * day_j + EX(i,j+1) * dax_i * day_jp1 + EX(i+1,j+1) * dax_ip1 * day_jp1
     E_Y = EY(i,j) * dax_i * day_j + EY(i+1,j) * dax_ip1 * day_j + EY(i,j+1) * dax_i * day_jp1 + EY(i+1,j+1) * dax_ip1 * day_jp1
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

! velocity advance: first half-acceleration due to electric field

     VX_minus = electron(k)%VX - 0.5_8 * E_X
     VY_minus = electron(k)%VY - 0.5_8 * E_Y
     VZ_minus = electron(k)%VZ - 0.5_8 * E_Z

! velocity advance: rotation in the magnetic field

     VX_plus = K11 * VX_minus + K12 * VY_minus + K13 * VZ_minus
     VY_plus = K21 * VX_minus + K22 * VY_minus + K23 * VZ_minus
     VZ_plus = K31 * VX_minus + K32 * VY_minus + K33 * VZ_minus

! velocity advance: second half-acceleration due to electric field

     updVX = VX_plus - 0.5_8 * E_X
     updVY = VY_plus - 0.5_8 * E_Y
     updVZ = VZ_plus - 0.5_8 * E_Z

! since we need the velocity advanced by half-a-time-step only, take the average of previous and updated velocities

     updVX = 0.5_8 * (electron(k)%VX + updVX)
     updVY = 0.5_8 * (electron(k)%VY + updVY)
     updVZ = 0.5_8 * (electron(k)%VZ + updVZ)

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
  
! adjust moments at the boundaries with material walls

  CALL ADJUST_DENSITY_AT_WALL_BOUNDARIES

END SUBROUTINE COLLECT_ELECTRON_MOMENTS

!--------------------------------------------------------
!
SUBROUTINE SYNCHRONIZE_MOMENTS_IN_OVERLAP_NODES

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE Snapshots

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER n1  ! number of nodes in the y-direction
  INTEGER n3  ! number of nodes in the x-direction

  INTEGER i, j

  INTEGER pos1, pos2

  INTEGER ALLOC_ERR
  REAL, ALLOCATABLE :: rbufer(:)

  n1 = c_indx_y_max - c_indx_y_min + 1
  n3 = c_indx_x_max - c_indx_x_min + 1

! include neighbor contributions in nodes which are one line away from the cluster boundary

  IF (WHITE_CLUSTER) THEN  
! "white processes"

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:(13*n1)), STAT=ALLOC_ERR)

     IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right moments in the right edge
        pos1=1
        pos2=n1
        rbufer(pos1:pos2) = cs_N(c_indx_x_max, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VX(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VY(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_max, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WX(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WY(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_max, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_max, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QX(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QY(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_max, c_indx_y_min:c_indx_y_max)

        CALL MPI_SEND(rbufer, 13*n1, MPI_REAL, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left moments in the left edge
        pos1=1
        pos2=n1
        rbufer(pos1:pos2) = cs_N(c_indx_x_min, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min, c_indx_y_min:c_indx_y_max)

        CALL MPI_SEND(rbufer, 13*n1, MPI_REAL, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left moments in the vertical line next to the left edge
        CALL MPI_RECV(rbufer, 13*n1, MPI_REAL, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)

        pos1 = -c_indx_y_min
        DO j = c_indx_y_min, c_indx_y_max
           cs_N(c_indx_x_min+1, j) = cs_N(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VX(c_indx_x_min+1, j) = cs_VX(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VY(c_indx_x_min+1, j) = cs_VY(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VZ(c_indx_x_min+1, j) = cs_VZ(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WX(c_indx_x_min+1, j) = cs_WX(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WY(c_indx_x_min+1, j) = cs_WY(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WZ(c_indx_x_min+1, j) = cs_WZ(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VXVY(c_indx_x_min+1, j) = cs_VXVY(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VXVZ(c_indx_x_min+1, j) = cs_VXVZ(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VYVZ(c_indx_x_min+1, j) = cs_VYVZ(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QX(c_indx_x_min+1, j) = cs_QX(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QY(c_indx_x_min+1, j) = cs_QY(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QZ(c_indx_x_min+1, j) = cs_QZ(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO

     END IF

     IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right moments in the vertical line next to the right edge
        CALL MPI_RECV(rbufer, 13*n1, MPI_REAL, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)

        pos1 = -c_indx_y_min
        DO j = c_indx_y_min, c_indx_y_max
           cs_N(c_indx_x_max-1, j) = cs_N(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VX(c_indx_x_max-1, j) = cs_VX(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VY(c_indx_x_max-1, j) = cs_VY(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VZ(c_indx_x_max-1, j) = cs_VZ(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WX(c_indx_x_max-1, j) = cs_WX(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WY(c_indx_x_max-1, j) = cs_WY(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WZ(c_indx_x_max-1, j) = cs_WZ(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VXVY(c_indx_x_max-1, j) = cs_VXVY(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VXVZ(c_indx_x_max-1, j) = cs_VXVZ(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VYVZ(c_indx_x_max-1, j) = cs_VYVZ(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QX(c_indx_x_max-1, j) = cs_QX(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QY(c_indx_x_max-1, j) = cs_QY(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QZ(c_indx_x_max-1, j) = cs_QZ(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO

     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:13*n3), STAT=ALLOC_ERR)

     IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up moments in the top edge
        pos1=1
        pos2=n3
        rbufer(pos1:pos2) = cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_max)

        CALL MPI_SEND(rbufer, 13*n3, MPI_REAL, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down moments in the bottom edge
        pos1=1
        pos2=n3
        rbufer(pos1:pos2) = cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_min)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_min)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_min)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_min)

        CALL MPI_SEND(rbufer, 13*n3, MPI_REAL, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below moments in the vertical line above the bottom line
        CALL MPI_RECV(rbufer, 13*n3, MPI_REAL, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)

        pos1 = -c_indx_x_min
        DO i = c_indx_x_min, c_indx_x_max
           cs_N(i, c_indx_y_min+1) = cs_N(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VX(i, c_indx_y_min+1) = cs_VX(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VY(i, c_indx_y_min+1) = cs_VY(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VZ(i, c_indx_y_min+1) = cs_VZ(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WX(i, c_indx_y_min+1) = cs_WX(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WY(i, c_indx_y_min+1) = cs_WY(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WZ(i, c_indx_y_min+1) = cs_WZ(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VXVY(i, c_indx_y_min+1) = cs_VXVY(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VXVZ(i, c_indx_y_min+1) = cs_VXVZ(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VYVZ(i, c_indx_y_min+1) = cs_VYVZ(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QX(i, c_indx_y_min+1) = cs_QX(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QY(i, c_indx_y_min+1) = cs_QY(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QZ(i, c_indx_y_min+1) = cs_QZ(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO

     END IF

     IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above moments in the vertical line under the top line
        CALL MPI_RECV(rbufer, 13*n3, MPI_REAL, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)

        pos1 = -c_indx_x_min
        DO i = c_indx_x_min, c_indx_x_max
           cs_N(i, c_indx_y_max-1) = cs_N(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VX(i, c_indx_y_max-1) = cs_VX(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VY(i, c_indx_y_max-1) = cs_VY(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VZ(i, c_indx_y_max-1) = cs_VZ(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WX(i, c_indx_y_max-1) = cs_WX(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WY(i, c_indx_y_max-1) = cs_WY(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WZ(i, c_indx_y_max-1) = cs_WZ(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VXVY(i, c_indx_y_max-1) = cs_VXVY(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VXVZ(i, c_indx_y_max-1) = cs_VXVZ(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VYVZ(i, c_indx_y_max-1) = cs_VYVZ(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QX(i, c_indx_y_max-1) = cs_QX(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QY(i, c_indx_y_max-1) = cs_QY(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QZ(i, c_indx_y_max-1) = cs_QZ(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO

     END IF

  ELSE
! "black" processes

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:13*n1), STAT=ALLOC_ERR)

     IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left moments in the vertical line next to the left edge
        CALL MPI_RECV(rbufer, 13*n1, MPI_REAL, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)

        pos1 = -c_indx_y_min
        DO j = c_indx_y_min, c_indx_y_max
           cs_N(c_indx_x_min+1, j) = cs_N(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VX(c_indx_x_min+1, j) = cs_VX(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VY(c_indx_x_min+1, j) = cs_VY(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VZ(c_indx_x_min+1, j) = cs_VZ(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WX(c_indx_x_min+1, j) = cs_WX(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WY(c_indx_x_min+1, j) = cs_WY(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WZ(c_indx_x_min+1, j) = cs_WZ(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VXVY(c_indx_x_min+1, j) = cs_VXVY(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VXVZ(c_indx_x_min+1, j) = cs_VXVZ(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VYVZ(c_indx_x_min+1, j) = cs_VYVZ(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QX(c_indx_x_min+1, j) = cs_QX(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QY(c_indx_x_min+1, j) = cs_QY(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QZ(c_indx_x_min+1, j) = cs_QZ(c_indx_x_min+1, j) + rbufer(j+pos1+1)
        END DO

     END IF

     IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right moments in the vertical line next to the right edge
        CALL MPI_RECV(rbufer, 13*n1, MPI_REAL, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)

        pos1 = -c_indx_y_min
        DO j = c_indx_y_min, c_indx_y_max
           cs_N(c_indx_x_max-1, j) = cs_N(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VX(c_indx_x_max-1, j) = cs_VX(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VY(c_indx_x_max-1, j) = cs_VY(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VZ(c_indx_x_max-1, j) = cs_VZ(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WX(c_indx_x_max-1, j) = cs_WX(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WY(c_indx_x_max-1, j) = cs_WY(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_WZ(c_indx_x_max-1, j) = cs_WZ(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VXVY(c_indx_x_max-1, j) = cs_VXVY(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VXVZ(c_indx_x_max-1, j) = cs_VXVZ(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_VYVZ(c_indx_x_max-1, j) = cs_VYVZ(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO

        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QX(c_indx_x_max-1, j) = cs_QX(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QY(c_indx_x_max-1, j) = cs_QY(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO
        pos1 = pos1+n1
        DO j = c_indx_y_min, c_indx_y_max
           cs_QZ(c_indx_x_max-1, j) = cs_QZ(c_indx_x_max-1, j) + rbufer(j+pos1+1)
        END DO

     END IF

     IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right moments in the right edge
        pos1=1
        pos2=n1
        rbufer(pos1:pos2) = cs_N(c_indx_x_max, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VX(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VY(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_max, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WX(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WY(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_max, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_max, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QX(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QY(c_indx_x_max, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_max, c_indx_y_min:c_indx_y_max)

        CALL MPI_SEND(rbufer, 13*n1, MPI_REAL, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left moments in the left edge
        pos1=1
        pos2=n1
        rbufer(pos1:pos2) = cs_N(c_indx_x_min, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min, c_indx_y_min:c_indx_y_max)

        CALL MPI_SEND(rbufer, 13*n1, MPI_REAL, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:13*n3), STAT=ALLOC_ERR)

     IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below moments in the vertical line above the bottom line
        CALL MPI_RECV(rbufer, 13*n3, MPI_REAL, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
        pos1 = -c_indx_x_min
        DO i = c_indx_x_min, c_indx_x_max
           cs_N(i, c_indx_y_min+1) = cs_N(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VX(i, c_indx_y_min+1) = cs_VX(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VY(i, c_indx_y_min+1) = cs_VY(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VZ(i, c_indx_y_min+1) = cs_VZ(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WX(i, c_indx_y_min+1) = cs_WX(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WY(i, c_indx_y_min+1) = cs_WY(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WZ(i, c_indx_y_min+1) = cs_WZ(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VXVY(i, c_indx_y_min+1) = cs_VXVY(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VXVZ(i, c_indx_y_min+1) = cs_VXVZ(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VYVZ(i, c_indx_y_min+1) = cs_VYVZ(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QX(i, c_indx_y_min+1) = cs_QX(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QY(i, c_indx_y_min+1) = cs_QY(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QZ(i, c_indx_y_min+1) = cs_QZ(i, c_indx_y_min+1) + rbufer(i+pos1+1)
        END DO

     END IF

     IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above moments in the vertical line under the top line
        CALL MPI_RECV(rbufer, 13*n3, MPI_REAL, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
        pos1 = -c_indx_x_min
        DO i = c_indx_x_min, c_indx_x_max
           cs_N(i, c_indx_y_max-1) = cs_N(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VX(i, c_indx_y_max-1) = cs_VX(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VY(i, c_indx_y_max-1) = cs_VY(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VZ(i, c_indx_y_max-1) = cs_VZ(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WX(i, c_indx_y_max-1) = cs_WX(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WY(i, c_indx_y_max-1) = cs_WY(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_WZ(i, c_indx_y_max-1) = cs_WZ(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VXVY(i, c_indx_y_max-1) = cs_VXVY(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VXVZ(i, c_indx_y_max-1) = cs_VXVZ(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_VYVZ(i, c_indx_y_max-1) = cs_VYVZ(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO

        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QX(i, c_indx_y_max-1) = cs_QX(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QY(i, c_indx_y_max-1) = cs_QY(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO
        pos1 = pos1+n3
        DO i = c_indx_x_min, c_indx_x_max
           cs_QZ(i, c_indx_y_max-1) = cs_QZ(i, c_indx_y_max-1) + rbufer(i+pos1+1)
        END DO

     END IF

     IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up moments in the top edge
        pos1=1
        pos2=n3
        rbufer(pos1:pos2) = cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_max)

        CALL MPI_SEND(rbufer, 13*n3, MPI_REAL, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down moments in the bottom edge
        pos1=1
        pos2=n3
        rbufer(pos1:pos2) = cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_min)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_min)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_min)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_min)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_min)

        CALL MPI_SEND(rbufer, 13*n3, MPI_REAL, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

  END IF   !###   IF (WHITE_CLUSTER) THEN

  IF (periodic_boundary_X_left) THEN
     
     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:(13*n1)), STAT=ALLOC_ERR)

! ## 1 ## send left moments in the left edge PLUS ONE
     pos1=1
     pos2=n1
     rbufer(pos1:pos2) = cs_N(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_WX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_WY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_WZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_QX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_QY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_QZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

     CALL MPI_SEND(rbufer, 13*n1, MPI_REAL, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 

! ## 2 ## receive from left moments in the vertical line AT the left edge
     CALL MPI_RECV(rbufer, 13*n1, MPI_REAL, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)

     pos1 = -c_indx_y_min
     DO j = c_indx_y_min, c_indx_y_max
        cs_N(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO

     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VX(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VY(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VZ(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO

     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_WX(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_WY(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_WZ(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO

     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VXVY(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VXVZ(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VYVZ(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO

     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_QX(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_QY(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_QZ(c_indx_x_min, j) = rbufer(j+pos1+1)
     END DO

  ELSE IF (periodic_boundary_X_right) THEN   !###   IF (periodic_boundary_X_left) THEN

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:(13*n1)), STAT=ALLOC_ERR)

! ## 1 ## receive from right moments in the vertical line AT the right edge
     CALL MPI_RECV(rbufer, 13*n1, MPI_REAL, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)

     pos1 = -c_indx_y_min
     DO j = c_indx_y_min, c_indx_y_max
        cs_N(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO

     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VX(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VY(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VZ(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO

     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_WX(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_WY(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_WZ(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO
     
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VXVY(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VXVZ(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_VYVZ(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO

     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_QX(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_QY(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO
     pos1 = pos1+n1
     DO j = c_indx_y_min, c_indx_y_max
        cs_QZ(c_indx_x_max, j) = rbufer(j+pos1+1)
     END DO
     
! ## 2 ## send right moments in the right edge MINUS ONE
     pos1=1
     pos2=n1
     rbufer(pos1:pos2) = cs_N(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_WX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_WY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_WZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VXVY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_QX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_QY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
     pos1=pos2+1
     pos2=pos2+n1
     rbufer(pos1:pos2) = cs_QZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

     CALL MPI_SEND(rbufer, 13*n1, MPI_REAL, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 

  END IF   !###   IF (periodic_boundary_X_left) THEN

  IF (periodic_boundary_Y_below) THEN

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:13*n3), STAT=ALLOC_ERR)

! ## 1 ## send down moments in the bottom edge PLUS ONE
     pos1=1
     pos2=n3
     rbufer(pos1:pos2) = cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

     CALL MPI_SEND(rbufer, 13*n3, MPI_REAL, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
 
! ## 2 ## receive from below moments in the bottom line
     CALL MPI_RECV(rbufer, 13*n3, MPI_REAL, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)

     pos1 = -c_indx_x_min
     DO i = c_indx_x_min, c_indx_x_max
        cs_N(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO

     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VX(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VY(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VZ(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO

     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_WX(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_WY(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_WZ(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO

     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VXVY(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VXVZ(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VYVZ(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO

     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_QX(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_QY(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_QZ(i, c_indx_y_min) = rbufer(i+pos1+1)
     END DO

  ELSE IF (periodic_boundary_Y_above) THEN   !###   IF (periodic_boundary_Y_below) THEN

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:13*n3), STAT=ALLOC_ERR)

! ## 1 ## receive from above moments in the top line
     CALL MPI_RECV(rbufer, 13*n3, MPI_REAL, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)

     pos1 = -c_indx_x_min
     DO i = c_indx_x_min, c_indx_x_max
        cs_N(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO

     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VX(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VY(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VZ(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO

     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_WX(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_WY(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_WZ(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO

     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VXVY(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VXVZ(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_VYVZ(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO

     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_QX(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_QY(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO
     pos1 = pos1+n3
     DO i = c_indx_x_min, c_indx_x_max
        cs_QZ(i, c_indx_y_max) = rbufer(i+pos1+1)
     END DO

! ## 2 ## send up moments in the top edge MINUS ONE
     pos1=1
     pos2=n3
     rbufer(pos1:pos2) = cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
     pos1=pos2+1
     pos2=pos2+n3
     rbufer(pos1:pos2) = cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

     CALL MPI_SEND(rbufer, 13*n3, MPI_REAL, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
 
  END IF   !###   IF (periodic_boundary_Y_below) THEN

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer,  STAT=ALLOC_ERR)

END SUBROUTINE SYNCHRONIZE_MOMENTS_IN_OVERLAP_NODES

!----------------------------------------------------
!
SUBROUTINE ADJUST_DENSITY_AT_WALL_BOUNDARIES

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE Snapshots

  IMPLICIT NONE

  INTEGER i, j

  IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! special case of self-connected X-periodic cluster
     IF (Rank_of_master_above.LT.0) THEN
        DO i = c_indx_x_min, c_indx_x_max
           cs_N(i, c_indx_y_max) = 2.0 * cs_N(i, c_indx_y_max)
        END DO
     END IF
  
     IF (Rank_of_master_below.LT.0) THEN
        DO i = c_indx_x_min, c_indx_x_max
           cs_N(i, c_indx_y_min) = 2.0 * cs_N(i, c_indx_y_min)
        END DO
     END IF

  ELSE
 
     IF (Rank_of_master_left.LT.0) THEN
        DO j = c_indx_y_min+1, c_indx_y_max-1
           cs_N(c_indx_x_min, j) = 2.0 * cs_N(c_indx_x_min, j)
        END DO
     END IF

     IF (Rank_of_master_right.LT.0) THEN
        DO j = c_indx_y_min+1, c_indx_y_max-1
           cs_N(c_indx_x_max, j) = 2.0 * cs_N(c_indx_x_max, j)
        END DO
     END IF
  
     IF (Rank_of_master_above.LT.0) THEN
        DO i = c_indx_x_min+1, c_indx_x_max-1
           cs_N(i, c_indx_y_max) = 2.0 * cs_N(i, c_indx_y_max)
        END DO
     END IF
  
     IF (Rank_of_master_below.LT.0) THEN
        DO i = c_indx_x_min+1, c_indx_x_max-1
           cs_N(i, c_indx_y_min) = 2.0 * cs_N(i, c_indx_y_min)
        END DO
     END IF

     SELECT CASE (c_left_bottom_corner_type)
        CASE (SURROUNDED_BY_WALL)
           cs_N(c_indx_x_min, c_indx_y_min) = 4.0 * cs_N(c_indx_x_min, c_indx_y_min) 
        CASE (FLAT_WALL_LEFT)
           cs_N(c_indx_x_min, c_indx_y_min) = 2.0 * cs_N(c_indx_x_min, c_indx_y_min) 
        CASE (FLAT_WALL_BELOW)
           cs_N(c_indx_x_min, c_indx_y_min) = 2.0 * cs_N(c_indx_x_min, c_indx_y_min) 
        CASE (EMPTY_CORNER_WALL_LEFT)
           cs_N(c_indx_x_min, c_indx_y_min+1) = 0.66666666666 * cs_N(c_indx_x_min, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
        CASE (EMPTY_CORNER_WALL_BELOW)
           cs_N(c_indx_x_min+1, c_indx_y_min) = 0.66666666666 * cs_N(c_indx_x_min+1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
     END SELECT
     
     SELECT CASE (c_left_top_corner_type)
        CASE (SURROUNDED_BY_WALL)
           cs_N(c_indx_x_min, c_indx_y_max) = 4.0 * cs_N(c_indx_x_min, c_indx_y_max) 
        CASE (FLAT_WALL_LEFT)
           cs_N(c_indx_x_min, c_indx_y_max) = 2.0 * cs_N(c_indx_x_min, c_indx_y_max) 
        CASE (FLAT_WALL_ABOVE)
           cs_N(c_indx_x_min, c_indx_y_max) = 2.0 * cs_N(c_indx_x_min, c_indx_y_max) 
        CASE (EMPTY_CORNER_WALL_LEFT)
           cs_N(c_indx_x_min, c_indx_y_max-1) = 0.66666666666 * cs_N(c_indx_x_min, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
        CASE (EMPTY_CORNER_WALL_ABOVE)
           cs_N(c_indx_x_min+1, c_indx_y_max) = 0.66666666666 * cs_N(c_indx_x_min+1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
     END SELECT

     SELECT CASE (c_right_bottom_corner_type)
        CASE (SURROUNDED_BY_WALL)
           cs_N(c_indx_x_max, c_indx_y_min) = 4.0 * cs_N(c_indx_x_max, c_indx_y_min) 
        CASE (FLAT_WALL_RIGHT)
           cs_N(c_indx_x_max, c_indx_y_min) = 2.0 * cs_N(c_indx_x_max, c_indx_y_min) 
        CASE (FLAT_WALL_BELOW)
           cs_N(c_indx_x_max, c_indx_y_min) = 2.0 * cs_N(c_indx_x_max, c_indx_y_min) 
        CASE (EMPTY_CORNER_WALL_RIGHT)
           cs_N(c_indx_x_max, c_indx_y_min+1) = 0.66666666666 * cs_N(c_indx_x_max, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
        CASE (EMPTY_CORNER_WALL_BELOW)
           cs_N(c_indx_x_max-1, c_indx_y_min) = 0.66666666666 * cs_N(c_indx_x_max-1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
     END SELECT
     
     SELECT CASE (c_right_top_corner_type)
        CASE (SURROUNDED_BY_WALL)
           cs_N(c_indx_x_max, c_indx_y_max) = 4.0 * cs_N(c_indx_x_max, c_indx_y_max) 
        CASE (FLAT_WALL_RIGHT)
           cs_N(c_indx_x_max, c_indx_y_max) = 2.0 * cs_N(c_indx_x_max, c_indx_y_max) 
        CASE (FLAT_WALL_ABOVE)
           cs_N(c_indx_x_max, c_indx_y_max) = 2.0 * cs_N(c_indx_x_max, c_indx_y_max) 
        CASE (EMPTY_CORNER_WALL_RIGHT)
           cs_N(c_indx_x_max, c_indx_y_max-1) = 0.66666666666 * cs_N(c_indx_x_max, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
        CASE (EMPTY_CORNER_WALL_ABOVE)
           cs_N(c_indx_x_max-1, c_indx_y_max) = 0.66666666666 * cs_N(c_indx_x_max-1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
     END SELECT

  END IF

END SUBROUTINE ADJUST_DENSITY_AT_WALL_BOUNDARIES

