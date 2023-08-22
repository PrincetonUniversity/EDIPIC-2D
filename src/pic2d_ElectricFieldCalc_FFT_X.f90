!------------------------------
!
SUBROUTINE SOLVE_POISSON_FFTX_LINSYSY

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries, ONLY : indx_x_min, indx_x_max, indx_y_min, indx_y_max
  USE ParallelFFTX
  USE Diagnostics, ONLY : Save_probes_e_data_T_cntr, N_of_probes_cluster, List_of_probes_cluster, Probe_position, probe_F_cluster
  USE SetupValues, ONLY : ht_soft_grid_requested, ht_grid_requested, grid_j, F_grid

  use mpi

  IMPLICIT NONE

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR
  INTEGER bufsize

!  REAL(8), ALLOCATABLE :: c_rho(:,:)
!  REAL(8), ALLOCATABLE :: c_phi(:,:)

  REAL(8) factor_rho   !
  REAL(8) factor_rho_2 !
  REAL(8) invfactor    !

  REAL(8), ALLOCATABLE :: xband(:,:)
  REAL(8), ALLOCATABLE :: yband(:,:)

  INTEGER i, j, k, n, n2, jj, jjstart
  INTEGER m, nwo
  INTEGER xsize, ysize
  INTEGER pos1, pos2, pos
  INTEGER shiftpos1, shiftpos2

  INTEGER indx_x_begin, indx_x_end
  INTEGER indx_y_begin, indx_y_end

  INTEGER jre, jim
  REAL(8), ALLOCATABLE :: rhs(:)

real(8) rho_vs_ij       ! temporary, prescribed charge density

character(18) proc_filename     ! abcdefghi_NNNN.dat
                                ! ----x----I----x---

  INTEGER npc, npa

  REAL(8), ALLOCATABLE :: corr_phi(:)     ! correction applied to get desired value of the average potential in the injection plane
                                          ! it is used only if ht_soft_grid_requested = .TRUE.

  REAL(8) myavg_phi, avg_phi

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface

character(10) myfilename ! A_NNNN.dat

  ALLOCATE(rbufer(1:maxbufsize), STAT = ALLOC_ERR)

!  factor_rho = 1.0_8 !### RESTORE ### 0.25_8 / DBLE(N_of_particles_cell)
  factor_rho   = 1.0_8 / DBLE(N_of_particles_cell)                                 !????????? was 1/4 here ???????????
!  factor_rho_2 = 2.0_8 / DBLE(N_of_particles_cell)                                 !????????? was 1/4 here ???????????

  IF (cluster_rank_key.EQ.0) THEN

! note that now c_rho and c_phi are global, 
! here c_rho contains electron density first, then it becomes the full charge density
! c_phi is  

! calculate the combined (electron+ion) charge density
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           c_rho(i,j) = factor_rho * (c_rho_i(i,j) - c_rho(i,j))   !### CHECK SIGN ###!
        END DO
     END DO

! account for the grid potential if it is included
     IF (ht_grid_requested) THEN
        IF ((grid_j.GE.c_indx_y_min).AND.(grid_j.LE.c_indx_y_max)) THEN
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i,j) = F_grid
           END DO
        END IF
     END IF

! note that in each cluster the surface charge at i=c_indx_x_min/max is nullified due to overlapping
! correct values of the surface charge density at i=0 will be implemented at later stage
! when the X-bands of charge density will be created 

     DO n = 1, c_N_of_local_object_parts_below
! if there is a wall below the cluster        
        m = c_index_of_local_object_part_below(n)
        nwo = c_local_object_part(m)%object_number
        IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
! account for surface charge density on dielectric
           DO i = c_local_object_part(m)%istart, c_local_object_part(m)%iend
              c_rho(i,c_indx_y_min) = 0.5_8 * c_rho(i,c_indx_y_min) + factor_rho * c_local_object_part(m)%surface_charge(i)
           END DO
        ELSE IF (whole_object(nwo)%object_type.EQ.VACUUM_GAP) THEN
! account for the potential in vacuum gap
           DO i = c_local_object_part(m)%istart, c_local_object_part(m)%iend
              c_rho(i,c_indx_y_min) = whole_object(nwo)%phi_profile(i)
           END DO
        ELSE IF (whole_object(nwo)%object_type.EQ.METAL_WALL) THEN
! account for the potential of metal wall
           DO i = c_local_object_part(m)%istart, c_local_object_part(m)%iend
              c_rho(i,c_indx_y_min) = whole_object(nwo)%phi
           END DO
!print '("proc ",i4,"sets bottom object ",i2," potential to ",f10.2," V")', Rank_of_process, nwo, whole_object(nwo)%phi * F_scale_V
        END IF
     END DO

     DO n = 1, c_N_of_local_object_parts_above
! if there is a wall above the cluster        
        m = c_index_of_local_object_part_above(n)
        nwo = c_local_object_part(m)%object_number
        IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
! account for surface charge density on dielectric
           DO i = c_local_object_part(m)%istart, c_local_object_part(m)%iend
              c_rho(i,c_indx_y_max) = 0.5_8 * c_rho(i,c_indx_y_max) + factor_rho * c_local_object_part(m)%surface_charge(i)
           END DO
        ELSE IF (whole_object(nwo)%object_type.EQ.VACUUM_GAP) THEN
! account for the potential in vacuum gap
            DO i = c_local_object_part(m)%istart, c_local_object_part(m)%iend
              c_rho(i,c_indx_y_max) = whole_object(nwo)%phi_profile(i)
           END DO
        ELSE IF (whole_object(nwo)%object_type.EQ.METAL_WALL) THEN
! account for the potential of metal wall
           DO i = c_local_object_part(m)%istart, c_local_object_part(m)%iend
              c_rho(i,c_indx_y_max) = whole_object(nwo)%phi
           END DO
!print '("proc ",i4,"sets top object ",i2," potential to ",f10.2," V")', Rank_of_process, nwo, whole_object(nwo)%phi * F_scale_V
        END IF
     END DO

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," passed A")', Rank_of_process

! master processes create X-bands of charge density (RHS of the Poisson's equation) ---------------------------------------------------------

  IF (cluster_rank_key.EQ.0) THEN

     ALLOCATE(xband(0:global_maximal_i, fftx_band_jmin:fftx_band_jmax), STAT=ALLOC_ERR)

! include its own contribution, exclude overlapping for now, see below
     indx_x_begin = c_indx_x_min+1
     indx_x_end   = c_indx_x_max-1
     DO j = fftx_band_jmin, fftx_band_jmax
        xband(indx_x_begin:indx_x_end, j) = c_rho(indx_x_begin:indx_x_end, j)
     END DO

     DO k = 1, N_of_comm_steps_X

        IF (comm_FFTX_step(k)%rank.LT.0) CYCLE
        
        IF (Rank_of_process.LT.comm_FFTX_step(k)%rank) THEN
! send
           indx_x_begin = c_indx_x_min+1                               ! exclude overlapping everywhere, including the end clusters
           indx_x_end   = c_indx_x_max-1                               ! because the density at i=0 may not be fully defined
                                                                       ! while the density at i=N+1=maximal_i is not necessary
                                                                       ! we set rho(0)=rho(N) later, already in the bands
                                                                       ! this also allows to apply correct surface charge density in i=0
                                                                       ! since it is stored not in the cluster with c_indx_x_min=0
                                                                       ! but in the cluster with c_indx_x_max=global_maximal_i in i=global_maximal_i-1
           xsize = indx_x_end - indx_x_begin + 1
           bufsize = (comm_FFTX_step(k)%fftx_band_jmax - comm_FFTX_step(k)%fftx_band_jmin + 1) * xsize
           
           pos2 = 0
           DO j = comm_FFTX_step(k)%fftx_band_jmin, comm_FFTX_step(k)%fftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              rbufer(pos1:pos2) = c_rho(indx_x_begin:indx_x_end, j)
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

! receive
           indx_x_begin = comm_FFTX_step(k)%c_indx_x_min+1
           indx_x_end   = comm_FFTX_step(k)%c_indx_x_max-1
           xsize = indx_x_end - indx_x_begin + 1
           bufsize = (fftx_band_jmax - fftx_band_jmin + 1) * xsize

           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, comm_FFTX_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)

           pos2 = 0
           DO j = fftx_band_jmin, fftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              xband(indx_x_begin:indx_x_end, j) = rbufer(pos1:pos2)
           END DO

        ELSE
! receive
           indx_x_begin = comm_FFTX_step(k)%c_indx_x_min+1
           indx_x_end   = comm_FFTX_step(k)%c_indx_x_max-1
           xsize = indx_x_end - indx_x_begin + 1
           bufsize = (fftx_band_jmax - fftx_band_jmin + 1) * xsize

           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, comm_FFTX_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)

           pos2 = 0
           DO j = fftx_band_jmin, fftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              xband(indx_x_begin:indx_x_end, j) = rbufer(pos1:pos2)
           END DO
! send
           indx_x_begin = c_indx_x_min+1
           indx_x_end   = c_indx_x_max-1      
           xsize = indx_x_end - indx_x_begin + 1
           bufsize = (comm_FFTX_step(k)%fftx_band_jmax - comm_FFTX_step(k)%fftx_band_jmin + 1) * xsize
           
           pos2 = 0
           DO j = comm_FFTX_step(k)%fftx_band_jmin, comm_FFTX_step(k)%fftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              rbufer(pos1:pos2) = c_rho(indx_x_begin:indx_x_end, j)
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 


        END IF

     END DO

  ELSE

     ALLOCATE(xband(0:global_maximal_i, fftx_strip_jmin:fftx_strip_jmax), STAT=ALLOC_ERR)

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," passed B")', Rank_of_process

! now each master distributes the charge density between associated field calculators ---------------------------------------------------------

  IF (cluster_rank_key.EQ.0) THEN
! include periodicity, this part must work even when the cluster has only one process 
     DO j = fftx_band_jmin, fftx_band_jmax
        xband(0, j) = xband(global_maximal_i-1, j)
     END DO

     DO k = 2, cluster_N_blocks
        bufsize = (global_maximal_i - 1) * (field_calculator(k)%fftx_strip_jmax - field_calculator(k)%fftx_strip_jmin + 1)
        pos2=0
        DO j = field_calculator(k)%fftx_strip_jmin, field_calculator(k)%fftx_strip_jmax
           pos1 = pos2 + 1
           pos2 = pos2 + global_maximal_i - 1
           rbufer(pos1:pos2) = xband(0:global_maximal_i-2, j)
        END DO
        CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 
     END DO

  ELSE

     bufsize = (global_maximal_i - 1) * (fftx_strip_jmax - fftx_strip_jmin + 1)
       
     CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_master, field_master, MPI_COMM_WORLD, stattus, ierr)
        
     pos2 = 0
     DO j = fftx_strip_jmin, fftx_strip_jmax
        pos1 = pos2 + 1
        pos2 = pos2 + global_maximal_i - 1
        xband(0:global_maximal_i-2, j) = rbufer(pos1:pos2)
     END DO

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," passed C")', Rank_of_process

! do fft transformation of charge density along X ----------------------------------------------------------------------------------------------

!proc_filename = 'xbdbeffft_NNNN.dat'
!               ! ----x----I----x---
!proc_filename(11:14) = convert_int_to_txt_string(Rank_of_process, 4)
!open (9, file = proc_filename)
!do j = fftx_strip_jmin, fftx_strip_jmax
!   do i = 0, global_maximal_i-2
!      write (9, '(2x,i5,2x,i5,2x,e14.7)') i, j, xband(i,j)
!   end do
!   write (9, '(" ")')
!end do
!print '("file ",A18," is ready")', proc_filename
!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 


  DO j = fftx_strip_jmin, fftx_strip_jmax

     CALL REAL8FT(xband(0:global_maximal_i-2,j), global_maximal_i-1, 1)

! this sub returns complex Fourier amplitudes with n = 0, 1, ..., N/2 
! each complex amplitude occupies two elements of aray xband(0:global_maximal_i-2,j), therefore total number of elements is (N/2+1)*2
! but amplitudes n=0 and n=N/2 are real and they are in the first and second elements of the array 
! so that total number of elements is (N/2+1)*2-2=N
! since array xband is longer than necessary (0:N-1) by two real numbers (0:N+1), here N is global_maximal_i-1
! one can place the complex value of the amplitude n=N/2 in the last two elements of the array

     xband(global_maximal_i-1,j) = xband(1,j)  ! real      ! parts of amplitude n=N/2
     xband(global_maximal_i,  j) = 0.0_8       ! imaginary !

     xband(1,j) = 0.0_8   ! set zero imaginary part of amplitude n=0

! the benefit of this approach is that all Fourier components will be treated in the same manner 
! when the linear system of equations along Y will be solved

  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!proc_filename = 'xbdaftfft_NNNN.dat'
!               ! ----x----I----x---
!proc_filename(11:14) = convert_int_to_txt_string(Rank_of_process, 4)
!open (9, file = proc_filename)
!do j = fftx_strip_jmin, fftx_strip_jmax
!   do i = 0, global_maximal_i-2
!      write (9, '(2x,i5,2x,i5,2x,e14.7)') i, j, xband(i,j)
!   end do
!   write (9, '(" ")')
!end do
!print '("file ",A18," is ready")', proc_filename
!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! each master collects Fourier components from the associated field calculators ---------------------------------------------------------

  IF (cluster_rank_key.EQ.0) THEN

     DO k = 2, cluster_N_blocks
        bufsize = (global_maximal_i + 1) * (field_calculator(k)%fftx_strip_jmax - field_calculator(k)%fftx_strip_jmin + 1)
        CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)
        pos2=0
        DO j = field_calculator(k)%fftx_strip_jmin, field_calculator(k)%fftx_strip_jmax
           pos1 = pos2 + 1
           pos2 = pos2 + global_maximal_i + 1
           xband(0:global_maximal_i, j) = rbufer(pos1:pos2)
        END DO
     END DO

  ELSE

     bufsize = (global_maximal_i + 1) * (fftx_strip_jmax - fftx_strip_jmin + 1)

     pos2 = 0
     DO j = fftx_strip_jmin, fftx_strip_jmax
        pos1 = pos2 + 1
        pos2 = pos2 + global_maximal_i + 1
        rbufer(pos1:pos2) = xband(0:global_maximal_i, j)
     END DO
        
     CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_master, Rank_of_process, MPI_COMM_WORLD, ierr) 

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," passed D")', Rank_of_process

! assemble density Fourier components in master processes ------------------------------------------------------------------------------

  IF (cluster_rank_key.EQ.0) THEN

! include the process itself
     indx_x_begin = c_indx_x_min
     indx_x_end   = c_indx_x_max-2
     IF (indx_x_end.EQ.global_maximal_i-2) indx_x_end = global_maximal_i
     DO j = fftx_band_jmin, fftx_band_jmax
        c_phi(indx_x_begin:indx_x_end, j) = xband(indx_x_begin:indx_x_end, j)
     END DO

     DO k = 1, N_of_comm_steps_X

        IF (comm_FFTX_step(k)%rank.LT.0) CYCLE
        
        IF (Rank_of_process.LT.comm_FFTX_step(k)%rank) THEN

! send
           indx_x_begin = comm_FFTX_step(k)%c_indx_x_min         ! to make sure that we have pairs corresponding to
           indx_x_end   = comm_FFTX_step(k)%c_indx_x_max - 2     ! real and imaginary parts of fourier components
           IF (indx_x_end.EQ.global_maximal_i-2) indx_x_end = global_maximal_i  ! to include component n=N/2 as a complex number
           xsize = indx_x_end - indx_x_begin + 1

           bufsize = (fftx_band_jmax - fftx_band_jmin + 1) * xsize

           pos2 = 0
           DO j = fftx_band_jmin, fftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              rbufer(pos1:pos2) = xband(indx_x_begin:indx_x_end, j)
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

! receive
           indx_x_begin = c_indx_x_min
           indx_x_end   = c_indx_x_max-2
           IF (indx_x_end.EQ.global_maximal_i-2) indx_x_end = global_maximal_i
           xsize = indx_x_end - indx_x_begin + 1

           bufsize = (comm_FFTX_step(k)%fftx_band_jmax - comm_FFTX_step(k)%fftx_band_jmin + 1) * xsize
 
           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, comm_FFTX_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)
          
           pos2 = 0
           DO j = comm_FFTX_step(k)%fftx_band_jmin, comm_FFTX_step(k)%fftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              c_phi(indx_x_begin:indx_x_end, j) = rbufer(pos1:pos2)
           END DO

        ELSE
! receive
           indx_x_begin = c_indx_x_min
           indx_x_end   = c_indx_x_max-2
           IF (indx_x_end.EQ.global_maximal_i-2) indx_x_end = global_maximal_i
           xsize = indx_x_end - indx_x_begin + 1
           bufsize = (comm_FFTX_step(k)%fftx_band_jmax - comm_FFTX_step(k)%fftx_band_jmin + 1) * xsize
 
           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, comm_FFTX_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)
          
           pos2 = 0
           DO j = comm_FFTX_step(k)%fftx_band_jmin, comm_FFTX_step(k)%fftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              c_phi(indx_x_begin:indx_x_end, j) = rbufer(pos1:pos2)
           END DO


! send
           indx_x_begin = comm_FFTX_step(k)%c_indx_x_min         ! to make sure that we have pairs corresponding to
           indx_x_end   = comm_FFTX_step(k)%c_indx_x_max - 2     ! real and imaginary parts of fourier components
           IF (indx_x_end.EQ.global_maximal_i-2) indx_x_end = global_maximal_i  ! to include component n=N/2 as a complex number
           xsize = indx_x_end - indx_x_begin + 1

           bufsize = (fftx_band_jmax - fftx_band_jmin + 1) * xsize

           pos2 = 0
           DO j = fftx_band_jmin, fftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              rbufer(pos1:pos2) = xband(indx_x_begin:indx_x_end, j)
           END DO

           CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

        END IF

     END DO

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," passed E")', Rank_of_process 

! re-distribute density Fourier components into bands along Y

! yband's will represent 0:global_maximal_j COMPLEX values for each value of n
! where n changes from sysy_band_nmin to sysy_band_nmax
! that is size of the line in yband is 2*(global_maximal_j+1)
! while size of the column is sysy_band_nmax - sysy_band_nmin + 1
! note:
! 0:global_maximal_j x2 => 0:2*global_maximal_j+1

  IF (cluster_rank_key.EQ.0) THEN

     ALLOCATE(yband(0:(2*global_maximal_j+1), sysy_band_nmin:sysy_band_nmax), STAT=ALLOC_ERR)   ! even/odd j's (first index) are real/imaginary parts, y-boundaries are excluded

! include its own contribution
     IF (c_indx_y_min.GT.0) THEN
        indx_y_begin = 2*(c_indx_y_min+1)       ! even, real part of the first element in a column
        jjstart = c_indx_y_min+1
     ELSE
        indx_y_begin = 0
        jjstart = 0
     END IF

     IF (c_indx_y_max.LT.global_maximal_j) THEN
        indx_y_end = 2*(c_indx_y_max-1)         ! even, real part of the last element in a column
     ELSE
        indx_y_end = 2*global_maximal_j
     END IF

     DO n = sysy_band_nmin, sysy_band_nmax
        n2 = n+n
        jj = jjstart 
        DO j = indx_y_begin, indx_y_end, 2
           yband(j, n)   = c_phi(n2,   jj)
           yband(j+1, n) = c_phi(n2+1, jj)
           jj = jj + 1
        END DO
     END DO

     DO k = 1, N_of_comm_steps_Y

        IF (comm_SYSY_step(k)%rank.LT.0) CYCLE
        
        IF (Rank_of_process.LT.comm_SYSY_step(k)%rank) THEN
! send
           IF (c_indx_y_min.GT.0) THEN
              indx_y_begin = c_indx_y_min+1
           ELSE
              indx_y_begin = 0
           END IF
           IF (c_indx_y_max.LT.global_maximal_j) THEN
              indx_y_end   = c_indx_y_max-1
           ELSE
              indx_y_end = global_maximal_j
           END IF
           ysize = 2 * (indx_y_end - indx_y_begin + 1)

           bufsize = (comm_SYSY_step(k)%sysy_band_nmax - comm_SYSY_step(k)%sysy_band_nmin + 1) * ysize
           
           pos2 = 0
           DO n = comm_SYSY_step(k)%sysy_band_nmin, comm_SYSY_step(k)%sysy_band_nmax
              pos1 = pos2 + 1
              pos2 = pos2 + ysize
              n2 = n+n
              jj = indx_y_begin
              DO j = pos1, pos2-1, 2
                 rbufer(j)   = c_phi(n2,  jj)
                 rbufer(j+1) = c_phi(n2+1,jj)
                 jj = jj + 1
              END DO
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_SYSY_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 
! receive
           IF (comm_SYSY_step(k)%c_indx_y_min.GT.0) THEN
              indx_y_begin = comm_SYSY_step(k)%c_indx_y_min+1
           ELSE
              indx_y_begin = 0
           END IF
           IF (comm_SYSY_step(k)%c_indx_y_max.LT.global_maximal_j) THEN
              indx_y_end   = comm_SYSY_step(k)%c_indx_y_max-1
           ELSE
              indx_y_end = global_maximal_j
           END IF

           bufsize = (sysy_band_nmax - sysy_band_nmin + 1) * 2 * (indx_y_end - indx_y_begin + 1)

           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_SYSY_step(k)%rank, comm_SYSY_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)

           shiftpos1 = indx_y_begin+indx_y_begin       ! if indx_y_begin=0 , shiftpos1=0 which is lower limit of first index of yband
           shiftpos2 = indx_y_end+indx_y_end+1         ! if indx_y_end=global_maximal_j , shiftpos2=2*global_maximal_j+1 which is upper limit of first index of yband
           ysize = shiftpos2-shiftpos1+1

           pos2 = 0
           DO n = sysy_band_nmin, sysy_band_nmax
              pos1 = pos2 + 1
              pos2 = pos2 + ysize
              yband(shiftpos1:shiftpos2,n) = rbufer(pos1:pos2)
           END DO

        ELSE
! receive
           IF (comm_SYSY_step(k)%c_indx_y_min.GT.0) THEN
              indx_y_begin = comm_SYSY_step(k)%c_indx_y_min+1
           ELSE
              indx_y_begin = 0
           END IF
           IF (comm_SYSY_step(k)%c_indx_y_max.LT.global_maximal_j) THEN
              indx_y_end   = comm_SYSY_step(k)%c_indx_y_max-1
           ELSE
              indx_y_end = global_maximal_j
           END IF

           bufsize = (sysy_band_nmax - sysy_band_nmin + 1) * 2 * (indx_y_end - indx_y_begin + 1)

           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_SYSY_step(k)%rank, comm_SYSY_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)

           shiftpos1 = indx_y_begin+indx_y_begin       ! if indx_y_begin=0 , shiftpos1=0 which is lower limit of first index of yband
           shiftpos2 = indx_y_end+indx_y_end+1         ! if indx_y_end=global_maximal_j , shiftpos2=2*global_maximal_j+1 which is upper limit of first index of yband
           ysize = shiftpos2-shiftpos1+1

           pos2 = 0
           DO n = sysy_band_nmin, sysy_band_nmax
              pos1 = pos2 + 1
              pos2 = pos2 + ysize
              yband(shiftpos1:shiftpos2,n) = rbufer(pos1:pos2)
           END DO
! send
           IF (c_indx_y_min.GT.0) THEN
              indx_y_begin = c_indx_y_min+1
           ELSE
              indx_y_begin = 0
           END IF
           IF (c_indx_y_max.LT.global_maximal_j) THEN
              indx_y_end   = c_indx_y_max-1
           ELSE
              indx_y_end = global_maximal_j
           END IF
           ysize = 2 * (indx_y_end - indx_y_begin + 1)

           bufsize = (comm_SYSY_step(k)%sysy_band_nmax - comm_SYSY_step(k)%sysy_band_nmin + 1) * ysize
           
           pos2 = 0
           DO n = comm_SYSY_step(k)%sysy_band_nmin, comm_SYSY_step(k)%sysy_band_nmax
              pos1 = pos2 + 1
              pos2 = pos2 + ysize
              n2 = n+n
              jj = indx_y_begin
              DO j = pos1, pos2-1, 2
                 rbufer(j)   = c_phi(n2,  jj)
                 rbufer(j+1) = c_phi(n2+1,jj)
                 jj = jj + 1
              END DO
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_SYSY_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

        END IF

     END DO

  ELSE

     ALLOCATE(yband(0:(2*global_maximal_j+1), sysy_strip_nmin:sysy_strip_nmax), STAT=ALLOC_ERR)

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)  

!print '("Process ",i4," passed F")', Rank_of_process

! now each master distributes the Fourier components between associated field calculators

  IF (cluster_rank_key.EQ.0) THEN

     DO k = 2, cluster_N_blocks
        ysize = 2 * (global_maximal_j + 1)
        bufsize = ysize * (field_calculator(k)%sysy_strip_nmax - field_calculator(k)%sysy_strip_nmin + 1)        ! try to use %indx_x_min/max instead
        pos2=0
        DO n = field_calculator(k)%sysy_strip_nmin, field_calculator(k)%sysy_strip_nmax
           pos1 = pos2 + 1
           pos2 = pos2 + ysize
           rbufer(pos1:pos2) = yband(0:(ysize-1), n)
        END DO
        CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 
     END DO

  ELSE

     ysize = 2 * (global_maximal_j + 1)
     bufsize = ysize * (sysy_strip_nmax - sysy_strip_nmin + 1)
        
     CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_master, field_master, MPI_COMM_WORLD, stattus, ierr)
        
     pos2 = 0
     DO n = sysy_strip_nmin, sysy_strip_nmax
        pos1 = pos2 + 1
        pos2 = pos2 + ysize
        yband(0:(ysize-1), n) = rbufer(pos1:pos2)
     END DO

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," passed G")', Rank_of_process

! solve linear system of equations for Fourier components of the potential for each column of nodes along Y

  ALLOCATE(rhs(0:(2*global_maximal_j+1)) , STAT = ALLOC_ERR)   !                                      0,1 are Re/Im parts of RHS at boundary y=0
                                                               ! 2*global_maximal_j, 2*global_maximal_j+1 are Re/Im parts of RHS at boundary y=global_maximal_j

  DO n = sysy_strip_nmin, sysy_strip_nmax

! set RHS 
     jre = 0
     jim = 1

     IF (a_eq(0,n).EQ.0.0_8) THEN
! the wall y=0 is a combination of METAL_WALL / VACUUM_GAP boundary objects, the potential is known
        rhs(jre) = yband(jre, n)
        rhs(jim) = yband(jim, n)
     ELSE
! the wall y=0 is DIELECTRIC
        rhs(jre) = -yband(jre, n) * a_eq(0, n)   ! if the wall y=0 is metal, a_eq(0,n)=0 :: this ensures zero rhs
        rhs(jim) = -yband(jim, n) * a_eq(0, n)   !
     END IF

     IF (two_dielectric_walls.AND.(n.EQ.0)) THEN
        rhs(jre) = 0.0_8
        rhs(jim) = 0.0_8
     END IF

     IF (ht_grid_requested) THEN
        DO j = 1, grid_j-1
           jre = jre + 2
           jim = jim + 2
           rhs(jre) = (-yband(jre, n) - rhs(jre-2)) * a_eq(j, n)
           rhs(jim) = (-yband(jim, n) - rhs(jim-2)) * a_eq(j, n)
        END DO
        j = grid_j
        jre = jre + 2
        jim = jim + 2
        rhs(jre) = yband(jre, n)
        rhs(jim) = yband(jim, n)
        DO j = grid_j+1, global_maximal_j-1
           jre = jre + 2
           jim = jim + 2
           rhs(jre) = (-yband(jre, n) - rhs(jre-2)) * a_eq(j, n)
           rhs(jim) = (-yband(jim, n) - rhs(jim-2)) * a_eq(j, n)
        END DO
     ELSE
        DO j = 1, global_maximal_j-1
           jre = jre + 2
           jim = jim + 2
           rhs(jre) = (-yband(jre, n) - rhs(jre-2)) * a_eq(j, n)
           rhs(jim) = (-yband(jim, n) - rhs(jim-2)) * a_eq(j, n)
        END DO
     END IF

     j = global_maximal_j
     jre = jre + 2
     jim = jim + 2

     IF (a_eq(j, n).EQ.0.0_8) THEN
! the wall y=global_maximal_j is a combination of METAL_WALL / VACUUM_GAP boundary objects, the potential is known
        rhs(jre) = yband(jre, n)
        rhs(jim) = yband(jim, n)
     ELSE
! the wall y=global_maximal_j is DIELECTRIC
        rhs(jre) = (-yband(jre, n) - rhs(jre-2)) * a_eq(j, n)
        rhs(jim) = (-yband(jim, n) - rhs(jim-2)) * a_eq(j, n)
     END IF

! solve the equation system

     yband(jre, n) = rhs(jre)
     yband(jim, n) = rhs(jim)

     DO j = global_maximal_j - 1, 0, -1  ! was global_maximal_j - 2, 1, -1 
        jre = jre - 2
        jim = jim - 2
        yband(jre, n) = rhs(jre) - a_eq(j, n) * yband(jre + 2, n)  ! for metal wall at j=0, a_eq=0, which gives yband=rhs
        yband(jim, n) = rhs(jim) - a_eq(j, n) * yband(jim + 2, n)  !
     END DO

  END DO

  DEALLOCATE(rhs, STAT = ALLOC_ERR)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!print '("Process ",i4," passed H")', Rank_of_process

! assemble potential Fourier components [as functions of Y] from field calculators in master processes

  IF (cluster_rank_key.EQ.0) THEN

     DO k = 2, cluster_N_blocks
        ysize = 2 * (global_maximal_j + 1)
        bufsize = ysize * (field_calculator(k)%sysy_strip_nmax - field_calculator(k)%sysy_strip_nmin + 1)
        CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)
        pos2=0
        DO n = field_calculator(k)%sysy_strip_nmin, field_calculator(k)%sysy_strip_nmax
           pos1 = pos2 + 1
           pos2 = pos2 + ysize
           yband(0:(ysize-1), n) = rbufer(pos1:pos2)
        END DO
     END DO

  ELSE

     ysize = 2 * (global_maximal_j + 1)
     bufsize = ysize * (sysy_strip_nmax - sysy_strip_nmin + 1)
     pos2 = 0
     DO n = sysy_strip_nmin, sysy_strip_nmax
        pos1 = pos2 + 1
        pos2 = pos2 + ysize
        rbufer(pos1:pos2) = yband(0:(ysize-1), n)
     END DO
     CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_master, Rank_of_process, MPI_COMM_WORLD, ierr) 

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," passed I")', Rank_of_process

! exchange potential FFT components between master processes in the same column to get full set of these in each master domain

  IF (cluster_rank_key.EQ.0) THEN

!if (T_cntr.eq.10) then
!myfilename = 'B_NNNN.dat'
!myfilename(3:6) = convert_int_to_txt_string(Rank_of_process,4)
!open (51, file = myfilename)
!do n=sysy_band_nmin,sysy_band_nmax
!   do j = 2, 2*global_maximal_j-1,2
!      write (51, '(2x,i4,2x,i4,2x,e14.7)') n+n, j/2, yband(j,n)
!      write (51, '(2x,i4,2x,i4,2x,e14.7)') n+n+1, j/2+1, yband(j+1,n)
!   end do
!   write (51, '(" ")')
!end do
!close (51, status = 'keep')
!print '("process ",i4," :: file ",A10," is ready")', Rank_of_process, myfilename
!end if

!     ALLOCATE(yband(2:(2*global_maximal_j-1), sysy_band_nmin:sysy_band_nmax), STAT=ALLOC_ERR) 

! include the process itself
     indx_y_begin = 2 * c_indx_y_min  !even, bottom line included, was MAX(2, 2 * c_indx_y_min)                         ! even, exclude bottom line j=0, numbering starts with 2 for real at j=1
     indx_y_end = 2 * c_indx_y_max    !even, top line included, was  MIN(2 * c_indx_y_max, 2 * global_maximal_j - 2)    ! even, exclude top line j=maxj, numbering corresponds to the real part of the last complex amplitude

     DO n = sysy_band_nmin, sysy_band_nmax
        n2 = n+n
        jj = c_indx_y_min  !MAX(c_indx_y_min, 1)
        DO j = indx_y_begin, indx_y_end, 2
           c_phi(n2, jj)   = yband(j, n)
           c_phi(n2+1, jj) = yband(j+1, n)
           jj = jj + 1
        END DO
     END DO

     DO k = 1, N_of_comm_steps_Y

        IF (comm_SYSY_step(k)%rank.LT.0) CYCLE
        
        IF (Rank_of_process.LT.comm_SYSY_step(k)%rank) THEN
! send
           indx_y_begin = 2 * comm_SYSY_step(k)%c_indx_y_min  !MAX(2, 2 * comm_SYSY_step(k)%c_indx_y_min)                           ! even, real part of the first complex amplitude
           indx_y_end = 2 * comm_SYSY_step(k)%c_indx_y_max + 1 !MIN(2 * comm_SYSY_step(k)%c_indx_y_max + 1, 2 * global_maximal_j - 1)  ! odd, includes the imaginary part of the last complex amplitude
           ysize = indx_y_end - indx_y_begin + 1

           bufsize = (sysy_band_nmax - sysy_band_nmin + 1) * ysize
           
           pos2 = 0
           DO n = sysy_band_nmin, sysy_band_nmax
              pos1 = pos2 + 1
              pos2 = pos2 + ysize
              rbufer(pos1:pos2) = yband(indx_y_begin:indx_y_end, n)
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_SYSY_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

! receive
           indx_y_begin = c_indx_y_min !MAX(c_indx_y_min, 1) 
           indx_y_end   = c_indx_y_max !MIN(c_indx_y_max, global_maximal_j-1)
           ysize = 2 * (indx_y_end - indx_y_begin + 1)

           bufsize = (comm_SYSY_step(k)%sysy_band_nmax - comm_SYSY_step(k)%sysy_band_nmin + 1) * ysize

           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_SYSY_step(k)%rank, comm_SYSY_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)

           pos2 = 0
           DO n = comm_SYSY_step(k)%sysy_band_nmin, comm_SYSY_step(k)%sysy_band_nmax
              pos1 = pos2 + 1
              pos2 = pos2 + ysize
              n2 = n+n
              jj = indx_y_begin
              DO j = pos1, pos2-1, 2
                 c_phi(n2, jj)   = rbufer(j)
                 c_phi(n2+1, jj) = rbufer(j+1)
                 jj = jj + 1
              END DO
           END DO

        ELSE
! receive
           indx_y_begin = c_indx_y_min !MAX(c_indx_y_min, 1) 
           indx_y_end   = c_indx_y_max !MIN(c_indx_y_max, global_maximal_j-1)
           ysize = 2 * (indx_y_end - indx_y_begin + 1)

           bufsize = (comm_SYSY_step(k)%sysy_band_nmax - comm_SYSY_step(k)%sysy_band_nmin + 1) * ysize

           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_SYSY_step(k)%rank, comm_SYSY_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)

           pos2 = 0
           DO n = comm_SYSY_step(k)%sysy_band_nmin, comm_SYSY_step(k)%sysy_band_nmax
              pos1 = pos2 + 1
              pos2 = pos2 + ysize
              n2 = n+n
              jj = indx_y_begin
              DO j = pos1, pos2-1, 2
                 c_phi(n2, jj)   = rbufer(j)
                 c_phi(n2+1, jj) = rbufer(j+1)
                 jj = jj + 1
              END DO
           END DO
! send
           indx_y_begin = 2 * comm_SYSY_step(k)%c_indx_y_min    !MAX(2, 2 * comm_SYSY_step(k)%c_indx_y_min)                             ! even, exclude bottom line j=0, numbering starts with 0 for real at j=1
           indx_y_end = 2 * comm_SYSY_step(k)%c_indx_y_max + 1  !MIN(2 * comm_SYSY_step(k)%c_indx_y_max + 1, 2 * global_maximal_j - 1)    ! odd, exclude top line j=maxj. numbering ends at 2*global_maximal_j-1 for imaginary at j=max-1
           ysize = indx_y_end - indx_y_begin + 1

           bufsize = (sysy_band_nmax - sysy_band_nmin + 1) * ysize
           
           pos2 = 0
           DO n = sysy_band_nmin, sysy_band_nmax
              pos1 = pos2 + 1
              pos2 = pos2 + ysize
              rbufer(pos1:pos2) = yband(indx_y_begin:indx_y_end, n)
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_SYSY_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

        END IF

     END DO

!if (T_cntr.eq.10) then
!myfilename = 'B_NNNN.dat'
!myfilename(3:6) = convert_int_to_txt_string(Rank_of_process,4)
!open (51, file = myfilename)
!do j = 2, 2*global_maximal_j-2, 2
!   do n=sysy_band_nmin,sysy_band_nmax
!      write (51, '(2x,i4,2x,i4,2x,e14.7)') n+n, j/2, yband(j,n)
!      write (51, '(2x,i4,2x,i4,2x,e14.7)') n+n+1, j/2, yband(j+1,n)
!   end do
!   write (51, '(" ")')
!end do
!close (51, status = 'keep')
!print '("process ",i4," :: file ",A10," is ready")', Rank_of_process, myfilename
!end if

  END IF

  DEALLOCATE(yband, STAT=ALLOC_ERR)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)           

!print '("Process ",i4," passed J")', Rank_of_process

! at this moment, each master process has complex inverse fft components =======================================================================================================================================================================
! for i=c_indx_x_min:c_indx_x_max-2       or (in the last column) i=c_indx_x_min:c_indx_x_max
! corresponding to 
!     n=c_indx_x_min/2:(c_indx_x_max-3)/2 or (in the last column) n=c_indx_x_min/2:(global_maximal_i-1)/2
! and for j=c_indx_y_min:c_indx_y_max everywhere except 
! at the bottom line for j=1:c_indx_y_max and 
! at the top    line for j=c_indx_y_min:global_maximal_j-1
! where the potential at j=0 and j=global_maximal_j is given

! re-distribute potential Fourier components into bands along X

  DEALLOCATE(xband, STAT = ALLOC_ERR)

  IF (cluster_rank_key.EQ.0) THEN

     ALLOCATE(xband(0:global_maximal_i, invfftx_band_jmin:invfftx_band_jmax), STAT=ALLOC_ERR)

! include its own contribution 
     indx_x_begin = c_indx_x_min
     indx_x_end   = c_indx_x_max-2
     IF (c_column.EQ.N_clusters_x) indx_x_end = c_indx_x_max !####################################??????????????????????????
     DO j = invfftx_band_jmin, invfftx_band_jmax
        xband(indx_x_begin:indx_x_end, j) = c_phi(indx_x_begin:indx_x_end, j)
     END DO

     DO k = 1, N_of_comm_steps_X

        IF (comm_FFTX_step(k)%rank.LT.0) CYCLE
        
        IF (Rank_of_process.LT.comm_FFTX_step(k)%rank) THEN
! send
           indx_x_begin = c_indx_x_min
           indx_x_end   = c_indx_x_max-2
           IF (c_column.EQ.N_clusters_x) indx_x_end = c_indx_x_max
           xsize = indx_x_end - indx_x_begin + 1

           bufsize = (comm_FFTX_step(k)%invfftx_band_jmax - comm_FFTX_step(k)%invfftx_band_jmin + 1) * xsize
           
           pos2 = 0
           DO j = comm_FFTX_step(k)%invfftx_band_jmin, comm_FFTX_step(k)%invfftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              rbufer(pos1:pos2) = c_phi(indx_x_begin:indx_x_end, j)
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

! receive
           indx_x_begin = comm_FFTX_step(k)%c_indx_x_min
           indx_x_end   = comm_FFTX_step(k)%c_indx_x_max-2
           IF (indx_x_end.EQ.global_maximal_i-2) indx_x_end = global_maximal_i
           xsize = indx_x_end - indx_x_begin + 1

           bufsize = (invfftx_band_jmax - invfftx_band_jmin + 1) * xsize

           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, comm_FFTX_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)

           pos2 = 0
           DO j = invfftx_band_jmin, invfftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              xband(indx_x_begin:indx_x_end, j) = rbufer(pos1:pos2)
           END DO

        ELSE
! receive
           indx_x_begin = comm_FFTX_step(k)%c_indx_x_min
           indx_x_end   = comm_FFTX_step(k)%c_indx_x_max-2
           IF (indx_x_end.EQ.global_maximal_i-2) indx_x_end = global_maximal_i
           xsize = indx_x_end - indx_x_begin + 1

           bufsize = (invfftx_band_jmax - invfftx_band_jmin + 1) * xsize

           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, comm_FFTX_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)

           pos2 = 0
           DO j = invfftx_band_jmin, invfftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              xband(indx_x_begin:indx_x_end, j) = rbufer(pos1:pos2)
           END DO
! send
           indx_x_begin = c_indx_x_min
           indx_x_end   = c_indx_x_max-2
           IF (c_column.EQ.N_clusters_x) indx_x_end = c_indx_x_max
           xsize = indx_x_end - indx_x_begin + 1

           bufsize = (comm_FFTX_step(k)%invfftx_band_jmax - comm_FFTX_step(k)%invfftx_band_jmin + 1) * xsize
          
           pos2 = 0
           DO j = comm_FFTX_step(k)%invfftx_band_jmin, comm_FFTX_step(k)%invfftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              rbufer(pos1:pos2) = c_phi(indx_x_begin:indx_x_end, j)
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

        END IF

     END DO

  ELSE

     ALLOCATE(xband(0:global_maximal_i, invfftx_strip_jmin:invfftx_strip_jmax), STAT=ALLOC_ERR)

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," passed K")', Rank_of_process

! now each master distributes the fft components of the potential between associated field calculators ---------------------------------------------------------

  IF (cluster_rank_key.EQ.0) THEN

     DO k = 2, cluster_N_blocks
        indx_y_begin = field_calculator(k)%invfftx_strip_jmin
        indx_y_end   = field_calculator(k)%invfftx_strip_jmax
! we expand the range of j in order to exclude Y-exchange of potentials between above/below neighbor clusters
! therefore we have to include endpoints c_indx_x_min and c_indx_x_max which were previously excluded to prevent overlapping
! the lower limit will be always included in the master process
! the upper limit will be included in the last field calculator (see the line below) or, if there are no field calculators, in the master process as well
!        IF (indx_y_end.EQ.(c_indx_y_max-1)) indx_y_end = MIN(c_indx_y_max, global_maximal_j-1)

        bufsize = (global_maximal_i + 1) * (indx_y_end - indx_y_begin + 1)

        pos2=0
        DO j = indx_y_begin, indx_y_end
           pos1 = pos2 + 1
           pos2 = pos2 + global_maximal_i + 1
           rbufer(pos1:pos2) = xband(0:global_maximal_i, j)
        END DO

        CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

     END DO

  ELSE

     bufsize = (global_maximal_i + 1) * (invfftx_strip_jmax - invfftx_strip_jmin + 1)
        
     CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_master, field_master, MPI_COMM_WORLD, stattus, ierr)
        
     pos2 = 0
     DO j = invfftx_strip_jmin, invfftx_strip_jmax
        pos1 = pos2 + 1
        pos2 = pos2 + global_maximal_i + 1
        xband(0:global_maximal_i, j) = rbufer(pos1:pos2)
     END DO

  END IF

!if (T_cntr.eq.10) then
!myfilename = 'C_NNNN.dat'
!myfilename(3:6) = convert_int_to_txt_string(Rank_of_process,4)
!open (51, file = myfilename)
!do j = invfftx_strip_jmin, invfftx_strip_jmax
!   do i = 0, global_maximal_i
!      write (51, '(2x,i4,2x,i4,2x,e14.7)') i, j, xband(i,j)
!   end do
!   write (51, '(" ")')
!end do
!close (51, status = 'keep')
!print '("process ",i4," :: file ",A10," is ready")', Rank_of_process, myfilename
!end if

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," passed L")', Rank_of_process

! perform inverse Fourier transformation

  invfactor = 2.0_8 / DBLE(global_maximal_i-1)
  DO j = invfftx_strip_jmin, invfftx_strip_jmax

! move the real part of the N/2 FFT component into the first pair
     xband(1, j) = xband(global_maximal_i-1,j)

     CALL REAL8FT(xband(0:global_maximal_i-2,j), global_maximal_i-1, -1)

     DO i = 0, global_maximal_i-2
        xband(i,j) = xband(i,j) * invfactor
     END DO

     xband(global_maximal_i-1,j) = xband(0,j)
     xband(global_maximal_i,  j) = xband(1,j)

  END DO

! assemble potential in master processes

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!print '("Process ",i4," passed M")', Rank_of_process

! each master collects the band of potential from the associated field calculators ---------------------------------------------------------

  IF (cluster_rank_key.EQ.0) THEN

     DO k = 2, cluster_N_blocks
        indx_y_begin = field_calculator(k)%invfftx_strip_jmin
        indx_y_end   = field_calculator(k)%invfftx_strip_jmax

        bufsize = (global_maximal_i + 1) * (indx_y_end - indx_y_begin + 1)

        CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)

        pos2=0
        DO j = indx_y_begin, indx_y_end
           pos1 = pos2 + 1
           pos2 = pos2 + global_maximal_i + 1
           xband(0:global_maximal_i, j) = rbufer(pos1:pos2)
        END DO
     END DO

  ELSE

     bufsize = (global_maximal_i + 1) * (invfftx_strip_jmax - invfftx_strip_jmin + 1)

     pos2 = 0
     DO j = invfftx_strip_jmin, invfftx_strip_jmax
        pos1 = pos2 + 1
        pos2 = pos2 + global_maximal_i + 1
        rbufer(pos1:pos2) = xband(0:global_maximal_i, j)
     END DO
        
     CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_master, Rank_of_process, MPI_COMM_WORLD, ierr) 

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! transform bands of potential into domains of master processes ------------------------------------------------------------------------------

!########## for each master, make sure you include contribution of the master itself :) ############ ???????????

  IF (cluster_rank_key.EQ.0) THEN

! include the own contribution of the master process
     DO j = invfftx_band_jmin, invfftx_band_jmax
        c_phi(c_indx_x_min:c_indx_x_max,j) = xband(c_indx_x_min:c_indx_x_max,j)
     END DO

     DO k = 1, N_of_comm_steps_X

        IF (comm_FFTX_step(k)%rank.LT.0) CYCLE
        
        IF (Rank_of_process.LT.comm_FFTX_step(k)%rank) THEN

! send
           indx_x_begin = comm_FFTX_step(k)%c_indx_x_min
           indx_x_end   = comm_FFTX_step(k)%c_indx_x_max
           xsize = indx_x_end - indx_x_begin + 1
           bufsize = (invfftx_band_jmax - invfftx_band_jmin + 1) * xsize

           pos2 = 0
           DO j = invfftx_band_jmin, invfftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              rbufer(pos1:pos2) = xband(indx_x_begin:indx_x_end, j)
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

! receive
           indx_x_begin = c_indx_x_min
           indx_x_end   = c_indx_x_max
           xsize = indx_x_end - indx_x_begin + 1
           bufsize = (comm_FFTX_step(k)%invfftx_band_jmax - comm_FFTX_step(k)%invfftx_band_jmin + 1) * xsize
 
           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, comm_FFTX_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)
          
           pos2 = 0
           DO j = comm_FFTX_step(k)%invfftx_band_jmin, comm_FFTX_step(k)%invfftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              c_phi(indx_x_begin:indx_x_end, j) = rbufer(pos1:pos2)
           END DO


        ELSE
! receive
           indx_x_begin = c_indx_x_min
           indx_x_end   = c_indx_x_max
           xsize = indx_x_end - indx_x_begin + 1
           bufsize = (comm_FFTX_step(k)%invfftx_band_jmax - comm_FFTX_step(k)%invfftx_band_jmin + 1) * xsize

           CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, comm_FFTX_step(k)%rank, MPI_COMM_WORLD, stattus, ierr)
          
           pos2 = 0
           DO j = comm_FFTX_step(k)%invfftx_band_jmin, comm_FFTX_step(k)%invfftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              c_phi(indx_x_begin:indx_x_end, j) = rbufer(pos1:pos2)
           END DO

! send
           indx_x_begin = comm_FFTX_step(k)%c_indx_x_min
           indx_x_end   = comm_FFTX_step(k)%c_indx_x_max
           xsize = indx_x_end - indx_x_begin + 1
           bufsize = (invfftx_band_jmax - invfftx_band_jmin + 1) * xsize

           pos2 = 0
           DO j = invfftx_band_jmin, invfftx_band_jmax
              pos1 = pos2 + 1
              pos2 = pos2 + xsize
              rbufer(pos1:pos2) = xband(indx_x_begin:indx_x_end, j)
           END DO

           CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, comm_FFTX_step(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 

        END IF

     END DO

! if requested, the average potential in the injection plane may be set to a given value
! this requires re-calculating the ext_phi array 
! note that here cluster_rank_key=0, that is only master processes execute the code below
! also note that the desired average value of the potential is set equal to the potential associated with the right boundary
! and once the corrected external field is applied, the potential of the right boundary is not necessary equal to 
! the initially associated value

     IF (ht_soft_grid_requested) THEN
! define average potential in the injection plane
! this option should be used only when the upper and lower boundaries are solid metal walls (as in LANDMARK tests)

        ALLOCATE(corr_phi(c_indx_y_min:c_indx_y_max), STAT = ALLOC_ERR)

        myavg_phi = 0.0_8
        IF ((grid_j.GT.c_indx_y_min).AND.(grid_j.LT.c_indx_y_max)) THEN
           DO i = c_indx_x_min+1, c_indx_x_max-1
              myavg_phi = myavg_phi + c_phi(i,grid_j)
           END DO
           IF (c_indx_x_min.EQ.0) myavg_phi = myavg_phi + c_phi(0,grid_j)
        END IF

        CALL MPI_REDUCE(myavg_phi, avg_phi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_HORIZONTAL, ierr)

        IF (Rank_horizontal.EQ.0) avg_phi = avg_phi / DBLE(global_maximal_i)

        CALL MPI_BCAST(avg_phi, 1, MPI_DOUBLE_PRECISION, 0, COMM_HORIZONTAL, ierr)

! recalculate the external potential array
        DO j = c_indx_y_min, c_indx_y_max
!           ext_phi(j) = whole_object(4)%phi + (whole_object(2)%phi - avg_phi - whole_object(4)%phi) * DBLE(j) / DBLE(grid_j)
           corr_phi(j) = (whole_object(2)%phi - avg_phi) * DBLE(j) / DBLE(grid_j)
        END DO

! add the external voltage
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              c_phi(i,j) = c_phi(i,j) + corr_phi(j)
           END DO
        END DO

        DEALLOCATE(corr_phi, STAT = ALLOC_ERR)

     END IF

! ################ diagnostics, electrostatic potential #################
!
     IF (T_cntr.EQ.Save_probes_e_data_T_cntr) THEN
        DO npc = 1, N_of_probes_cluster
           npa = List_of_probes_cluster(npc)
           i = Probe_position(1,npa)
           j = Probe_position(2,npa)
           probe_F_cluster(npc) = c_phi(i,j)
        END DO
     END IF
     
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)                       !------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>

  DEALLOCATE(xband, STAT = ALLOC_ERR)
  DEALLOCATE(rbufer, STAT=ALLOC_ERR)

END SUBROUTINE SOLVE_POISSON_FFTX_LINSYSY

!-------------------------------------------------------
!
SUBROUTINE CALCULATE_ELECTRIC_FIELD_FFTX_LINSYSY

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries
  USE ClusterAndItsBoundaries

  use mpi

  IMPLICIT NONE

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  REAL(8) factor2_E_from_F
  REAL(8) factor1_E_from_F

  INTEGER i, j, k, pos, n1, n3, pos1, pos2
  INTEGER shift, bufsize 

  REAL(8), ALLOCATABLE :: rbufer(:)
  REAL(8), ALLOCATABLE :: loc_EX(:,:)
  REAL(8), ALLOCATABLE :: loc_EY(:,:)
  INTEGER ALLOC_ERR

  factor2_E_from_F = F_scale_V / (E_scale_Vm * 2.0_8 * delta_x_m)
  factor1_E_from_F = F_scale_V / (E_scale_Vm * delta_x_m)

  IF (cluster_rank_key.EQ.0) THEN

! volume
     DO j = c_indx_y_min+1, c_indx_y_max-1
        DO i = c_indx_x_min+1, c_indx_x_max-1
           EX(i,j) = factor2_E_from_F * (c_phi(i-1,j)-c_phi(i+1,j))
           EY(i,j) = factor2_E_from_F * (c_phi(i,j-1)-c_phi(i,j+1))
        END DO
     END DO

! edge above if there is no neighbour
     IF (Rank_of_master_above.LT.0) THEN
        j = c_indx_y_max
        DO i = c_indx_x_min+1, c_indx_x_max-1
           EX(i,j) = factor2_E_from_F * (c_phi(i-1,j)-c_phi(i+1,j))
           EY(i,j) = factor1_E_from_F * (c_phi(i,j-1)-c_phi(i,j))           
        END DO
     END IF
     
! edge below if there is no neighbour
     IF (Rank_of_master_below.LT.0) THEN
        j = c_indx_y_min
        DO i = c_indx_x_min+1, c_indx_x_max-1
           EX(i,j) = factor2_E_from_F * (c_phi(i-1,j)-c_phi(i+1,j))
           EY(i,j) = factor1_E_from_F * (c_phi(i,j)-c_phi(i,j+1))           
        END DO
     END IF

! exchange boundary field values with neighbor masters

     n1 = c_indx_y_max - c_indx_y_min + 1
     n3 = c_indx_x_max - c_indx_x_min + 1

     IF (WHITE_CLUSTER) THEN  
! "white processes"

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:(n1+n1)), STAT=ALLOC_ERR)

        IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right electric fields in column left of the right edge
           rbufer(1:n1)           = EX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer((n1+1):(n1+n1)) = EY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left electric field in column right of the left edge
           rbufer(1:n1)           = EX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer((n1+1):(n1+n1)) = EY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left electric field along the left edge
           CALL MPI_RECV(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(1:n1)   
           EY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer((n1+1):(n1+n1))
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right electric field along the right edge
           CALL MPI_RECV(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(1:n1)
           EY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer((n1+1):(n1+n1))
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:(n3+n3)), STAT=ALLOC_ERR)

        IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up electric fields in the row below the top edge
           rbufer(1:n3)           = EX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer((n3+1):(n3+n3)) = EY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           CALL MPI_SEND(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down electric fields in the row above the bottom edge
           rbufer(1:n3)           = EX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer((n3+1):(n3+n3)) = EY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           CALL MPI_SEND(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below electric fields in the bottom edge
           CALL MPI_RECV(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(1:n3)
           EY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer((n3+1):(n3+n3))
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above electric fields in the top edge
           CALL MPI_RECV(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(1:n3)
           EY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer((n3+1):(n3+n3))
        END IF
     
     ELSE
! "black" processes

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:(n1+n1)), STAT=ALLOC_ERR)

        IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left electric field along the left edge
           CALL MPI_RECV(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(1:n1)   
           EY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer((n1+1):(n1+n1))
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right electric field along the right edge
           CALL MPI_RECV(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(1:n1)
           EY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer((n1+1):(n1+n1))
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right electric fields in column left of the right edge
           rbufer(1:n1)           = EX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer((n1+1):(n1+n1)) = EY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left electric field in column right of the left edge
           rbufer(1:n1)           = EX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer((n1+1):(n1+n1)) = EY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, n1+n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF
     
        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:(n3+n3)), STAT=ALLOC_ERR)

        IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below electric fields in the bottom edge
           CALL MPI_RECV(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(1:n3)
           EY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer((n3+1):(n3+n3))
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above electric fields in the top edge
           CALL MPI_RECV(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           EX(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(1:n3)
           EY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer((n3+1):(n3+n3))
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up electric fields in the row below the top edge
           rbufer(1:n3)           = EX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer((n3+1):(n3+n3)) = EY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           CALL MPI_SEND(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down electric fields in the row above the bottom edge
           rbufer(1:n3)           = EX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer((n3+1):(n3+n3)) = EY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           CALL MPI_SEND(rbufer, n3+n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)

! send complete field array to all members of the cluster
     bufsize = (c_indx_x_max - c_indx_x_min + 1) * (c_indx_y_max - c_indx_y_min + 1)
     CALL MPI_BCAST(EX, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)
     CALL MPI_BCAST(EY, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)

  ELSE

! receive complete field array from the master of the cluster
     bufsize = (c_indx_x_max - c_indx_x_min + 1) * (c_indx_y_max - c_indx_y_min + 1)
     CALL MPI_BCAST(EX, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)
     CALL MPI_BCAST(EY, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)

  END IF

! each member (master and non-master) in each cluster accumulates fields for ions in the whoole cluster domain 
! (in this case there is no need for communications at all)
  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min, c_indx_x_max
        acc_EX(i,j) = acc_EX(i,j) + EX(i,j)
     END DO
  END DO
  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min, c_indx_x_max
        acc_EY(i,j) = acc_EY(i,j) + EY(i,j)
     END DO
  END DO
  
END SUBROUTINE CALCULATE_ELECTRIC_FIELD_FFTX_LINSYSY

!-------------------------------------------------------------
!
! BORROWED FROM NUMERICAL RECIPIES
!
! Calculates the Fourier transform of a set of n real-valued data points.
!
! Replaces this data (which is stored in array data(1:n)) by the positive (!!!) frequency half
! of its complex Fourier transform.
!
! The real-valued first and last components of the complex transform are returned as elements
! data(1) and data(2), respectively.
!
! n must be a power of 2.
!
! The routine also calculates the inverse transform of a complex data array if it is the transform of real data.
! Result in this case must be multiplied by 2/n.
!
! The routine uses FOUR1 
! 
SUBROUTINE REAL8FT(data, n, isign)

  IMPLICIT NONE

  INTEGER isign, n
  REAL(8) data(n)

  INTEGER i, i1, i2, i3, i4, n2p3
  REAL(8) c1, c2, h1i, h1r, h2i, h2r, wis, wrs
  REAL(8) theta, wi, wpi, wpr, wr, wtemp  ! double precision for the trigonometric recurrences

  theta = 3.141592653589793_8 / DBLE(n/2) ! initialize the recurrence

  c1 = 0.5_8
  IF (isign.EQ.1) THEN
     c2 = -0.5_8
     CALL FOUR1_8(data, n/2, +1)   ! the forward transform is here
  ELSE
     c2 = 0.5_8                    ! set up for inverse transform
     theta = -theta
  END IF

  wpr = -2.0_8 * SIN(0.5_8 * theta)**2
  wpi = SIN(theta)
  wr = 1.0_8 + wpr
  wi = wpi
  n2p3 = n+3

  DO i = 2, n/4                        ! case i=1 done separately below

     i1 = 2 * i - 1
     i2 = i1 + 1
     i3 = n2p3 - i2
     i4 = i3 + 1
!     wrs = REAL(wr)
!     wis = REAL(wi)
     wrs = wr
     wis = wi
! the two separate transforms are separated out of data
     h1r = c1 * (data(i1) + data(i3))   
     h1i = c1 * (data(i2) - data(i4)) 
     h2r =-c2 * (data(i2) + data(i4))
     h2i = c2 * (data(i1) - data(i3))
! here they are recombined to form the true transform of the original real data
     data(i1) = h1r + wrs * h2r - wis * h2i   
     data(i2) = h1i + wrs * h2i + wis * h2r
     data(i3) = h1r - wrs * h2r + wis * h2i
     data(i4) =-h1i + wrs * h2i + wis * h2r
! the recurrence
     wtemp = wr
     wr = wr * wpr -    wi * wpi + wr
     wi = wi * wpr + wtemp * wpi + wi

  END DO

  IF (isign.EQ.1) THEN
     h1r = data(1)
     data(1) = h1r + data(2)  ! squeeze the first and last data together to get
     data(2) = h1r - data(2)  ! them all within the original array
  ELSE
     h1r = data(1)
     data(1) = c1 * (h1r + data(2))
     data(2) = c1 * (h1r - data(2))
     CALL FOUR1_8(data, n/2, -1)       ! this is the inverse transform for the case isign=-1
  END IF

  RETURN

END SUBROUTINE REAL8FT

!--------------------------------------
!
! Replaces data (1:2*nn) by its discrete Fourier transform, if isign is input as 1.
!
! Replaces data (1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as -1.
!
! data is a complex array of length nn or, equivalently, a real array of length 2*nn.
!
! nn must be an integer power of 2.
!
SUBROUTINE FOUR1_8(data, nn, isign)

  INTEGER isign, nn
  REAL(8) data(2*nn)

  INTEGER i, istep, j, m, mmax, n
  REAL(8) tempi, tempr
  REAL(8) theta, wi, wpi, wpr, wr, wtemp  ! double precision for the trigonometric recurrences

  n = 2*nn
  j = 1

! the bit-reversal section of the routine
  DO i = 1, n, 2
     IF (j.GT.i) THEN         ! exchange the two complex numbers 
        tempr = data(j)
        tempi = data(j+1)
        data(j)   = data(i)
        data(j+1) = data(i+1)
        data(i)   = tempr
        data(i+1) = tempi
     END IF
     m = nn
     DO WHILE ((m.GE.2).AND.(j.GT.m))
        j = j - m
        m = m / 2
     END DO
     j = j + m
  END DO
  mmax = 2

! Danielson-Lanczos section of the routine
  DO WHILE (n.GT.mmax)
     istep = 2 * mmax
     theta = 6.28318530717959_8 / (isign * mmax)   ! initializa the trigonometric recurrence
     wpr = -2.0_8 * SIN(0.5_8 * theta)**2
     wpi = SIN(theta)
     wr = 1.0_8
     wi = 0.0_8
     DO m = 1, mmax, 2
        DO i = m, n, istep
! Danielson - Lanczos formula:
           j = i + mmax
!           tempr = REAL(wr) * data(j)   - REAL(wi) * data(j+1)
!           tempi = REAL(wr) * data(j+1) + REAL(wi) * data(j)
           tempr = wr * data(j)   - wi * data(j+1)
           tempi = wr * data(j+1) + wi * data(j)
           data(j)   = data(i)   - tempr
           data(j+1) = data(i+1) - tempi
           data(i)   = data(i)   + tempr
           data(i+1) = data(i+1) + tempi
        END DO
! trigonometric recurrence
        wtemp = wr
        wr = wr * wpr -    wi * wpi + wr
        wi = wi * wpr + wtemp * wpi + wi
     END DO
     mmax = istep
  END DO
  
  RETURN

END SUBROUTINE FOUR1_8

!-------------------------------------------------------------
!

!-------------------------------------------------------------
!
! BORROWED FROM NUMERICAL RECIPIES
!
! Calculates the Fourier transform of a set of n real-valued data points.
!
! Replaces this data (which is stored in array data(1:n)) by the positive (!!!) frequency half
! of its complex Fourier transform.
!
! The real-valued first and last components of the complex transform are returned as elements
! data(1) and data(2), respectively.
!
! n must be a power of 2.
!
! The routine also calculates the inverse transform of a complex data array if it is the transform of real data.
! Result in this case must be multiplied by 2/n.
!
! The routine uses FOUR1 
! 
SUBROUTINE REALFT(data, n, isign)

  IMPLICIT NONE

  INTEGER isign, n
  REAL data(n)

  INTEGER i, i1, i2, i3, i4, n2p3
  REAL c1, c2, h1i, h1r, h2i, h2r, wis, wrs
  REAL(8) theta, wi, wpi, wpr, wr, wtemp  ! double precision for the trigonometric recurrences

  theta = 3.141592653589793_8 / DBLE(n/2) ! initialize the recurrence

  c1 = 0.5
  IF (isign.EQ.1) THEN
     c2 = -0.5
     CALL FOUR1(data, n/2, +1)   ! the forward transform is here
  ELSE
     c2 = 0.5                    ! set up for inverse transform
     theta = -theta
  END IF

  wpr = -2.0_8 * SIN(0.5_8 * theta)**2
  wpi = SIN(theta)
  wr = 1.0_8 + wpr
  wi = wpi
  n2p3 = n+3

  DO i = 2, n/4                        ! case i=1 done separately below

     i1 = 2 * i - 1
     i2 = i1 + 1
     i3 = n2p3 - i2
     i4 = i3 + 1
     wrs = REAL(wr)
     wis = REAL(wi)
! the two separate transforms are separated out of data
     h1r = c1 * (data(i1) + data(i3))   
     h1i = c1 * (data(i2) - data(i4)) 
     h2r =-c2 * (data(i2) + data(i4))
     h2i = c2 * (data(i1) - data(i3))
! here they are recombined to form the true transform of the original real data
     data(i1) = h1r + wrs * h2r - wis * h2i   
     data(i2) = h1i + wrs * h2i + wis * h2r
     data(i3) = h1r - wrs * h2r + wis * h2i
     data(i4) =-h1i + wrs * h2i + wis * h2r
! the recurrence
     wtemp = wr
     wr = wr * wpr -    wi * wpi + wr
     wi = wi * wpr + wtemp * wpi + wi

  END DO

  IF (isign.EQ.1) THEN
     h1r = data(1)
     data(1) = h1r + data(2)  ! squeeze the first and last data together to get
     data(2) = h1r - data(2)  ! them all within the original array
  ELSE
     h1r = data(1)
     data(1) = c1 * (h1r + data(2))
     data(2) = c1 * (h1r - data(2))
     CALL FOUR1(data, n/2, -1)       ! this is the inverse transform for the case isign=-1
  END IF

  RETURN

END SUBROUTINE REALFT

!--------------------------------------
!
! Replaces data (1:2*nn) by its discrete Fourier transform, if isign is input as 1.
!
! Replaces data (1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as -1.
!
! data is a complex array of length nn or, equivalently, a real array of length 2*nn.
!
! nn must be an integer power of 2.
!
SUBROUTINE FOUR1(data, nn, isign)

  INTEGER isign, nn
  REAL data(2*nn)

  INTEGER i, istep, j, m, mmax, n
  REAL tempi, tempr
  REAL(8) theta, wi, wpi, wpr, wr, wtemp  ! double precision for the trigonometric recurrences

  n = 2*nn
  j = 1

! the bit-reversal section of the routine
  DO i = 1, n, 2
     IF (j.GT.i) THEN         ! exchange the two complex numbers 
        tempr = data(j)
        tempi = data(j+1)
        data(j)   = data(i)
        data(j+1) = data(i+1)
        data(i)   = tempr
        data(i+1) = tempi
     END IF
     m = nn
     DO WHILE ((m.GE.2).AND.(j.GT.m))
        j = j - m
        m = m / 2
     END DO
     j = j + m
  END DO
  mmax = 2

! Danielson-Lanczos section of the routine
  DO WHILE (n.GT.mmax)
     istep = 2 * mmax
     theta = 6.28318530717959_8 / (isign * mmax)   ! initializa the trigonometric recurrence
     wpr = -2.0_8 * SIN(0.5_8 * theta)**2
     wpi = SIN(theta)
     wr = 1.0_8
     wi = 0.0_8
     DO m = 1, mmax, 2
        DO i = m, n, istep
! Danielson - Lanczos formula:
           j = i + mmax
           tempr = REAL(wr) * data(j)   - REAL(wi) * data(j+1)
           tempi = REAL(wr) * data(j+1) + REAL(wi) * data(j)
           data(j)   = data(i)   - tempr
           data(j+1) = data(i+1) - tempi
           data(i)   = data(i)   + tempr
           data(i+1) = data(i+1) + tempi
        END DO
! trigonometric recurrence
        wtemp = wr
        wr = wr * wpr -    wi * wpi + wr
        wi = wi * wpr + wtemp * wpi + wi
     END DO
     mmax = istep
  END DO
  
  RETURN

END SUBROUTINE FOUR1



!-----------------------------------------
!
real(8) function rho_vs_ij(i,j)

  use CurrentProblemValues, ONLY : global_maximal_i, global_maximal_j

  implicit none

  real(8), parameter :: pi = 3.141592653589793_8

  integer N_grid_x, N_grid_y

  integer i, j

  N_grid_x = global_maximal_i-1
  N_grid_y = global_maximal_j

  rho_vs_ij = pi * pi &
            & * &
            & (( 4.0_8 / (dble(N_grid_x))**2 + 1.0_8 / (dble(N_grid_y))**2) * sin(2.0_8 * pi * dble(i) / dble(N_grid_x)) + &
            &  (16.0_8 / (dble(N_grid_x))**2 + 1.0_8 / (dble(N_grid_y))**2) * cos(4.0_8 * pi * dble(i) / dble(N_grid_x))) &
            & * &
            & sin(pi * dble(j) / dble(N_grid_y))

end function rho_vs_ij

!-----------------------------------------
!
real(8) function phi_vs_ij(i,j)

  use CurrentProblemValues, ONLY : global_maximal_i, global_maximal_j

  implicit none

  real(8), parameter :: pi = 3.141592653589793_8

  integer N_grid_x, N_grid_y

  integer i, j

  N_grid_x = global_maximal_i-1
  N_grid_y = global_maximal_j

  phi_vs_ij =  (sin(2.0_8 * pi * dble(i) / dble(N_grid_x)) + cos(4.0_8 * pi * dble(i) / dble(N_grid_x))) &
            & * &
            & sin(pi * dble(j) / dble(N_grid_y))

end function phi_vs_ij
