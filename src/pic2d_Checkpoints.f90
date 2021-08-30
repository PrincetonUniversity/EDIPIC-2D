
!-------------------------------------------
! this subroutine is called by every process
!
SUBROUTINE SAVE_CHECKPOINT_MPIIO_2(n_sub)

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE LoadBalancing, ONLY : T_cntr_global_load_balance, T_cntr_cluster_load_balance
  USE Diagnostics, ONLY : Save_probes_data_T_cntr, N_of_saved_records
  USE Snapshots, ONLY : current_snap
  USE ElectronParticles
  USE IonParticles
  USE Checkpoints
  USE ClusterAndItsBoundaries, ONLY : c_local_object_part, c_N_of_local_object_parts
  USE SetupValues, ONLY : total_cathode_N_e_to_inject

  USE rng_wrapper

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER, INTENT(IN) :: n_sub

  INTEGER*4 state_ind, func, state_var(624)

  INTEGER ibufer_length
  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  INTEGER s, pos, pos2

  INTEGER len_surf_charge
  INTEGER n, nwo, jstart, jend, istart, iend

  REAL(8), ALLOCATABLE :: rbufer(:)

  INTEGER file_handle

  CHARACTER(20) filename_check         ! Tcntr_TTTTTTTT.check
                                       ! ----x----I----x----I

  INTEGER, PARAMETER :: maxNsend=100000
  INTEGER Nsend

  CHARACTER(21) filename_report        ! report_TTTTTTTT.write
                                       ! ----x----I----x----I-

  INTEGER myibufer(1)

  INTEGER, ALLOCATABLE :: jbufer(:)
  REAL(8), ALLOCATABLE :: dbufer(:)
  INTEGER k, bufsize

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

! fool proof
  IF (n_sub.NE.0) THEN
     PRINT '("Process :: ",i5," ERROR :: SAVE_CHECKPOINT called at wrong time, n_sub = ",i4)', Rank_of_process, n_sub
     STOP
  END IF

  CALL get_rng_state(state_var, state_ind, func)

  ibufer_length = 636+2*N_spec+1      ! 636=624+1+1+7+2+1

  ALLOCATE(ibufer(1:ibufer_length), STAT = ALLOC_ERR)

  ibufer(1:624) = state_var(1:624)
  ibufer(625) = state_ind
  ibufer(626) = func

  ibufer(627) = T_cntr
  ibufer(628) = T_cntr_global_load_balance
  ibufer(629) = T_cntr_cluster_load_balance
  ibufer(630) = Save_probes_data_T_cntr
  ibufer(631) = N_of_saved_records
  ibufer(632) = current_snap
  ibufer(633) = particle_master

  ibufer(634) = max_N_electrons
  ibufer(635) = N_electrons

  ibufer(636) = total_cathode_N_e_to_inject
  
  pos=637

  DO s = 1, N_spec
     ibufer(pos) = max_N_ions(s)
     ibufer(pos+1) = N_ions(s)
     pos = pos+2
  END DO
  
  len_surf_charge = 0

  IF (cluster_rank_key.EQ.0) THEN
! save surface charge arrays if applicable     
     DO n = 1, c_N_of_local_object_parts
        nwo = c_local_object_part(n)%object_number
        IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
           jstart = c_local_object_part(n)%jstart
           jend = c_local_object_part(n)%jend
           istart = c_local_object_part(n)%istart
           iend = c_local_object_part(n)%iend
           IF (jend.GT.jstart) THEN
! vertical boundary
              len_surf_charge = len_surf_charge + jend - jstart + 1
!              WRITE (1) c_local_object_part(n)%surface_charge(jstart:jend)
           ELSE IF (iend.GT.istart) THEN
! horizontal boundary
              len_surf_charge = len_surf_charge + iend - istart + 1
!              WRITE(1) c_local_object_part(n)%surface_charge(istart:iend)
           END IF
        END IF
     END DO
! include dielectric inner objects, only for left bottom cluster master (all masters have same copy)
     IF (Rank_of_process.EQ.Rank_of_bottom_left_cluster_master) THEN
        DO nwo = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
              len_surf_charge = len_surf_charge + whole_object(nwo)%N_boundary_nodes
           END IF
        END DO
     END IF
  END IF

  ibufer(pos) = len_surf_charge

  ALLOCATE(rbufer(1:MAX(1,len_surf_charge)), STAT = ALLOC_ERR)

  pos2=0

  IF (cluster_rank_key.EQ.0) THEN
! save surface charge arrays if applicable     
     DO n = 1, c_N_of_local_object_parts
        nwo = c_local_object_part(n)%object_number
        IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
           jstart = c_local_object_part(n)%jstart
           jend = c_local_object_part(n)%jend
           istart = c_local_object_part(n)%istart
           iend = c_local_object_part(n)%iend
           IF (jend.GT.jstart) THEN
! vertical boundary
              pos = pos2 + 1
              pos2 = pos2 + jend - jstart + 1
              rbufer(pos:pos2) = c_local_object_part(n)%surface_charge(jstart:jend)
           ELSE IF (iend.GT.istart) THEN
! horizontal boundary
              pos = pos2 + 1
              pos2 = pos2 + iend - istart + 1
              rbufer(pos:pos2) = c_local_object_part(n)%surface_charge(istart:iend)
           END IF
        END IF
     END DO
! include dielectric inner objects, only for left bottom cluster master (all masters have same copy)
     IF (Rank_of_process.EQ.Rank_of_bottom_left_cluster_master) THEN
        DO nwo = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
              pos = pos2 + 1
              pos2 = pos2 + whole_object(nwo)%N_boundary_nodes
              rbufer(pos:pos2) = whole_object(nwo)%surface_charge(1:whole_object(nwo)%N_boundary_nodes)
           END IF
        END DO
     END IF
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! save everything into a single binary file 

! create filename
  filename_check = 'Tcntr_TTTTTTTT.check'
  filename_check(7:14) = convert_int_to_txt_string(T_cntr, 8)

  CALL MPI_FILE_OPEN( MPI_COMM_WORLD, &
                    & filename_check,  &
                    & MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                    & MPI_INFO_NULL, &
                    & file_handle, &
                    & ierr )
  
  CALL MPI_FILE_WRITE_ORDERED( file_handle, ibufer, ibufer_length, MPI_INTEGER, stattus, ierr )
  
  CALL MPI_FILE_WRITE_ORDERED( file_handle, rbufer(1:len_surf_charge), len_surf_charge, MPI_DOUBLE_PRECISION, stattus, ierr )

  bufsize = 5 * N_electrons
  ALLOCATE(dbufer(1:bufsize), STAT = ALLOC_ERR)
  pos = 1
  DO k = 1, N_electrons
     dbufer(pos)   = electron(k)%X
     dbufer(pos+1) = electron(k)%Y
     dbufer(pos+2) = electron(k)%VX
     dbufer(pos+3) = electron(k)%VY
     dbufer(pos+4) = electron(k)%VZ
     pos = pos + 5
  END DO

  CALL MPI_FILE_WRITE_ORDERED( file_handle, dbufer, bufsize, MPI_DOUBLE_PRECISION, stattus, ierr )

  bufsize = N_electrons
  ALLOCATE(jbufer(1:bufsize), STAT = ALLOC_ERR)
  DO k = 1, N_electrons
     jbufer(k) = electron(k)%tag
  END DO

  CALL MPI_FILE_WRITE_ORDERED( file_handle, jbufer, bufsize, MPI_INTEGER, stattus, ierr )

  DO s = 1, N_spec

     bufsize = 5 * N_ions(s)
     IF (ALLOCATED(dbufer)) DEALLOCATE(dbufer, STAT = ALLOC_ERR)
     ALLOCATE(dbufer(1:bufsize), STAT = ALLOC_ERR)

     pos = 1
     DO k = 1, N_ions(s)
        dbufer(pos)   = ion(s)%part(k)%X
        dbufer(pos+1) = ion(s)%part(k)%Y
        dbufer(pos+2) = ion(s)%part(k)%VX
        dbufer(pos+3) = ion(s)%part(k)%VY
        dbufer(pos+4) = ion(s)%part(k)%VZ
        pos = pos + 5
     END DO

     CALL MPI_FILE_WRITE_ORDERED( file_handle, dbufer, bufsize, MPI_DOUBLE_PRECISION, stattus, ierr )

     bufsize = N_ions(s)
     IF (ALLOCATED(jbufer)) DEALLOCATE(jbufer, STAT = ALLOC_ERR)
     ALLOCATE(jbufer(1:bufsize), STAT = ALLOC_ERR)

     DO k = 1, N_ions(s)
        jbufer(k)   = ion(s)%part(k)%tag
     END DO

     CALL MPI_FILE_WRITE_ORDERED( file_handle, jbufer, bufsize, MPI_INTEGER, stattus, ierr )

  END DO

  CALL MPI_FILE_CLOSE(file_handle, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! REPORT==============================================================================================================

  IF (Rank_of_process.LT.(N_of_processes-1)) THEN
     CALL MPI_RECV(myibufer(1), 1, MPI_INTEGER, Rank_of_process+1, Rank_of_process+1, MPI_COMM_WORLD, stattus, ierr)
  END IF

  filename_report = 'report_TTTTTTTT.write'
  filename_report(8:15) = convert_int_to_txt_string(T_cntr, 8)
  IF (Rank_of_process.EQ.(N_of_processes-1)) THEN
     OPEN (2, FILE = filename_report)
  ELSE
     OPEN (2, FILE = filename_report, POSITION = 'APPEND')
  END IF
  WRITE (2, '(2x,i5,"-----------------------------------------------------")') Rank_of_process

  WRITE (2, '(8(2x,i12))') ibufer(1:624) 
  WRITE (2, '(2(2x,i10))') ibufer(625:626)
  WRITE (2, '(7(2x,i10))') ibufer(627:633) 
  WRITE (2, '(2(2x,i10))') ibufer(634:635)
  WRITE (2, '(2x,i10)') ibufer(636)
  WRITE (2, '(10(2x,i10))') ibufer(637:ibufer_length-1)
  WRITE (2, '(2x,i10)') ibufer(ibufer_length)

  IF (len_surf_charge.GT.0) WRITE (2, '("surface charge first/last :: ",2(2x,e21.14))') rbufer(1), rbufer(len_surf_charge)

  Nsend = MIN(maxNsend,N_electrons)
  pos2=0
  DO WHILE (pos2.LT.N_electrons)
     pos  = pos2 + 1
     pos2 = MIN(pos2 + Nsend, N_electrons)

     WRITE (2, '("electrons # ",i8," / ",i8," :: ",5(2x,e21.14),2x,i3," / ",5(2x,e21.14),2x,i3)') pos, pos2, &
       & electron(pos)%X,  electron(pos)%Y,  electron(pos)%VX,  electron(pos)%VY,  electron(pos)%VZ,  electron(pos)%tag, &
       & electron(pos2)%X, electron(pos2)%Y, electron(pos2)%VX, electron(pos2)%VY, electron(pos2)%VZ, electron(pos2)%tag

  END DO

  DO s = 1, N_spec

     Nsend = MIN(maxNsend,N_ions(s))
     pos2=0
     DO WHILE (pos2.LT.N_ions(s))
        pos  = pos2 + 1
        pos2 = MIN(pos2 + Nsend, N_ions(s))

        WRITE (2, '("ions ",i2," # ",i8," / ",i8,"  :: "5(2x,e21.14),2x,i3," / ",5(2x,e21.14),2x,i3)') s, pos, pos2, &
          & ion(s)%part(pos)%X,  ion(s)%part(pos)%Y,  ion(s)%part(pos)%VX,  ion(s)%part(pos)%VY,  ion(s)%part(pos)%VZ,  ion(s)%part(pos)%tag, &
          & ion(s)%part(pos2)%X, ion(s)%part(pos2)%Y, ion(s)%part(pos2)%VX, ion(s)%part(pos2)%VY, ion(s)%part(pos2)%VZ, ion(s)%part(pos2)%tag

     END DO
  END DO

  CLOSE (2, STATUS = 'KEEP')

  IF (Rank_of_process.GT.0) THEN  
     myibufer(1) = Rank_of_process
     CALL MPI_SEND(myibufer(1), 1, MPI_INTEGER, Rank_of_process-1, Rank_of_process, MPI_COMM_WORLD, request, ierr)
  END IF

! END OF REPORT=======================================================================================================

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (ALLOCATED(ibufer)) DEALLOCATE(ibufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(dbufer)) DEALLOCATE(dbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(jbufer)) DEALLOCATE(jbufer, STAT = ALLOC_ERR)

END SUBROUTINE SAVE_CHECKPOINT_MPIIO_2

!-------------------------------------------
! this subroutine is called by every process
!
SUBROUTINE READ_CHECKPOINT_MPIIO_2

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE LoadBalancing, ONLY : T_cntr_global_load_balance, T_cntr_cluster_load_balance
  USE ElectronParticles
  USE IonParticles
  USE Checkpoints
  USE ClusterAndItsBoundaries, ONLY : c_local_object_part, c_N_of_local_object_parts
  USE SetupValues, ONLY : total_cathode_N_e_to_inject

  USE rng_wrapper

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  LOGICAL exists

  INTEGER*4 state_ind, func, state_var(624)

  INTEGER ibufer_length
  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  INTEGER file_handle

  CHARACTER(20) filename_check         ! Tcntr_TTTTTTTT.check
                                       ! ----x----I----x----I

  INTEGER len_surf_charge
  INTEGER n, nwo, jstart, jend, istart, iend

  INTEGER s, pos, pos2

  REAL(8), ALLOCATABLE :: rbufer(:)

  INTEGER, PARAMETER :: maxNsend=100000
  INTEGER Nsend

  CHARACTER(20) filename_report        ! report_TTTTTTTT.read
                                       ! ----x----I----x----I

  INTEGER myibufer(1)

  INTEGER, ALLOCATABLE :: jbufer(:)
  REAL(8), ALLOCATABLE :: dbufer(:)
  INTEGER k, bufsize

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

! foolproof: make sure that the unpacked file has proper name

  filename_check = 'Tcntr_TTTTTTTT.check'
  filename_check(7:14) = convert_int_to_txt_string(T_cntr_to_continue, 8)

  INQUIRE (FILE = filename_check, EXIST = exists)
  IF (.NOT.exists) THEN
     PRINT '("Process :: ",i5,", ERROR in READ_CHECKPOINT_MPIIO :: file ",A20," not found program terminated")', &
          & Rank_of_process, filename_check
     STOP
  END IF

! read the binary file 

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  call MPI_FILE_OPEN( MPI_COMM_WORLD, &
                    & filename_check,  &
                    & MPI_MODE_RDONLY, & 
                    & MPI_INFO_NULL, &
                    & file_handle, &
                    & ierr )

  ibufer_length = 636+2*N_spec+1      ! 636=624+1+1+7+2+1

  ALLOCATE(ibufer(1:ibufer_length), STAT = ALLOC_ERR)

  CALL MPI_FILE_READ_ORDERED( file_handle, ibufer, ibufer_length, MPI_INTEGER, stattus, ierr)

  state_var(1:624) = ibufer(1:624)
  state_ind = ibufer(625)
  func = ibufer(626)

  CALL set_rng_state(state_var, state_ind, func)

  Start_T_cntr                  = ibufer(627) != T_cntr
  T_cntr_global_load_balance    = ibufer(628) != T_cntr_global_load_balance
  T_cntr_cluster_load_balance   = ibufer(629) != T_cntr_cluster_load_balance
  Save_probes_data_T_cntr_check = ibufer(630) != Save_probes_data_T_cntr
  N_of_saved_records_check      = ibufer(631) != N_of_saved_records
  current_snap_check            = ibufer(632) != current_snap
  particle_master               = ibufer(633) != particle_master

  max_N_electrons = ibufer(634)  != max_N_electrons
  N_electrons     = ibufer(635)  != N_electrons

  total_cathode_N_e_to_inject = ibufer(636)

  pos=637

  DO s = 1, N_spec
     max_N_ions(s) = ibufer(pos) != max_N_ions(s)
     N_ions(s)     = ibufer(pos+1) != N_ions(s)
     pos = pos+2
  END DO
  
  len_surf_charge = ibufer(pos) != len_surf_charge

  ALLOCATE(electron(1:max_N_electrons), STAT=ALLOC_ERR)

  ALLOCATE(ion(1:N_spec), STAT = ALLOC_ERR)
  DO s = 1, N_spec
     ALLOCATE(ion(s)%part(1:max_N_ions(s)), STAT = ALLOC_ERR)
  END DO

  ALLOCATE(rbufer(1:MAX(1,len_surf_charge)), STAT = ALLOC_ERR)

  CALL MPI_FILE_READ_ORDERED( file_handle, rbufer(1:len_surf_charge), len_surf_charge, MPI_DOUBLE_PRECISION, stattus, ierr )

  pos2 = 0

  IF (cluster_rank_key.EQ.0) THEN
! read surface charge arrays if applicable     
     DO n = 1, c_N_of_local_object_parts
        nwo = c_local_object_part(n)%object_number
        IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
           jstart = c_local_object_part(n)%jstart
           jend = c_local_object_part(n)%jend
           istart = c_local_object_part(n)%istart
           iend = c_local_object_part(n)%iend
           IF (jend.GT.jstart) THEN
! vertical boundary
              pos = pos2 + 1
              pos2 = pos2 + jend - jstart + 1
              c_local_object_part(n)%surface_charge(jstart:jend) = rbufer(pos:pos2)
           ELSE IF (iend.GT.istart) THEN
! horizontal boundary
              pos = pos2 + 1
              pos2 = pos2 + iend - istart + 1
              c_local_object_part(n)%surface_charge(istart:iend) = rbufer(pos:pos2)
           END IF
        END IF
     END DO
! get surface charge on dielectric inner objects, only for left bottom cluster master (all other masters will get the copy later, after SET_COMMUNICATIONS is called)
     IF (Rank_of_process.EQ.Rank_of_bottom_left_cluster_master) THEN
        DO nwo = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (whole_object(nwo)%object_type.EQ.DIELECTRIC) THEN
              pos = pos2 + 1
              pos2 = pos2 + whole_object(nwo)%N_boundary_nodes
              whole_object(nwo)%surface_charge(1:whole_object(nwo)%N_boundary_nodes) = rbufer(pos:pos2)
           END IF
        END DO
     END IF
 END IF

  bufsize = 5 * N_electrons
  ALLOCATE(dbufer(1:bufsize), STAT = ALLOC_ERR)

  CALL MPI_FILE_READ_ORDERED( file_handle, dbufer, bufsize, MPI_DOUBLE_PRECISION, stattus, ierr )

  pos = 1
  DO k = 1, N_electrons
     electron(k)%X  = dbufer(pos)
     electron(k)%Y  = dbufer(pos+1) 
     electron(k)%VX = dbufer(pos+2) 
     electron(k)%VY = dbufer(pos+3) 
     electron(k)%VZ = dbufer(pos+4) 
     pos = pos + 5
  END DO

  bufsize = N_electrons
  ALLOCATE(jbufer(1:bufsize), STAT = ALLOC_ERR)

  CALL MPI_FILE_READ_ORDERED( file_handle, jbufer, bufsize, MPI_INTEGER, stattus, ierr )

  DO k = 1, N_electrons
     electron(k)%tag = jbufer(k)
  END DO

  DO s = 1, N_spec

     bufsize = 5 * N_ions(s)
     IF (ALLOCATED(dbufer)) DEALLOCATE(dbufer, STAT = ALLOC_ERR)
     ALLOCATE(dbufer(1:bufsize), STAT = ALLOC_ERR)

     CALL MPI_FILE_READ_ORDERED( file_handle, dbufer, bufsize, MPI_DOUBLE_PRECISION, stattus, ierr )

     pos = 1
     DO k = 1, N_ions(s)
        ion(s)%part(k)%X  = dbufer(pos)
        ion(s)%part(k)%Y  = dbufer(pos+1)
        ion(s)%part(k)%VX = dbufer(pos+2)
        ion(s)%part(k)%VY = dbufer(pos+3)
        ion(s)%part(k)%VZ = dbufer(pos+4)
        pos = pos + 5
     END DO

     bufsize = N_ions(s)
     IF (ALLOCATED(jbufer)) DEALLOCATE(jbufer, STAT = ALLOC_ERR)
     ALLOCATE(jbufer(1:bufsize), STAT = ALLOC_ERR)

     CALL MPI_FILE_READ_ORDERED( file_handle, jbufer, bufsize, MPI_INTEGER, stattus, ierr )

     DO k = 1, N_ions(s)
        ion(s)%part(k)%tag = jbufer(k)
     END DO

  END DO

  CALL MPI_FILE_CLOSE(file_handle, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! REPORT==============================================================================================================

  IF (Rank_of_process.LT.(N_of_processes-1)) THEN
     CALL MPI_RECV(myibufer(1), 1, MPI_INTEGER, Rank_of_process+1, Rank_of_process+1, MPI_COMM_WORLD, stattus, ierr)
  END IF

  filename_report = 'report_TTTTTTTT.read'
  filename_report(8:15) = convert_int_to_txt_string(Start_T_cntr, 8)
  IF (Rank_of_process.EQ.(N_of_processes-1)) THEN
     OPEN (2, FILE = filename_report)
  ELSE
     OPEN (2, FILE = filename_report, POSITION = 'APPEND')
  END IF
  WRITE (2, '(2x,i5,"-----------------------------------------------------")') Rank_of_process

  WRITE (2, '(8(2x,i12))') ibufer(1:624) 
  WRITE (2, '(2(2x,i10))') ibufer(625:626)
  WRITE (2, '(7(2x,i10))') ibufer(627:633) 
  WRITE (2, '(2(2x,i10))') ibufer(634:635)
  WRITE (2, '(2x,i10)') ibufer(636)
  WRITE (2, '(10(2x,i10))') ibufer(637:ibufer_length-1)
  WRITE (2, '(2x,i10)') ibufer(ibufer_length)

  IF (len_surf_charge.GT.0) WRITE (2, '("surface charge first/last :: ",2(2x,e21.14))') rbufer(1), rbufer(len_surf_charge)

  Nsend = MIN(maxNsend,N_electrons)
  pos2=0
  DO WHILE (pos2.LT.N_electrons)
     pos  = pos2 + 1
     pos2 = MIN(pos2 + Nsend, N_electrons)

     WRITE (2, '("electrons # ",i8," / ",i8," :: ",5(2x,e21.14),2x,i3," / ",5(2x,e21.14),2x,i3)') pos, pos2, &
       & electron(pos)%X,  electron(pos)%Y,  electron(pos)%VX,  electron(pos)%VY,  electron(pos)%VZ,  electron(pos)%tag, &
       & electron(pos2)%X, electron(pos2)%Y, electron(pos2)%VX, electron(pos2)%VY, electron(pos2)%VZ, electron(pos2)%tag

  END DO

  DO s = 1, N_spec

     Nsend = MIN(maxNsend,N_ions(s))
     pos2=0
     DO WHILE (pos2.LT.N_ions(s))
        pos  = pos2 + 1
        pos2 = MIN(pos2 + Nsend, N_ions(s))

        WRITE (2, '("ions ",i2," # ",i8," / ",i8,"  :: "5(2x,e21.14),2x,i3," / ",5(2x,e21.14),2x,i3)') s, pos, pos2, &
          & ion(s)%part(pos)%X,  ion(s)%part(pos)%Y,  ion(s)%part(pos)%VX,  ion(s)%part(pos)%VY,  ion(s)%part(pos)%VZ,  ion(s)%part(pos)%tag, &
          & ion(s)%part(pos2)%X, ion(s)%part(pos2)%Y, ion(s)%part(pos2)%VX, ion(s)%part(pos2)%VY, ion(s)%part(pos2)%VZ, ion(s)%part(pos2)%tag

     END DO
  END DO

  CLOSE (2, STATUS = 'KEEP')

  IF (Rank_of_process.GT.0) THEN  
     myibufer(1) = Rank_of_process
     CALL MPI_SEND(myibufer(1), 1, MPI_INTEGER, Rank_of_process-1, Rank_of_process, MPI_COMM_WORLD, request, ierr)
  END IF

! END OF REPORT=======================================================================================================

  IF (ALLOCATED(ibufer)) DEALLOCATE(ibufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(dbufer)) DEALLOCATE(dbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(jbufer)) DEALLOCATE(jbufer, STAT = ALLOC_ERR)

END SUBROUTINE READ_CHECKPOINT_MPIIO_2

!--------------------------------------------------
!
SUBROUTINE ADJUST_T_CNTR_SAVE_CHECKPOINT

  USE CurrentProblemValues, ONLY : N_subcycles
  USE ParallelOperationValues, ONLY : Rank_of_process
  USE Snapshots, ONLY : current_snap, N_of_all_snaps, Tcntr_snapshot, save_ionization_rates_2d
  USE Checkpoints, ONLY : T_cntr_save_checkpoint
  USE MCCollisions, ONLY : en_collisions_turned_off, no_ionization_collisions

  IMPLICIT NONE

  INTEGER T_cntr_save_checkpoint_old
  INTEGER n

  IF (en_collisions_turned_off) RETURN
  IF (no_ionization_collisions) RETURN

  T_cntr_save_checkpoint_old = T_cntr_save_checkpoint

! find first snapshot after the desired checkpoint time

  DO n = current_snap, N_of_all_snaps
     IF (Tcntr_snapshot(n).GE.T_cntr_save_checkpoint) THEN
        IF (.NOT.save_ionization_rates_2d(n)) RETURN
! if the snapshot saves ionization rates, move checkpoint save to the timestep following the last ion advance before the snapshot save
        T_cntr_save_checkpoint = N_subcycles * (Tcntr_snapshot(n) / N_subcycles)
        IF (T_cntr_save_checkpoint.NE.T_cntr_save_checkpoint_old) THEN
           IF (Rank_of_process.EQ.0) PRINT '("### ADJUST_T_CNTR_SAVE_CHECKPOINT :: changed T_cntr_save_checkpoint from ",i9," : to ",i9," ###")', &
                & T_cntr_save_checkpoint_old, T_cntr_save_checkpoint
        END IF
        RETURN
     END IF
  END DO

END SUBROUTINE ADJUST_T_CNTR_SAVE_CHECKPOINT
