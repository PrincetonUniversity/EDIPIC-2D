!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER jbelow, jabove
  REAL(8) dqbelow, dqabove

  particle_not_processed = .TRUE.

  y = MIN(MAX(y,c_Y_area_min),c_Y_area_max)

  DO n = 1, c_N_of_local_object_parts_left

     m = c_index_of_local_object_part_left(n)

     IF ( (y.GE.c_local_object_part(m)%jstart).AND. &
        & (y.LE.c_local_object_part(m)%jend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = y - jbelow
              dqbelow = 1.0_8 - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) - dqbelow
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) - dqabove
        END SELECT

!###        IF (whole_object(nwo)%SEE_enabled) THEN
!###           CALL PROCESS_SEE_EVENT(x, y, vx, vy, vz, 1, m)   ! "1" is for a left wall 
!###        END IF

!>>>> processing begins >>>>
! the simplest case - specular reflection everywhere
!###        CALL ADD_ELECTRON_TO_ADD_LIST(c_X_area_min, y, -vx, vy, vz, tag)
! here the long min/max operator is necessary to process a particle hitting a corner
!<<<< processing ends <<<<
        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     STOP
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT

!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER jbelow, jabove
  REAL(8) dqbelow, dqabove

  particle_not_processed = .TRUE.

  y = MIN(MAX(y,c_Y_area_min),c_Y_area_max)

  DO n = 1, c_N_of_local_object_parts_right

     m = c_index_of_local_object_part_right(n)

     IF ( (y.GE.c_local_object_part(m)%jstart).AND. &
        & (y.LE.c_local_object_part(m)%jend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = y - jbelow
              dqbelow = 1.0_8 - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) - dqbelow
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) - dqabove
        END SELECT

!>>>> processing begins >>>>
! the simplest case - specular reflection everywhere
!###         CALL ADD_ELECTRON_TO_ADD_LIST(c_X_area_max, y, -vx, vy, vz, tag)
!<<<< processing ends <<<<
        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     STOP
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT

!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER ileft, iright
  REAL(8) dqleft, dqright

  particle_not_processed = .TRUE.

  x = MIN(MAX(x,c_X_area_min),c_X_area_max)

  DO n = 1, c_N_of_local_object_parts_below

     m = c_index_of_local_object_part_below(n)

     IF ( (x.GE.c_local_object_part(m)%istart).AND. &
        & (x.LE.c_local_object_part(m)%iend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = x - ileft
              dqleft = 1.0_8 - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft) - dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) - dqright
        END SELECT

!>>>> processing begins >>>>
! the simplest case - specular reflection everywhere
!###         CALL ADD_ELECTRON_TO_ADD_LIST(x, c_Y_area_min, vx, -vy, vz, tag)
!<<<< processing ends <<<<
        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     STOP
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW

!-----------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(x, y, vx, vy, vz, tag)

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object, VACUUM_GAP, METAL_WALL, DIELECTRIC

  IMPLICIT NONE

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  INTEGER n, m, nwo  ! nwo stands for number of the whole object
  LOGICAL particle_not_processed

  INTEGER ileft, iright
  REAL(8) dqleft, dqright

  particle_not_processed = .TRUE.

  x = MIN(MAX(x,c_X_area_min),c_X_area_max)

  DO n = 1, c_N_of_local_object_parts_above

     m = c_index_of_local_object_part_above(n)

     IF ( (x.GE.c_local_object_part(m)%istart).AND. &
        & (x.LE.c_local_object_part(m)%iend) ) THEN

        nwo = c_local_object_part(m)%object_number

        whole_object(nwo)%electron_hit_count = whole_object(nwo)%electron_hit_count + 1

        SELECT CASE (whole_object(nwo)%object_type)
           CASE (VACUUM_GAP)
           CASE (METAL_WALL)
           CASE (DIELECTRIC)
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = x - ileft
              dqleft = 1.0_8 - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft) - dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) - dqright
        END SELECT

!>>>> processing begins >>>>
! the simplest case - specular reflection everywhere
!###        CALL ADD_ELECTRON_TO_ADD_LIST(x, c_Y_area_max, vx, -vy, vz, tag)
!<<<< processing ends <<<<
        particle_not_processed = .FALSE.
        EXIT

     END IF
     CYCLE
  END DO

  IF (particle_not_processed) THEN
     PRINT '("Process ",i4,": ERROR in PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE")', Rank_of_process
     PRINT '("particle x= ",e14.7," y= ",e14.7)', x, y
     STOP
  END IF

END SUBROUTINE PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE

!------------------------------------------------------
!
SUBROUTINE COLLECT_ELECTRON_BOUNDARY_HITS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : N_spec

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
!  INTEGER stattus(MPI_STATUS_SIZE)
!  INTEGER request

  INTEGER, ALLOCATABLE :: ibuf_send(:)
  INTEGER, ALLOCATABLE :: ibuf_receive(:)
  INTEGER ALLOC_ERR

  INTEGER k

  IF (c_N_of_local_object_parts.GT.0) THEN

     ALLOCATE(ibuf_send(1:N_of_boundary_objects), STAT = ALLOC_ERR)
     ALLOCATE(ibuf_receive(1:N_of_boundary_objects), STAT = ALLOC_ERR)

! each cluster adjacent to a boundary assembles electron-boundary hit counters from all cluster members in the master of the cluster

     ibuf_send(1:N_of_boundary_objects) = whole_object(1:N_of_boundary_objects)%electron_hit_count
     ibuf_receive = 0

     CALL MPI_REDUCE(ibuf_send, ibuf_receive, N_of_boundary_objects, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)

     IF ((cluster_rank_key.EQ.0).AND.(c_N_of_local_object_parts.GT.0)) THEN
! now counters from all boundary cluster masters are assembled in the process with global rank zero

        whole_object(1:N_of_boundary_objects)%electron_hit_count = ibuf_receive(1:N_of_boundary_objects)

        ibuf_send(1:N_of_boundary_objects) = whole_object(1:N_of_boundary_objects)%electron_hit_count
        ibuf_receive = 0
! 
        CALL MPI_REDUCE(ibuf_send, ibuf_receive, N_of_boundary_objects, MPI_INTEGER, MPI_SUM, 0, COMM_BOUNDARY, ierr)

        IF (Rank_of_process.EQ.0) THEN
           whole_object(1:N_of_boundary_objects)%electron_hit_count = ibuf_receive(1:N_of_boundary_objects)
print '("electrons hit boundaries :: ",10(2x,i8))', whole_object(1:N_of_boundary_objects)%electron_hit_count  

! set the ion hit counters here to zero because when this subroutine is called the ions do not move
           DO k = 1, N_of_boundary_objects
              whole_object(k)%ion_hit_count(1:N_spec) = 0
           END DO

        END IF

     END IF

     DEALLOCATE(ibuf_send, STAT = ALLOC_ERR)
     DEALLOCATE(ibuf_receive, STAT = ALLOC_ERR)

  END IF

END SUBROUTINE COLLECT_ELECTRON_BOUNDARY_HITS

!-------------------------------------------------------------------------------------------
!
SUBROUTINE INITIATE_WALL_DIAGNOSTICS

  USE ParallelOperationValues
  USE CurrentProblemValues, ONLY : N_of_boundary_objects, Start_T_cntr
  USE Checkpoints, ONLY : use_checkpoint
!  USE Diagnostics, ONLY : N_of_saved_records
  USE SetupValues, ONLY : ht_use_e_emission_from_cathode, ht_use_e_emission_from_cathode_zerogradf, ht_emission_constant

  IMPLICIT NONE

                                    ! ----x----I----x--
  CHARACTER(17) historybo_filename  ! history_bo_NN.dat

  LOGICAL exists
  INTEGER i, k
  INTEGER i_dummy

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (Rank_of_process.NE.0) RETURN

  IF (ht_use_e_emission_from_cathode.OR.ht_use_e_emission_from_cathode_zerogradf.OR.ht_emission_constant) RETURN

! hardwired for objects #2 (cathode) and #4 (anode)

  IF (use_checkpoint.EQ.1) THEN
! start from checkpoint, must trim the time dependences

     DO k = 1, N_of_boundary_objects

        historybo_filename = 'history_bo_NN.dat'
        historybo_filename(12:13) = convert_int_to_txt_string(k, 2)

        INQUIRE (FILE = historybo_filename, EXIST = exists)
        IF (exists) THEN                                                       
           OPEN (21, FILE = historybo_filename, STATUS = 'OLD')          
           DO i = 1, Start_T_cntr   !N_of_saved_records             ! these files are updated at every electron timestep
              READ (21, '(2x,i8,10(2x,i8))') i_dummy
           END DO
           ENDFILE 21       
           CLOSE (21, STATUS = 'KEEP')        
        END IF

     END DO

  ELSE
! fresh start, empty files, clean up whatever garbage there might be

     DO k = 1, N_of_boundary_objects

        historybo_filename = 'history_bo_NN.dat'
        historybo_filename(12:13) = convert_int_to_txt_string(k, 2)

        OPEN  (21, FILE = historybo_filename, STATUS = 'REPLACE')          
        CLOSE (21, STATUS = 'KEEP')

     END DO

  END IF

END SUBROUTINE INITIATE_WALL_DIAGNOSTICS

!-------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE SetupValues, ONLY : ht_use_e_emission_from_cathode, ht_use_e_emission_from_cathode_zerogradf, ht_emission_constant
  USE IonParticles, ONLY : N_spec

  IMPLICIT NONE

  INTEGER k
                                    ! ----x----I----x--
  CHARACTER(17) historybo_filename  ! history_bo_NN.dat

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  IF (Rank_of_process.NE.0) RETURN

  IF (ht_use_e_emission_from_cathode.OR.ht_use_e_emission_from_cathode_zerogradf.OR.ht_emission_constant) RETURN

  DO k = 1, N_of_boundary_objects

     historybo_filename = 'history_bo_NN.dat'
     historybo_filename(12:13) = convert_int_to_txt_string(k, 2)

     OPEN (21, FILE = historybo_filename, POSITION = 'APPEND')
     WRITE (21, '(2x,i8,10(2x,i8))') &
          & T_cntr, &
          & whole_object(k)%electron_hit_count , &
          & whole_object(k)%ion_hit_count(1:N_spec), &
          & whole_object(k)%electron_emit_count

     CLOSE (21, STATUS = 'KEEP')
  
  END DO

END SUBROUTINE SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS

!-------------------------------------------------------------------------------------------
! Prepares the tabulated values of integral of the maxwell distribution function
!  
SUBROUTINE PrepareMaxwellDistribIntegral

  USE ParallelOperationValues
  USE MaxwellVelocity
!  USE CurrentProblemValues, ONLY : N_box_vel
  IMPLICIT NONE

  INTEGER i
  INTEGER N_pnts
  INTEGER count
  REAL(8) V_min, V_max
  REAL(8) F(0:180003)      ! to be sure that we overcome V_max
  REAL(8) temp
  REAL(8) dV

  LOGICAL check1, check2

  check1 = .FALSE.
  check2 = .FALSE.

! ------- for symmetrical maxwellian
  N_pnts = 180000  !30000
  V_min  = -U_max
  V_max  =  U_max
  dV = (V_max - V_min) / N_pnts

  F = 0.0_8
  DO i = 1, N_pnts + 3
     F(i) = F(i-1) + EXP( - (V_min + (DBLE(i)-0.5_8) * dV)**2 )
  END DO

  temp = F(N_pnts)
  F = F * R_max / temp   ! normalize integral such that F(N_pnts) = R_max

  v(0) = V_min
  count = 0
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        v(count) = V_min + i * dV
        IF (count.EQ.R_max) THEN
           check1 = .TRUE.
           EXIT
        END IF
     END IF
  END DO

!  v = v * N_box_vel

!--------- for asymmetrical maxwellian * v (used for injection, v > 0)

  N_pnts = 90000   !15000
  V_min  = 0.0_8
  V_max  = U_max
  dV = (V_max - V_min) / N_pnts

  F = 0.0_8
  DO i = 1, N_pnts + 3
     temp = V_min + (REAL(i)-0.5_8) * dV
     F(i) = F(i-1) + EXP( - temp**2 ) * temp
  END DO

  temp = F(N_pnts)
  F(1:(N_pnts+3)) = F(1:(N_pnts+3)) * R_max_inj / temp   ! normalize integral such that F(N_pnts) = R_max_inj

  v_inj(0) = V_min
  count = 0
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        v_inj(count) = V_min + i * dV
        IF (count.EQ.R_max_inj) THEN
           check2 = .TRUE.
           EXIT
        END IF
     END IF
  END DO

!  v_inj = v_inj * N_box_vel

  IF (check1.AND.check2) THEN
!     PRINT '(2x,"Process ",i3," : Integrals for producing maxwell distributions are successfully obtained ...")', &
!                                                                                                  & Rank_of_process
  ELSE
     PRINT '(2x,"Process ",i3," : ERROR in PrepareMaxwellDistribIntegral !!!")', Rank_of_process
     PRINT '(2x,"The initialization in PrepareMaxwellDistribIntegral is not performed !!!")'
     PRINT '(2x,"The program will be terminated now :(")'
     STOP
  END IF

END SUBROUTINE PrepareMaxwellDistribIntegral

!-------------------------------------------------------------------------------------------
!  
SUBROUTINE GetInjMaxwellVelocity(U) 

  USE MaxwellVelocity
 
  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) U

  REAL(8) R
  INTEGER indx
  
  R = R_max_inj * well_random_number()

  indx = INT(R)

  IF (indx.LT.R_max_inj) THEN
     U = v_inj(indx) + (R - indx) * (v_inj(indx+1) - v_inj(indx))
  ELSE
     U = v_inj(R_max_inj)
  END IF
  RETURN
  
END SUBROUTINE GetInjMaxwellVelocity

!-------------------------------------------------------------------------------------------
!  
SUBROUTINE GetMaxwellVelocity(U) 

  USE MaxwellVelocity

  USE rng_wrapper

  IMPLICIT NONE

  REAL(8) U

  REAL(8) R
  INTEGER indx
  
  R = R_max * well_random_number()

  indx = INT(R)

  IF (indx.LT.R_max) THEN
     U = v(indx) + (R - indx) * (v(indx+1) - v(indx))
  ELSE
     U = v(R_max)
  END IF
  RETURN
  
END SUBROUTINE GetMaxwellVelocity
