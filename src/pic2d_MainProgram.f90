!===========================================
PROGRAM MainProg

  USE CurrentProblemValues
  USE ParallelOperationValues
  USE LoadBalancing
  USE ClusterAndItsBoundaries
  USE Checkpoints

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  REAL(8) start, finish
  INTEGER n_sub

  REAL(8) t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, Rank_of_process, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_of_processes, ierr)

  CALL PrepareMaxwellDistribIntegral

  Start_T_cntr = 0
  T_cntr_global_load_balance = Start_T_cntr
  T_cntr_cluster_load_balance = Start_T_cntr

  CALL INITIATE_PARAMETERS
print *, "did INITIATE_PARAMETERS"

  CALL INITIATE_ELECTRON_NEUTRAL_COLLISIONS

print *, "did INITIATE_ELECTRON_NEUTRAL_COLLISIONS"
CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
!stop

  CALL INITIATE_PROBE_DIAGNOSTICS

  CALL INITIATE_WALL_DIAGNOSTICS_HT_SETUP   ! only one of the two actually works
  CALL INITIATE_WALL_DIAGNOSTICS            !

  CALL INITIATE_SNAPSHOTS

  start = MPI_WTIME()

  n_sub = 0

  DO T_cntr = Start_T_cntr, Max_T_cntr

if (Rank_of_process.eq.0) print *, T_cntr

     t0 = MPI_WTIME()

     IF (T_cntr.EQ.T_cntr_save_checkpoint) THEN
        CALL SAVE_CHECKPOINT_MPIIO_2(n_sub)
        T_cntr_save_checkpoint = T_cntr_save_checkpoint + dT_save_checkpoint
     END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     t1 = MPI_WTIME()

call report_total_number_of_particles

     IF (T_cntr.EQ.T_cntr_global_load_balance) THEN
        IF (n_sub.NE.0) THEN
           PRINT '("Process ",i5," :: ERROR-1 in MainProg :: GLOBAL_LOAD_BALANCE is about to be called at wrong time :: T_cntr = ",i8," n_sub = ",i8)', Rank_of_process, T_cntr, n_sub
           STOP
        END IF
        CALL GLOBAL_LOAD_BALANCE  ! includes calls to SET_COMMUNICATIONS 
                                  !                   DISTRIBUTE_CLUSTER_PARAMETERS
        T_cntr_global_load_balance = T_cntr_global_load_balance + dT_global_load_balance
     END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     t2 = MPI_WTIME()

     IF (T_cntr.EQ.T_cntr_cluster_load_balance) THEN
        CALL BALANCE_LOAD_WITHIN_CLUSTER
        T_cntr_cluster_load_balance = T_cntr_cluster_load_balance + dT_cluster_load_balance
     END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     t3 = MPI_WTIME()

     IF (n_sub.EQ.0) CALL GATHER_ION_CHARGE_DENSITY 

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     t4 = MPI_WTIME()

     CALL GATHER_ELECTRON_CHARGE_DENSITY      ! 

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL PERFORM_ELECTRON_EMISSION_HT_SETUP_ZERO_GRAD_F   ! this procedure is used only when axial-azimuthal periodic model of a Hall thruster is simulated

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL UPDATE_WALL_POTENTIALS(T_cntr)

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     t5 = MPI_WTIME()

     IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

!        CALL SOR_GIVEN_POTENTIAL_4_WALLS                ! n
        CALL SOLVE_POTENTIAL_WITH_PETSC
        
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        t6 = MPI_WTIME()

        CALL CALCULATE_ELECTRIC_FIELD               ! n

     ELSE IF (periodicity_flag.EQ.PERIODICITY_X) THEN

        CALL SOLVE_POISSON_FFTX_LINSYSY

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        t6 = MPI_WTIME()

        CALL CALCULATE_ELECTRIC_FIELD_FFTX_LINSYSY
     
     END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     t7 = MPI_WTIME()

     CALL DO_PROBE_DIAGNOSTICS(n_sub)    ! n_sub notifies the sub that the ion densities were refreshed, this avoids accumulation

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     t8 = MPI_WTIME()

     CALL CREATE_SNAPSHOT

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     t9 = MPI_WTIME()

     CALL ADVANCE_ELECTRONS                      !   velocity: n-1/2 ---> n+1/2
                                                 ! coordinate: n     ---> n+1

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     t10 = MPI_WTIME()

     n_sub = n_sub + 1
     IF (n_sub.EQ.N_subcycles) THEN            ! N_subcycles is odd

        if (Rank_of_process.eq.0) print '("----- doing ions at step ",i6," ------")', T_cntr

        CALL ADVANCE_IONS                      !   velocity: n-N_e_subcycles+1/2 ---> n+1/2
                                               ! coordinate: n-int(N_e_subcycles/2) ---> n-int(N_e_subcycles/2)+N_e_subcycles
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        t11 = MPI_WTIME()

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN  ! send/receive BOTH electrons and ions crossing the borders
           CALL EXCHANGE_PARTICLES_WITH_ABOVE_BELOW_NEIGHBOURS            ! need only one X-pass for self-connected X-periodic clusters
        ELSE                                                              !
           CALL EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS             ! in general, three passes X-Y-X are needed
           CALL EXCHANGE_PARTICLES_WITH_ABOVE_BELOW_NEIGHBOURS            !
           CALL EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS             !
        END IF

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL PROCESS_ADDED_ELECTRONS                ! add the new electrons to the main array   !### NEW       
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        t12 = MPI_WTIME()

        CALL COLLECT_PARTICLE_BOUNDARY_HITS

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        t13 = MPI_WTIME()

        CALL PERFORM_ELECTRON_NEUTRAL_COLLISIONS

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        CALL PERFORM_IONIZATION_HT_SETUP

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        t14 = MPI_WTIME()

        CALL PROCESS_ADDED_IONS                  ! add the new ions to the main array

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        t15 = MPI_WTIME()

        CALL CLEAR_ACCUMULATED_FIELDS
        n_sub = 0

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        t16 = MPI_WTIME()

     ELSE

        t11 = t10

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN  ! send/receive ONLY electrons crossing the borders
           CALL EXCHANGE_ELECTRONS_WITH_ABOVE_BELOW_NEIGHBOURS            ! need only one X-pass for self-connected X-periodic clusters
        ELSE                                                              !
           CALL EXCHANGE_ELECTRONS_WITH_LEFT_RIGHT_NEIGHBOURS             ! in general, three passes X-Y-X are needed
           CALL EXCHANGE_ELECTRONS_WITH_ABOVE_BELOW_NEIGHBOURS            !
           CALL EXCHANGE_ELECTRONS_WITH_LEFT_RIGHT_NEIGHBOURS             !
        END IF

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        t12 = MPI_WTIME()

        CALL COLLECT_ELECTRON_BOUNDARY_HITS

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        t13 = MPI_WTIME()
        t14 = t13
        t15 = t13
        t16 = t13

     END IF

     CALL PERFORM_ELECTRON_EMISSION_HT_SETUP        ! either this or
                                                    ! PERFORM_ELECTRON_EMISSION_HT_SETUP_ZERO_GRAD_F (called above) works
                                                    ! not both
     CALL PERFORM_ELECTRON_EMISSION_SETUP           ! this works for non-HT setup

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     t17 = MPI_WTIME()

     CALL PROCESS_ADDED_ELECTRONS                ! add the new electrons to the main array

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     t18 = MPI_WTIME()

     CALL SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS_HT_SETUP  ! only one of the two will work
     CALL SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS           ! 

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     t19 = MPI_WTIME()

     CALL GATHER_SURFACE_CHARGE_DENSITY      ! 

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     t20 = MPI_WTIME()

     IF (Rank_of_process.EQ.0) PRINT '(2x,i7,2x,f8.3,2x,20(1x,f5.1))', &
          & T_cntr, &
          & REAL(t20 - t0), &
          & 100.0 * REAL((t1  - t0 ) / (t20 - t0)), &
          & 100.0 * REAL((t2  - t1 ) / (t20 - t0)), &
          & 100.0 * REAL((t3  - t2 ) / (t20 - t0)), &
          & 100.0 * REAL((t4  - t3 ) / (t20 - t0)), &
          & 100.0 * REAL((t5  - t4 ) / (t20 - t0)), &
          & 100.0 * REAL((t6  - t5 ) / (t20 - t0)), &
          & 100.0 * REAL((t7  - t6 ) / (t20 - t0)), &
          & 100.0 * REAL((t8  - t7 ) / (t20 - t0)), &
          & 100.0 * REAL((t9  - t8 ) / (t20 - t0)), &
          & 100.0 * REAL((t10 - t9 ) / (t20 - t0)), &
          & 100.0 * REAL((t11 - t10) / (t20 - t0)), &
          & 100.0 * REAL((t12 - t11) / (t20 - t0)), &
          & 100.0 * REAL((t13 - t12) / (t20 - t0)), &
          & 100.0 * REAL((t14 - t13) / (t20 - t0)), &
          & 100.0 * REAL((t15 - t14) / (t20 - t0)), &
          & 100.0 * REAL((t16 - t15) / (t20 - t0)), &
          & 100.0 * REAL((t17 - t16) / (t20 - t0)), &
          & 100.0 * REAL((t18 - t17) / (t20 - t0)), &
          & 100.0 * REAL((t19 - t18) / (t20 - t0)), &
          & 100.0 * REAL((t20 - t19) / (t20 - t0))

  END DO

  finish = MPI_WTIME()

  PRINT '(2x,"**** Process ",i3" : Simulation time is  : ", f12.3," sec")', Rank_of_process, finish - start

  CALL FINISH_SNAPSHOTS

  CALL MPI_FINALIZE(ierr)

END PROGRAM MainProg

!------------------------------------
subroutine try_notexchange_arrays

  USE CurrentProblemValues
  USE ParallelOperationValues
!  USE SendToNeighbor
  USE BlockAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR

  INTEGER i, m, n

  ALLOCATE(rbufer(1:((N_grid_block_x+2)*(N_grid_block_y+2))), STAT=ALLOC_ERR)

  IF (WHITE) THEN  ! WHITE='TRUE' ############################################################################################
! "white processes"

! send right array of charge density
     IF (Rank_of_process_right.GE.0) THEN
! prepare array to send to right neighbor
!        i=0
!        DO n = indx_y_min, indx_y_max
!           do m = indx_x_min, indx_x_max
!              i=i+1
!              rbufer(i) = rho_e(m,n)
!           end do
!        END DO
! send particles to right neighbor
        CALL MPI_SEND(rho_e, (N_grid_block_x+2)*(N_grid_block_y+2), MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process, MPI_COMM_WORLD, request, ierr)     
     END IF

! receive array of charge density from left neighbor
     IF (Rank_of_process_left.GE.0) THEN
! receive particles from left neighbor
        CALL MPI_RECV(rbufer, (N_grid_block_x+2)*(N_grid_block_y+2), MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process_left, MPI_COMM_WORLD, stattus, ierr)     
! process the received particles
        i=0
        DO n = indx_y_min, indx_y_max
           do m = indx_x_min, indx_x_max
              i=i+1
              rho_e(m,n)=rho_e(m,n)+rbufer(i)
           end do
        END DO
     END IF

  ELSE             ! WHITE='FALSE' ############################################################################################
! "black" processes

! receive array of charge density from left neighbor
     IF (Rank_of_process_left.GE.0) THEN
! receive particles from left neighbor
        CALL MPI_RECV(rbufer, (N_grid_block_x+2)*(N_grid_block_y+2), MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process_left, MPI_COMM_WORLD, stattus, ierr)     
! process the received particles
        i=0
        DO n = indx_y_min, indx_y_max
           do m = indx_x_min, indx_x_max
              i=i+1
              rho_e(m,n)=rho_e(m,n)+rbufer(i)
           end do
        END DO
     END IF

! send right array of charge density
     IF (Rank_of_process_right.GE.0) THEN
! prepare array to send to right neighbor
!        i=0
!        DO n = indx_y_min, indx_y_max
!           do m = indx_x_min, indx_x_max
!              i=i+1
!              rbufer(i) = rho_e(m,n)
!           end do
!        END DO
! send particles to right neighbor
        CALL MPI_SEND(rho_e, (N_grid_block_x+2)*(N_grid_block_y+2), MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process, MPI_COMM_WORLD, request, ierr)     
     END IF

  END IF

  DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        

end subroutine try_notexchange_arrays


