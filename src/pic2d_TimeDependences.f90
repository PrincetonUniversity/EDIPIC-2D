!------------------------------
! 
SUBROUTINE INITIATE_PROBE_DIAGNOSTICS
    
  USE ParallelOperationValues
  USE Diagnostics
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries
  USE IonParticles, ONLY : N_spec
  USE Checkpoints, ONLY : use_checkpoint

  use mpi

  IMPLICIT NONE

  INTEGER errcode,ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  LOGICAL exists

  CHARACTER (1) buf

  INTEGER Save_probes_data_step_rff

  INTEGER N_of_given_probes
  INTEGER n, s, k, i
  INTEGER IOS
  INTEGER temp_pos(1:2,1:Max_N_of_probes)         ! 1 = i_x, 2 = j_y 

  INTEGER count_skip

  INTEGER temp_N_of_probes_cluster

  INTEGER npa, npc, npb  ! used to denote index of probe global (a=all), in cluster (c), and in block (b)

  REAL(8) start_time_ns
  REAL(8) time_ns

  CHARACTER(34) clusterprobes_filename    ! _probelocs_cluster_NNNNN_NNNNN.dat
                                          ! ----x----I----x----I----x----I----

  CHARACTER(16) Nis_filename
  CHARACTER(17) NNis_filename

  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  INTEGER bufsize

  INTERFACE
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  END INTERFACE

  INQUIRE (FILE = 'init_probes.dat', EXIST = exists)

  IF (.NOT.exists) THEN
     IF (Rank_of_process.EQ.0) PRINT '(2x,"### init_probes.dat not found, time dependencies in probes will not be created ###")'
     N_of_probes = 0
     N_of_probes_cluster = 0
     N_of_probes_block = 0
     Save_probes_e_data_T_cntr = 0
     Save_probes_i_data_T_cntr = N_subcycles-1
     Save_probes_e_data_step = 1
     Save_probes_i_data_step = N_subcycles
     Save_probes_data_step = 1   ! have to keep this for INITIATE_SNAPSHOT
     text_output_counter = 0  ! obsolete values, keep them to have something to use in checkpoints
     N_of_saved_records = 0
     RETURN
  END IF

  IF (Rank_of_process.EQ.0) THEN
     PRINT '(2x,"Process ",i5," : init_probes.dat is found. Reading the data file...")', Rank_of_process
  END IF

  OPEN (9, FILE = 'init_probes.dat')

  READ (9, '(A1)') buf !Step for saving (ion timesteps if >=1, 0 means save as frequent as possible), type below
  READ (9, *) Save_probes_data_step_rff
  READ (9, '(A1)') buf !Start saving data at (timesteps, >=0), type below
  READ (9, *) Save_probes_data_T_cntr_rff
  READ (9, '(A1)') buf !Skip periods of writing between text outputs (>=0), type below
  READ (9, *) TextOut_skip
  READ (9, '(A1)') buf !Number of probes (off if <=0), type below
  READ (9, *) N_of_given_probes
  READ (9, '(A1)') buf !Probe coordinates (x/y, node numbers), type below

  IF (N_of_given_probes.GT.Max_N_of_probes) THEN
     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i5," : requested number of probes ",i4," was reduced to maximal allowed value of ",i4)', Rank_of_process, Max_N_of_probes
     END IF
     N_of_given_probes = Max_N_of_probes
  END IF
     
  DO n = 1, N_of_given_probes
     READ (9, *, IOSTAT=IOS) temp_pos(1,n), temp_pos(2,n)
     IF (IOS.NE.0) THEN
        PRINT '(2x,"Process ",i5," :: INITIATE_PROBE_DIAGNOSTICS : ERROR-1 : while reading file init_probes.dat : wrong coordinates of probe ",i4," program terminated.")', Rank_of_process, n
        errcode=420
        CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
     END IF
  END DO

  CLOSE (9, STATUS = 'KEEP')

  IF (Save_probes_data_step_rff.LE.0) THEN
     Save_probes_e_data_step = 1
     Save_probes_i_data_step = N_subcycles

     Save_probes_e_data_T_cntr = 0
     Save_probes_i_data_T_cntr = N_subcycles-1
  ELSE
     Save_probes_e_data_step = Save_probes_data_step_rff * N_subcycles
     Save_probes_i_data_step = Save_probes_e_data_step

     Save_probes_e_data_T_cntr = N_subcycles/2
     Save_probes_i_data_T_cntr = N_subcycles-1
  END IF

  Save_probes_data_step = Save_probes_e_data_step   ! have to keep this for INITIATE_SNAPSHOT

  text_output_counter = 0  ! obsolete values, keep them to have something to use in checkpoints
  N_of_saved_records = 0

  Save_probes_data_T_cntr_rff = (MAX(0, Save_probes_data_T_cntr_rff) / N_subcycles) * N_subcycles    ! have to keep this for INITIATE_SNAPSHOT

  Save_probes_e_data_T_cntr = Save_probes_e_data_T_cntr + Save_probes_data_T_cntr_rff
  Save_probes_i_data_T_cntr = Save_probes_i_data_T_cntr + Save_probes_data_T_cntr_rff

  DO WHILE (Save_probes_e_data_T_cntr.LT.Start_T_cntr)
     Save_probes_e_data_T_cntr = Save_probes_e_data_T_cntr + Save_probes_e_data_step
  END DO

  DO WHILE (Save_probes_i_data_T_cntr.LT.Start_T_cntr)
     Save_probes_i_data_T_cntr = Save_probes_i_data_T_cntr + Save_probes_i_data_step
  END DO

! there are two situations:
! N_of_given_probes<=0 : no any time dependencies will be created at all
! N_of_given_probes>0  : time dependencies will be created in probes if there is at least one probe with valid coordinates

  N_of_probes = 0
  count_skip = 0

  IF (N_of_given_probes.GT.0) THEN
! define number of legal probes
     DO n = 1, N_of_given_probes
        IF ( (temp_pos(1,n).LT.0).OR. &
           & (temp_pos(1,n).GT.global_maximal_i).OR. &
           & (temp_pos(2,n).LT.0).OR. &
           & (temp_pos(2,n).GT.global_maximal_j) ) THEN
! skip probes with coordinates out of the simulation domain
           count_skip = count_skip + 1
           CYCLE
        END IF
! coordinates of the requested probe are valid
        N_of_probes = N_of_probes + 1
     END DO
  END IF

  IF (count_skip.GT.0) THEN
     IF (Rank_of_process.EQ.0) PRINT '("### While reading init_probes.dat skipped ",i3," invalid probe locations... ###")', count_skip
  END IF

  IF (N_of_probes.EQ.0) THEN
     N_of_probes_cluster = 0
     N_of_probes_block = 0
     IF (Rank_of_process.EQ.0) PRINT '("###  Temporal dependence files will NOT be created. ###")'
     RETURN
  END IF

! since we are here, there is at least one valid probe location

! each process stores information about all probes
! this makes it easier to assemble/transfer information about probes inside some particular cluster or block

  ALLOCATE(Probe_position(1:2,1:N_of_probes), STAT = ALLOC_ERR)
  npa = 0                      ! npa = number of probes all
  DO n = 1, N_of_given_probes

     IF ( (temp_pos(1,n).LT.0).OR. &
        & (temp_pos(1,n).GT.global_maximal_i).OR. &
        & (temp_pos(2,n).LT.0).OR. &
        & (temp_pos(2,n).GT.global_maximal_j) ) THEN
! skip probes with coordinates out of the simulation domain
        CYCLE
     END IF

     npa = npa+1
     Probe_position(1:2, npa) = temp_pos(1:2,n)

  END DO

! report
  IF (Rank_of_process.EQ.0) THEN

     PRINT '("### In total, there is ",i3," valid probe locations ###")', N_of_probes

! save probe positions
     OPEN (20, FILE = '_probelocs.dat')
     WRITE (20, '("# col  1 is the global probe number")')
     WRITE (20, '("# col  2 is index i of probe position in computational grid (along x)")')
     WRITE (20, '("# col  3 is index j of probe position in computational grid (along y)")')
     WRITE (20, '("# col  4 is probe cartesian coordinate x [mm]")')
     WRITE (20, '("# col  5 is probe cartesian coordinate y [mm]")')
     WRITE (20, '("# ----------")')
     DO n = 1, N_of_probes
        WRITE (20, '(2x,i3,2(2x,i5),2(2x,e14.7))') &
             & n, &
             & Probe_position(1,n), &
             & Probe_position(2,n), &
             & Probe_position(1,n) * delta_x_m * 1000.0_8, &
             & Probe_position(2,n) * delta_x_m * 1000.0_8
     END DO
     CLOSE (20, STATUS = 'KEEP')
     PRINT '("### file _probelocs.dat is ready ###")'

  END IF
        
! master processes will keep their own probes only
! note, N_of_probes_cluster = 0 is in INITIATE_PARAMETERS to ensure correct call of DISTRIBUTE_CLUSTER_PARAMETERS 
! note, if a probe is at the boundary between two clusters, it will be included in both clusters
!
! in this case, at the final assembly stage 
! for electric field or potential the data from one cluster will be overwritten by the data from the other cluster
! while for the electron and ion densities a sum of data from both (or as many as 4 at the corners) clusters will be taken
! this eliminates the need for additional communications to account for particles of the same cell split between multiple processes
! 
  IF (cluster_rank_key.EQ.0) THEN
! master process checks whether the probe belongs to its domain

     DO npa = 1, N_of_probes
        IF ( (Probe_position(1,npa).GE.c_indx_x_min).AND. &
           & (Probe_position(1,npa).LE.c_indx_x_max).AND. &
           & (Probe_position(2,npa).GE.c_indx_y_min).AND. &
           & (Probe_position(2,npa).LE.c_indx_y_max) ) THEN
           N_of_probes_cluster = N_of_probes_cluster + 1
        END IF
     END DO

     PRINT '("Cluster ",i5," : has ",i3," probes")', Rank_of_process, N_of_probes_cluster

     IF (Rank_horizontal.EQ.0) THEN
        ALLOCATE(probe_JXsum(1:N_of_probes), STAT = ALLOC_ERR)
        ALLOCATE(probe_JYsum(1:N_of_probes), STAT = ALLOC_ERR)
        ALLOCATE(probe_JZsum(1:N_of_probes), STAT = ALLOC_ERR)
     END IF

     IF (N_of_probes_cluster.GT.0) THEN

        ALLOCATE(List_of_probes_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)

        ALLOCATE(probe_F_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)

        ALLOCATE(probe_Ne_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)

        ALLOCATE(probe_VXe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        ALLOCATE(probe_VYe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        ALLOCATE(probe_VZe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)

        ALLOCATE(probe_VXVYe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        ALLOCATE(probe_VXVZe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        ALLOCATE(probe_VYVZe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)

        ALLOCATE(probe_WXe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        ALLOCATE(probe_WYe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        ALLOCATE(probe_WZe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)

        ALLOCATE(probe_QXe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        ALLOCATE(probe_QYe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        ALLOCATE(probe_QZe_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)

        ALLOCATE(probe_Ni_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)

        ALLOCATE(probe_VXi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)
        ALLOCATE(probe_VYi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)
        ALLOCATE(probe_VZi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)

        ALLOCATE(probe_VXVYi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)
        ALLOCATE(probe_VXVZi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)
        ALLOCATE(probe_VYVZi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)

        ALLOCATE(probe_WXi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)
        ALLOCATE(probe_WYi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)
        ALLOCATE(probe_WZi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)

        ALLOCATE(probe_QXi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)
        ALLOCATE(probe_QYi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)
        ALLOCATE(probe_QZi_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)

        npc = 0
        DO npa = 1, N_of_probes

           IF ( (Probe_position(1,npa).GE.c_indx_x_min).AND. &
              & (Probe_position(1,npa).LE.c_indx_x_max).AND. &
              & (Probe_position(2,npa).GE.c_indx_y_min).AND. &
              & (Probe_position(2,npa).LE.c_indx_y_max) ) THEN

              npc = npc+1
              List_of_probes_cluster(npc) = npa
              
           END IF
        END DO

! save probe positions
        clusterprobes_filename = '_probelocs_cluster_NNNNN_NNNNN.dat'
!                                 ----x----I----x----I----x----I----
        clusterprobes_filename(20:24) = convert_int_to_txt_string(Rank_of_process, 5)
        clusterprobes_filename(26:30) = convert_int_to_txt_string(Rank_horizontal, 5)

        OPEN (20, FILE = clusterprobes_filename)
        WRITE (20, '("# col  1 is the global probe number")')
        WRITE (20, '("# col  2 is index i of probe position in computational grid (along x)")')
        WRITE (20, '("# col  3 is index j of probe position in computational grid (along y)")')
        WRITE (20, '("# col  4 is probe cartesian coordinate x [mm]")')
        WRITE (20, '("# col  5 is probe cartesian coordinate y [mm]")')
        WRITE (20, '("# ----------")')
        DO npc = 1, N_of_probes_cluster
           npa = List_of_probes_cluster(npc)
           WRITE (20, '(2x,i3,2(2x,i5),2(2x,e14.7))') &
                & npa, &
                & Probe_position(1,npa), &
                & Probe_position(2,npa), &
                & Probe_position(1,npa) * delta_x_m * 1000.0_8, &
                & Probe_position(2,npa) * delta_x_m * 1000.0_8
        END DO
        CLOSE (20, STATUS = 'KEEP')
        PRINT '("Process ",i5," :: file ",A34," is ready")', Rank_of_process, clusterprobes_filename

     END IF

! particle calculators are already assigned, so we must tell them about the probes in the master process

     CALL MPI_BCAST(N_of_probes_cluster, 1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

     IF (N_of_probes_cluster.GT.0) THEN
        CALL MPI_BCAST(List_of_probes_cluster, N_of_probes_cluster, MPI_INTEGER, 0, COMM_CLUSTER, ierr)        
     END IF
     

  ELSE   !###   IF (cluster_rank_key.EQ.0) THEN
! receive number of probes in cluster from the master
     CALL MPI_BCAST(N_of_probes_cluster, 1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

     IF (N_of_probes_cluster.GT.0) THEN

        ALLOCATE(List_of_probes_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        CALL MPI_BCAST(List_of_probes_cluster, N_of_probes_cluster, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
     END IF

  END IF   !###   IF (cluster_rank_key.EQ.0) THEN

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

! field calculators will be providing electrostatic potential (unless FFT-based solver is on) and also must identify which probes belong to their domains
  
     ALLOCATE(ibufer(1:N_of_probes), STAT = ALLOC_ERR)

     IF (cluster_rank_key.EQ.0) THEN

! account for that a master is a field calculator as well
! identify how many probes from the list are inside the field calculator domain
        N_of_probes_block = 0
        DO npc = 1, N_of_probes_cluster
           npa = List_of_probes_cluster(npc)
           IF ( (Probe_position(1,npa).GE.indx_x_min).AND. &
              & (Probe_position(1,npa).LE.indx_x_max).AND. &
              & (Probe_position(2,npa).GE.indx_y_min).AND. &
              & (Probe_position(2,npa).LE.indx_y_max) ) THEN
              N_of_probes_block = N_of_probes_block + 1
           END IF
        END DO
! save probes from the list which are inside the field calculator domain
        IF (N_of_probes_block.GT.0) THEN
           ALLOCATE(Probe_params_block_list(1:3,1:N_of_probes_block), STAT = ALLOC_ERR)
           npb = 0
           DO npc = 1, N_of_probes_cluster
              npa = List_of_probes_cluster(npc)
              IF ( (Probe_position(1,npa).GE.indx_x_min).AND. &
                 & (Probe_position(1,npa).LE.indx_x_max).AND. &
                 & (Probe_position(2,npa).GE.indx_y_min).AND. &
                 & (Probe_position(2,npa).LE.indx_y_max) ) THEN
                 npb = npb + 1
! field calculators report the potential to the master so it is better to actually store probe coordinates 
                 Probe_params_block_list(1,npb) = Probe_position(1,npa)
                 Probe_params_block_list(2,npb) = Probe_position(2,npa)
                 Probe_params_block_list(3,npb) = npc
              END IF
           END DO
! allocate array to save the potential
           ALLOCATE(probe_F_block(1:N_of_probes_block), STAT = ALLOC_ERR)
        END IF

! master sends to the associated field calculators number of its probes
        ibufer(1) = N_of_probes_cluster
        DO k = 2, cluster_N_blocks
           CALL MPI_SEND(ibufer(1), 1, MPI_INTEGER, field_calculator(k)%rank, Rank_of_process, MPI_COMM_WORLD, ierr) 
        END DO

        IF (N_of_probes_cluster.GT.0) THEN
! master sends to the associated field calculators a list of its probes
           DO k = 2, cluster_N_blocks
              CALL MPI_SEND(List_of_probes_cluster(1:N_of_probes_cluster), N_of_probes_cluster, MPI_INTEGER, field_calculator(k)%rank, Rank_of_process+100000, MPI_COMM_WORLD, ierr) 
           END DO
! receive number of probes in each block
           ALLOCATE(diag_block(2:cluster_N_blocks), STAT = ALLOC_ERR)
           DO k = 2, cluster_N_blocks
              CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)
              diag_block(k)%N_probes = ibufer(1)
           END DO
! receive lists of probe positions in the cluster's probe list
           DO k = 2, cluster_N_blocks
              IF (diag_block(k)%N_probes.GT.0) THEN
                 bufsize = diag_block(k)%N_probes
                 ALLOCATE(diag_block(k)%probe_number(1:bufsize), STAT = ALLOC_ERR)
                 CALL MPI_RECV(ibufer(1:bufsize), bufsize, MPI_INTEGER, field_calculator(k)%rank, field_calculator(k)%rank+100000, MPI_COMM_WORLD, stattus, ierr)
                 diag_block(k)%probe_number(1:bufsize) = ibufer(1:bufsize)
              END IF
           END DO
        END IF

     ELSE   !###   IF (cluster_rank_key.EQ.0) THEN

! field calculators receive number of probes from the master
        CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, field_master, field_master, MPI_COMM_WORLD, stattus, ierr)
        temp_N_of_probes_cluster = ibufer(1)  ! here temp_N_of_probes_cluster is a temporary value which will be replaced
                                         ! once the process becomes associated with a cluster as a particle calculator

        N_of_probes_block = 0   ! default assumption

        IF (temp_N_of_probes_cluster.GT.0) THEN
! field calculators receive list of probes from the master
           CALL MPI_RECV(ibufer(1:temp_N_of_probes_cluster), temp_N_of_probes_cluster, MPI_INTEGER, field_master, field_master+100000, MPI_COMM_WORLD, stattus, ierr)
! identify how many probes from the list are inside the field calculator domain
           DO npc = 1, temp_N_of_probes_cluster
              npa = ibufer(npc)
              IF ( (Probe_position(1,npa).GE.indx_x_min).AND. &
                 & (Probe_position(1,npa).LE.indx_x_max).AND. &
                 & (Probe_position(2,npa).GE.indx_y_min).AND. &
                 & (Probe_position(2,npa).LE.indx_y_max) ) THEN
                 N_of_probes_block = N_of_probes_block + 1
              END IF
           END DO
! save those probes from the list which are inside the field calculator domain
           IF (N_of_probes_block.GT.0) THEN
              ALLOCATE(Probe_params_block_list(1:3,1:N_of_probes_block), STAT = ALLOC_ERR)
              npb = 0
              DO npc = 1, temp_N_of_probes_cluster
                 npa = ibufer(npc)
                 IF ( (Probe_position(1,npa).GE.indx_x_min).AND. &
                    & (Probe_position(1,npa).LE.indx_x_max).AND. &
                    & (Probe_position(2,npa).GE.indx_y_min).AND. &
                    & (Probe_position(2,npa).LE.indx_y_max) ) THEN
                    npb = npb + 1
                    Probe_params_block_list(1,npb) = Probe_position(1,npa)
                    Probe_params_block_list(2,npb) = Probe_position(2,npa)
                    Probe_params_block_list(3,npb) = npc
                 END IF
              END DO
           END IF
! field calculators tell their master about the probes that they have
! report the probe number
           ibufer(1) = N_of_probes_block
           CALL MPI_SEND(ibufer(1), 1, MPI_INTEGER, field_master, Rank_of_process, MPI_COMM_WORLD, ierr) 
! report the vector of probe order numbers in the cluster's probe list 
           IF (N_of_probes_block.GT.0) THEN
              ibufer(1:N_of_probes_block) = Probe_params_block_list(3,1:N_of_probes_block)
              CALL MPI_SEND(ibufer(1:N_of_probes_block), N_of_probes_block, MPI_INTEGER, field_master, Rank_of_process+100000, MPI_COMM_WORLD, ierr) 
! allocate array to save the potential
              ALLOCATE(probe_F_block(1:N_of_probes_block), STAT = ALLOC_ERR)
           END IF
        END IF

     END IF   !###   IF (cluster_rank_key.EQ.0) THEN

     DEALLOCATE(ibufer, STAT = ALLOC_ERR)

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  END IF   !###   IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

! now all clusters must tell process 0 about the probes they have 

  IF (cluster_rank_key.EQ.0) THEN

     ALLOCATE(ibufer(1:N_of_probes), STAT = ALLOC_ERR)

     IF (Rank_horizontal.EQ.0) THEN    ! this process has Rank_of_process=0 as well

        ALLOCATE(diag_cluster(1:N_processes_horizontal-1), STAT=ALLOC_ERR)
! receive number of probes from each cluster
        DO n = 1, N_processes_horizontal-1
           CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, n, n, COMM_HORIZONTAL, stattus, ierr)
           diag_cluster(n)%N_probes = ibufer(1)
        END DO

! receive sets of indices [in the global array of probes] of probes from each cluster
        DO n = 1, N_processes_horizontal-1
           IF (diag_cluster(n)%N_probes.LE.0) CYCLE
           bufsize = diag_cluster(n)%N_probes
           CALL MPI_RECV(ibufer(1:bufsize), bufsize, MPI_INTEGER, n, n+100000, COMM_HORIZONTAL, stattus, ierr)
           ALLOCATE(diag_cluster(n)%probe_number(1:bufsize), STAT=ALLOC_ERR)
           diag_cluster(n)%probe_number(1:bufsize) = ibufer(1:bufsize)
        END DO

        IF (use_checkpoint.EQ.1) THEN
! start from checkpoint, must trim the time dependences

           start_time_ns = (DBLE(Start_T_cntr) - 0.5_8) * (delta_t_s * 1.0d9)  ! -0.5 to avoid uncertainty due to possible round off error

! electric field EX 
           INQUIRE (FILE = 'dim_EX_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_EX_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electric field EY
           INQUIRE (FILE = 'dim_EY_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_EY_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electrostatic potential F 
           INQUIRE (FILE = 'dim_F_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_F_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! full electric current density JXsum
           INQUIRE (FILE = 'dim_JXsum_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_JXsum_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! full electric current density JYsum
           INQUIRE (FILE = 'dim_JYsum_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_JYsum_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! full electric current density JZsum
           INQUIRE (FILE = 'dim_JZsum_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_JZsum_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron density Ne 
           INQUIRE (FILE = 'dim_Ne_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_Ne_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron flow velocity VXe 
           INQUIRE (FILE = 'dim_VXe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_VXe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron flow velocity VYe 
           INQUIRE (FILE = 'dim_VYe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_VYe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron flow velocity VZe 
           INQUIRE (FILE = 'dim_VZe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_VZe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electric current density due to electrons JXe 
           INQUIRE (FILE = 'dim_JXe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_JXe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electric current density due to electrons JYe 
           INQUIRE (FILE = 'dim_JYe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_JYe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electric current density due to electrons JZe 
           INQUIRE (FILE = 'dim_JZe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_JZe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! average electron energy WXe 
           INQUIRE (FILE = 'dim_WXe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_WXe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! average electron energy WYe 
           INQUIRE (FILE = 'dim_WYe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_WYe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! average electron energy WZe 
           INQUIRE (FILE = 'dim_WZe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_WZe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron temperature TXe
           INQUIRE (FILE = 'dim_TXe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_TXe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron temperature TYe 
           INQUIRE (FILE = 'dim_TYe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_TYe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron temperature TZe 
           INQUIRE (FILE = 'dim_TZe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_TZe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron heat flow QXe
           INQUIRE (FILE = 'dim_QXe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_QXe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron heat flow QYe
           INQUIRE (FILE = 'dim_QYe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_QYe_vst.dat', STATUS = 'OLD')          
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron heat flow QZe
           INQUIRE (FILE = 'dim_QZe_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_QZe_vst.dat', STATUS = 'OLD')          
              DO
                 READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                 IF (ios.NE.0) EXIT
                 IF (time_ns.GE.start_time_ns) EXIT
              END DO
              BACKSPACE(21)
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! ion species data
           DO s = 1, N_spec

! density
              Nis_filename = 'dim_Ni_S_vst.dat'
              Nis_filename(8:8) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = Nis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = Nis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! flow velocity VX
              NNis_filename = 'dim_VXi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! flow velocity VY
              NNis_filename = 'dim_VYi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! flow velocity VZ
              NNis_filename = 'dim_VZi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! electric curent density due to ions JX
              NNis_filename = 'dim_JXi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! electric curent density due to ions JY
              NNis_filename = 'dim_JYi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! electric curent density due to ions JZ
              NNis_filename = 'dim_JZi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! average ion energy WX
              NNis_filename = 'dim_WXi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! average ion energy WY
              NNis_filename = 'dim_WYi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! average ion energy WZ
              NNis_filename = 'dim_WZi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! ion temperature TX
              NNis_filename = 'dim_TXi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! ion temperature TY
              NNis_filename = 'dim_TYi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! ion temperature TZ
              NNis_filename = 'dim_TZi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! ion heat flow QX
              NNis_filename = 'dim_QXi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! ion heat flow QY
              NNis_filename = 'dim_QYi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

! ion heat flow QZ
              NNis_filename = 'dim_QZi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = NNis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = NNis_filename, STATUS = 'OLD')          
                 DO
                    READ (21, '(1x,f15.6,100(1x,e14.7))', iostat = ios) time_ns
                    IF (ios.NE.0) EXIT
                    IF (time_ns.GE.start_time_ns) EXIT
                 END DO
                 BACKSPACE(21)
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF

           END DO
 
        ELSE
! fresh start, empty files, remove old files if there are any

! electric field EX 
           OPEN  (21, FILE = 'dim_EX_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electric field EY
           OPEN  (21, FILE = 'dim_EY_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electrostatic potential F 
           OPEN  (21, FILE = 'dim_F_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! full electric current density JXsum
           OPEN  (21, FILE = 'dim_JXsum_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! full electric current density JYsum
           OPEN  (21, FILE = 'dim_JYsum_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! full electric current density JZsum
           OPEN  (21, FILE = 'dim_JZsum_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron density Ne 
           OPEN  (21, FILE = 'dim_Ne_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron flow velocity VXe 
           OPEN  (21, FILE = 'dim_VXe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron flow velocity VYe 
           OPEN  (21, FILE = 'dim_VYe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron flow velocity VZe 
           OPEN  (21, FILE = 'dim_VZe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electric current density due to electrons JXe 
           OPEN  (21, FILE = 'dim_JXe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electric current density due to electrons JYe 
           OPEN  (21, FILE = 'dim_JYe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electric current density due to electrons JZe 
           OPEN  (21, FILE = 'dim_JZe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! average electron energy WXe 
           OPEN  (21, FILE = 'dim_WXe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! average electron energy WYe 
           OPEN  (21, FILE = 'dim_WYe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! average electron energy WZe 
           OPEN  (21, FILE = 'dim_WZe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron temperature TXe 
           OPEN  (21, FILE = 'dim_TXe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron temperature TYe 
           OPEN  (21, FILE = 'dim_TYe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron temperature TZe 
           OPEN  (21, FILE = 'dim_TZe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron heat flow QXe 
           OPEN  (21, FILE = 'dim_QXe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron heat flow QYe 
           OPEN  (21, FILE = 'dim_QYe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron heat flow QZe 
           OPEN  (21, FILE = 'dim_QZe_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! ion species data
           DO s = 1, N_spec

! density
              Nis_filename = 'dim_Ni_S_vst.dat'
              Nis_filename(8:8) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = Nis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! flow velocity VX
              NNis_filename = 'dim_VXi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! flow velocity VY
              NNis_filename = 'dim_VYi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! flow velocity VZ
              NNis_filename = 'dim_VZi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! electric curent density due to ions JX
              NNis_filename = 'dim_JXi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! electric curent density due to ions JY
              NNis_filename = 'dim_JYi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! electric curent density due to ions JZ
              NNis_filename = 'dim_JZi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! average ion energy WX
              NNis_filename = 'dim_WXi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! average ion energy WY
              NNis_filename = 'dim_WYi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! average ion energy WZ
              NNis_filename = 'dim_WZi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! ion temperature TX
              NNis_filename = 'dim_TXi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! ion temperature TY
              NNis_filename = 'dim_TYi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! ion temperature TZ
              NNis_filename = 'dim_TZi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! ion heat flow QX
              NNis_filename = 'dim_QXi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! ion heat flow QY
              NNis_filename = 'dim_QYi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

! ion heat flow QZ
              NNis_filename = 'dim_QZi_S_vst.dat'
              NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = NNis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')

           END DO
        END IF

     ELSE

        ibufer(1) = N_of_probes_cluster
! send number of probes to the master with zero rank     
        CALL MPI_SEND(ibufer(1:1), 1, MPI_INTEGER, 0, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        
! send global indices of probes to the master with zero rank
        IF (N_of_probes_cluster.GT.0) THEN
           CALL MPI_SEND(List_of_probes_cluster, N_of_probes_cluster, MPI_INTEGER, 0, Rank_horizontal+100000, COMM_HORIZONTAL, ierr) 
        END IF

     END IF   !###   IF (Rank_horizontal.EQ.0) THEN 

     DEALLOCATE(ibufer, STAT = ALLOC_ERR)

  END IF   !###   IF (cluster_rank_key.EQ.0) THEN

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!if (Rank_of_process.eq.0 )print *, "Save_probes_data_T_cntr is ", Save_probes_data_T_cntr

END SUBROUTINE INITIATE_PROBE_DIAGNOSTICS

!-------------------------------------------------------
!
SUBROUTINE DO_PROBE_DIAGNOSTICS

  USE CurrentProblemValues, ONLY : T_cntr
  USE Diagnostics, ONLY : N_of_probes, Save_probes_e_data_T_cntr, Save_probes_e_data_step, Save_probes_i_data_T_cntr, Save_probes_i_data_step

  use mpi

  IMPLICIT NONE

  INTEGER ierr

  IF (N_of_probes.LE.0) RETURN

  IF (T_cntr.EQ.Save_probes_e_data_T_cntr) THEN
     CALL DO_PROBE_DIAGNOSTICS_e_DATA
     Save_probes_e_data_T_cntr = Save_probes_e_data_T_cntr + Save_probes_e_data_step 
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  IF (T_cntr.EQ.Save_probes_i_data_T_cntr) THEN
     CALL DO_PROBE_DIAGNOSTICS_i_DATA
     Save_probes_i_data_T_cntr = Save_probes_i_data_T_cntr + Save_probes_i_data_step 
  END IF

END SUBROUTINE DO_PROBE_DIAGNOSTICS

!------------------------------
!
SUBROUTINE DO_PROBE_DIAGNOSTICS_e_DATA
    
  USE ParallelOperationValues
  USE Diagnostics
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  use mpi

  IMPLICIT NONE  

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER bufsize
  REAL, ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR
  INTEGER pos1, pos2

  INTEGER k, n, i, j
  INTEGER npa, npb, npc

  REAL, ALLOCATABLE :: probe_EX(:)
  REAL, ALLOCATABLE :: probe_EY(:)
  REAL, ALLOCATABLE :: probe_F(:)

  REAL, ALLOCATABLE :: probe_Ne(:)

  REAL, ALLOCATABLE :: probe_JXe(:)
  REAL, ALLOCATABLE :: probe_JYe(:)
  REAL, ALLOCATABLE :: probe_JZe(:)

  REAL, ALLOCATABLE :: probe_WXe(:)
  REAL, ALLOCATABLE :: probe_WYe(:)
  REAL, ALLOCATABLE :: probe_WZe(:)

  REAL, ALLOCATABLE :: probe_VXe(:)
  REAL, ALLOCATABLE :: probe_VYe(:)
  REAL, ALLOCATABLE :: probe_VZe(:)

  REAL, ALLOCATABLE :: probe_VXVYe(:)
  REAL, ALLOCATABLE :: probe_VXVZe(:)
  REAL, ALLOCATABLE :: probe_VYVZe(:)

  REAL, ALLOCATABLE :: probe_TXe(:)
  REAL, ALLOCATABLE :: probe_TYe(:)
  REAL, ALLOCATABLE :: probe_TZe(:)

  REAL, ALLOCATABLE :: probe_QXe(:)
  REAL, ALLOCATABLE :: probe_QYe(:)
  REAL, ALLOCATABLE :: probe_QZe(:)

  INTEGER recsize

  REAL rvx, rvy, rvz
  REAL rvx2, rvy2, rvz2
  REAL rtemp, inv_N

  REAL time_ns

! collect potential from field calculators (for systems which do not use FFT solver) ---------------------------------------------

  IF ((periodicity_flag.EQ.PERIODICITY_NONE).OR.(periodicity_flag.EQ.PERIODICITY_X_PETSC).OR.(periodicity_flag.EQ.PERIODICITY_X_Y)) THEN

     IF (cluster_rank_key.EQ.0) THEN
        IF (N_of_probes_cluster.GT.0) THEN
! account for that a master process is a field calculator as well
           DO npb = 1, N_of_probes_block
              npc = Probe_params_block_list(3, npb)
              probe_F_cluster(npc) = probe_F_block(npb)
           END DO
! now consider other field calculators
           DO k = 2, cluster_N_blocks
              IF (diag_block(k)%N_probes.LE.0) CYCLE
              bufsize = diag_block(k)%N_probes
              ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)
              CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_REAL, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)
              DO npb = 1, diag_block(k)%N_probes
                 npc = diag_block(k)%probe_number(npb)
                 probe_F_cluster(npc) = rbufer(npb)
              END DO
              DEALLOCATE(rbufer, STAT = ALLOC_ERR)
           END DO
        END IF
     ELSE
        IF (N_of_probes_block.GT.0) THEN
           CALL MPI_SEND(probe_F_block, N_of_probes_block, MPI_REAL, field_master, Rank_of_process, MPI_COMM_WORLD, ierr) 
        END IF
     END IF
  
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  IF (cluster_rank_key.NE.0) RETURN

! the rest is carried out by masters only

  IF (Rank_horizontal.EQ.0) THEN
! the zero-rank master process assembles all probe data 

     ALLOCATE(probe_EX(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_EY(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_F( 1:N_of_probes), STAT = ALLOC_ERR)

     ALLOCATE(probe_Ne(1:N_of_probes), STAT = ALLOC_ERR)

     ALLOCATE(probe_JXe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_JYe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_JZe(1:N_of_probes), STAT = ALLOC_ERR)

     ALLOCATE(probe_VXe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_VYe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_VZe(1:N_of_probes), STAT = ALLOC_ERR)

     ALLOCATE(probe_VXVYe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_VXVZe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_VYVZe(1:N_of_probes), STAT = ALLOC_ERR)

     ALLOCATE(probe_WXe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_WYe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_WZe(1:N_of_probes), STAT = ALLOC_ERR)

     ALLOCATE(probe_TXe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_TYe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_TZe(1:N_of_probes), STAT = ALLOC_ERR)

     ALLOCATE(probe_QXe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_QYe(1:N_of_probes), STAT = ALLOC_ERR)
     ALLOCATE(probe_QZe(1:N_of_probes), STAT = ALLOC_ERR)

! values of moments N, V, W, V*V*, Q in probes are already calculated but are NOT synchronized in overlapping nodes

     probe_Ne = 0.0

     probe_VXe = 0.0
     probe_VYe = 0.0
     probe_VZe = 0.0

     probe_WXe = 0.0
     probe_WYe = 0.0
     probe_WZe = 0.0

     probe_VXVYe = 0.0
     probe_VXVZe = 0.0
     probe_VYVZe = 0.0

     probe_QXe = 0.0
     probe_QYe = 0.0
     probe_QZe = 0.0

! includes its own probes (if there are any)
     DO npc = 1, N_of_probes_cluster

        npa = List_of_probes_cluster(npc)    ! npa = n-umber of p-robe, a-ll   npc = n-umber of p-robe in c-luster

        i = Probe_position(1, npa)
        j = Probe_position(2, npa)

        probe_EX(npa) = EX(i, j)
        probe_EY(npa) = EY(i, j)
        probe_F(npa) = probe_F_cluster(npc)

        probe_Ne(npa) = probe_Ne_cluster(npc)

        probe_VXe(npa) = probe_VXe_cluster(npc)
        probe_VYe(npa) = probe_VYe_cluster(npc)
        probe_VZe(npa) = probe_VZe_cluster(npc)

        probe_WXe(npa) = probe_WXe_cluster(npc)
        probe_WYe(npa) = probe_WYe_cluster(npc)
        probe_WZe(npa) = probe_WZe_cluster(npc)

        probe_VXVYe(npa) = probe_VXVYe_cluster(npc)
        probe_VXVZe(npa) = probe_VXVZe_cluster(npc)
        probe_VYVZe(npa) = probe_VYVZe_cluster(npc)

        probe_QXe(npa) = probe_QXe_cluster(npc)
        probe_QYe(npa) = probe_QYe_cluster(npc)
        probe_QZe(npa) = probe_QZe_cluster(npc)

     END DO

! assemble probe data from all other clusters
     DO n = 1, N_processes_horizontal-1

        IF (diag_cluster(n)%N_probes.LE.0) CYCLE

        recsize = diag_cluster(n)%N_probes
        bufsize = (3 + 13) * recsize   !{EX,EY,F,  e::N,VX,VY,VZ,WX,WY,WZ,QX,QY,QZ,VXVY,VXVZ,VYVZ}
        ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)

        CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_REAL, n, n, COMM_HORIZONTAL, stattus, ierr)

! since moments in overlapping areas are not synchronized, values in probes belonging to more than one clusters are accumulated
        DO npc = 1, recsize
           npa = diag_cluster(n)%probe_number(npc)

           probe_EX(npa) = rbufer(npc)
           probe_EY(npa) = rbufer(npc+recsize)
           probe_F(npa)  = rbufer(npc+2*recsize)

           probe_Ne(npa) = probe_Ne(npa) + rbufer(npc+3*recsize) 

           probe_VXe(npa) = probe_VXe(npa) + rbufer(npc+4*recsize)
           probe_VYe(npa) = probe_VYe(npa) + rbufer(npc+5*recsize)
           probe_VZe(npa) = probe_VZe(npa) + rbufer(npc+6*recsize)

           probe_WXe(npa) = probe_WXe(npa) + rbufer(npc+7*recsize)
           probe_WYe(npa) = probe_WYe(npa) + rbufer(npc+8*recsize)
           probe_WZe(npa) = probe_WZe(npa) + rbufer(npc+9*recsize)

           probe_QXe(npa) = probe_QXe(npa) + rbufer(npc+10*recsize)
           probe_QYe(npa) = probe_QYe(npa) + rbufer(npc+11*recsize)
           probe_QZe(npa) = probe_QZe(npa) + rbufer(npc+12*recsize)

           probe_VXVYe(npa) = probe_VXVYe(npa) + rbufer(npc+13*recsize)
           probe_VXVZe(npa) = probe_VXVZe(npa) + rbufer(npc+14*recsize)
           probe_VYVZe(npa) = probe_VYVZe(npa) + rbufer(npc+15*recsize)
        END DO   !###  DO npc = 1, recsize 
           
        DEALLOCATE(rbufer, STAT = ALLOC_ERR)

     END DO   !### DO n = 1, N_processes_horizontal-1

! finalize moment calculation
     DO npa = 1, N_of_probes

        IF (probe_Ne(npa).GT.1.0e-9) THEN    ! note this is small but not zero

           inv_N = 1.0 / probe_Ne(npa)

           probe_VXe(npa) = probe_VXe(npa) * inv_N
           probe_VYe(npa) = probe_VYe(npa) * inv_N
           probe_VZe(npa) = probe_VZe(npa) * inv_N

           probe_WXe(npa) = probe_WXe(npa) * inv_N
           probe_WYe(npa) = probe_WYe(npa) * inv_N
           probe_WZe(npa) = probe_WZe(npa) * inv_N

           probe_VXVYe(npa) = 2.0 * probe_VXVYe(npa) * inv_N
           probe_VXVZe(npa) = 2.0 * probe_VXVZe(npa) * inv_N
           probe_VYVZe(npa) = 2.0 * probe_VYVZe(npa) * inv_N

           rvx = probe_VXe(npa)
           rvy = probe_VYe(npa)
           rvz = probe_VZe(npa)

           rvx2 = probe_WXe(npa)
           rvy2 = probe_WYe(npa)
           rvz2 = probe_WZe(npa)

!           rtemp = 2.0 * (rvx * rvx + rvy * rvy + rvz * rvz) - rvx2 - rvy2 - rvz2
           rtemp = rvx * rvx + rvy * rvy + rvz * rvz
           rtemp = rtemp + rtemp - rvx2 - rvy2 - rvz2

           probe_QXe(npa) = probe_QXe(npa) * inv_N + &
                          & (rtemp - rvx2 - rvx2) * rvx - &
                          & probe_VXVYe(npa) * rvy - & 
                          & probe_VXVZe(npa) * rvz

           probe_QYe(npa) = probe_QYe(npa) * inv_N - &
                          & probe_VXVYe(npa) * rvx + &
                          & (rtemp - rvy2 - rvy2) * rvy - &
                          & probe_VYVZe(npa) * rvz

           probe_QZe(npa) = probe_QZe(npa) * inv_N - &
                          & probe_VXVZe(npa) * rvx - &
                          & probe_VYVZe(npa) * rvy + &
                          & (rtemp - rvz2 - rvz2) * rvz

        ELSE
           probe_Ne(npa) = 0.0

           probe_VXe(npa) = 0.0
           probe_VYe(npa) = 0.0
           probe_VZe(npa) = 0.0

           probe_WXe(npa) = 0.0
           probe_WYe(npa) = 0.0
           probe_WZe(npa) = 0.0

           probe_QXe(npa) = 0.0
           probe_QYe(npa) = 0.0
           probe_QZe(npa) = 0.0
        END IF   !###   IF (probe_Ne(npa).GT.1.0e-9) THEN

     END DO   !###   DO npa = 1, N_of_probes

!########## TO DO :: ADJUST DENSITIES at walls #########################

! calculate derived values and convert to dimensional units

     DO npa = 1, N_of_probes

        probe_EX(npa) = REAL(E_scale_Vm) * probe_EX(npa)
        probe_EY(npa) = REAL(E_scale_Vm) * probe_EY(npa)

        probe_F(npa) = REAL(F_scale_V) * probe_F(npa)

        probe_JXe(npa) = -probe_VXe(npa) * probe_Ne(npa) * REAL(current_factor_Am2)
        probe_JYe(npa) = -probe_VYe(npa) * probe_Ne(npa) * REAL(current_factor_Am2)
        probe_JZe(npa) = -probe_VZe(npa) * probe_Ne(npa) * REAL(current_factor_Am2)

        probe_TXe(npa) = MAX(0.0, probe_WXe(npa) - probe_VXe(npa)**2) * REAL(temperature_factor_eV)
        probe_TYe(npa) = MAX(0.0, probe_WYe(npa) - probe_VYe(npa)**2) * REAL(temperature_factor_eV)
        probe_TZe(npa) = MAX(0.0, probe_WZe(npa) - probe_VZe(npa)**2) * REAL(temperature_factor_eV)

        probe_QXe(npa) = probe_QXe(npa) * probe_Ne(npa) * REAL(heat_flow_factor_Wm2)
        probe_QYe(npa) = probe_QYe(npa) * probe_Ne(npa) * REAL(heat_flow_factor_Wm2)
        probe_QZe(npa) = probe_QZe(npa) * probe_Ne(npa) * REAL(heat_flow_factor_Wm2)

        probe_Ne(npa) = probe_Ne(npa) * REAL(N_scale_part_m3)

        probe_VXe(npa) = probe_VXe(npa) * REAL(V_scale_ms)
        probe_VYe(npa) = probe_VYe(npa) * REAL(V_scale_ms)
        probe_VZe(npa) = probe_VZe(npa) * REAL(V_scale_ms)

        probe_WXe(npa) = probe_WXe(npa) * REAL(energy_factor_eV)
        probe_WYe(npa) = probe_WYe(npa) * REAL(energy_factor_eV)
        probe_WZe(npa) = probe_WZe(npa) * REAL(energy_factor_eV)

        IF (T_cntr.EQ.(Save_probes_i_data_T_cntr - N_subcycles/2)) THEN
! full current is saved only at the time when the ion current is available
! full current is saved with the ion data
           probe_JXsum(npa) = probe_JXe(npa)
           probe_JYsum(npa) = probe_JYe(npa)
           probe_JZsum(npa) = probe_JZe(npa)
        END IF

     END DO   !### DO npa = 1, N_of_probes

! ############ 
! TODO
! correct densities and currents in the corners and near walls 
! ###########

!  write to files

     time_ns = REAL(T_cntr * delta_t_s * 1.0d9)

! electric field EX 
     OPEN  (21, FILE = 'dim_EX_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_EX
     CLOSE (21, STATUS = 'KEEP')

! electric field EY 
     OPEN  (21, FILE = 'dim_EY_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_EY
     CLOSE (21, STATUS = 'KEEP')

! electrostatic potential
     OPEN  (21, FILE = 'dim_F_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_F
     CLOSE (21, STATUS = 'KEEP')

! electron density Ne 
     OPEN  (21, FILE = 'dim_Ne_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_Ne
     CLOSE (21, STATUS = 'KEEP')

! electron flow velocity VXe 
     OPEN  (21, FILE = 'dim_VXe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_VXe
     CLOSE (21, STATUS = 'KEEP')

! electron flow velocity VYe 
     OPEN  (21, FILE = 'dim_VYe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_VYe
     CLOSE (21, STATUS = 'KEEP')

! electron flow velocity VZe 
     OPEN  (21, FILE = 'dim_VZe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_VZe
     CLOSE (21, STATUS = 'KEEP')

! electric current density due to electrons JXe 
     OPEN  (21, FILE = 'dim_JXe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_JXe
     CLOSE (21, STATUS = 'KEEP')

! electric current density due to electrons JYe 
     OPEN  (21, FILE = 'dim_JYe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_JYe
     CLOSE (21, STATUS = 'KEEP')

! electric current density due to electrons JZe 
     OPEN  (21, FILE = 'dim_JZe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_JZe
     CLOSE (21, STATUS = 'KEEP')

! average electron energy WXe 
     OPEN  (21, FILE = 'dim_WXe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_WXe
     CLOSE (21, STATUS = 'KEEP')

! average electron energy WYe 
     OPEN  (21, FILE = 'dim_WYe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_WYe
     CLOSE (21, STATUS = 'KEEP')

! average electron energy WZe 
     OPEN  (21, FILE = 'dim_WZe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_WZe
     CLOSE (21, STATUS = 'KEEP')

! electron temperature TXe 
     OPEN  (21, FILE = 'dim_TXe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_TXe
     CLOSE (21, STATUS = 'KEEP')

! electron temperature TYe 
     OPEN  (21, FILE = 'dim_TYe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_TYe
     CLOSE (21, STATUS = 'KEEP')

! electron temperature TZe 
     OPEN  (21, FILE = 'dim_TZe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_TZe
     CLOSE (21, STATUS = 'KEEP')

! electron heat flow QXe 
     OPEN  (21, FILE = 'dim_QXe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_QXe
     CLOSE (21, STATUS = 'KEEP')

! electron heat flow QYe 
     OPEN  (21, FILE = 'dim_QYe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_QYe
     CLOSE (21, STATUS = 'KEEP')

! electron heat flow QZe 
     OPEN  (21, FILE = 'dim_QZe_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_QZe
     CLOSE (21, STATUS = 'KEEP')

     DEALLOCATE(probe_EX, STAT = ALLOC_ERR)
     DEALLOCATE(probe_EY, STAT = ALLOC_ERR)
     DEALLOCATE(probe_F, STAT = ALLOC_ERR)

     DEALLOCATE(probe_Ne, STAT = ALLOC_ERR)

     DEALLOCATE(probe_VXe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_VYe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_VZe, STAT = ALLOC_ERR)

     DEALLOCATE(probe_JXe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_JYe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_JZe, STAT = ALLOC_ERR)

     DEALLOCATE(probe_WXe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_WYe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_WZe, STAT = ALLOC_ERR)

     DEALLOCATE(probe_TXe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_TYe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_TZe, STAT = ALLOC_ERR)

     DEALLOCATE(probe_QXe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_QYe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_QZe, STAT = ALLOC_ERR)

     DEALLOCATE(probe_VXVYe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_VXVZe, STAT = ALLOC_ERR)
     DEALLOCATE(probe_VYVZe, STAT = ALLOC_ERR)

  ELSE

! send data to the zero cluster master
     IF (N_of_probes_cluster.GT.0) THEN

        bufsize = (3 + 13) * N_of_probes_cluster
        ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)

        DO npc = 1, N_of_probes_cluster
           npa = List_of_probes_cluster(npc)
           i = Probe_position(1,npa)
           j = Probe_position(2,npa)
           rbufer(npc)                     = EX(i, j)
           rbufer(npc+N_of_probes_cluster) = EY(i, j)
        END DO

        pos1 = 2 * N_of_probes_cluster+1
        pos2 = 3 * N_of_probes_cluster
        rbufer(pos1:pos2) = probe_F_cluster(1:N_of_probes_cluster)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_Ne_cluster(1:N_of_probes_cluster)   ! +3*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VXe_cluster(1:N_of_probes_cluster)   ! +4*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VYe_cluster(1:N_of_probes_cluster)   ! +5*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VZe_cluster(1:N_of_probes_cluster)   ! +6*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_WXe_cluster(1:N_of_probes_cluster)   ! +7*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_WYe_cluster(1:N_of_probes_cluster)   ! +8*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_WZe_cluster(1:N_of_probes_cluster)   ! +9*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_QXe_cluster(1:N_of_probes_cluster)   ! +10*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_QYe_cluster(1:N_of_probes_cluster)   ! +11*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_QZe_cluster(1:N_of_probes_cluster)   ! +12*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VXVYe_cluster(1:N_of_probes_cluster)   ! +13*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VXVZe_cluster(1:N_of_probes_cluster)   ! +14*recsize

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VYVZe_cluster(1:N_of_probes_cluster)   ! +15*recsize

        CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_REAL, 0, Rank_horizontal, COMM_HORIZONTAL, ierr) 
 
        DEALLOCATE(rbufer)

     END IF   !### IF (N_of_probes_cluster.GT.0) THEN

  END IF   !### IF (Rank_horizontal.EQ.0) THEN

END SUBROUTINE DO_PROBE_DIAGNOSTICS_e_DATA


!------------------------------
!
SUBROUTINE DO_PROBE_DIAGNOSTICS_i_DATA
    
  USE ParallelOperationValues
  USE Diagnostics
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : N_spec, Qs, Ms

  use mpi

  IMPLICIT NONE  

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER bufsize
  REAL, ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR
  INTEGER pos1, pos2

  INTEGER s, k, n, i, j
  INTEGER npa, npb, npc

  REAL, ALLOCATABLE :: probe_Ni(:,:)

  REAL, ALLOCATABLE :: probe_VXi(:,:)
  REAL, ALLOCATABLE :: probe_VYi(:,:)
  REAL, ALLOCATABLE :: probe_VZi(:,:)

  REAL, ALLOCATABLE :: probe_WXi(:,:)
  REAL, ALLOCATABLE :: probe_WYi(:,:)
  REAL, ALLOCATABLE :: probe_WZi(:,:)

  REAL, ALLOCATABLE :: probe_VXVYi(:,:)
  REAL, ALLOCATABLE :: probe_VXVZi(:,:)
  REAL, ALLOCATABLE :: probe_VYVZi(:,:)

  REAL, ALLOCATABLE :: probe_QXi(:,:)
  REAL, ALLOCATABLE :: probe_QYi(:,:)
  REAL, ALLOCATABLE :: probe_QZi(:,:)

  REAL, ALLOCATABLE :: probe_JXi(:,:)
  REAL, ALLOCATABLE :: probe_JYi(:,:)
  REAL, ALLOCATABLE :: probe_JZi(:,:)

  REAL, ALLOCATABLE :: probe_TXi(:,:)
  REAL, ALLOCATABLE :: probe_TYi(:,:)
  REAL, ALLOCATABLE :: probe_TZi(:,:)

  INTEGER pos
  INTEGER recsize

  REAL rvx, rvy, rvz
  REAL rvx2, rvy2, rvz2
  REAL rtemp, inv_N

  REAL time_ns

  CHARACTER(16) Nis_filename   ! dim_Ni_S_vst.dat
                               ! ----*----I----*--
  CHARACTER(17) NNis_filename  ! dim_VXi_S_vst.dat

  INTERFACE
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  END INTERFACE

  IF (cluster_rank_key.NE.0) RETURN

! the rest is carried out by masters only

  IF (Rank_horizontal.EQ.0) THEN
! the zero-rank master process assembles all probe data 

     ALLOCATE(probe_Ni(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)

     ALLOCATE(probe_VXi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_VYi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_VZi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)

     ALLOCATE(probe_WXi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_WYi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_WZi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)

     ALLOCATE(probe_VXVYi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_VXVZi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_VYVZi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)

     ALLOCATE(probe_QXi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_QYi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_QZi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)

     ALLOCATE(probe_JXi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_JYi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_JZi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)

     ALLOCATE(probe_TXi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_TYi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)
     ALLOCATE(probe_TZi(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)

! values of moments N, V, W, V*V*, Q in probes are already calculated but are NOT synchronized in overlapping nodes

     probe_Ni = 0.0

     probe_VXi = 0.0
     probe_VYi = 0.0
     probe_VZi = 0.0

     probe_WXi = 0.0
     probe_WYi = 0.0
     probe_WZi = 0.0

     probe_VXVYi = 0.0
     probe_VXVZi = 0.0
     probe_VYVZi = 0.0

     probe_QXi = 0.0
     probe_QYi = 0.0
     probe_QZi = 0.0

! includes its own probes (if there are any)
     DO npc = 1, N_of_probes_cluster

        npa = List_of_probes_cluster(npc)    ! npa = n-umber of p-robe, a-ll   npc = n-umber of p-robe in c-luster

        i = Probe_position(1, npa)
        j = Probe_position(2, npa)
        DO s = 1, N_spec
           probe_Ni(npa, s) = probe_Ni_cluster(npc, s)

           probe_VXi(npa, s) = probe_VXi_cluster(npc, s)
           probe_VYi(npa, s) = probe_VYi_cluster(npc, s)
           probe_VZi(npa, s) = probe_VZi_cluster(npc, s)

           probe_WXi(npa, s) = probe_WXi_cluster(npc, s)
           probe_WYi(npa, s) = probe_WYi_cluster(npc, s)
           probe_WZi(npa, s) = probe_WZi_cluster(npc, s)

           probe_VXVYi(npa, s) = probe_VXVYi_cluster(npc, s)
           probe_VXVZi(npa, s) = probe_VXVZi_cluster(npc, s)
           probe_VYVZi(npa, s) = probe_VYVZi_cluster(npc, s)

           probe_QXi(npa, s) = probe_QXi_cluster(npc, s)
           probe_QYi(npa, s) = probe_QYi_cluster(npc, s)
           probe_QZi(npa, s) = probe_QZi_cluster(npc, s)
        END DO
     END DO

! assemble probe data from all other clusters
     DO n = 1, N_processes_horizontal-1

        IF (diag_cluster(n)%N_probes.LE.0) CYCLE

        recsize = diag_cluster(n)%N_probes
        bufsize = 13 * N_spec * recsize   ! i::N,VX,VY,VZ,WX,WY,WZ,QX,QY,QZ,VXVY,VXVZ,VYVZ
        ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)

        CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_REAL, n, n, COMM_HORIZONTAL, stattus, ierr)

! since moments in overlapping areas are not synchronized, values in probes belonging to more than one clusters are accumulated
        DO s = 1, N_spec
           pos = (s-1) * 13 * recsize
           DO npc = 1, recsize
              npa = diag_cluster(n)%probe_number(npc)

              probe_Ni(npa, s)  = probe_Ni(npa, s) + rbufer(npc+pos)

!print *, npc, npa, rbufer(npc+pos)

              probe_VXi(npa, s) = probe_VXi(npa, s) + rbufer(npc+pos+  recsize)
              probe_VYi(npa, s) = probe_VYi(npa, s) + rbufer(npc+pos+2*recsize)
              probe_VZi(npa, s) = probe_VZi(npa, s) + rbufer(npc+pos+3*recsize)

              probe_WXi(npa, s) = probe_WXi(npa, s) + rbufer(npc+pos+4*recsize)
              probe_WYi(npa, s) = probe_WYi(npa, s) + rbufer(npc+pos+5*recsize)
              probe_WZi(npa, s) = probe_WZi(npa, s) + rbufer(npc+pos+6*recsize)

              probe_QXi(npa, s) = probe_QXi(npa, s) + rbufer(npc+pos+7*recsize)
              probe_QYi(npa, s) = probe_QYi(npa, s) + rbufer(npc+pos+8*recsize)
              probe_QZi(npa, s) = probe_QZi(npa, s) + rbufer(npc+pos+9*recsize)

              probe_VXVYi(npa, s) = probe_VXVYi(npa, s) + rbufer(npc+pos+10*recsize)
              probe_VXVZi(npa, s) = probe_VXVZi(npa, s) + rbufer(npc+pos+11*recsize)
              probe_VYVZi(npa, s) = probe_VYVZi(npa, s) + rbufer(npc+pos+12*recsize)
           END DO
        END DO
           
        DEALLOCATE(rbufer, STAT = ALLOC_ERR)

     END DO   !### DO n = 1, N_processes_horizontal-1

! finalize moment calculation
     DO s = 1, N_spec
        DO npa = 1, N_of_probes

           IF (probe_Ni(npa, s).GT.1.0e-9) THEN    ! note this is small but not zero

              inv_N = 1.0 / probe_Ni(npa, s)

              probe_VXi(npa, s) = probe_VXi(npa, s) * inv_N
              probe_VYi(npa, s) = probe_VYi(npa, s) * inv_N
              probe_VZi(npa, s) = probe_VZi(npa, s) * inv_N

              probe_WXi(npa, s) = probe_WXi(npa, s) * inv_N
              probe_WYi(npa, s) = probe_WYi(npa, s) * inv_N
              probe_WZi(npa, s) = probe_WZi(npa, s) * inv_N

              probe_VXVYi(npa, s) = 2.0 * probe_VXVYi(npa, s) * inv_N
              probe_VXVZi(npa, s) = 2.0 * probe_VXVZi(npa, s) * inv_N
              probe_VYVZi(npa, s) = 2.0 * probe_VYVZi(npa, s) * inv_N

              rvx = probe_VXi(npa, s)
              rvy = probe_VYi(npa, s)
              rvz = probe_VZi(npa, s)

              rvx2 = probe_WXi(npa, s)
              rvy2 = probe_WYi(npa, s)
              rvz2 = probe_WZi(npa, s)

!           rtemp = 2.0 * (rvx * rvx + rvy * rvy + rvz * rvz) - rvx2 - rvy2 - rvz2
              rtemp = rvx * rvx + rvy * rvy + rvz * rvz
              rtemp = rtemp + rtemp - rvx2 - rvy2 - rvz2

              probe_QXi(npa, s) = probe_QXi(npa, s) * inv_N + &
                             & (rtemp - rvx2 - rvx2) * rvx - &
                             & probe_VXVYi(npa, s) * rvy - & 
                             & probe_VXVZi(npa, s) * rvz

              probe_QYi(npa, s) = probe_QYi(npa, s) * inv_N - &
                             & probe_VXVYi(npa, s) * rvx + &
                             & (rtemp - rvy2 - rvy2) * rvy - &
                             & probe_VYVZi(npa, s) * rvz

              probe_QZi(npa, s) = probe_QZi(npa, s) * inv_N - &
                             & probe_VXVZi(npa, s) * rvx - &
                             & probe_VYVZi(npa, s) * rvy + &
                             & (rtemp - rvz2 - rvz2) * rvz

           ELSE
              probe_Ni(npa, s) = 0.0

              probe_VXi(npa, s) = 0.0
              probe_VYi(npa, s) = 0.0
              probe_VZi(npa, s) = 0.0

              probe_WXi(npa, s) = 0.0
              probe_WYi(npa, s) = 0.0
              probe_WZi(npa, s) = 0.0

              probe_QXi(npa, s) = 0.0
              probe_QYi(npa, s) = 0.0
              probe_QZi(npa, s) = 0.0
           END IF   !###   IF (probe_Ni(npa, s).GT.1.0e-9) THEN

        END DO   !###   DO npa = 1, N_of_probes
     END DO   !###   DO s = 1, N_spec

!### TO DO :: adjust densities near walls #####################

! calculate derived values and convert to dimensional units

     DO npa = 1, N_of_probes

        DO s = 1, N_spec

           probe_JXi(npa, s) = probe_VXi(npa, s) * probe_Ni(npa, s) * REAL(Qs(s) * current_factor_Am2)
           probe_JYi(npa, s) = probe_VYi(npa, s) * probe_Ni(npa, s) * REAL(Qs(s) * current_factor_Am2)
           probe_JZi(npa, s) = probe_VZi(npa, s) * probe_Ni(npa, s) * REAL(Qs(s) * current_factor_Am2)
          
           probe_TXi(npa, s) = MAX(0.0, probe_WXi(npa, s) - probe_VXi(npa, s)**2) * REAL(Ms(s) * temperature_factor_eV)
           probe_TYi(npa, s) = MAX(0.0, probe_WYi(npa, s) - probe_VYi(npa, s)**2) * REAL(Ms(s) * temperature_factor_eV)
           probe_TZi(npa, s) = MAX(0.0, probe_WZi(npa, s) - probe_VZi(npa, s)**2) * REAL(Ms(s) * temperature_factor_eV)

           probe_QXi(npa, s) = probe_Ni(npa, s) * probe_QXi(npa, s) * REAL(Ms(s) * heat_flow_factor_Wm2)
           probe_QYi(npa, s) = probe_Ni(npa, s) * probe_QYi(npa, s) * REAL(Ms(s) * heat_flow_factor_Wm2)
           probe_QZi(npa, s) = probe_Ni(npa, s) * probe_QZi(npa, s) * REAL(Ms(s) * heat_flow_factor_Wm2)

           probe_Ni(npa, s) = probe_Ni(npa, s) * REAL(N_scale_part_m3)

           probe_VXi(npa, s) = probe_VXi(npa, s) * REAL(V_scale_ms)
           probe_VYi(npa, s) = probe_VYi(npa, s) * REAL(V_scale_ms)
           probe_VZi(npa, s) = probe_VZi(npa, s) * REAL(V_scale_ms)

           probe_WXi(npa, s) = probe_WXi(npa, s) * REAL(Ms(s) * energy_factor_eV)
           probe_WYi(npa, s) = probe_WYi(npa, s) * REAL(Ms(s) * energy_factor_eV)
           probe_WZi(npa, s) = probe_WZi(npa, s) * REAL(Ms(s) * energy_factor_eV)

           probe_JXsum(npa) = probe_JXsum(npa) + probe_JXi(npa, s)
           probe_JYsum(npa) = probe_JYsum(npa) + probe_JYi(npa, s)
           probe_JZsum(npa) = probe_JZsum(npa) + probe_JZi(npa, s)
        END DO

     END DO   !### DO npa = 1, N_of_probes

!  write to files

     time_ns = REAL((T_cntr - N_subcycles/2) * delta_t_s * 1.0d9)

! full electric current density JXsum 
     OPEN  (21, FILE = 'dim_JXsum_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_JXsum
     CLOSE (21, STATUS = 'KEEP')

! full electric current density JYsum 
     OPEN  (21, FILE = 'dim_JYsum_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_JYsum
     CLOSE (21, STATUS = 'KEEP')

! full electric current density JZsum 
     OPEN  (21, FILE = 'dim_JZsum_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_JZsum
     CLOSE (21, STATUS = 'KEEP')

! ion species data
     DO s = 1, N_spec

! density
        Nis_filename = 'dim_Ni_S_vst.dat'
        Nis_filename(8:8) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = Nis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_Ni(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! flow velocity VX
        NNis_filename = 'dim_VXi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_VXi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! flow velocity VY
        NNis_filename = 'dim_VYi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_VYi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! flow velocity VZ
        NNis_filename = 'dim_VZi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_VZi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! electric curent density due to ions JX
        NNis_filename = 'dim_JXi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_JXi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! electric curent density due to ions JY
        NNis_filename = 'dim_JYi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_JYi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! electric curent density due to ions JZ
        NNis_filename = 'dim_JZi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_JZi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! average ion energy WX
        NNis_filename = 'dim_WXi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_WXi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! average ion energy WY
        NNis_filename = 'dim_WYi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_WYi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! average ion energy WZ
        NNis_filename = 'dim_WZi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_WZi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! ion temperature TX
        NNis_filename = 'dim_TXi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_TXi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! ion temperature TY
        NNis_filename = 'dim_TYi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_TYi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! ion temperature TY
        NNis_filename = 'dim_TZi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_TZi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! ion heat flow QX
        NNis_filename = 'dim_QXi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_QXi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! ion heat flow QY
        NNis_filename = 'dim_QYi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_QYi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

! ion heat flow QZ
        NNis_filename = 'dim_QZi_S_vst.dat'
        NNis_filename(9:9) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = NNis_filename, POSITION = 'APPEND')
        WRITE (21, '(1x,f15.6,100(1x,e14.7))') time_ns, probe_QZi(1:N_of_probes, s)
        CLOSE (21, STATUS = 'KEEP')

     END DO

     DEALLOCATE(probe_Ni, STAT = ALLOC_ERR)

     DEALLOCATE(probe_VXi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_VYi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_VZi, STAT = ALLOC_ERR)

     DEALLOCATE(probe_VXVYi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_VXVZi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_VYVZi, STAT = ALLOC_ERR)

     DEALLOCATE(probe_JXi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_JYi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_JZi, STAT = ALLOC_ERR)

     DEALLOCATE(probe_WXi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_WYi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_WZi, STAT = ALLOC_ERR)

     DEALLOCATE(probe_TXi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_TYi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_TZi, STAT = ALLOC_ERR)

     DEALLOCATE(probe_QXi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_QYi, STAT = ALLOC_ERR)
     DEALLOCATE(probe_QZi, STAT = ALLOC_ERR)

  ELSE   !### IF (Rank_horizontal.EQ.0) THEN

! send data to the zero cluster master
     IF (N_of_probes_cluster.LE.0) RETURN

     bufsize = 13 * N_spec * N_of_probes_cluster
     ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)

     pos1 = 1 - N_of_probes_cluster
     pos2 = 0

     DO s = 1, N_spec
        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_Ni_cluster(1:N_of_probes_cluster, s)

!print *, rbufer(pos1:pos2)
        
        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VXi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VYi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VZi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_WXi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_WYi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_WZi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_QXi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_QYi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_QZi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VXVYi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VXVZi_cluster(1:N_of_probes_cluster, s)

        pos1 = pos1 + N_of_probes_cluster
        pos2 = pos2 + N_of_probes_cluster
        rbufer(pos1:pos2) = probe_VYVZi_cluster(1:N_of_probes_cluster, s)
     END DO
     
     CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_REAL, 0, Rank_horizontal, COMM_HORIZONTAL, ierr) 
 
     DEALLOCATE(rbufer)

  END IF   !### IF (Rank_horizontal.EQ.0) THEN

END SUBROUTINE DO_PROBE_DIAGNOSTICS_i_DATA

