!------------------------------
! 
SUBROUTINE INITIATE_PROBE_DIAGNOSTICS
    
  USE ParallelOperationValues
  USE Diagnostics
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries
  USE IonParticles, ONLY : N_spec
  USE Checkpoints, ONLY : use_checkpoint, Save_probes_data_T_cntr_check, N_of_saved_records_check

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  LOGICAL exists

  CHARACTER (1) buf

  INTEGER N_of_given_probes
  INTEGER n, s, k, i
  INTEGER IOS
  INTEGER temp_pos(1:2,1:Max_N_of_probes)         ! 1 = i_x, 2 = j_y 

  INTEGER count_skip

  INTEGER npa, npc, npb  ! used to denote index of probe global (a=all), in cluster (c), and in block (b)

  REAL r_dummy

  CHARACTER(34) clusterprobes_filename    ! _probelocs_cluster_NNNNN_NNNNN.dat
                                          ! ----x----I----x----I----x----I----

  CHARACTER(16) Nis_filename

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

  IF (exists) THEN

     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i5," : init_probes.dat is found. Reading the data file...")', Rank_of_process
     END IF

     OPEN (9, FILE = 'init_probes.dat')

     READ (9, '(A1)') buf !----dddddd--- Step for saving (timesteps, >=1)
     READ (9, '(4x,i6)') Save_probes_data_step
     READ (9, '(A1)') buf !---ddddddd--- Start saving data at (timesteps, >=0)
     READ (9, '(3x,i7)') Save_probes_data_T_cntr_rff                                          !######## was Save_probes_data_T_cntr
     READ (9, '(A1)') buf !------dddd--- Skip periods of writing between text outputs (>=0)
     READ (9, '(6x,i4)') TextOut_skip
     READ (9, '(A1)') buf !------dddd--- Number of probes (off if <=0)
     READ (9, '(6x,i4)') N_of_given_probes
     READ (9, '(A1)') buf !--ddddd--ddddd-- Probe coordinates (x/y, node numbers)

     IF (N_of_given_probes.GT.Max_N_of_probes) THEN
        IF (Rank_of_process.EQ.0) THEN
           PRINT '(2x,"Process ",i5," : requested number of probes ",i4," was reduced to maximal allowed value of ",i4)', Rank_of_process, Max_N_of_probes
        END IF
        N_of_given_probes = Max_N_of_probes
     END IF
     
     DO n = 1, N_of_given_probes
        READ (9, '(2x,i5,2x,i5)', IOSTAT=IOS) temp_pos(1,n), temp_pos(2,n)
        IF (IOS.NE.0) THEN
           PRINT '(2x,"Process ",i5," :: INITIATE_PROBE_DIAGNOSTICS : ERROR-1 : while reading file init_probes.dat : wrong coordinates of probe ",i4," program terminated.")', Rank_of_process, n
           STOP
        END IF
     END DO

     CLOSE (9, STATUS = 'KEEP')

  ELSE
     
     PRINT '(2x,"Process ",i5," :: INITIATE_PROBE_DIAGNOSTICS : ERROR-2 : init_probes.dat not found, program terminated.")', Rank_of_process
     STOP

  END IF

  text_output_counter = 0
  N_of_saved_records = 0
  Save_probes_data_T_cntr = Save_probes_data_T_cntr_rff   

! overwrite if checkpoint is used to initialize the system
  IF (use_checkpoint.EQ.1) THEN
     Save_probes_data_T_cntr = Save_probes_data_T_cntr_check
     N_of_saved_records = N_of_saved_records_check
  END IF

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

     IF (N_of_probes_cluster.GT.0) THEN

        ALLOCATE(List_of_probes_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)

        ALLOCATE(probe_Ne_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        ALLOCATE(probe_F_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        ALLOCATE(probe_Ni_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)

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
     

  ELSE
! receive number of probes in cluster from the master
     CALL MPI_BCAST(N_of_probes_cluster, 1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

     IF (N_of_probes_cluster.GT.0) THEN
        ALLOCATE(probe_Ni_cluster(1:N_of_probes_cluster,1:N_spec), STAT = ALLOC_ERR)
        ALLOCATE(List_of_probes_cluster(1:N_of_probes_cluster), STAT = ALLOC_ERR)
        CALL MPI_BCAST(List_of_probes_cluster, N_of_probes_cluster, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
     END IF

  END IF

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
           CALL MPI_SEND(ibufer(1), 1, MPI_INTEGER, field_calculator(k)%rank, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
        END DO

        IF (N_of_probes_cluster.GT.0) THEN
! master sends to the associated field calculators a list of its probes
           DO k = 2, cluster_N_blocks
              CALL MPI_SEND(List_of_probes_cluster(1:N_of_probes_cluster), N_of_probes_cluster, MPI_INTEGER, field_calculator(k)%rank, Rank_of_process+100000, MPI_COMM_WORLD, request, ierr) 
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

     ELSE

! field calculators receive number of probes from the master
        CALL MPI_RECV(ibufer(1), 1, MPI_INTEGER, field_master, field_master, MPI_COMM_WORLD, stattus, ierr)
        N_of_probes_cluster = ibufer(1)  ! here N_of_probes_cluster is a temporary value which will be replaced
                                         ! once the process becomes associated with a cluster as a particle calculator

        N_of_probes_block = 0   ! default assumption

        IF (N_of_probes_cluster.GT.0) THEN
! field calculators receive list of probes from the master
           CALL MPI_RECV(ibufer(1:N_of_probes_cluster), N_of_probes_cluster, MPI_INTEGER, field_master, field_master+100000, MPI_COMM_WORLD, stattus, ierr)
! identify how many probes from the list are inside the field calculator domain
           DO npc = 1, N_of_probes_cluster
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
              DO npc = 1, N_of_probes_cluster
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
           CALL MPI_SEND(ibufer(1), 1, MPI_INTEGER, field_master, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
! report the vector of probe order numbers in the cluster's probe list 
           IF (N_of_probes_block.GT.0) THEN
              ibufer(1:N_of_probes_block) = Probe_params_block_list(3,1:N_of_probes_block)
              CALL MPI_SEND(ibufer(1:N_of_probes_block), N_of_probes_block, MPI_INTEGER, field_master, Rank_of_process+100000, MPI_COMM_WORLD, request, ierr) 
! allocate array to save the potential
              ALLOCATE(probe_F_block(1:N_of_probes_block), STAT = ALLOC_ERR)
           END IF
        END IF

     END IF

     DEALLOCATE(ibufer, STAT = ALLOC_ERR)

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  END IF

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

! electric field EX 
           INQUIRE (FILE = 'dim_EX_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_EX_vst.dat', STATUS = 'OLD')          
              DO i = 1, N_of_saved_records
                 READ (21, '(2x,f14.6,100(1x,e14.7))') r_dummy
              END DO
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electric field EY
           INQUIRE (FILE = 'dim_EY_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_EY_vst.dat', STATUS = 'OLD')          
              DO i = 1, N_of_saved_records
                 READ (21, '(2x,f14.6,100(1x,e14.7))') r_dummy
              END DO
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electrostatic potential F 
           INQUIRE (FILE = 'dim_F_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_F_vst.dat', STATUS = 'OLD')          
              DO i = 1, N_of_saved_records
                 READ (21, '(2x,f14.6,100(1x,e14.7))') r_dummy
              END DO
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! electron density Ne 
           INQUIRE (FILE = 'dim_Ne_vst.dat', EXIST = exists)
           IF (exists) THEN                                                       
              OPEN (21, FILE = 'dim_Ne_vst.dat', STATUS = 'OLD')          
              DO i = 1, N_of_saved_records
                 READ (21, '(2x,f14.6,100(1x,e14.7))') r_dummy
              END DO
              ENDFILE 21       
              CLOSE (21, STATUS = 'KEEP')        
           END IF

! ion density
           DO s = 1, N_spec
              Nis_filename = 'dim_Ni_S_vst.dat'
              Nis_filename(8:8) = convert_int_to_txt_string(s, 1)
              INQUIRE (FILE = Nis_filename, EXIST = exists)
              IF (exists) THEN                                                       
                 OPEN (21, FILE = Nis_filename, STATUS = 'OLD')          
                 DO i = 1, N_of_saved_records
                    READ (21, '(2x,f14.6,100(1x,e14.7))') r_dummy
                 END DO
                 ENDFILE 21       
                 CLOSE (21, STATUS = 'KEEP')        
              END IF
           END DO
 
        ELSE
! fresh start, empty files, clean up whatever garbage there might be

! electric field EX 
           OPEN  (21, FILE = 'dim_EX_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electric field EY
           OPEN  (21, FILE = 'dim_EY_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electrostatic potential F 
           OPEN  (21, FILE = 'dim_F_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! electron density Ne 
           OPEN  (21, FILE = 'dim_Ne_vst.dat', STATUS = 'REPLACE')          
           CLOSE (21, STATUS = 'KEEP')

! ion density
           DO s = 1, N_spec
              Nis_filename = 'dim_Ni_S_vst.dat'
              Nis_filename(8:8) = convert_int_to_txt_string(s, 1)
              OPEN  (21, FILE = Nis_filename, STATUS = 'REPLACE')          
              CLOSE (21, STATUS = 'KEEP')
           END DO
        END IF

     ELSE

        ibufer(1) = N_of_probes_cluster
! send number of probes to the master with zero rank     
        CALL MPI_SEND(ibufer(1:1), 1, MPI_INTEGER, 0, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        
! send global indices of probes to the master with zero rank
        IF (N_of_probes_cluster.GT.0) THEN
           CALL MPI_SEND(List_of_probes_cluster, N_of_probes_cluster, MPI_INTEGER, 0, Rank_horizontal+100000, COMM_HORIZONTAL, request, ierr) 
        END IF

     END IF

     DEALLOCATE(ibufer, STAT = ALLOC_ERR)

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!if (Rank_of_process.eq.0 )print *, "Save_probes_data_T_cntr is ", Save_probes_data_T_cntr

END SUBROUTINE INITIATE_PROBE_DIAGNOSTICS

!------------------------------
!
SUBROUTINE DO_PROBE_DIAGNOSTICS(n_sub)
    
  USE ParallelOperationValues
  USE Diagnostics
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : N_spec

  IMPLICIT NONE  

  INCLUDE 'mpif.h'

  INTEGER n_sub

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  LOGICAL, SAVE :: probe_Ni_cluster_was_updated = .FALSE.

  INTEGER bufsize
  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR
  INTEGER pos1, pos2

  INTEGER s, k, n, i, j
  INTEGER npa, npb, npc

  REAL(8), ALLOCATABLE :: probe_EX(:)
  REAL(8), ALLOCATABLE :: probe_EY(:)
  REAL(8), ALLOCATABLE :: probe_F(:)
  REAL(8), ALLOCATABLE :: probe_Ne(:)
  REAL(8), ALLOCATABLE :: probe_Ni(:,:)

  INTEGER recsize

  REAL(8) time_ns

  CHARACTER(16) Nis_filename

  INTERFACE
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  END INTERFACE

! check whether the GATHER_ION_CHARGE_DENSITY subroutine was called
  IF (n_sub.EQ.0) THEN
     probe_Ni_cluster_was_updated = .TRUE.
  END IF

! exit if the current time layer is not assigned for the diagnostic output 
  IF (T_cntr.NE.Save_probes_data_T_cntr) RETURN

! modify counter to skip extra diagnostics periods between two text outputs
  text_output_counter = text_output_counter + 1

! this counter is saved in checkpoints to be able to trim the time dependence data files
  N_of_saved_records = N_of_saved_records + 1

! perform text output
  IF (text_output_counter.GT.TextOut_skip) THEN
!     CALL DoTextOutput
     text_output_counter = 0
  END IF

! modify the diagnostics time for future use
  Save_probes_data_T_cntr = Save_probes_data_T_cntr + Save_probes_data_step 

  IF (N_of_probes.LE.0) RETURN

! collect ion densities from particle calculators -------------------------------------

  IF (probe_Ni_cluster_was_updated) THEN
! note, if the probe data are saved with interval less than the number of electron subcycles
! then it is possible that probe_Ni_cluster will not be flushed/recalculated 
! and the use of MPI_REDUCE will lead to probe densities accumulating 
! this is why we have to use the "switch" variable probe_Ni_cluster_was_updated
     bufsize = N_spec * N_of_probes_cluster
     IF (bufsize.GT.0) THEN
        ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)
        pos2 = 0
        DO s = 1, N_spec
           pos1 = pos2+1
           pos2 = pos2 + N_of_probes_cluster
           rbufer(pos1:pos2) = probe_Ni_cluster(pos1:pos2, s)
        END DO
        CALL MPI_REDUCE(rbufer, probe_Ni_cluster, bufsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_CLUSTER, ierr)
        DEALLOCATE(rbufer, STAT = ALLOC_ERR)
     END IF
     probe_Ni_cluster_was_updated = .FALSE.
  END IF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

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
              CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)
              DO npb = 1, diag_block(k)%N_probes
                 npc = diag_block(k)%probe_number(npb)
                 probe_F_cluster(npc) = rbufer(npb)
              END DO
              DEALLOCATE(rbufer, STAT = ALLOC_ERR)
           END DO
        END IF
     ELSE
        IF (N_of_probes_block.GT.0) THEN
           CALL MPI_SEND(probe_F_block, N_of_probes_block, MPI_DOUBLE_PRECISION, field_master, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
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
     ALLOCATE(probe_Ni(1:N_of_probes,1:N_spec), STAT = ALLOC_ERR)

     probe_Ne = 0.0_8   ! set zeros because these values will be accumulated (in case of probes at boundaries between clusters)
     probe_Ni = 0.0_8   !

! includes its own probes (if there are any)
     DO npc = 1, N_of_probes_cluster
        npa = List_of_probes_cluster(npc)    ! npa = n-umber of p-robe, a-ll   npc = n-umber of p-robe in c-luster
        i = Probe_position(1, npa)
        j = Probe_position(2, npa)
        probe_EX(npa) = EX(i, j)
        probe_EY(npa) = EY(i, j)
        probe_F(npa) = probe_F_cluster(npc)
        probe_Ne(npa) = probe_Ne(npa) + probe_Ne_cluster(npc)
        DO s = 1, N_spec
           probe_Ni(npa, s) = probe_Ni(npa, s) + probe_Ni_cluster(npc, s)
        END DO
     END DO

! assemble probe data from all other clusters
     DO n = 1, N_processes_horizontal-1

        IF (diag_cluster(n)%N_probes.LE.0) CYCLE

        recsize = diag_cluster(n)%N_probes
        bufsize = (4 + N_spec) * recsize
        ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)

        CALL MPI_RECV(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, n, n, COMM_HORIZONTAL, stattus, ierr)

        DO npc = 1, recsize
           npa = diag_cluster(n)%probe_number(npc)
           probe_EX(npa) = rbufer(npc)                    ! if a probe is at a boundary between clusters, it is included in both clusters
           probe_EY(npa) = rbufer(npc+recsize)            ! for such probes, data like EX, EY, and F are simply overwritten
           probe_F(npa)  = rbufer(npc+2*recsize)          !
           probe_Ne(npa) = probe_Ne(npa) + rbufer(npc+3*recsize)                 !
           DO s = 1, N_spec                                                      !  while Ne and Ni are accumulated
              probe_Ni(npa, s) = probe_Ni(npa, s) + rbufer(npc+(3+s)*recsize)    ! 
           END DO                                                                !
        END DO
           
        DEALLOCATE(rbufer, STAT = ALLOC_ERR)

     END DO

! ############ 
! correct densities in the corners and near walls 
! (the simplest version which works only for a box or 
! for a system periodic in one direction with flat walls in the other direction) 
! note that correct processing of plasma density in probes at periodic boundaries
! is yet to be implemented 
! ###########

     IF (periodicity_flag.EQ.PERIODICITY_NONE) THEN

        DO npa = 1, N_of_probes
           i = Probe_position(1, npa)
           j = Probe_position(2, npa)
           IF ((i.EQ.0).OR.(i.EQ.global_maximal_i)) THEN
              probe_Ne(npa) = 2.0_8 * probe_Ne(npa)
              DO s = 1, N_spec
                 probe_Ni(npa, s) = 2.0_8 * probe_Ni(npa, s)
              END DO
           END IF
           IF ((j.EQ.0).OR.(j.EQ.global_maximal_j)) THEN
              probe_Ne(npa) = 2.0_8 * probe_Ne(npa)
              DO s = 1, N_spec
                 probe_Ni(npa, s) = 2.0_8 * probe_Ni(npa, s)
              END DO
           END IF
        END DO

     ELSE IF (periodicity_flag.EQ.PERIODICITY_X) THEN

        DO npa = 1, N_of_probes
           i = Probe_position(1, npa)
           j = Probe_position(2, npa)
           IF ((j.EQ.0).OR.(j.EQ.global_maximal_j)) THEN
              probe_Ne(npa) = 2.0_8 * probe_Ne(npa)
              DO s = 1, N_spec
                 probe_Ni(npa, s) = 2.0_8 * probe_Ni(npa, s)
              END DO
           END IF
        END DO

     END IF

!  write to files

     time_ns = T_cntr * delta_t_s * 1.0d9

! electric field EX 
     OPEN  (21, FILE = 'dim_EX_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(2x,f14.6,100(1x,e14.7))') time_ns, (REAL(probe_EX(n) * E_scale_Vm), n=1, N_of_probes, 1)
     CLOSE (21, STATUS = 'KEEP')

! electric field EY 
     OPEN  (21, FILE = 'dim_EY_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(2x,f14.6,100(1x,e14.7))') time_ns, (REAL(probe_EY(n) * E_scale_Vm), n=1, N_of_probes, 1)
     CLOSE (21, STATUS = 'KEEP')

! electrostatic potential
     OPEN  (21, FILE = 'dim_F_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(2x,f14.6,100(1x,e14.7))') time_ns, (REAL(probe_F(n) * F_scale_V), n=1, N_of_probes, 1)
     CLOSE (21, STATUS = 'KEEP')

! electron density Ne 
     OPEN  (21, FILE = 'dim_Ne_vst.dat', POSITION = 'APPEND')
     WRITE (21, '(2x,f14.6,100(1x,e14.7))') time_ns, (REAL(probe_Ne(n) * N_scale_part_m3), n=1, N_of_probes, 1)  !####### don't forget to multiply by converting factor 1/n_part_per_cell
     CLOSE (21, STATUS = 'KEEP')

! ion density
     DO s = 1, N_spec
        Nis_filename = 'dim_Ni_S_vst.dat'
        Nis_filename(8:8) = convert_int_to_txt_string(s, 1)
        OPEN  (21, FILE = Nis_filename, POSITION = 'APPEND')
        WRITE (21, '(2x,f14.6,100(1x,e14.7))') time_ns, (REAL(probe_Ni(n, s) * N_scale_part_m3), n=1, N_of_probes, 1)  !####### don't forget to multiply by converting factor 1/n_part_per_cell
        CLOSE (21, STATUS = 'KEEP')
     END DO

     DEALLOCATE(probe_EX, STAT = ALLOC_ERR)
     DEALLOCATE(probe_EY, STAT = ALLOC_ERR)
     DEALLOCATE(probe_F, STAT = ALLOC_ERR)
     DEALLOCATE(probe_Ne, STAT = ALLOC_ERR)
     DEALLOCATE(probe_Ni, STAT = ALLOC_ERR)

  ELSE

! send data to the zero cluster master
     IF (N_of_probes_cluster.GT.0) THEN

        recsize = N_of_probes_cluster
        bufsize = (4 + N_spec) * recsize
        ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)

        DO npc = 1, recsize
           npa = List_of_probes_cluster(npc)
           i = Probe_position(1,npa)
           j = Probe_position(2,npa)
           rbufer(npc) = EX(i, j)
           rbufer(npc+recsize) = EY(i, j)
           rbufer(npc+2*recsize) = probe_F_cluster(npc)
           rbufer(npc+3*recsize) = probe_Ne_cluster(npc)
           DO s = 1, N_spec
              rbufer(npc+(3+s)*recsize) = probe_Ni_cluster(npc, s)
           END DO
        END DO

        CALL MPI_SEND(rbufer(1:bufsize), bufsize, MPI_DOUBLE_PRECISION, 0, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
 
        DEALLOCATE(rbufer)

     END IF

  END IF

END SUBROUTINE DO_PROBE_DIAGNOSTICS
