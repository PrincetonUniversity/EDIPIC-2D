!-----------------------------------------
!
SUBROUTINE GLOBAL_LOAD_BALANCE

  USE ParallelOperationValues
  USE LoadBalancing
  USE ElectronParticles
  USE IonParticles
  USE CurrentProblemValues, ONLY : N_subcycles

  USE ClusterAndItsBoundaries
  use mpi
 
  IMPLICIT NONE

! these variables are used to print ranks of processes in cluster communicator
  integer clustergroup, worldgroup
  integer, allocatable :: clusterranks(:), worldranks(:)
  integer clustergroup_size

  INTEGER errcode,ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER eff_N_particles
  INTEGER eff_N_particles_cluster

  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  INTEGER N_of_all_free_processes

  INTEGER, ALLOCATABLE :: new_particle_masters(:)

  INTEGER i, m, n, k, s

  INTEGER, ALLOCATABLE :: free_process_new_master(:)

  INTEGER N_processes_to_free_in_cluster, delta_N
  INTEGER pos_begin, pos_end

  INTEGER total_part_number_in_bufer
  INTEGER       N_part_to_receive_from_proc(0:N_spec, 1:N_processes_cluster)
  INTEGER total_N_part_to_receive_from_proc(1:N_processes_cluster)

  INTEGER total_new_part_number(0:N_spec)
  INTEGER new_size

  REAL(8), ALLOCATABLE :: rbufer(:)

!  INTEGER ibufer(0:N_spec)
  INTEGER sum_N_part_to_send

! collect number of particles from all processes in a cluster
  eff_N_particles = 0
  DO s = 1, N_spec
     eff_N_particles = eff_N_particles + N_ions(s)
  END DO
  eff_N_particles = eff_N_particles / N_subcycles + N_electrons

  CALL MPI_REDUCE(eff_N_particles, eff_N_particles_cluster, 1, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)

  IF (cluster_rank_key.EQ.0) THEN
     ALLOCATE(ibufer(0:N_processes_horizontal-1), STAT=ALLOC_ERR)
! collect total numbers of particles [i.e. total loads] from all clusters
     CALL MPI_GATHER(eff_N_particles_cluster, 1, MPI_INTEGER, ibufer, 1, MPI_INTEGER, 0, COMM_HORIZONTAL, ierr)
     IF (Rank_horizontal.EQ.0) THEN
        cluster%N_particles = ibufer
        CALL CALCULATE_N_PROCESSES_BALANCED(N_of_all_free_processes)
     END IF
     DEALLOCATE(ibufer, STAT=ALLOC_ERR)
! let all masters know whether to proceed with the balancing (cancel if N_of_all_free_processes is zero)
     CALL MPI_BCAST(N_of_all_free_processes, 1, MPI_INTEGER, 0, COMM_HORIZONTAL, ierr)
  END IF

! let all cluster members know whether to proceed with balancing (cancel if N_of_all_free_processes is zero)
  CALL MPI_BCAST(N_of_all_free_processes, 1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

  IF (N_of_all_free_processes.EQ.0) RETURN   !########### cancel balancing ###########

  ALLOCATE(new_particle_masters(1:N_processes_cluster-1), STAT=ALLOC_ERR)

  IF (cluster_rank_key.EQ.0) THEN

! prepare THE MASTERS MESSAGE with future values of particle_master for all members of the cluster (except the master)
     DO i = 1, N_processes_cluster-1
        new_particle_masters(i) = particle_master
     END DO
     IF (Rank_horizontal.EQ.0) THEN
        ALLOCATE(free_process_new_master(1:N_of_all_free_processes), STAT=ALLOC_ERR)
! save ranks of masters for processes to be assigned
        i = 0
        DO m = 0, N_processes_horizontal-1
           DO n = 1, cluster(m)%N_processes_balanced-cluster(m)%N_processes
              i = i + 1
              free_process_new_master(i) = cluster(m)%particle_master
           END DO
        END DO
! send messages to other masters about the changes in the number of processes
        DO m = 1, N_processes_horizontal-1
           CALL MPI_SEND(cluster(m)%N_processes-cluster(m)%N_processes_balanced, 1, MPI_INTEGER, m, 0, COMM_HORIZONTAL, ierr)            
        END DO
! prepare messages to masters of clusters that will be releasing processes containing values of particle_master for the new cluster of these processes
        i = 0
! message to himself, if necessary [will not be transmitted]
        N_processes_to_free_in_cluster = cluster(0)%N_processes - cluster(0)%N_processes_balanced
        pos_end = 0
        IF (N_processes_to_free_in_cluster.GT.0) THEN
           pos_begin = pos_end + 1
           pos_end = pos_begin + N_processes_to_free_in_cluster - 1
! cluster will lose N_processes_to_free_in_cluster processes, the message with new_particle_masters gets new particle_master values for processes to be released
           new_particle_masters(N_processes_cluster-N_processes_to_free_in_cluster:N_processes_cluster-1) = free_process_new_master(pos_begin:pos_end)
        END IF
! messages to other masters
        DO m = 1, N_processes_horizontal-1
           delta_N = cluster(m)%N_processes-cluster(m)%N_processes_balanced  ! use separate variable delta_N because the process already defined N_processes_to_free_in_cluster for itself
           IF (delta_N.GT.0) THEN
              pos_begin = pos_end + 1
              pos_end = pos_begin + delta_N - 1
              CALL MPI_SEND(free_process_new_master(pos_begin:pos_end), delta_N, MPI_INTEGER, m, SHIFT1, COMM_HORIZONTAL, ierr)
           END IF
        END DO
! update the cluster information
        DO m = 0, N_processes_horizontal-1
           cluster(m)%N_processes = cluster(m)%N_processes_balanced
        END DO
        DEALLOCATE(free_process_new_master, STAT=ALLOC_ERR)
     ELSE
! receive message about the changes in the number of processes
        CALL MPI_RECV(N_processes_to_free_in_cluster, 1, MPI_INTEGER, 0, 0, COMM_HORIZONTAL, stattus, ierr)
        IF (N_processes_to_free_in_cluster.GT.0) THEN
! cluster will lose N_processes_to_free_in_cluster processes, the message with new_particle_masters gets new particle_master values for processes to be released
           CALL MPI_RECV(new_particle_masters(N_processes_cluster-N_processes_to_free_in_cluster:N_processes_cluster-1), N_processes_to_free_in_cluster, MPI_INTEGER, 0, SHIFT1, COMM_HORIZONTAL, stattus, ierr)
        END IF
     END IF      ! IF (Rank_horizontal.EQ.0) THEN/ELSE

! send THE MASTERS MESSAGE to all cluster members
     CALL MPI_BCAST(new_particle_masters, N_processes_cluster-1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

! if the cluster loses processes, get particles from the processes that will be freed
     IF (N_processes_to_free_in_cluster.GT.0) THEN
! receive numbers of particles 
        total_part_number_in_bufer = 0
        DO n = N_processes_cluster-N_processes_to_free_in_cluster, N_processes_cluster-1
           CALL MPI_RECV(N_part_to_receive_from_proc(0:N_spec,n), N_spec+1, MPI_INTEGER, n, 0, COMM_CLUSTER, stattus, ierr)
           total_N_part_to_receive_from_proc(n) = N_part_to_receive_from_proc(0,n)
           DO s = 1, N_spec
              total_N_part_to_receive_from_proc(n) = total_N_part_to_receive_from_proc(n) + N_part_to_receive_from_proc(s,n)
           END DO
           total_part_number_in_bufer = total_part_number_in_bufer + total_N_part_to_receive_from_proc(n)
        END DO
        ALLOCATE(rbufer(1:total_part_number_in_bufer*6), STAT=ALLOC_ERR)
        pos_end = 0
! receive particles
        DO n = N_processes_cluster-N_processes_to_free_in_cluster, N_processes_cluster-1
           IF (total_N_part_to_receive_from_proc(n).GT.0) THEN
              pos_begin = pos_end + 1
              pos_end = pos_begin + total_N_part_to_receive_from_proc(n)*6 - 1
              CALL MPI_RECV(rbufer(pos_begin:pos_end), total_N_part_to_receive_from_proc(n)*6, MPI_DOUBLE_PRECISION, n, n, COMM_CLUSTER, stattus, ierr)
           END IF
        END DO
! add the newly received particles 
! resize arrays
        total_new_part_number(0:N_spec) = 0
        DO s = 0, N_spec
           DO n = N_processes_cluster-N_processes_to_free_in_cluster, N_processes_cluster-1
              total_new_part_number(s) = total_new_part_number(s) + N_part_to_receive_from_proc(s,n)
           END DO
        END DO
! electrons 
        IF ((total_new_part_number(0)+N_electrons).GT.max_N_electrons) THEN
           new_size = INT(MAX(1.05*REAL(N_electrons+total_new_part_number(0)), REAL(N_electrons+total_new_part_number(0)+50)))
           CALL RESIZE_ELECTRON_ARRAY(new_size)   ! sets max_N_electrons = new_size
        END IF
! ions
        DO s = 1, N_spec
           IF ((total_new_part_number(s)+N_ions(s)).GT.max_N_ions(s)) THEN
              new_size = INT(MAX(1.05*REAL(N_ions(s)+total_new_part_number(s)), REAL(N_ions(s)+total_new_part_number(s)+50)))
              CALL RESIZE_ION_ARRAY(s, new_size)   ! sets max_N_ions(s) = new_size
           END IF
        END DO
! now we don't need to worry about the size of the arrays
! retrieve particles from buffer
        i = 1       
        DO n = N_processes_cluster-N_processes_to_free_in_cluster, N_processes_cluster-1
! electrons
           DO k = N_electrons+1, N_electrons + N_part_to_receive_from_proc(0,n)
              electron(k)%X   =     rbufer(i)
              electron(k)%Y   =     rbufer(i+1)
              electron(k)%VX  =     rbufer(i+2)
              electron(k)%VY  =     rbufer(i+3)
              electron(k)%VZ  =     rbufer(i+4)
              electron(k)%tag = INT(rbufer(i+5))
              i = i + 6

              if ( (electron(k)%X.lt.c_X_area_min).or. &
                   & (electron(k)%X.gt.c_X_area_max).or. &
                   & (electron(k)%Y.lt.c_Y_area_min).or. &
                   & (electron(k)%Y.gt.c_Y_area_max) ) then
                 print '("Process ",i4," : Error-1 in GLOBAL_LOAD_BALANCE : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
                 print '("Process ",i4," : k/N_electrons : ",i8,2x,i8)', Rank_of_process, k, N_electrons
                 print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag
                 errcode=330
                 CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
              end if

           END DO
           N_electrons = N_electrons+N_part_to_receive_from_proc(0,n)
! ions
           DO s = 1, N_spec
              DO k = N_ions(s)+1, N_ions(s) + N_part_to_receive_from_proc(s,n)
                 ion(s)%part(k)%X   =     rbufer(i)
                 ion(s)%part(k)%Y   =     rbufer(i+1)
                 ion(s)%part(k)%VX  =     rbufer(i+2)
                 ion(s)%part(k)%VY  =     rbufer(i+3)
                 ion(s)%part(k)%VZ  =     rbufer(i+4)
                 ion(s)%part(k)%tag = INT(rbufer(i+5))
                 i = i + 6
                 
                 if ( (ion(s)%part(k)%X.lt.c_X_area_min).or. &
                      & (ion(s)%part(k)%X.gt.c_X_area_max).or. &
                      & (ion(s)%part(k)%Y.lt.c_Y_area_min).or. &
                      & (ion(s)%part(k)%Y.gt.c_Y_area_max) ) then
                    print '("Process ",i4," : Error-2 in GLOBAL_LOAD_BALANCE : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
                    print '("Process ",i4," : s/k/N_ions : ",i8,2x,i8,2x,i8)', Rank_of_process, s, k, N_ions(s)
                    print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag
                    errcode=331
                    CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
                 end if

              END DO
              N_ions(s) = N_ions(s)+N_part_to_receive_from_proc(s,n)
           END DO
        END DO

        DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     END IF

  ELSE  ! IF (cluster_rank_key.EQ.0) THEN

! receive THE MASTERS MESSAGE from the master
     CALL MPI_BCAST(new_particle_masters, N_processes_cluster-1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
! check part of the message related to this process

print '("MASTER MESSAGE :: Rank_of_process ",i4," Rank_cluster ",i4," particle_master old/new ",i4," / ",i4)', Rank_of_process, Rank_cluster, particle_master, new_particle_masters(Rank_cluster)  !#######

     IF (particle_master.NE.new_particle_masters(Rank_cluster)) THEN
! process will be released, send all particles to the master
! send number of particles
        ALLOCATE(ibufer(0:N_spec), STAT=ALLOC_ERR)
        ibufer(0) = N_electrons
        ibufer(1:N_spec) = N_ions(1:N_spec)
        CALL MPI_SEND(ibufer(0:N_spec), N_spec+1, MPI_INTEGER, 0, 0, COMM_CLUSTER, ierr) 
        DEALLOCATE(ibufer, STAT=ALLOC_ERR)
        sum_N_part_to_send = N_electrons
        DO s = 1, N_spec
           sum_N_part_to_send = sum_N_part_to_send + N_ions(s)
        END DO
! send particles
        IF (sum_N_part_to_send.GT.0) THEN
           ALLOCATE(rbufer(1:sum_N_part_to_send*6), STAT=ALLOC_ERR)
! place electrons first
           i=1
           DO k = 1, N_electrons
              rbufer(i) = electron(k)%X
              rbufer(i+1) = electron(k)%Y
              rbufer(i+2) = electron(k)%VX
              rbufer(i+3) = electron(k)%VY
              rbufer(i+4) = electron(k)%VZ
              rbufer(i+5) = DBLE(electron(k)%tag)
              i = i + 6
           END DO
! place ions
           DO s = 1, N_spec
              DO k = 1, N_ions(s)
                 rbufer(i)   = ion(s)%part(k)%X
                 rbufer(i+1) = ion(s)%part(k)%Y
                 rbufer(i+2) = ion(s)%part(k)%VX
                 rbufer(i+3) = ion(s)%part(k)%VY
                 rbufer(i+4) = ion(s)%part(k)%VZ
                 rbufer(i+5) = DBLE(ion(s)%part(k)%tag)
                 i = i+6
              END DO
           END DO
! send particles
           CALL MPI_SEND(rbufer, sum_N_part_to_send*6, MPI_DOUBLE_PRECISION, 0, Rank_cluster, COMM_CLUSTER, ierr)     
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
! clear the counters
        N_electrons = 0
        N_ions = 0
!print '("GLOBAL_LOAD_BALANCE-1 : ",i4,2x,i4,2x,i8,2x,i4,2x,i4)', Rank_of_process, Rank_cluster, N_electrons, particle_master, new_particle_masters(Rank_cluster)
! change master code
        particle_master = new_particle_masters(Rank_cluster)
     END IF

  END IF !IF (cluster_rank_key.EQ.0) THEN /ELSE
  
  DEALLOCATE(new_particle_masters, STAT=ALLOC_ERR)

!  CALL MPI_COMM_FREE(COMM_HORIZONTAL, ierr)
  CALL MPI_COMM_DISCONNECT(COMM_HORIZONTAL, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!  CALL MPI_COMM_FREE(COMM_CLUSTER, ierr)
  CALL MPI_COMM_DISCONNECT(COMM_CLUSTER, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  CALL SET_COMMUNICATIONS

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  CALL DISTRIBUTE_CLUSTER_PARAMETERS
!  CALL BALANCE_LOAD_WITHIN_CLUSTER

!print '("GLOBAL_LOAD_BALANCE-2 : ",i4,2x,i4,2x,i8,2x,i4,2x,i8)', Rank_of_process, Rank_cluster, N_electrons, particle_master, max_N_electrons

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  if (Rank_cluster.ne.0) return

  call mpi_comm_group(COMM_CLUSTER, clustergroup, ierr)
  call mpi_comm_group(MPI_COMM_WORLD, worldgroup, ierr)
  call mpi_group_size(clustergroup, clustergroup_size, ierr)
  allocate(clusterranks(clustergroup_size))
  do i = 1, clustergroup_size
     clusterranks(i) = i-1
  end do
  allocate(worldranks(clustergroup_size))
  call mpi_group_translate_ranks(clustergroup, clustergroup_size, clusterranks, worldgroup, worldranks, ierr)
  print '(">>>>>>>>>>>",2x,i4,2x,i4,10(2x,i4))', Rank_of_process, clustergroup_size, worldranks
  if (allocated(clusterranks)) deallocate(clusterranks)
  if (allocated(worldranks)) deallocate(worldranks)
  call mpi_group_free(clustergroup, ierr)
  call mpi_group_free(worldgroup, ierr)

END SUBROUTINE GLOBAL_LOAD_BALANCE

!------------------------------------------
! this subroutine is called by the master process with zero rank in the COMM_HORIZONTAL communicator for masters
! it is expected that the values of cluster%N_particles are already known
!
SUBROUTINE CALCULATE_N_PROCESSES_BALANCED(N_of_all_free_processes)

  USE ParallelOperationValues  
  USE LoadBalancing

  use mpi

  IMPLICIT NONE

  INTEGER errcode,ierr

  INTEGER N_of_all_free_processes

  REAL, PARAMETER :: max_allowed_avg_load_factor = 1.2   ! factor shows how much the fluctuation of the load 
                                                         ! can exceed the ideal average load value

  REAL maximal_process_load, new_maximal_process_load
  REAL max_allowed_avg_load

  INTEGER m, n

  REAL avg_load
  REAL global_ideal_avg_load

  INTEGER max_N_proc_per_cluster
  INTEGER N_proc_to_release
  INTEGER N_proc_to_assign
  INTEGER m_min_load, m_start
  REAL avg_load_min

  INTEGER m_max_load
  LOGICAL error_detected

! find maximal process load
  maximal_process_load = 0.0
  DO m = 0, N_processes_horizontal - 1
     avg_load = REAL(cluster(m)%N_particles) / REAL(cluster(m)%N_processes)
     IF (avg_load.GT.maximal_process_load) maximal_process_load = avg_load
  END DO

! calculate ideal average load
  global_ideal_avg_load = 0.0
  DO m = 0, N_processes_horizontal-1
     global_ideal_avg_load = global_ideal_avg_load + cluster(m)%N_particles
  END DO
  global_ideal_avg_load = global_ideal_avg_load / REAL(N_of_processes)
  max_allowed_avg_load = max_allowed_avg_load_factor * global_ideal_avg_load

  max_N_proc_per_cluster = N_of_processes - N_processes_horizontal + 1  ! don't need this elsewhere so "recalculate, don't store"

! find the smallest number of processes which ensures that max_allowed_avg_load is not exceeded
  DO m = 0, N_processes_horizontal-1
     DO n = 1, max_N_proc_per_cluster
        avg_load = REAL(cluster(m)%N_particles) / REAL(n)
        IF (avg_load.LT.max_allowed_avg_load) EXIT
     END DO
! the way it is done now is expected to work both in ordinary situation (the exit condition goes off) and 
! if the whole loop above finished but the exit condition was not achieved
     cluster(m)%N_processes_balanced = MIN(n, max_N_proc_per_cluster)
     cluster(m)%avg_load_balanced = avg_load
  END DO

! find numbers of processes to be released and processes to be assigned
  N_proc_to_release=0
  N_proc_to_assign=0
  DO m = 0, N_processes_horizontal-1
     IF (cluster(m)%N_processes.GT.cluster(m)%N_processes_balanced) THEN
        N_proc_to_release = N_proc_to_release + cluster(m)%N_processes - cluster(m)%N_processes_balanced
     ELSE IF (cluster(m)%N_processes.LT.cluster(m)%N_processes_balanced) THEN
        N_proc_to_assign = N_proc_to_assign - cluster(m)%N_processes + cluster(m)%N_processes_balanced
     END IF
  END DO

  IF (N_proc_to_release.LT.N_proc_to_assign) THEN

! a potentially dangerous situation we want to avoid is that N_proc_to_assign can exceed maximal number of assigned processes in a cluster
! at the same time N_proc_to_release will never exceed this value because it cannot exceed the current number of assigned processes 
! for each cluster so the sum cannot exceed the current total number of assigned processes which is expected to be valid
! this is an extreme case of another dangerous situation when the algorithm defines the number of assigned processes larger than 
! the number of released processes
! a universal solution is to reconsider the NEW "BALANCED" system and release some additional processes where appropriate

! release additional processes
     DO n = 1, N_proc_to_assign-N_proc_to_release
! find process with lowest average balanced load after an additional process is released 
! among processes with balanced number of processes exceeding 1 (cannot use clusters with only 1 process)
        DO m = 0, N_processes_horizontal-1
! here we simply find first cluster with more than one process to start
           IF (cluster(m)%N_processes_balanced.GT.1) THEN 
              m_min_load = m
              avg_load_min = REAL(cluster(m_min_load)%N_particles) / REAL(cluster(m_min_load)%N_processes_balanced-1)
              EXIT
           END IF
        END DO
        m_start = m_min_load+1  ! allows to immediately skip first clusters if they have one process only
        DO m = m_start, N_processes_horizontal-1
! here we find a cluster with more than one process where removing one process will create minimal increase of the average load
           IF (cluster(m)%N_processes_balanced.EQ.1) CYCLE
           avg_load = REAL(cluster(m)%N_particles) / REAL(cluster(m)%N_processes_balanced-1)
           IF (avg_load.LT.avg_load_min) THEN
              m_min_load = m
              avg_load_min = avg_load
           END IF
        END DO
! release one process
        cluster(m_min_load)%N_processes_balanced = cluster(m_min_load)%N_processes_balanced-1
        cluster(m_min_load)%avg_load_balanced    = REAL(cluster(m_min_load)%N_particles) / REAL(cluster(m_min_load)%N_processes_balanced)
     END DO

  ELSE IF (N_proc_to_release.GT.N_proc_to_assign) THEN
  
! another possible situation is that the number of released processes exceeds the number of assigned processes
! then we can assign some additional processes to clusters with the highest load
     
     DO n = 1, N_proc_to_release-N_proc_to_assign
! find cluster with maximal average load
        m_max_load = 0
        DO m = 1, N_processes_horizontal-1
           IF (cluster(m)%avg_load_balanced.GT.cluster(m_max_load)%avg_load_balanced) THEN
              m_max_load=m
           END IF
        END DO
! assign an additional process to this cluster
        cluster(m_max_load)%N_processes_balanced = cluster(m_max_load)%N_processes_balanced + 1
        cluster(m_max_load)%avg_load_balanced    = REAL(cluster(m_max_load)%N_particles) / REAL(cluster(m_max_load)%N_processes_balanced)
     END DO

  END IF

! error check >>>>>>>>>>>>>>>>>>>>>>>>>>>>

  error_detected = .FALSE.

! recalculate numbers of processes to be released and processes to be assigned
  N_proc_to_release=0
  N_proc_to_assign=0
  DO m = 0, N_processes_horizontal-1
     IF (cluster(m)%N_processes.GT.cluster(m)%N_processes_balanced) THEN
        N_proc_to_release = N_proc_to_release + cluster(m)%N_processes - cluster(m)%N_processes_balanced
     ELSE IF (cluster(m)%N_processes.LT.cluster(m)%N_processes_balanced) THEN
        N_proc_to_assign = N_proc_to_assign - cluster(m)%N_processes + cluster(m)%N_processes_balanced
     END IF
  END DO

! compare the numbers, they must be equal at this stage     
  IF (N_proc_to_assign.NE.N_proc_to_release) THEN
     PRINT '("### Error-1 in CALCULATE_N_PROCESSES_BALANCED: number of assigned processes ",i9," is not equal to number of released processes ",i9)', &
          & N_proc_to_assign, N_proc_to_release
     error_detected = .TRUE.
  END IF

! check that the total number of processes after balancing did not change
! and that there is no clusters with zero or negative number of processes
  n=0
  DO m = 0, N_processes_horizontal-1
     IF (cluster(m)%N_processes_balanced.LT.1) THEN
        PRINT '("### Error-2 in CALCULATE_N_PROCESSES_BALANCED: cluster ",i8," has wrong balanced number of processes ",i8)', &
          & m, cluster(m)%N_processes_balanced
        error_detected = .TRUE.
     END IF
     n = n + cluster(m)%N_processes_balanced
  END DO
  IF (n.NE.N_of_processes) THEN
     PRINT '("### Error-3 in CALCULATE_N_PROCESSES_BALANCED: the total number of processes after balancing ",i8," is not the initial number ",i8)', &
          & n, N_of_processes
     error_detected = .TRUE.
  END IF

  IF (error_detected) THEN
     OPEN (40, FILE = 'Error_in_CALCULATE_N_PROCESSES_BALANCED.dat')
     DO m = 0, N_processes_horizontal-1
        WRITE (40, '(4(2x,i9),2x,f10.1)') &
             & m, &
             & cluster(m)%N_particles, &
             & cluster(m)%N_processes, &
             & cluster(m)%N_processes_balanced, &
             & cluster(m)%avg_load_balanced 
     END DO
     CLOSE (40, STATUS = 'KEEP')
     PRINT '("### Check file Error_in_CALCULATE_N_PROCESSES_BALANCED.dat")'
     PRINT '("### Terminating the program...")'
     errcode=332
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

! find new maximal process load
  new_maximal_process_load = 0.0
  DO m = 0, N_processes_horizontal-1
     IF (cluster(m)%avg_load_balanced.GT.new_maximal_process_load) new_maximal_process_load = cluster(m)%avg_load_balanced
  END DO

! apply the balancing only if it actually reduces the maximal average load
  IF (new_maximal_process_load.LT.maximal_process_load) THEN
     PRINT '("CALCULATE_N_PROCESSES_BALANCED: maximal load before ",f10.1," is larger than after ",f10.1," : balancing APPROVED")', maximal_process_load, new_maximal_process_load
     N_of_all_free_processes = N_proc_to_release
  ELSE
     PRINT '("CALCULATE_N_PROCESSES_BALANCED: maximal load before ",f10.1," is not larger than after ",f10.1," : balancing CANCELLED")', maximal_process_load, new_maximal_process_load
     DO m = 0, N_processes_horizontal - 1
        cluster(m)%N_processes_balanced = cluster(m)%N_processes
     END DO  
     N_of_all_free_processes = 0
  END IF

END SUBROUTINE CALCULATE_N_PROCESSES_BALANCED

!---------------------------------------------
!
SUBROUTINE BALANCE_LOAD_WITHIN_CLUSTER

  USE ParallelOperationValues
  USE ElectronParticles
  USE IonParticles

  use mpi

  IMPLICIT NONE

  INTEGER errcode,ierr
  INTEGER stattus(MPI_STATUS_SIZE)

  INTEGER initial_N_particles(0:N_spec)        ! used for diagnostics only

  INTEGER ALLOC_ERR
  INTEGER, ALLOCATABLE :: N_particles_spec_proc(:,:)
  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ibufer(0:N_spec)
  
  INTEGER m, k, s
  INTEGER avg_N_particles(0:N_spec)
  INTEGER N_to_distribute(0:N_spec)
  INTEGER total_N_to_distribute
  INTEGER distr_pos(0:N_spec)
  INTEGER max_distr_pos(0:N_spec)

  INTEGER message_size, probed_message_size
  INTEGER pos1, pos2

  INTEGER delta_N, ibuf
  INTEGER new_size

  IF (N_processes_cluster.LE.1) RETURN

  initial_N_particles(0) = N_electrons
  initial_N_particles(1:N_spec) = N_ions(1:N_spec)

  IF (Rank_cluster.EQ.0) THEN
     
     ALLOCATE(N_particles_spec_proc(0:N_spec,1:N_processes_cluster-1), STAT=ALLOC_ERR)

! receive number of particles from each particle calculator   
     DO m = 1, N_processes_cluster-1
        CALL MPI_RECV(N_particles_spec_proc(0:N_spec,m), N_spec+1, MPI_INTEGER, m, m, COMM_CLUSTER, stattus, ierr)
     END DO

! calculate average load for each species
! electrons
     avg_N_particles(0) = N_electrons
     DO m = 1, N_processes_cluster-1
        avg_N_particles(0) = avg_N_particles(0) + N_particles_spec_proc(0,m)
     END DO
     avg_N_particles(0) = avg_N_particles(0) / N_processes_cluster
! ions
     DO s = 1, N_spec
        avg_N_particles(s) = N_ions(s)
        DO m = 1, N_processes_cluster-1
           avg_N_particles(s) = avg_N_particles(s) + N_particles_spec_proc(s,m)
        END DO
        avg_N_particles(s) = avg_N_particles(s) / N_processes_cluster
     END DO     

! calculate the number of particles of each species to be acquired from all processes where N_particles > avg_N_particles
! electrons
     N_to_distribute(0) = MAX(N_electrons - avg_N_particles(0), 0)
     DO m = 1, N_processes_cluster-1
        N_to_distribute(0)  = N_to_distribute(0) + MAX(N_particles_spec_proc(0,m) - avg_N_particles(0), 0)
     END DO
! ions
     DO s = 1, N_spec
        N_to_distribute(s) = MAX(N_ions(s) - avg_N_particles(s), 0)
        DO m = 1, N_processes_cluster-1
           N_to_distribute(s)  = N_to_distribute(s) + MAX(N_particles_spec_proc(s,m) - avg_N_particles(s), 0)
        END DO
     END DO

! calaculate total number of particles to be acquired
     total_N_to_distribute = N_to_distribute(0)
     DO s = 1, N_spec
        total_N_to_distribute = total_N_to_distribute + N_to_distribute(s)
     END DO

! create a buffer of proper size where all particles to distribute will be temporary stored
     ALLOCATE(rbufer(1:6*total_N_to_distribute), STAT=ALLOC_ERR)

! if the master process has excessive number of particles, move them to the buffer
! prepare start positions for each species
     distr_pos(0) = 1
     DO s = 1, N_spec
        distr_pos(s) = distr_pos(s-1) + 6*N_to_distribute(s-1)
     END DO
! electrons
     IF (N_electrons.GT.avg_N_particles(0)) THEN
        DO k = avg_N_particles(0)+1, N_electrons
           rbufer(distr_pos(0))   = electron(k)%X
           rbufer(distr_pos(0)+1) = electron(k)%Y
           rbufer(distr_pos(0)+2) = electron(k)%VX
           rbufer(distr_pos(0)+3) = electron(k)%VY
           rbufer(distr_pos(0)+4) = electron(k)%VZ
           rbufer(distr_pos(0)+5) = electron(k)%tag
           distr_pos(0) = distr_pos(0)+6
        END DO
        N_electrons = avg_N_particles(0)
     END IF
! ions
     DO s = 1, N_spec
        IF (N_ions(s).GT.avg_N_particles(s)) THEN
           DO k = avg_N_particles(s)+1, N_ions(s)
              rbufer(distr_pos(s))   = ion(s)%part(k)%X
              rbufer(distr_pos(s)+1) = ion(s)%part(k)%Y
              rbufer(distr_pos(s)+2) = ion(s)%part(k)%VX
              rbufer(distr_pos(s)+3) = ion(s)%part(k)%VY
              rbufer(distr_pos(s)+4) = ion(s)%part(k)%VZ
              rbufer(distr_pos(s)+5) = ion(s)%part(k)%tag
              distr_pos(s) = distr_pos(s)+6
           END DO
           N_ions(s) = avg_N_particles(s)
        END IF
     END DO

! send the balanced [average] number of electrons back to particle calculators
     CALL MPI_BCAST(avg_N_particles(0:N_spec), N_spec+1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

! receive particles from calculators which have excessive particles 
     DO m = 1, N_processes_cluster-1
        DO s = 0, N_spec
           IF (N_particles_spec_proc(s,m).LE.avg_N_particles(s)) CYCLE
           message_size = 6 * (N_particles_spec_proc(s,m) - avg_N_particles(s))
           pos1 = distr_pos(s)
           pos2 = distr_pos(s)+message_size-1

           CALL MPI_PROBE(m, m+SHIFT1+s, COMM_CLUSTER, stattus, ierr)
           CALL MPI_GET_COUNT(stattus, MPI_DOUBLE_PRECISION, probed_message_size, ierr)
           IF (message_size.NE.probed_message_size) THEN
              PRINT '("Proc ",i4," :: Error-1 in BALANCE_LOAD_WITHIN_CLUSTER :: ",2x,i4,2x,i8,2x,i8,2x,i2)', &
                   & Rank_of_process, m, message_size, probed_message_size, s
              errcode=333
              CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
           END IF

           CALL MPI_RECV(rbufer(pos1:pos2), message_size, MPI_DOUBLE_PRECISION, m, m+SHIFT1+s, COMM_CLUSTER, stattus, ierr)
           distr_pos(s) = distr_pos(s) + message_size
        END DO
     END DO

! send particles to calculators where the number of particles is lower than the balanced value
! at this time, the rbufer is full, so we need to reset the starting positions
     distr_pos(0) = 1
     DO s = 1, N_spec
        distr_pos(s) = distr_pos(s-1) + 6*N_to_distribute(s-1)
     END DO
     DO m = 1, N_processes_cluster-1
        DO s = 0, N_spec
           IF (avg_N_particles(s).LE.N_particles_spec_proc(s,m)) CYCLE
           message_size = 6 * (avg_N_particles(s) - N_particles_spec_proc(s,m))
           pos1 = distr_pos(s)
           pos2 = distr_pos(s)+message_size-1
           CALL MPI_SEND(rbufer(pos1:pos2), message_size, MPI_DOUBLE_PRECISION, m, s, COMM_CLUSTER, ierr)     
           distr_pos(s) = distr_pos(s) + message_size
        END DO
     END DO

! if at this stage there are unprocessed particles left in rbufer,
! these particles must be assigned to the master process

! prepare maximal values of index for each species
     max_distr_pos(0) = 6 * N_to_distribute(0)
     DO s = 1, N_spec
        max_distr_pos(s) = max_distr_pos(s-1) + 6 * N_to_distribute(s)
     END DO

! electrons
     IF (distr_pos(0).LT.max_distr_pos(0)) THEN
        delta_N = (max_distr_pos(0)-distr_pos(0)+1)/6
! check and resize if necessary the size of the electron array
        IF ((N_electrons+delta_N).GT.max_N_electrons) THEN
           new_size = INT(MAX(1.05*REAL(N_electrons+delta_N), REAL(N_electrons+delta_N)+50))
           CALL RESIZE_ELECTRON_ARRAY(new_size)   ! sets max_N_electrons = new_size
        END IF
        DO k = N_electrons+1, N_electrons+delta_N
           electron(k)%X   =     rbufer(distr_pos(0))
           electron(k)%Y   =     rbufer(distr_pos(0)+1)
           electron(k)%VX  =     rbufer(distr_pos(0)+2)
           electron(k)%VY  =     rbufer(distr_pos(0)+3)
           electron(k)%VZ  =     rbufer(distr_pos(0)+4)
           electron(k)%tag = INT(rbufer(distr_pos(0)+5))
           distr_pos(0) = distr_pos(0) + 6
        END DO
        N_electrons = N_electrons+delta_N
     END IF

! ions 
     DO s = 1, N_spec
        IF (distr_pos(s).LT.max_distr_pos(s)) THEN
           delta_N = (max_distr_pos(s)-distr_pos(s)+1)/6
! check and resize if necessary the size of the electron array
           IF ((N_ions(s)+delta_N).GT.max_N_ions(s)) THEN
              new_size = INT(MAX(1.05*REAL(N_ions(s)+delta_N), REAL(N_ions(s)+delta_N+50)))
              CALL RESIZE_ION_ARRAY(s, new_size)   ! sets max_N_ions(s) = new_size                     !### bug found, was RESIZE_ION_ARRAY(new_size)
           END IF
           DO k = N_ions(s)+1, N_ions(s)+delta_N
              ion(s)%part(k)%X   =     rbufer(distr_pos(s))
              ion(s)%part(k)%Y   =     rbufer(distr_pos(s)+1)
              ion(s)%part(k)%VX  =     rbufer(distr_pos(s)+2)
              ion(s)%part(k)%VY  =     rbufer(distr_pos(s)+3)
              ion(s)%part(k)%VZ  =     rbufer(distr_pos(s)+4)
              ion(s)%part(k)%tag = INT(rbufer(distr_pos(s)+5))
              distr_pos(s) = distr_pos(s) + 6
           END DO
           N_ions(s) = N_ions(s)+delta_N
        END IF
     END DO

     DEALLOCATE(N_particles_spec_proc, STAT=ALLOC_ERR)

  ELSE

! send number of particles to the master
     ibufer(0) = N_electrons
     ibufer(1:N_spec) = N_ions(1:N_spec)
     CALL MPI_SEND(ibufer(0:N_spec), N_spec+1, MPI_INTEGER, 0, Rank_cluster, COMM_CLUSTER, ierr)

! receive the balanced number of particles
     CALL MPI_BCAST(avg_N_particles(0:N_spec), N_spec+1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)

! send extra electrons, if applicable
     IF (N_electrons.GT.avg_N_particles(0)) THEN
        message_size = 6 * (N_electrons - avg_N_particles(0))
        ALLOCATE(rbufer(1:message_size), STAT=ALLOC_ERR)
        ibuf = 1
        DO k = avg_N_particles(0)+1, N_electrons
           rbufer(ibuf)   =      electron(k)%X
           rbufer(ibuf+1) =      electron(k)%Y
           rbufer(ibuf+2) =      electron(k)%VX
           rbufer(ibuf+3) =      electron(k)%VY
           rbufer(ibuf+4) =      electron(k)%VZ
           rbufer(ibuf+5) = DBLE(electron(k)%tag)
           ibuf = ibuf + 6
        END DO
        N_electrons = avg_N_particles(0)
        CALL MPI_SEND(rbufer, message_size, MPI_DOUBLE_PRECISION, 0, Rank_cluster+SHIFT1+0, COMM_CLUSTER, ierr)
        DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     END IF

! send extra ions if applicable
     DO s = 1, N_spec
        IF (N_ions(s).GT.avg_N_particles(s)) THEN
           message_size = 6 * (N_ions(s) - avg_N_particles(s))
           ALLOCATE(rbufer(1:message_size), STAT=ALLOC_ERR)
           ibuf = 1
           DO k = avg_N_particles(s)+1, N_ions(s)
              rbufer(ibuf)   =      ion(s)%part(k)%X
              rbufer(ibuf+1) =      ion(s)%part(k)%Y
              rbufer(ibuf+2) =      ion(s)%part(k)%VX
              rbufer(ibuf+3) =      ion(s)%part(k)%VY
              rbufer(ibuf+4) =      ion(s)%part(k)%VZ
              rbufer(ibuf+5) = DBLE(ion(s)%part(k)%tag)
              ibuf = ibuf + 6
           END DO
           N_ions(s) = avg_N_particles(s)
           CALL MPI_SEND(rbufer, message_size, MPI_DOUBLE_PRECISION, 0, Rank_cluster+SHIFT1+s, COMM_CLUSTER, ierr)
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        END IF
     END DO

! receive additional electrons, if applicable
     IF (N_electrons.LT.avg_N_particles(0)) THEN
! make sure that particle arrays have size enough to accommodate new electrons
        IF (avg_N_particles(0).GT.max_N_electrons) THEN
           new_size = INT(MAX(1.05*REAL(avg_N_particles(0)), REAL(avg_N_particles(0)+50)))
           CALL RESIZE_ELECTRON_ARRAY(new_size)
        END IF
        message_size = 6 * (avg_N_particles(0) - N_electrons)

        CALL MPI_PROBE(0, 0, COMM_CLUSTER, stattus, ierr)
        CALL MPI_GET_COUNT(stattus, MPI_DOUBLE_PRECISION, probed_message_size, ierr)
        IF (message_size.NE.probed_message_size) THEN
           PRINT '("Proc ",i4," :: Error-2 in BALANCE_LOAD_WITHIN_CLUSTER :: ",2x,i4,2x,i8,2x,i8)', &
                & Rank_of_process, m, message_size, probed_message_size
           errcode=334
           CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
        END IF

        ALLOCATE(rbufer(1:message_size), STAT = ALLOC_ERR)
        CALL MPI_RECV(rbufer, message_size, MPI_DOUBLE_PRECISION, 0, 0, COMM_CLUSTER, stattus, ierr)
        ibuf = 1
        DO k = N_electrons+1, avg_N_particles(0)
           electron(k)%X   =     rbufer(ibuf)
           electron(k)%Y   =     rbufer(ibuf+1)
           electron(k)%VX  =     rbufer(ibuf+2)
           electron(k)%VY  =     rbufer(ibuf+3)
           electron(k)%VZ  =     rbufer(ibuf+4)
           electron(k)%tag = INT(rbufer(ibuf+5))
           ibuf = ibuf + 6
        END DO
        N_electrons = avg_N_particles(0)
        DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     END IF

! receive additional ions, if applicable
     DO s = 1, N_spec
        IF (N_ions(s).LT.avg_N_particles(s)) THEN
! make sure that particle arrays have size enough to accommodate new electrons
           IF (avg_N_particles(s).GT.max_N_ions(s)) THEN
              new_size = INT(MAX(1.05*REAL(avg_N_particles(s)), REAL(avg_N_particles(s)+50)))
              CALL RESIZE_ION_ARRAY(s, new_size)
           END IF
           message_size = 6 * (avg_N_particles(s) - N_ions(s))

           CALL MPI_PROBE(0, s, COMM_CLUSTER, stattus, ierr)
           CALL MPI_GET_COUNT(stattus, MPI_DOUBLE_PRECISION, probed_message_size, ierr)
           IF (message_size.NE.probed_message_size) THEN
              PRINT '("Proc ",i4," :: Error-3 in BALANCE_LOAD_WITHIN_CLUSTER :: ",2x,i4,2x,i8,2x,i8,2x,i2)', &
                   & Rank_of_process, m, message_size, probed_message_size, s
              errcode=335
              CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
           END IF

           ALLOCATE(rbufer(1:message_size), STAT = ALLOC_ERR)
           CALL MPI_RECV(rbufer, message_size, MPI_DOUBLE_PRECISION, 0, s, COMM_CLUSTER, stattus, ierr)
           ibuf = 1
           DO k = N_ions(s)+1, avg_N_particles(s)
              ion(s)%part(k)%X   =     rbufer(ibuf)
              ion(s)%part(k)%Y   =     rbufer(ibuf+1)
              ion(s)%part(k)%VX  =     rbufer(ibuf+2)
              ion(s)%part(k)%VY  =     rbufer(ibuf+3)
              ion(s)%part(k)%VZ  =     rbufer(ibuf+4)
              ion(s)%part(k)%tag = INT(rbufer(ibuf+5))
              ibuf = ibuf + 6
           END DO
           N_ions(s) = avg_N_particles(s)
           DEALLOCATE(rbufer, STAT=ALLOC_ERR)
         END IF
     END DO

  END IF

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)

! if the number of particles in a process was shrinking, it may be useful to reduce the size of the array
  new_size = INT( MAX(1.05*REAL(N_electrons), REAL(N_electrons+50)) )
  IF (new_size.LT.max_N_electrons) CALL RESIZE_ELECTRON_ARRAY(new_size)

  DO s = 1, N_spec
     new_size = INT( MAX(1.05*REAL(N_ions(s)), REAL(N_ions(s)+50)) )
     IF (new_size.LT.max_N_ions(s)) CALL RESIZE_ION_ARRAY(s, new_size)
  END DO

print '("Rank_of_process ",i4," Rank_cluster ",i4," particle_master ",i4," N_electrons before/after/max ",i8,"/",i8,"/",i8," N_ions(1) b/a/m ",i8,"/",i8,"/",i8)',  &
     & Rank_of_process, Rank_cluster, particle_master, &
     & initial_N_particles(0), N_electrons, max_N_electrons, &
     & initial_N_particles(1), N_ions(1), max_N_ions(1)

END SUBROUTINE BALANCE_LOAD_WITHIN_CLUSTER

!--------------------------------------------------------------
!
! Note that this subroutine changes the value of max_N_electrons.
! It is expected that before the subroutine is called, 
! max_N_electrons is the valid current size of the array.
!
SUBROUTINE RESIZE_ELECTRON_ARRAY(new_size)

  USE ElectronParticles, ONLY : N_electrons, max_N_electrons, electron

  IMPLICIT NONE

  INTEGER new_size

  INTEGER N_part_to_save

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k

  N_part_to_save = MIN(N_electrons, new_size, max_N_electrons)

  ALLOCATE(bufer(1:N_part_to_save), STAT=ALLOC_ERR)

  DO k = 1, N_part_to_save
     bufer(k)%X   = electron(k)%X
     bufer(k)%Y   = electron(k)%Y
     bufer(k)%VX  = electron(k)%VX
     bufer(k)%VY  = electron(k)%VY
     bufer(k)%VZ  = electron(k)%VZ
     bufer(k)%tag = electron(k)%tag
  END DO

  DEALLOCATE(electron, STAT=DEALLOC_ERR)

  ALLOCATE(electron(1:new_size), STAT=ALLOC_ERR)

  max_N_electrons = new_size  !###

  DO k = 1, N_part_to_save
     electron(k)%X   = bufer(k)%X
     electron(k)%Y   = bufer(k)%Y
     electron(k)%VX  = bufer(k)%VX
     electron(k)%VY  = bufer(k)%VY
     electron(k)%VZ  = bufer(k)%VZ
     electron(k)%tag = bufer(k)%tag
  END DO

  DEALLOCATE(bufer, STAT=DEALLOC_ERR)

END SUBROUTINE RESIZE_ELECTRON_ARRAY

!--------------------------------------------------------------
!
! Note that this subroutine changes the value of max_N_electrons.
! It is expected that before the subroutine is called, 
! max_N_electrons is the valid current size of the array.
!
SUBROUTINE RESIZE_ION_ARRAY(s, new_size)

  USE IonParticles, ONLY : N_spec, N_ions, max_N_ions, ion

  IMPLICIT NONE

  INTEGER s
  INTEGER new_size

  INTEGER N_part_to_save

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k

  N_part_to_save = MIN(N_ions(s), new_size, max_N_ions(s))

  ALLOCATE(bufer(1:N_part_to_save), STAT=ALLOC_ERR)

  DO k = 1, N_part_to_save
     bufer(k)%X   = ion(s)%part(k)%X
     bufer(k)%Y   = ion(s)%part(k)%Y
     bufer(k)%VX  = ion(s)%part(k)%VX
     bufer(k)%VY  = ion(s)%part(k)%VY
     bufer(k)%VZ  = ion(s)%part(k)%VZ
     bufer(k)%tag = ion(s)%part(k)%tag
  END DO

  DEALLOCATE(ion(s)%part, STAT=DEALLOC_ERR)
  !NULLIFY(ion(s)%part)

  ALLOCATE(ion(s)%part(1:new_size), STAT=ALLOC_ERR)

  max_N_ions(s) = new_size  !###

  DO k = 1, N_part_to_save
     ion(s)%part(k)%X   = bufer(k)%X
     ion(s)%part(k)%Y   = bufer(k)%Y
     ion(s)%part(k)%VX  = bufer(k)%VX
     ion(s)%part(k)%VY  = bufer(k)%VY
     ion(s)%part(k)%VZ  = bufer(k)%VZ
     ion(s)%part(k)%tag = bufer(k)%tag
  END DO

  DEALLOCATE(bufer, STAT=DEALLOC_ERR)

END SUBROUTINE RESIZE_ION_ARRAY
