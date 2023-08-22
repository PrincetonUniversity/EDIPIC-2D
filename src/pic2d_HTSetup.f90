!-----------------------------------------
!
SUBROUTINE PREPARE_HT_SETUP_VALUES

  USE SetupValues
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles
  USE ClusterAndItsBoundaries

  use mpi

  IMPLICIT NONE

  LOGICAl exists
  
  CHARACTER(1) buf 

  INTEGER use_e_emission_from_cathode_flag
  INTEGER grid_requested_flag
  INTEGER use_ionization_source_flag

  REAL(8) T_inj_eV
  REAL(8) Y_ion_anode
  REAL(8) Y_ion_cathode
  REAL(8) T_ion_e_eV
  REAL(8) T_ion_i_eV
  REAL(8) N_i_outflow_m3

  INTEGER ALLOC_ERR

  INTEGER s, j
  INTEGER N_part_outflow

! function
  REAL(8) RateDistrFun

  ht_emission_constant = .FALSE.
  ht_use_e_emission_from_cathode = .FALSE.
  ht_use_e_emission_from_cathode_zerogradf = .FALSE.
  ht_use_ionization_source = .FALSE.
  ht_grid_requested = .FALSE.
  ht_soft_grid_requested = .FALSE.

  total_cathode_N_e_to_inject = 0

  INQUIRE (FILE = 'init_setup.dat', EXIST = exists)

  IF (exists) THEN

     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i5," : init_setup.dat is found. Reading the data file...")', Rank_of_process
     END IF

     OPEN (9, FILE = 'init_setup.dat')

     READ (9, '(A1)') buf !"NOTE::")')
     READ (9, '(A1)') buf !"electrons supporting the discharge current can be injected through the cathode boundary")')
     READ (9, '(A1)') buf !"or at the emission plane inside the domain")')
     READ (9, '(A1)') buf !"or at the neighbor nodes of the cathode [one cell from the cathode]")')

     READ (9, '(A1)') buf !"the number of injected electrons is either equal to the number of electrons escaping through the anode boundary")')
     READ (9, '(A1)') buf !"or is a given constant")')
     READ (9, '(A1)') buf !"or makes the difference of potential averaged over x between the cathode and the neighbor node equal zero")')

     READ (9, '(A1)') buf !"a metal grid electrode with given potential can be placed at the injection plane")')

     READ (9, '(A1)') buf !"the ionization area limited along Y is introduced which produces electron-ion pairs")')
     READ (9, '(A1)') buf !"the rate of ionization is set to balance the neutral plasma outflow through the cathode")')
     READ (9, '(A1)') buf !"the density of this flow is an input parameter, see below")')
     READ (9, '(A1)') buf !"the flow velocity is that of ions with the energy corresponding to anode-cathode voltage")')

     READ (9, '(A1)') buf !"------d--------- electron emission mode (0/1/2/3 = OFF/ON,linked to anode current/ON,constant/ON,makes zero potential gradient at cathode )")')
     READ (9, '(6x,i1)') use_e_emission_from_cathode_flag
     READ (9, '(A1)') buf !"---dddd--------- number of macroparticles to be injected each timestep for constant injection [dim-less]")')
     READ (9, '(3x,i4)') N_macro_constant_injection
     READ (9, '(A1)') buf !"---dddd.ddd----- temperature of electron emission [eV]")')
     READ (9, '(3x,f8.3)') T_inj_eV
     READ (9, '(A1)') buf !"---dddd.ddd----- y-coordinate of the emission plane [mm]")')
     READ (9, '(3x,f8.3)') injection_y
     READ (9, '(A1)') buf !"------d--------- place a grid electrode at the emission plane? (1/0/2= Yes/No/Use soft grid)")')
     READ (9, '(6x,i1)') grid_requested_flag
     READ (9, '(A1)') buf !"---dddd.ddd----- potential of the grid electrode [V]")')
     READ (9, '(3x,f8.3)') F_grid
 
     READ (9, '(A1)') buf !"------d--------- use ionization source? (1/0 = Yes/No)")')
     READ (9, '(6x,i1)') use_ionization_source_flag
     READ (9, '(A1)') buf !"----ddd.ddd----- position of anode end of ionization source [%, 0:100]")')
     READ (9, '(3x,f8.3)') Y_ion_anode
     READ (9, '(A1)') buf !"----ddd.ddd----- position of cathode end of ionization source [%, 0:100]")')
     READ (9, '(3x,f8.3)') Y_ion_cathode
     READ (9, '(A1)') buf !"---dddd.ddd----- temperature of electrons produced by the ionization source [eV]")')
     READ (9, '(3x,f8.3)') T_ion_e_eV
     READ (9, '(A1)') buf !"---dddd.ddd----- temperature of ions produced by the ionization source [eV]")')
     READ (9, '(3x,f8.3)') T_ion_i_eV
     READ (9, '(A1)') buf !"-+d.dddE+dd----- number density of neutral plasma flowing out through the cathode boundary [m^-3]")')
     READ (9, '(1x,e10.3)') N_i_outflow_m3
 
     IF (use_e_emission_from_cathode_flag.EQ.1) THEN
        ht_use_e_emission_from_cathode = .TRUE.
     ELSE IF (use_e_emission_from_cathode_flag.EQ.2) THEN
        ht_use_e_emission_from_cathode = .TRUE.
        ht_emission_constant = .TRUE.
     ELSE IF (use_e_emission_from_cathode_flag.EQ.3) THEN
        ht_use_e_emission_from_cathode_zerogradf = .TRUE.
     END IF

     IF (grid_requested_flag.EQ.1) ht_grid_requested = .TRUE.
     IF (grid_requested_flag.EQ.2) ht_soft_grid_requested = .TRUE.
     
     IF (use_ionization_source_flag.NE.0) ht_use_ionization_source = .TRUE.

     CLOSE (9, STATUS = 'KEEP')

  ELSE
     
     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i5," : WARNING : init_setup.dat not found. Default setup assumed")', Rank_of_process
     END IF

     RETURN

  END IF

! place grid in such a way that it is not exactly on any boundary
  grid_j = MAX(1, MIN(INT(injection_y * 0.001_8 / delta_x_m), global_maximal_j-1))

  IF ((injection_y * 0.001_8).LT.(global_maximal_j * delta_x_m)) THEN
     injection_y = DBLE(grid_j)
     ht_injection_inside = .TRUE.
  ELSE
     injection_y = DBLE(global_maximal_j)-1.0d-6
     ht_injection_inside = .FALSE.
  END IF

  F_grid = F_grid / F_scale_V

! allocate arrays
  ALLOCATE(N_to_ionize_total(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(N_to_ionize_cluster(1:N_spec), STAT = ALLOC_ERR)
  ALLOCATE(factor_convert_vion_i(1:N_spec), STAT = ALLOC_ERR)

! calculate node index boundaries for the ionization source

  j_ion_source_1 = INT(REAL(global_maximal_j) * (Y_ion_anode / 100.0))
  j_ion_source_2 = INT(REAL(global_maximal_j) * (Y_ion_cathode / 100.0))

  IF (cluster_rank_key.EQ.0) THEN

! note that the algorithm of calculation of c_j_ion_source_[12] shown below
! sets c_j_ion_source_1>c_j_ion_source_2 if the cluster is out of the ionization region
     c_j_ion_source_1 = MAX(j_ion_source_1, c_indx_y_min)
     IF (c_indx_y_max.LT.global_maximal_j) THEN
        c_j_ion_source_2 = MIN(j_ion_source_2, c_indx_y_max-1)
     ELSE
        c_j_ion_source_2 = MIN(j_ion_source_2, global_maximal_j)
     END IF

! calculate total number of particles to be produced each [ion] time step  {the ion timestep is N_subcycles times the electron timestep)

     s=1  ! do it for a single ion species only, account for multiple ion species later (need to know composition) ########

     N_part_outflow = INT( DBLE(N_subcycles) * &
                         & (N_i_outflow_m3 / N_plasma_m3) * &
                         & (SQRT(ABS(whole_object(2)%phi) * F_scale_V  / (Ms(s) * T_e_eV)) / DBLE(N_max_vel)) * &
                         & DBLE(N_of_particles_cell) * &
                         & DBLE(whole_object(2)%L))

     IF (Rank_of_process.EQ.0) THEN
        PRINT '(" ### Expected ion outflow through the whole cathode :: ",e12.5," [1/m^2s] or ",i6," macroparticles per ",i3," electron time steps [to be balanced by ionization] ###")', &
             & N_i_outflow_m3 * SQRT(2.0_8 * e_Cl * ABS(whole_object(2)%phi) * F_scale_V  / (Ms(s) * m_e_kg)), &
             & N_part_outflow, &
             & N_subcycles
     END IF

     N_to_ionize_total(s) = N_part_outflow

     CALL PrepareYIonizRateDistribIntegral

  END IF

! prepare velocity conversion factors
  factor_convert_vinj = SQRT(T_inj_eV / T_e_eV) / N_max_vel

  factor_convert_vion_e = SQRT(T_ion_e_eV / T_e_eV) / N_max_vel
  factor_convert_vion_i(1:N_spec) = SQRT(T_ion_i_eV / T_e_eV) / (N_max_vel * SQRT(Ms(1:N_spec))) 

  IF (Rank_of_process.EQ.0) THEN
     OPEN (10, FILE = 'setup_ionization_vs_y.dat')
     WRITE (10, '("# column 1 is the y-node number [dim-less]")')
     WRITE (10, '("# column 2 is the y-node coordinate [cm]")')
     WRITE (10, '("# column 3 is the intensity of ionization source [arb.un.]")')
     DO j = 0, global_maximal_j
        WRITE (10, '(2x,i6,2x,f12.9,2x,e12.5)') j, j*delta_x_m*100.0_8, RateDistrFun(DBLE(j))
     END DO
     CLOSE (10, STATUS = 'KEEP')
     PRINT '("Process 0 created file setup_ionization_vs_y.dat")'
  END IF

END SUBROUTINE PREPARE_HT_SETUP_VALUES

!-------------------------------------------------------------------------------------------
! Prepares the tabulated values of integral of the ionization rate distribution function
! is called only by cluster masters
!  
SUBROUTINE PrepareYIonizRateDistribIntegral

  USE ParallelOperationValues, ONLY : Rank_of_process
  USE SetupValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues, ONLY : whole_object
  USE IonParticles, ONLY : N_spec

  use mpi

  IMPLICIT NONE

  INTEGER N_pnts
  REAL(8) y_min, y_max, dy

  INTEGER i, s
  REAL(8) int_all, int_cluster

  LOGICAL IncreaseNpnts
  INTEGER count
  REAL(8), ALLOCATABLE :: F(:)
  INTEGER ALLOC_ERR
  REAL(8) temp

! function
  REAL(8) RateDistrFun

! exclude clusters which do not have the ionization area
  IF (c_j_ion_source_1.GE.c_j_ion_source_2) THEN
     N_to_ionize_cluster = 0
     yi = 0.0_8
     RETURN
  END IF

! ------- integrate the whole ionization area and the part which belongs to this cluster
!
  N_pnts = 20000
  y_min  = DBLE(j_ion_source_1)
  y_max  = DBLE(j_ion_source_2)
  dy = (y_max - y_min) / N_pnts

  int_all = 0.0_8
  DO i = 1, N_pnts-1
     int_all = int_all + RateDistrFun(y_min + (DBLE(i)-0.5_8) * dy)
  END DO
  int_all = dy * (int_all + RateDistrFun(y_max - 0.5_8 * dy))

! ------- integrate part of the ionization area which belongs to this cluster

  y_min = DBLE(c_j_ion_source_1)
  y_max = DBLE(c_j_ion_source_2)
  N_pnts = (y_max - y_min) / dy
  dy = (y_max - y_min) / N_pnts

  int_cluster = 0.0_8
  DO i = 1, N_pnts-1
     int_cluster = int_cluster + RateDistrFun(y_min + (DBLE(i)-0.5_8) * dy)
  END DO
  int_cluster = dy * (int_cluster + RateDistrFun(y_max - 0.5_8 * dy))

! calculate number of e-i(s) pairs to be produced by this cluster

  DO s = 1, N_spec
     N_to_ionize_cluster(s) = N_to_ionize_total(s) * (int_cluster / int_all) * (DBLE(c_indx_x_max-1 - c_indx_x_min) / DBLE(whole_object(2)%L) )
  END DO

  PRINT '("Cluster ",i5," each ion time step will produce ion(s)-electron pairs :: ",5(2x,i7,"(",i1,")"))', Rank_of_process, N_to_ionize_cluster(1), 1

! ----- prepare the tabulated integral for the part of the ionization area belonging to this cluster

  y_min = DBLE(c_j_ion_source_1)
  y_max = DBLE(c_j_ion_source_2)

  IncreaseNpnts = .TRUE.
  N_pnts = 20000
  count = 0
  DO WHILE (IncreaseNpnts)

     ALLOCATE(F(0:N_pnts), STAT = ALLOC_ERR)

     dy = (y_max - y_min) / N_pnts

     F(0) = 0.0_8
     DO i = 1, N_pnts-1
        F(i) = F(i-1) + RateDistrFun(y_min + (DBLE(i)-0.5_8) * dy)
     END DO
     F(N_pnts) = F(N_pnts-1) + RateDistrFun(y_max - 0.5_8 * dy)

     temp = F(N_pnts)
     F = F * c_R_max / temp   ! normalize integral such that F(N_pnts) = c_R_max
     F(N_pnts) = c_R_max

! check that in the normalized integral the difference between values in neighbor points does not exceed 1

     IncreaseNpnts = .FALSE.

     DO i = 1, N_pnts
        IF ((INT(F(i))-INT(F(i-1))).GT.1) THEN
           IncreaseNpnts = .TRUE.
           N_pnts = N_pnts * 2
           count = count + 1
           DEALLOCATE(F, STAT = ALLOC_ERR)

           IF (count.GT.4) THEN
              PRINT '(2x,"Process ",i3," : ERROR-1 in PrepareYIonizRateDistribIntegral !!!")', Rank_of_process
              STOP
           END IF

           EXIT
        END IF
     END DO

  END DO

  yi(0) = y_min
  count = 0
  DO i = 1, N_pnts-1
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        yi(count) = y_min + i * dy
     END IF
  END DO

  IF (count.NE.c_R_max-1) THEN
     PRINT '(2x,"Process ",i3," : ERROR-2 in PrepareYIonizRateDistribIntegral !!!")', Rank_of_process
     STOP
  END IF

  yi(c_R_max) = y_max

  DEALLOCATE(F, STAT = ALLOC_ERR)

END SUBROUTINE PrepareYIonizRateDistribIntegral

!--------------------------------------------
!
SUBROUTINE PERFORM_IONIZATION_HT_SETUP

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE SetupValues
  USE IonParticles, ONLY : N_spec

  USE rng_wrapper

  use mpi

  IMPLICIT NONE


  INTEGER ierr

  INTEGER add_N_i_ionize(1:N_spec)
  INTEGER s, n
  INTEGER tag
  REAL(8) x, y, vx, vy, vz

  IF (.NOT.ht_use_ionization_source) RETURN

! let the master process of each cluster notify its members (particle calculators) about the amount of e-i pairs to be produced
! this allows to easily include time-varying intensity of ionization source

  add_N_i_ionize = 0

!  ibufer(1:N_spec) = N_to_ionize_cluster(1:N_spec)
  CALL MPI_BCAST(N_to_ionize_cluster(1:N_spec), N_spec, MPI_INTEGER, 0, COMM_CLUSTER, ierr)


  DO s = 1, N_spec
! identify number of e-i pairs to be produced in each process  for each ion species
     add_N_i_ionize(s) = N_to_ionize_cluster(s) / N_processes_cluster
     IF (Rank_cluster.EQ.(N_processes_cluster-1)) THEN
        add_N_i_ionize(s) = N_to_ionize_cluster(s) - (N_processes_cluster-1) * (N_to_ionize_cluster(s) / N_processes_cluster)
     END IF

     DO n = 1, add_N_i_ionize(s)

        tag = 0

        x = DBLE(c_indx_x_min) + well_random_number() * DBLE(c_indx_x_max-1-c_indx_x_min)
        x = MIN(MAX(x, DBLE(c_indx_x_min)), DBLE(c_indx_x_max-1))

        CALL GetYCoordIoniz(y)

        CALL GetMaxwellVelocity(vx)
        CALL GetMaxwellVelocity(vy)
        CALL GetMaxwellVelocity(vz)
     
        vx = vx * factor_convert_vion_e
        vy = vy * factor_convert_vion_e
        vz = vz * factor_convert_vion_e

        CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)

        CALL GetMaxwellVelocity(vx)
        CALL GetMaxwellVelocity(vy)
        CALL GetMaxwellVelocity(vz)
     
        vx = vx * factor_convert_vion_i(s)
        vy = vy * factor_convert_vion_i(s)
        vz = vz * factor_convert_vion_i(s)

        CALL ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, tag)
     END DO
  END DO

END SUBROUTINE PERFORM_IONIZATION_HT_SETUP

!-------------------------------------
!
SUBROUTINE GetYCoordIoniz(y)

  USE SetupValues

  USE rng_wrapper

  use mpi

  IMPLICIT NONE

  REAL(8) y

  REAL(8) R
  INTEGER indx
  
  R = c_R_max * well_random_number()

  indx = INT(R)

  IF (indx.LT.c_R_max) THEN
     y = yi(indx) + DBLE(R - indx) * (yi(indx+1) - yi(indx))
  ELSE
     y = yi(c_R_max)
  END IF

  y = MIN(MAX(y, DBLE(c_j_ion_source_1)), DBLE(c_j_ion_source_2))

  RETURN
  
END SUBROUTINE GetYCoordIoniz

!---------------------------------------
!
REAL(8) FUNCTION RateDistrFun(y)

  USE SetupValues, ONLY : j_ion_source_1, j_ion_source_2
  

  use mpi

  IMPLICIT NONE

  REAL(8), PARAMETER :: pi = 3.141592653589793_8

  REAL(8) y, y1, y2, ymax

  y1 = DBLE(j_ion_source_1)
  y2 = DBLE(j_ion_source_2)
  ymax = 0.5_8 * (y1 + y2)

  RateDistrFun = 0.0_8

  IF ((y.LE.y1).OR.(y.GE.y2)) RETURN

!  RateDistrFun = 0.25_8 * (y - y1) * (y2 - y) / (y2 - y1)**2
  RateDistrFun = COS(pi * (y - ymax) / (y2 - y1))

END FUNCTION RateDistrFun

!--------------------------------------------
!
SUBROUTINE PERFORM_ELECTRON_EMISSION_HT_SETUP

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE SetupValues

  USE rng_wrapper

  use mpi

  IMPLICIT NONE


  INTEGER ierr

  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  INTEGER add_N_e_to_emit
  INTEGER n, m, nwo, N_plus, N_minus

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  IF (.NOT.ht_use_e_emission_from_cathode) RETURN

  add_N_e_to_emit = 0

  ALLOCATE(ibufer(1:N_of_boundary_objects), STAT = ALLOC_ERR)

  IF ((cluster_rank_key.EQ.0).AND.(c_N_of_local_object_parts.GT.0)) THEN

! define number of additionally injected particles (zero-rank process of COMM_BOUNDARY communicator)
! note that everywhere it is assuemed that this process also has zero rank in MPI_COMM_WORLD
! this may become a potential source of errors if the configuration will be such that this process will not be on a boundary 

     IF (Rank_of_process.EQ.0) THEN
        ibufer = 0

        IF (ht_emission_constant) THEN

           ibufer(2) = N_macro_constant_injection

        ELSE

           total_cathode_N_e_to_inject = total_cathode_N_e_to_inject + whole_object(4)%electron_hit_count - whole_object(4)%ion_hit_count(1)       !######## hardwired now, must be self-organizing in future

           ibufer(2) = MAX(0, total_cathode_N_e_to_inject)                      ! inject only if total_cathode_N_e_to_inject>0
           
           total_cathode_N_e_to_inject = MIN(total_cathode_N_e_to_inject, 0)    ! flush counter only if emission is permitted

        END IF

        total_cathode_N_e_injected = ibufer(2)                               ! save for diagnostics

     END IF

! send all boundary cluster master processes info about the number of additionally injected particles for each boundary object
     CALL MPI_BCAST(ibufer, N_of_boundary_objects, MPI_INTEGER, 0, COMM_BOUNDARY, ierr)

! the master process determines how many additional particles it has to inject
     DO n = 1, c_N_of_local_object_parts_above
        m = c_index_of_local_object_part_above(n)

        nwo = c_local_object_part(m)%object_number

        IF (nwo.EQ.2) THEN
! we are here if the cluster has a segment of the boundary object that will perform emission (hardwared object #2 now)
! below it is assumed that object #2 has a shape of a straight line
           N_plus = INT(ibufer(2) * REAL( MIN(c_local_object_part(m)%iend, c_indx_x_max-1) - whole_object(nwo)%segment(1)%istart ) / REAL(whole_object(nwo)%L))
           IF (c_local_object_part(m)%iend.EQ.global_maximal_i) N_plus = ibufer(2)

           N_minus = INT(ibufer(2) * REAL( MIN(c_local_object_part(m)%istart, c_indx_x_min) - whole_object(nwo)%segment(1)%istart ) / REAL(whole_object(nwo)%L))

! we use N_plus and N_minus to guarantee that the sum of all emitted particles is always ibufer(2)
           add_N_e_to_emit = N_plus - N_minus

!print '("Cluster ",i4," will emit ",i5," electrons")', Rank_of_process, add_N_e_to_emit

        END IF
     END DO

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! all processes in each cluster receive proper value of add_N_e_to_emit

  ibufer(1) = add_N_e_to_emit
  CALL MPI_BCAST(ibufer(1:1), 1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
  add_N_e_to_emit = ibufer(1)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! clusters that do not have to inject any additional particles may leave

  IF (add_N_e_to_emit.EQ.0) RETURN

! produce particles to be emitted and place them into the add list

! calculate number of particles to be emitted in each process of the cluster
!     temp = add_N_e_to_emit / N_processes_cluster
  IF (Rank_cluster.EQ.N_processes_cluster-1) THEN
     add_N_e_to_emit = add_N_e_to_emit - (N_processes_cluster-1) * (add_N_e_to_emit / N_processes_cluster) !temp
  ELSE
     add_N_e_to_emit = add_N_e_to_emit / N_processes_cluster !temp
  END IF

  DO n = 1, add_N_e_to_emit

     x = DBLE(c_indx_x_min) + well_random_number() * DBLE(c_indx_x_max-1-c_indx_x_min)
     x = MIN(MAX(x, DBLE(c_indx_x_min)), DBLE(c_indx_x_max-1))
!###     y = DBLE(c_indx_y_max)-1.0d-6

     y = injection_y

!     y = DBLE(c_indx_y_min + 1*(c_indx_y_max-c_indx_y_min)/3)   ! 3/4 does not work

     CALL GetMaxwellVelocity(vx)
     vx = vx * factor_convert_vinj  !###0.0_8
     vz = 0.0_8
!###     CALL GetInjMaxwellVelocity(vy)
!###     vy = -vy * factor_convert_vinj !SQRT(T_inj_eV / T_e_eV)   ?????????????
     IF (ht_injection_inside) THEN 
        CALL GetMaxwellVelocity(vy)
        vy = vy * factor_convert_vinj
     ELSE
        CALL GetInjMaxwellVelocity(vy)
        vy = -vy * factor_convert_vinj
     END IF

!???###      vy = vy * factor_convert_vinj !SQRT(T_inj_eV / T_e_eV)   ?????????????

     tag = 0

     CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
  END DO

!print '("process ",i4," of cluster ",i4," emitted ",i5," electrons from top boundary")', Rank_of_process, particle_master, add_N_e_to_emit

END SUBROUTINE PERFORM_ELECTRON_EMISSION_HT_SETUP

!-------------------------------------------------
!
SUBROUTINE PERFORM_ELECTRON_EMISSION_HT_SETUP_ZERO_GRAD_F

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE SetupValues

  USE rng_wrapper

  use mpi

  IMPLICIT NONE


  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER add_N_e_to_emit
 
  INTEGER, ALLOCATABLE :: ibufer(:)
  INTEGER ALLOC_ERR

  REAL(8), ALLOCATABLE :: rbufer(:)
  REAL(8), ALLOCATABLE :: rho_avg_x(:)

  INTEGER i, j

  REAL(8), ALLOCATABLE ::  rhs(:)
  REAL(8), ALLOCATABLE :: myeq(:)

  REAL(8) density_corr_avg_x

  INTEGER n, m, nwo, N_plus, N_minus

  REAL(8), ALLOCATABLE :: temp_c_rho_band(:)
  INTEGER n3

  REAL(8) x, y, vx, vy, vz
  INTEGER tag

  REAL(8) ax_ip1, ax_i

  IF (.NOT.ht_use_e_emission_from_cathode_zerogradf) RETURN

  IF (periodicity_flag.NE.PERIODICITY_X) RETURN

  add_N_e_to_emit = 0

  ALLOCATE(ibufer(1:N_of_boundary_objects), STAT = ALLOC_ERR)

  IF (cluster_rank_key.EQ.0) THEN

! master processes assemble array with charge densities integrated over x in the master process with rank zero

     ALLOCATE(   rbufer(1:global_maximal_j-1), STAT=ALLOC_ERR)
     ALLOCATE(rho_avg_x(1:global_maximal_j-1), STAT=ALLOC_ERR)
     rho_avg_x = 0.0_8
     rbufer = 0.0_8

     DO j = c_indx_y_min+1, c_indx_y_max-1
        DO i = c_indx_x_min+1, c_indx_x_max-1
           rbufer(j) = rbufer(j) - c_rho(i,j) + c_rho_i(i,j)
        END DO
     END DO

     CALL MPI_REDUCE(rbufer, rho_avg_x, global_maximal_j-1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_HORIZONTAL, ierr)

!print '("process ",i4," did CALL MPI_REDUCE(rbufer, rho_avg_x, global_maximal_j-1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_HORIZONTAL, ierr)")', Rank_of_process

! master process with the zero rank solves system of linear equations to get the charge density correction in node j=global_maximal_j-1
! which would make the x-averaged potential at this point equal to the potential at the boundary j=global_maximal_j

     IF (Rank_horizontal.EQ.0) THEN

        ALLOCATE( rhs(1:global_maximal_j-2), STAT=ALLOC_ERR)
        ALLOCATE(myeq(1:global_maximal_j-2), STAT=ALLOC_ERR)

        DO j = 1, global_maximal_j-2
           rhs(j) = -rho_avg_x(j) / DBLE(N_of_particles_cell) 
        END DO
        rhs(1)                  = rhs(1)                  - (global_maximal_i-1) * whole_object(4)%phi !potential_of_bottom_y_boundary_avg_x   ! ??????
        rhs(global_maximal_j-2) = rhs(global_maximal_j-2) - (global_maximal_i-1) * whole_object(2)%phi !potential_of_top_y_boundary_avg_x      ! ??????

        myeq(1) = 1.0_8 / (-2.0_8)
        rhs(1) = rhs(1) * myeq(1)

!        myeq(2) = 1.0_8 / (-2.0_8 - myeq(1))
!        rhs(2) = (rhs(2) - rhs(1)) * myeq(2)
!        myeq(3) = 1.0_8 / (-2.0_8 - myeq(2))
!        rhs(3) = (rhs(3) - rhs(2)) * myeq(3)

        DO j = 2, global_maximal_j-2
           myeq(j) = 1.0_8 / (-2.0_8 - myeq(j-1))
           rhs(j) = (rhs(j) - rhs(j-1)) * myeq(j)
        END DO

!        potenial_at_global_maximal_j_m_2_avg_x = rhs(global_maximal_j-2)

!       density_corr_avg_x = -rho_avg_x(global_maximal_j-1) + potential_of_top_y_boundary_avg_x - potenial_at_global_maximal_j_m_2_avg_x

        density_corr_avg_x = -rho_avg_x(global_maximal_j-1) / DBLE(N_of_particles_cell) + (global_maximal_i-1) * whole_object(2)%phi - rhs(global_maximal_j-2)

        total_cathode_N_e_to_inject = 0
        IF (density_corr_avg_x.LT.0.0_8) THEN
           total_cathode_N_e_to_inject = ABS(density_corr_avg_x * DBLE(N_of_particles_cell))
        END IF

print '(" total number of electron macroparticles to inject from cathode is ",i8,2x,i8)', total_cathode_N_e_to_inject, INT(density_corr_avg_x * DBLE(N_of_particles_cell))

     END IF    !### IF (Rank_horizontal.EQ.0) THEN

  END IF    !### IF (cluster_rank_key.EQ.0) THEN

  IF ((cluster_rank_key.EQ.0).AND.(c_N_of_local_object_parts.GT.0)) THEN

! define number of additionally injected particles (zero-rank process of COMM_BOUNDARY communicator)
! note that everywhere it is assumed that this process also has zero rank in MPI_COMM_WORLD
! this may become a potential source of errors if the configuration will be such that this process will not be on a boundary 

     IF (Rank_of_process.EQ.0) THEN
        ibufer = 0
        ibufer(2) = MAX(0, total_cathode_N_e_to_inject)                      ! inject only if total_cathode_N_e_to_inject>0
        total_cathode_N_e_injected = ibufer(2)                               ! save for diagnostics
     END IF

! send all boundary cluster master processes info about the number of additionally injected particles for each boundary object
     CALL MPI_BCAST(ibufer, N_of_boundary_objects, MPI_INTEGER, 0, COMM_BOUNDARY, ierr)

! a master process determines how many additional particles it has to inject
     DO n = 1, c_N_of_local_object_parts_above
        m = c_index_of_local_object_part_above(n)

        nwo = c_local_object_part(m)%object_number

        IF (nwo.EQ.2) THEN
! we are here if the cluster has a segment of the boundary object that will perform emission (hardwared object #2 now)
! below it is assumed that object #2 has a shape of a straight line
           N_plus = INT(ibufer(2) * REAL( MIN(c_local_object_part(m)%iend, c_indx_x_max-1) - whole_object(nwo)%segment(1)%istart ) / REAL(whole_object(nwo)%L))
           IF (c_local_object_part(m)%iend.EQ.global_maximal_i) N_plus = ibufer(2)

           N_minus = INT(ibufer(2) * REAL( MIN(c_local_object_part(m)%istart, c_indx_x_min) - whole_object(nwo)%segment(1)%istart ) / REAL(whole_object(nwo)%L))

! we use N_plus and N_minus to guarantee that the sum of all emitted particles is always ibufer(2)
           add_N_e_to_emit = N_plus - N_minus

print '("Cluster ",i4," will emit ",i7," electron macroparticles")', Rank_of_process, add_N_e_to_emit

        END IF
     END DO

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! all processes in each cluster receive proper value of add_N_e_to_emit

  ibufer(1) = add_N_e_to_emit
  CALL MPI_BCAST(ibufer(1:1), 1, MPI_INTEGER, 0, COMM_CLUSTER, ierr)
  add_N_e_to_emit = ibufer(1)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! cleanup
  IF (ALLOCATED(ibufer))    DEALLOCATE(ibufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(rbufer))    DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  IF (ALLOCATED(rho_avg_x)) DEALLOCATE(rho_avg_x, STAT = ALLOC_ERR)
  IF (ALLOCATED(rhs))       DEALLOCATE(rhs, STAT = ALLOC_ERR)
  IF (ALLOCATED(myeq))      DEALLOCATE(myeq, STAT = ALLOC_ERR)

! clusters that do not have to inject any additional particles may leave

  IF (add_N_e_to_emit.EQ.0) RETURN

! produce particles to be emitted and place them into the add list

! calculate number of particles to be emitted in each process of the cluster
!     temp = add_N_e_to_emit / N_processes_cluster
  IF (Rank_cluster.EQ.N_processes_cluster-1) THEN
     add_N_e_to_emit = add_N_e_to_emit - (N_processes_cluster-1) * (add_N_e_to_emit / N_processes_cluster) !temp
  ELSE
     add_N_e_to_emit = add_N_e_to_emit / N_processes_cluster !temp
  END IF

print '("      Process ",i4," will emit ",i7," electron macroparticles")', Rank_of_process, add_N_e_to_emit

  ALLOCATE(         rbufer(c_indx_x_min:c_indx_x_max), STAT=ALLOC_ERR)
  ALLOCATE(temp_c_rho_band(c_indx_x_min:c_indx_x_max), STAT=ALLOC_ERR)
  rbufer = 0.0_8
  temp_c_rho_band = 0.0_8

  n3 = c_indx_x_max - c_indx_x_min + 1

  DO n = 1, add_N_e_to_emit

     x = DBLE(c_indx_x_min) + well_random_number() * DBLE(c_indx_x_max-1-c_indx_x_min)
     x = MIN(MAX(x, DBLE(c_indx_x_min)), DBLE(c_indx_x_max-1))
     y = DBLE(c_indx_y_max-1)

     i = INT(x)
     IF (x.EQ.c_X_area_max) i = c_indx_x_max-1

     ax_ip1 = x - DBLE(i)
     ax_i   = 1.0_8 - ax_ip1
     rbufer(i)   = rbufer(i)   + ax_i
     rbufer(i+1) = rbufer(i+1) + ax_ip1

     CALL GetMaxwellVelocity(vx)
     vx = vx * factor_convert_vinj

!     CALL GetInjMaxwellVelocity(vy)
!     vy = -vy * factor_convert_vinj
     CALL GetMaxwellVelocity(vy)
     vy = vy * factor_convert_vinj

     vz = 0.0_8

     tag = 0

     CALL ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, tag)
  END DO

! collect densities from all processes in a cluster
  CALL MPI_REDUCE(rbufer, temp_c_rho_band, n3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_CLUSTER, ierr)

! now cluster masters exchange information about densities in overlapping nodes
  IF (cluster_rank_key.EQ.0) THEN

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:2), STAT=ALLOC_ERR)

     IF (WHITE_CLUSTER) THEN  
! "white processes"

        IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right densities in the right edge
           rbufer(1:2) = temp_c_rho_band(c_indx_x_max-1:c_indx_x_max)
           CALL MPI_SEND(rbufer, 2, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left densities in the left edge
           rbufer(1:2) = temp_c_rho_band(c_indx_x_min:c_indx_x_min+1)
           CALL MPI_SEND(rbufer, 2, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left densities in the vertical line next to the left edge
           CALL MPI_RECV(rbufer, 2, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           temp_c_rho_band(c_indx_x_min:c_indx_x_min+1) = rbufer(1:2)
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right densities in the vertical line next to the right edge
           CALL MPI_RECV(rbufer, 2, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           temp_c_rho_band(c_indx_x_max-1:c_indx_x_max) = rbufer(1:2)
        END IF

     ELSE
! "black" processes

        IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left densities in the vertical line next to the left edge
           CALL MPI_RECV(rbufer, 2, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           temp_c_rho_band(c_indx_x_min:c_indx_x_min+1) = temp_c_rho_band(c_indx_x_min:c_indx_x_min+1) + rbufer(1:2)
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right densities in the vertical line next to the right edge
           CALL MPI_RECV(rbufer, 2, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           temp_c_rho_band(c_indx_x_max-1:c_indx_x_max) = temp_c_rho_band(c_indx_x_max-1:c_indx_x_max) + rbufer(1:2)
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right densities in the right edge
           rbufer(1:2) = temp_c_rho_band(c_indx_x_max-1:c_indx_x_max)
           CALL MPI_SEND(rbufer, 2, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left densities in the left edge
           rbufer(1:2) = temp_c_rho_band(c_indx_x_min:c_indx_x_min+1)
           CALL MPI_SEND(rbufer, 2, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, ierr) 
        END IF

     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)

! adjust densities 
     c_rho(c_indx_x_min:c_indx_x_max,c_indx_y_max-1) = c_rho(c_indx_x_min:c_indx_x_max,c_indx_y_max-1) + temp_c_rho_band(c_indx_x_min:c_indx_x_max) 

  END IF   !### IF (cluster_rank_key.EQ.0) THEN

! cleanup

  CALL PROCESS_ADDED_ELECTRONS

END SUBROUTINE PERFORM_ELECTRON_EMISSION_HT_SETUP_ZERO_GRAD_F




!-------------------------------------------------------------------------------------------
!
SUBROUTINE INITIATE_WALL_DIAGNOSTICS_HT_SETUP

  USE ParallelOperationValues
  USE Checkpoints
  USE Diagnostics
  USE SetupValues, ONLY : ht_use_e_emission_from_cathode, ht_use_e_emission_from_cathode_zerogradf, ht_emission_constant

  use mpi

  IMPLICIT NONE

  LOGICAL exists
  INTEGER i
  INTEGER i_dummy

  IF (Rank_of_process.NE.0) RETURN

  IF ((.NOT.ht_use_e_emission_from_cathode).AND.(.NOT.ht_use_e_emission_from_cathode_zerogradf).AND.(.NOT.ht_emission_constant)) RETURN

! hardwired for objects #2 (cathode) and #4 (anode)

  IF (use_checkpoint.EQ.1) THEN
! start from checkpoint, must trim the time dependences

! upper wall 
     INQUIRE (FILE = 'history_bo_02.dat', EXIST = exists)
     IF (exists) THEN                                                       
        OPEN (21, FILE = 'history_bo_02.dat', STATUS = 'OLD')          
        DO i = 1, N_of_saved_records
           READ (21, '(2x,i8,3(2x,i8))') i_dummy
        END DO
        ENDFILE 21       
        CLOSE (21, STATUS = 'KEEP')        
     END IF

! lower wall 
     INQUIRE (FILE = 'history_bo_04.dat', EXIST = exists)
     IF (exists) THEN                                                       
        OPEN (21, FILE = 'history_bo_04.dat', STATUS = 'OLD')          
        DO i = 1, N_of_saved_records
           READ (21, '(2x,i8,3(2x,i8))') i_dummy
        END DO
        ENDFILE 21       
        CLOSE (21, STATUS = 'KEEP')        
     END IF

  ELSE
! fresh start, empty files, clean up whatever garbage there might be

     OPEN  (21, FILE = 'history_bo_02.dat', STATUS = 'REPLACE')          
     CLOSE (21, STATUS = 'KEEP')

     OPEN  (21, FILE = 'history_bo_04.dat', STATUS = 'REPLACE')          
     CLOSE (21, STATUS = 'KEEP')

  END IF

END SUBROUTINE INITIATE_WALL_DIAGNOSTICS_HT_SETUP

!-------------------------------------------------------------------------------------------
!
SUBROUTINE SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS_HT_SETUP

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE SetupValues, ONLY : total_cathode_N_e_to_inject, total_cathode_N_e_injected, &
                        & ht_use_e_emission_from_cathode, ht_use_e_emission_from_cathode_zerogradf, ht_emission_constant

  use mpi

  IMPLICIT NONE

  IF (Rank_of_process.NE.0) RETURN

  IF ((.NOT.ht_use_e_emission_from_cathode).AND.(.NOT.ht_use_e_emission_from_cathode_zerogradf).AND.(.NOT.ht_emission_constant)) RETURN

! hardwired for objects #2 (cathode) and #4 (anode) and 1 ion species
! object #2 (cathode) emits electrons  

  OPEN (21, FILE = 'history_bo_02.dat', POSITION = 'APPEND')

  WRITE (21, '(2x,i8,3(2x,i8))') &
       & T_cntr, &
       & whole_object(2)%electron_hit_count , &
       & whole_object(2)%ion_hit_count(1), &
       & total_cathode_N_e_injected

  CLOSE (21, STATUS = 'KEEP')
  
! object #4 (anode) does not emit anything at this time         

  OPEN (21, FILE = 'history_bo_04.dat', POSITION = 'APPEND')

  WRITE (21, '(2x,i8,3(2x,i8))') &
       & T_cntr, &
       & whole_object(4)%electron_hit_count , &
       & whole_object(4)%ion_hit_count(1), &
       & total_cathode_N_e_to_inject
  
  CLOSE (21, STATUS = 'KEEP')
  
END SUBROUTINE SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS_HT_SETUP
