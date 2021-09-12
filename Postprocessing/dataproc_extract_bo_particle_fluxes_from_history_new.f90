!------------------------------------
! 
! This program is for the newer version of the code which includes inner material objects.
! Running this program with data produced by the older version of the code (which did not have the inner objects) 
! may cause an error while reading the input file init_configuration.dat, 
! especially if there are characters in the line with the number of boundary objects after this number
! (for example some comments).
!
! To process data created by the older code, better use program
!       dataproc_extract_bo_particle_fluxes_from_history.f90
! also available on github.
!
! This program reads files history_bo_NN.dat, where NN is the boundary object number. 
! Such a file contains records of how many electrons and ions of different species 
! collided with the boundary object at each time step. 
! The program calculates average fluxes of the collided and emitted particles in units of the electric current [Amperes, A]
! and saves those as functions of time in files history_bo_NN_avg.dat, here NN is again the boundary object number.
! The averaging is performed over the requested number of ion cycles (the program will ask for it).
! The program also needs the following input data files from the simulation:
! init_configuration.dat
! init_simcontrol.dat
! init_particles.dat 
!
! In simulation output files history_bo_NN.dat 
! column 1 is the time step counter, 
! column 2 is the number of electron macroparticles collided with boundary object NN,
! columns 3:3+N_spec-1 are the numbers of ion macroparticles of species 1:N_spec collided with boundary object NN, 
! column 3+N_spec is the number of electrons emitted by boundary object NN, N_spec is the number of ion species.
!
! In processed files history_bo_NN_avg.dat
! column 1 is the time [ns],
! column 2 is the average flux of electron collided with the boundary object NN [A],
! columns 3:3+N_spec are the average fluxes of ions of species 1:N_spec collided with the boundary object NN [A],
! column 3+N_spec is the average flux of electron emitted the boundary object NN [A], N_spec is the number of ion species.
!
! The fluxes/currents in processed files history_bo_NN_avg.dat are all positive, 
! it is user's responsibility to use the proper sign.
!
! Also, in the top of the processed files there is a line with the length of the corresponding boundary object, 
! which can be used to calculate the current density, if necessary.
!
!----------------------------
!
module physical_constants

  REAL(8), PARAMETER :: e_Cl     = 1.602189d-19      ! Charge of single electron [Cl]
  REAL(8), PARAMETER :: m_e_kg   = 9.109534d-31      ! Mass of single electron [kg]
  REAL(8), PARAMETER :: eps_0_Fm = 8.854188d-12      ! The dielectric constant [F/m]

end module physical_constants

!------------------------------------
!
program get_particle_wall_fluxes

  use physical_constants

  implicit none

  real delta_x_m
  real delta_t_s
  real N_plasma_m3

  integer N_of_particles_cell
  integer N_of_boundary_objects
  integer N_subcycles
  integer N_spec

  real bo_L_m(100)
  integer, allocatable :: intflux(:), intfluxavg(:)
  integer alloc_err

  real dt_ns
  real factor_avgJ_A

  INTEGER k
                                       ! ----x----I----x----I-
  CHARACTER(17) historybo_filename     ! history_bo_NN.dat
  CHARACTER(21) historyboavg_filename  ! history_bo_NN_avg.dat

  logical exists
  integer iostatus

  integer nsub, navgi
  integer N_of_records
  integer idummy

  integer pos, n

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  call read_init_configuration(delta_x_m, delta_t_s, N_plasma_m3, N_of_particles_cell, N_of_boundary_objects, N_subcycles, N_spec, bo_L_m)

  allocate(intflux(0:N_spec+1), stat = alloc_err)
  allocate(intfluxavg(0:N_spec+1), stat = alloc_err)

  print *, "average over how many ion records?"
  read *, navgi

  dt_ns = delta_t_s * 1.0e9

  factor_avgJ_A = (e_Cl * N_plasma_m3 * delta_x_m * delta_x_m) / (delta_t_s * N_of_particles_cell * navgi * N_subcycles)

!-------------------------

  intflux = 0 

  do k = 1, N_of_boundary_objects

     historybo_filename = 'history_bo_NN.dat'
     historybo_filename(12:13) = convert_int_to_txt_string(k, 2)

     historyboavg_filename = 'history_bo_NN_avg.dat'
     historyboavg_filename(12:13) = convert_int_to_txt_string(k, 2)

     inquire (file = historybo_filename, exist = exists)
     if (.not.exists) then
        print '("File ",A17," not found, skip")', historybo_filename
        cycle
     end if

     print '("File ",A17," found, reading the data file ...")', historybo_filename

     open (9, file = historybo_filename)
     N_of_records = 0
     do
        read (9, '(2x,i8,10(2x,i8))', iostat = iostatus) idummy, intflux  !(0:N_spec+1)
        if (iostatus.ne.0) exit
        N_of_records = N_of_records + 1
     end do
     close (9, status = 'keep')

     print '("### The total number of records in the data file is ",i8," (first record has number 1)")', N_of_records 

! re-read data file
     open ( 9, file = historybo_filename)
     open (10, file = historyboavg_filename)
     write (10, '("# whole boundary object length ",f10.8," m")') bo_L_m(k)
     pos=0
     do while (pos.lt.(N_of_records-navgi*N_subcycles))

        intfluxavg = 0
        do n = 1, navgi*N_subcycles
           read (9, '(2x,i8,10(2x,i8))', iostat = iostatus) idummy, intflux  !(0:N_spec+1)
           intfluxavg = intfluxavg + intflux
        end do

        pos=pos+navgi*N_subcycles

        write (10, '(2x,f14.6,10(2x,e12.5))') &
             & (pos-1) * dt_ns, &
             & factor_avgJ_A * intfluxavg

     end do
     close (9, status = 'keep')
     close (10, status = 'keep')

     print '("Reading from file ",A17," done")', historybo_filename
     print '("Created file ",A21," with time and currents in dimensional units...")', historyboavg_filename

  end do

  deallocate(intflux, stat = alloc_err)
  deallocate(intfluxavg, stat = alloc_err)

end program get_particle_wall_fluxes

!--------------------------------
!
subroutine read_init_configuration(delta_x_m, delta_t_s, N_plasma_m3, N_of_particles_cell, N_of_boundary_and_inner_objects, N_subcycles, N_spec, bo_L_m)

  use physical_constants

  implicit none

  real, intent(out) :: delta_x_m
  real, intent(out) :: delta_t_s
  real, intent(out) :: N_plasma_m3

  integer, intent(out) :: N_of_particles_cell
  integer, intent(out) :: N_of_boundary_and_inner_objects
  integer, intent(out) :: N_subcycles
  integer, intent(out) :: N_spec

  integer N_of_boundary_objects
  integer N_of_inner_objects

  real, intent(out) :: bo_L_m(100) ! array with full lengthes of boundary objects [m], 100 should be more than enough 

  logical exists

  character(1) buf

  real T_e_eV
  integer N_of_cells_debye
  integer N_max_vel
  real v_Te_ms
  real W_plasma_s1
  real L_debye_m

  integer idummy
  real rdummy

  integer number_of_segments
  integer n, m, wo_L
  integer istart, jstart, iend, jend
  integer ileft, jbottom, iright, jtop

  INQUIRE (FILE = 'init_configuration.dat', EXIST = exists)
  IF (.NOT.exists) THEN
     PRINT '(2x,"ERROR : init_configuration.dat not found. Program terminated")'
     STOP
  END IF
 
  INQUIRE (FILE = 'init_simcontrol.dat', EXIST = exists)
  IF (.NOT.exists) THEN
     PRINT '(2x,"ERROR : init_simcontrol.dat not found. Program terminated")'
     STOP
  END IF

  INQUIRE (FILE = 'init_particles.dat', EXIST = exists)
  IF (.NOT.exists) THEN
     PRINT '(2x,"ERROR : init_particles.dat not found. Program terminated")'
     STOP
  END IF

  PRINT '(2x,"init_configuration.dat is found. Reading the data file...")'

  OPEN (9, FILE = 'init_configuration.dat')

  READ (9, '(A1)') buf !"---dddd.ddddddd----- scale electron temperature [eV]")')
  READ (9, '(3x,f12.7)') T_e_eV
  READ (9, '(A1)') buf !"---+d.dddddddE+dd--- scale electron density [m^-3]")')
  READ (9, '(3x,e14.7)') N_plasma_m3
  READ (9, '(A1)') buf !"---ddd---------- number of cells per scale electron Debye length")')
  READ (9, '(3x,i3)') N_of_cells_debye
  READ (9, '(A1)') buf !"---ddd---------- maximal expected velocity [units of scale thermal electron velocity]")')
  READ (9, '(3x,i3)') N_max_vel

  v_Te_ms     = SQRT(2.0_8 * DBLE(T_e_eV) * e_Cl / m_e_kg)
  W_plasma_s1 = SQRT(DBLE(N_plasma_m3) * e_Cl**2 / (eps_0_Fm * m_e_kg))
  L_debye_m   = v_Te_ms / W_plasma_s1
  delta_x_m   = L_debye_m / N_of_cells_debye    
  delta_t_s   = delta_x_m / (N_max_vel * v_Te_ms) 

  READ (9, '(A1)') buf !"---ddd---------- number of blocks (processes) along the X (horizontal) direction")')
  READ (9, '(3x,i3)') idummy ! N_blocks_x
  READ (9, '(A1)') buf !"---ddd---------- number of blocks (processes) along the Y (vertical) direction")')
  READ (9, '(3x,i3)') idummy ! N_blocks_y
  READ (9, '(A1)') buf !"---ddd---------- number of cells along the X-direction in a block")')
  READ (9, '(3x,i3)') idummy ! N_grid_block_x
  READ (9, '(A1)') buf !"---ddd---------- number of cells along the Y-direction in a block")')
  READ (9, '(3x,i3)') idummy ! N_grid_block_y

  READ (9, '(A1)') buf !"--dddd---------- number of macroparticles per cell for the scale density")')
  READ (9, '(2x,i4)') N_of_particles_cell

  READ (9, '(A1)') buf !"-----d---------- number of blocks in a cluster along the X-direction")')
  READ (9, '(5x,i1)') idummy ! cluster_N_blocks_x
  READ (9, '(A1)') buf !"-----d---------- number of blocks in a cluster along the Y-direction")')
  READ (9, '(5x,i1)') idummy ! cluster_N_blocks_y

  N_of_inner_objects = 0

  READ (9, '(A1)') buf !"---ddd---ddd---- number of objects along domain boundary // number of material inner objects (>=0), each inner objects is a rectangle")')
  READ (9, '(3x,i3,3x,i3)') N_of_boundary_objects, N_of_inner_objects

  N_of_boundary_and_inner_objects = N_of_boundary_objects + N_of_inner_objects

! calculate total length of boundary objects
  bo_L_m = 0.0
  DO n = 1, N_of_boundary_objects
     READ (9, '(A1)') buf !"===dd===dd=== object type, number of segments")')
     READ (9, '(3x,i2,3x,i2)') idummy, number_of_segments
     READ (9, '(A1)') buf !"---dddddd---dddddd---dddddd---dddddd---segment start X/Y end X/Y [global node index]")')
     wo_L = 0
     DO m = 1, number_of_segments
        READ (9, '(3x,i6,3x,i6,3x,i6,3x,i6)') istart, jstart, iend, jend
        wo_L = wo_L + ABS(jend - jstart) + ABS(iend - istart)
     END DO
     bo_L_m(n) = real(wo_L * delta_x_m)
  END DO

  DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects  !N_of_inner_objects
     READ (9, '(A1)') buf !"===dd=== object type")')
     READ (9, '(3x,i2)') idummy
     READ (9, '(A1)') buf !"---dddddd---dddddd---dddddd---dddddd--- coordinates of left bottom X/Y corner and right top X/Y corners [global node index]")')
     READ (9, '(3x,i6,3x,i6,3x,i6,3x,i6)') ileft, jbottom, iright, jtop  ! these are left bottom corner and right top corner coordinates of the whole object
     bo_L_m(n) = real((jtop - jbottom + iright - ileft) * 2 * delta_x_m)   ! this is unattenuated full length which does not account that the object may be partially covered by other objects
  END DO

  CLOSE (9, STATUS = 'KEEP')

  PRINT '(2x,"Process ",i5," : init_simcontrol.dat is found. Reading the data file...")'
  
  OPEN (9, FILE = 'init_simcontrol.dat')

  READ (9, '(A1)') buf !"---ddddddd.ddd----- simulation time [ns]")')
  READ (9, '(3x,f11.3)') rdummy ! t_sim_ns
  READ (9, '(A1)') buf !"--------dd--------- number of electron sub-cycles per ion cycle (odd)")')
  READ (9, '(8x,i2)') N_subcycles
!  READ (9, '(A1)') buf !"------dddd--------- number of ion cycles between internal cluster load balancing events")')
!  READ (9, '(6x,i4)') dT_cluster_load_balance
!  READ (9, '(A1)') buf !"------dddd--------- number of internal cluster load balancing events between global load balancing events")')
!  READ (9, '(6x,i4)') dT_global_load_balance
!  READ (9, '(A1)') buf !"---ddddddd--------- number of ion cycles between checkpoints (no checkpoints if <=0)")')
!  READ (9, '(3x,i7)') dT_save_checkpoint
!  READ (9, '(A1)') buf !"---------d--------- use checkpoint (2/1/0 = Yes, to start/Yes, to continue/No)")')
!  READ (9, '(9x,i1)') use_checkpoint
!  READ (9, '(A1)') buf !"--dddddddd--------- time step when the checkpoint to be used was saved (dim-less)")')
!  READ (9, '(2x,i8)') T_cntr_to_continue

  CLOSE (9, STATUS = 'KEEP')

  print '("N_subcycles = ",i4)', N_subcycles

  PRINT '(2x,"Process ",i3," : init_particles.dat is found. Reading the data file...")'
  
  OPEN (9, FILE = 'init_particles.dat')

  READ (9, '(A1)') buf !"---dddd.ddd----- initial electron temperature [eV]")')
  READ (9, '(3x,f8.3)') rdummy ! init_Te_eV
  READ (9, '(A1)') buf !"---+d.dddE+dd--- initial electron density [m^-3]")')
  READ (9, '(3x,e10.3)') rdummy ! init_Ne_m3
  READ (9, '(A1)') buf !"------d--------- number of ion species")')
  READ (9, '(6x,i1)') N_spec

  CLOSE (9, STATUS = 'KEEP')

  print '("N_spec = ",i4)', N_spec

  return

end subroutine read_init_configuration

!-----------------------------------
! creates a string of length "length_of_string" out of an integer number "int_number"
!
function convert_int_to_txt_string(int_number, length_of_string)

  implicit none

  integer length_of_string
  character*(length_of_string) convert_int_to_txt_string
  integer int_number
!  integer length_of_string
  character(5) format_string
  character(2) length2_txt
  character(1) length1_txt

  character*(length_of_string) number_txt
  
  integer blanks_number
  integer i

! produce format string
  if ((length_of_string.gt.0).and.(length_of_string.lt.10)) then
     write (length1_txt, '(i1)') length_of_string
     format_string = '(iN) '
     format_string(3:3) = length1_txt
  else if ((length_of_string.ge.10).and.(length_of_string.lt.100)) then
     write (length2_txt, '(i2)') length_of_string
     format_string = '(iNN)'
     format_string(3:4) = length2_txt
  else
     print *, "ERROR in CONVERT_INT_TO_TXT_STRING:"
     print *, "incorrect string length requested: ", length_of_string
     stop
  end if

  WRITE (number_txt, format_string) int_number
  number_txt = ADJUSTL(TRIM(number_txt))
  blanks_number = length_of_string - LEN_TRIM(number_txt)
  number_txt = ADJUSTR(number_txt)
  do i = 1, blanks_number
     number_txt(i:i) = '0'
  end do

  convert_int_to_txt_string = number_txt

end function convert_int_to_txt_string
