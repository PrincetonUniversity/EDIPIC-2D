!------------------------------------------------------------------
! This program reads files with parameters of ion particles that collided with a given boundary object.
! Then this program saves for each segment of the boundary object and for each ion species: 
! a file with the particle flux and particle energy flux as function of coordinate along the segment, and
! a file with the energy spectrum of collided ions (number of particles in energy bins) as function of coordinate along the segment.
!
! During a simulation, the EDPIC-2D code will save, if requested, all ions collided with a given boundary object into a set of files.
! Each file contains particles collided with the boundary object during time interval between two consecutive snapshots.
!
! For example, file _0003_ions_collided_with_bo_01.bin, contains parameters of all ion particles that collided with 
! boundary object #1 between times when snapshots 2 and 3 were created.
!
! The binary file with particles has all the information necessary to run this program, 
! no additional input files from the original simulation are necessary.
!
! General parameters required for processing are saved in the beginning of the binary file. The particles are saved in blocks called records.
! Each record begins with the time of the record and the number of particles in the record, then follows the particle data.
! For each particle there are 5 numbers: 
! a token which contains the ion species index, the number of segment of boundary object that the particle collided with, and the tag of the particle
! coordinate of collision, which is either x or y of the point of intersection with the horizontal or vertical segment of the boundary object, respectively
! three velocity vector components.
!
! When the program runs for the very first time, it asks the user to specify a number of parameters. 
! The input of the user is saved into a control text file dtp_ion_wall_coll.dat. The user may edit this file before running the program next time.
! if the program does not find the control file, it will ask to enter the required parameters manually again.
! 
! There is an option to create a file where particle flux integrated over the surface of the object is saved as a function of time.
! This is a convenient opportunity to check whether the particle fluxes calculated by this program are consistent with other diagnostics,
! specifically, with the averaged fluxes created by program dataproc_extract_bo_particle_fluxes_from_history_new.f90
! The best way to use this option is when the code created a sequence of snapshots over some time interval, 
! and this program processes these snapshots all at once.
!
program get_energy_spectrum_of_ions_collided_with_walls

  implicit none

  REAL(8), PARAMETER :: e_Cl     = 1.602189d-19      ! Charge of single electron [Cl]
  REAL(8), PARAMETER :: amu_kg   = 1.660565d-27      ! atomic mass unit [kg]

  type cell_data
     real x_m
     real y_m
     real, allocatable :: particle_number_flux(:)
     real, allocatable :: particle_energy_flux(:)
     integer, allocatable :: N_of_parts_in_energy_box(:,:)
   end type cell_data

  type bo_segment
     integer istart
     integer jstart
     integer iend
     integer jend

     integer start_indx
     integer end_indx
     type (cell_data), allocatable :: spatial_box(:)
  end type bo_segment

  logical exists
  character(1) buf

  integer first_snap, last_snap
  integer nbo

  integer print_record_details_flag
  integer print_each_particle_flag
  integer save_avg_currents_vs_time_flag

  real w_min_eV, w_max_eV, dw_eV
  integer N_of_energy_boxes
!  real, allocatable :: energy_box_middle_eV(:)

  character(47) Javg_filename ! avg_current_of_ions_collided_with_bo_NN_vst.dat
!                               ----*----I----*----I----*----I----*----I----*--

  integer snap
!                                   ! ----x----I----x----I----x----I----
  character(34) filename_ioncoll    ! _NNNN_ions_collided_with_bo_NN.bin

  integer N_of_segments
  type (bo_segment), allocatable :: segment(:)
  integer alloc_err

  integer n

  integer N_spec
  real t_begin_ns, t_end_ns, t_record_ns
  real delta_x_m
  real N_scale_part_m3
  real V_scale_ms

  real(8) energy_factor_eV
  real, allocatable :: M_i_amu(:)
  integer, allocatable ::N_part_of_spec(:)   ! number of particles of each species
  real, allocatable :: avgJ_A_of_ion_spec(:) ! average current of ion species (average flux integrated over the surface, in Amperes)
                                             ! calculated assuming that the ion charge is +1e
  integer s

  integer rec_count     ! counter of different time moments (records)
  integer n_parts
  integer skip_count    ! counter of particles with energies beyond limits assignet to teh energy distribution function

  real minimal_ion_energy_eV
  real maximal_ion_energy_eV

  integer ios
  real rdummy
  real r_N_parts_in_record
  integer N_parts_in_record

  integer k

  real(8) number_flux_factor_mm2ns

  real rtoken           ! particle values read from file
  real coll_coord       !
  real vx, vy, vz       !

  integer token
  integer tag    ! particle tag, not used here, may be necessary in future

  integer isb  ! index of spatial box
  real w_eV
  integer iwb  ! index of energy box

!                                           ! ----*----I----*----I----*----I----*
  character(35) filename_segionfluxes       ! _NNNN_bo_NN_seg_NN_ion_S_fluxes.asc
!                                           ! ----*----I----*----I----*----I----*----I----
  character(44) filename_segionspectr       ! _NNNN_bo_NN_seg_NN_ion_S_energy_spectrum.asc

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE
 
  inquire (file = 'dtp_ion_wall_coll.dat', exist = exists)

  if (exists) then

     open (9,file = 'dtp_ion_wall_coll.dat')
     read (9, '(A1)') buf  ! # the first snapshot number
     read (9, *) first_snap
     read (9, '(A1)') buf  ! # the last snapshot number
     read (9, *) last_snap
     read (9, '(A1)') buf  ! # the number of the boundary object
     read (9, *) nbo
     read (9, '(A1)') buf  ! # do you want to print time and number of particles in each record? (1/0=Yes/No)
     read (9, *) print_record_details_flag
     read (9, '(A1)') buf  ! # do you want to print each particle? (1/0=Yes/No)
     read (9, *) print_each_particle_flag
     read (9, '(A1)') buf  ! # do you want to save average total particle flux integrated over object surface as a function of time? (1/0=Yes/No)
     read (9, *) save_avg_currents_vs_time_flag
     read (9, '(A1)') buf  ! # the minimal energy for ion energy spectrum (eV)
     read (9, *) w_min_eV
     read (9, '(A1)') buf  ! # the maximal energy for ion energy spectrum (eV)
     read (9, *) w_max_eV
     read (9, '(A1)') buf  ! # the energy bin size for ion energy spectrum (eV)
     read (9, *) dw_eV
    close (9, status = 'keep')

  else

     print '("file dtp_ion_wall_coll.dat not found, do manual input...")'
     print '("### inputs will be saved into file dtp_ion_wall_coll.dat, edit this file as necessary for future use or remove ###")'

     print *, "enter the first snapshot number"
     read *, first_snap

     print *, "enter the last snapshot number"
     read *, last_snap

     print *, "enter the number of the boundary object"
     read *, nbo

     print *, "do you want to print time and number of particles in each record? (1/0=Yes/No)"
     read *, print_record_details_flag

     print *, "do you want to print each particle? (1/0=Yes/No)"
     read *, print_each_particle_flag

     print *, "do you want to save average total particle flux integrated over object surface as a function of time? (1/0=Yes/No)"
     read *, save_avg_currents_vs_time_flag

     print *, "enter the minimal energy for ion energy spectrum (eV)"
     read *, w_min_eV

     print *, "enter the maximal energy for ion energy spectrum (eV)"
     read *, w_max_eV

     print *, "enter the energy bin size for ion energy spectrum (eV)"
     read *, dw_eV

     open (9, file = 'dtp_ion_wall_coll.dat')
     write (9, '("# the first snapshot number")')
     write (9, '(2x,i4)') first_snap 
     write (9, '("# the last snapshot number")')
     write (9, '(2x,i4)') last_snap
     write (9, '("# the number of the boundary object")')
     write (9, '(2x,i4)') nbo
     write (9, '("# do you want to print time and number of particles in each record? (1/0=Yes/No)")')
     write (9, '(2x,i4)') print_record_details_flag
     write (9, '("# do you want to print each particle? (1/0=Yes/No)")')
     write (9, '(2x,i4)') print_each_particle_flag
     write (9, '("# do you want to save average total particle flux integrated over object surface as a function of time? (1/0=Yes/No)")')
     write (9, '(2x,i4)') save_avg_currents_vs_time_flag
     write (9, '("# the minimal energy for ion energy spectrum (eV)")')
     write (9, '(2x,f10.3)') w_min_eV
     write (9, '("# the maximal energy for ion energy spectrum (eV)")')
     write (9, '(2x,f10.3)') w_max_eV
     write (9, '("# the energy bin size for ion energy spectrum (eV)")')
     write (9, '(2x,f10.3)') dw_eV
     close (9, status = 'keep')

     print '("created file dtp_ion_wall_coll.dat")'

  end if

  N_of_energy_boxes = (w_max_eV - w_min_eV) / dw_eV
  w_max_eV = w_min_eV + N_of_energy_boxes * dw_eV

  print '("Number of energy boxes ",i6)', N_of_energy_boxes
  print '("Adjusted maximal energy ",f10.3," eV")', w_max_eV

!  allocate(energy_box_middle_eV(1:N_of_energy_boxes), stat = alloc_err)
!  do iw = 1, N_of_energy_boxes
!     energy_box_middle_eV(iw) = w_min_eV + (real(iw) - 0.5) * dw_eV
!  end do

  if (save_avg_currents_vs_time_flag.gt.0) then
     Javg_filename = 'avg_current_of_ions_collided_with_bo_NN_vst.dat'
!                     ----*----I----*----I----*----I----*----I----*--
     Javg_filename(38:39) = convert_int_to_txt_string(nbo, 2)
     open (10, file = Javg_filename)
  end if

  do snap = first_snap, last_snap

     filename_ioncoll = '_NNNN_ions_collided_with_bo_NN.bin'
     filename_ioncoll(2:5) = convert_int_to_txt_string(snap, 4)
     filename_ioncoll(29:30) = convert_int_to_txt_string(nbo, 2)
     inquire (file = filename_ioncoll, exist = exists)

     if (.not.exists) then
        print '("file ",A34," not found, skip")', filename_ioncoll
        cycle
     end if

     open (19, file = filename_ioncoll, access = 'stream', form = 'unformatted')
     print '("processing file ",A34)', filename_ioncoll

     read (19) N_of_segments
     print '("N_of_segments = ",i3)', N_of_segments
     if (allocated(segment)) deallocate(segment, stat = alloc_err)
     allocate(segment(1:N_of_segments), stat = alloc_err)

     do n = 1, N_of_segments
        read (19) segment(n)%istart
        read (19) segment(n)%jstart
        read (19) segment(n)%iend
        read (19) segment(n)%jend
        print '("segment ",i2," istart/jstart/iend/jend :: ",4(2x,i4))', n, segment(n)%istart, segment(n)%jstart, segment(n)%iend, segment(n)%jend
     end do

     read (19) N_spec
     print '("N_spec = ",i2)', N_spec

     read (19) t_begin_ns
     print '("start time = ",f12.5," ns")', t_begin_ns

     read (19) delta_x_m
     print '("mesh size = ",f10.6," mm")', delta_x_m * 1.0e3

     read (19) N_scale_part_m3
     print '("density of a single macroparticle ",e12.5," m-3")', N_scale_part_m3

     read (19) V_scale_ms
     print '("scale velocity ",f12.2," m/s")', V_scale_ms

     energy_factor_eV = 0.5_8 * amu_kg * V_scale_ms * V_scale_ms / e_Cl

     allocate(M_i_amu(1:N_spec), stat=alloc_err)
     do s = 1, N_spec
        read (19) M_i_amu(s)
        print '("ion species ",i2," mass ",f7.3," amu")', s, M_i_amu(s)
     end do
     allocate(N_part_of_spec(1:N_spec), stat=alloc_err)
     allocate(avgJ_A_of_ion_spec(1:N_spec), stat=alloc_err)

     do n = 1, N_of_segments
        if (segment(n)%istart.eq.segment(n)%iend) then
! vertical segment
           segment(n)%start_indx = min(segment(n)%jend, segment(n)%jstart)
           segment(n)%end_indx   = max(segment(n)%jend, segment(n)%jstart)-1
           allocate(segment(n)%spatial_box(segment(n)%start_indx:segment(n)%end_indx), stat = alloc_err)

           do isb = segment(n)%start_indx, segment(n)%end_indx
              segment(n)%spatial_box(isb)%x_m = segment(n)%istart * delta_x_m
              segment(n)%spatial_box(isb)%y_m = (real(isb) + 0.5) * delta_x_m

              allocate(segment(n)%spatial_box(isb)%particle_number_flux(1:N_spec), stat = alloc_err)
              allocate(segment(n)%spatial_box(isb)%particle_energy_flux(1:N_spec), stat = alloc_err)
              allocate(segment(n)%spatial_box(isb)%N_of_parts_in_energy_box(1:N_of_energy_boxes,1:N_spec), stat = alloc_err)
           end do

        else if (segment(n)%jstart.eq.segment(n)%jend) then
! horizontal segment
           segment(n)%start_indx = min(segment(n)%iend, segment(n)%istart)
           segment(n)%end_indx   = max(segment(n)%iend, segment(n)%istart)-1
           allocate(segment(n)%spatial_box(segment(n)%start_indx:segment(n)%end_indx), stat = alloc_err)

           do isb = segment(n)%start_indx, segment(n)%end_indx
              segment(n)%spatial_box(isb)%y_m = segment(n)%jstart * delta_x_m
              segment(n)%spatial_box(isb)%x_m = (real(isb) + 0.5) * delta_x_m

              allocate(segment(n)%spatial_box(isb)%particle_number_flux(1:N_spec), stat = alloc_err)
              allocate(segment(n)%spatial_box(isb)%particle_energy_flux(1:N_spec), stat = alloc_err)
              allocate(segment(n)%spatial_box(isb)%N_of_parts_in_energy_box(1:N_of_energy_boxes,1:N_spec), stat = alloc_err)
           end do
        else
! error
           print '("error-1")'
           stop
        end if
     end do   !###   do n = 1, N_of_segments

! clear all arrays
     do n = 1, N_of_segments
        do isb = segment(n)%start_indx, segment(n)%end_indx
           segment(n)%spatial_box(isb)%particle_number_flux = 0.0
           segment(n)%spatial_box(isb)%particle_energy_flux = 0.0
           segment(n)%spatial_box(isb)%N_of_parts_in_energy_box = 0
        end do
     end do

! read particle data

     rec_count = 0
     n_parts = 0
     skip_count = 0
     N_part_of_spec = 0
     minimal_ion_energy_eV =  1.0e7  ! initialize with values that make no sense
     maximal_ion_energy_eV = -1.0e7  !

     do  ! cycle over records, a record includes record time, number of particles, and all particles saved at this time
        
        read (19, iostat = ios) rdummy
        if (ios.ne.0) then 
! EOF
           exit
        end if

        read (19, iostat = ios) r_N_parts_in_record
        if (ios.ne.0) then
! this is an error
           print '("error-2")'
           stop
        end if

        rec_count = rec_count + 1
        t_record_ns = rdummy
        N_parts_in_record = r_N_parts_in_record
        n_parts = n_parts + N_parts_in_record

        if (print_record_details_flag.gt.0) print '(" record ",i6," time ",f12.5," ns ",i6," ion particles")', rec_count, t_record_ns, N_parts_in_record

        do k = 1, N_parts_in_record

           read (19) rtoken
           read (19) coll_coord
           read (19) vx
           read (19) vy
           read (19) vz

           token = rtoken

           s = token / 10000
           N_part_of_spec(s) = N_part_of_spec(s) + 1
           n = (token - s * 10000) / 100
           tag = token - s * 10000 - n * 100

           isb = max(segment(n)%start_indx,min(int(coll_coord),segment(n)%end_indx))

           segment(n)%spatial_box(isb)%particle_number_flux(s) = segment(n)%spatial_box(isb)%particle_number_flux(s) + 1.0

           w_eV = M_i_amu(s) * (vx * vx + vy * vy + vz * vz)  ! here it is not in eV yet

           segment(n)%spatial_box(isb)%particle_energy_flux(s) = segment(n)%spatial_box(isb)%particle_energy_flux(s) + w_eV

           w_eV = w_eV * energy_factor_eV                     ! now it is in eV

           if (w_eV.lt.minimal_ion_energy_eV) minimal_ion_energy_eV = w_eV
           if (w_eV.gt.maximal_ion_energy_eV) maximal_ion_energy_eV = w_eV

           if ((w_eV.ge.w_min_eV).and.(w_eV.le.w_max_eV)) then
              iwb = 1 + (w_eV - w_min_eV) / dw_eV
              iwb=max(1, min(iwb, N_of_energy_boxes))
              segment(n)%spatial_box(isb)%N_of_parts_in_energy_box(iwb,s) = segment(n)%spatial_box(isb)%N_of_parts_in_energy_box(iwb,s) + 1
              if (print_each_particle_flag.ne.0) then
                 print '("part ",i6," token ",i6," spec ",i1," seg ",i2," tag ",i3," coll.coord ",f10.4," coll.indx ",i5," vx/vy/vz ",e12.5," ",e12.5," ",e12.5,"   W ",f10.3," eV, energy box ",i6)', &
                      & k, token, s, n, tag, coll_coord, isb, vx, vy, vz, w_eV, iwb
              end if
           else          
              skip_count = skip_count + 1
              if (print_each_particle_flag.ne.0) then
                 print '("part ",i6," token ",i6," spec ",i1," seg ",i2," tag ",i3," coll.coord ",f10.4," coll.indx ",i5," vx/vy/vz ",e12.5," ",e12.5," ",e12.5,"   W ",f10.3," eV")', &
                      & k, token, s, n, tag, coll_coord, isb, vx, vy, vz, w_eV
              end if
           end if

        end do   !###   do k = 1, N_parts_in_record
     end do   !### do

     print '("actual energy range from ",f10.3," eV to ",f10.3," eV")', minimal_ion_energy_eV, maximal_ion_energy_eV
     print '("total number of saved particles ",i6)', n_parts

     if (skip_count.gt.0) print '("### warning, there were ",i7," particles with energies outside the specified spectrum energy range ###")', skip_count

!     read (19) t_end_ns
!     print '("end time = ",f12.5," ns")', t_end_ns
     t_end_ns = t_record_ns

     if (t_end_ns.le.t_begin_ns) then
        print '("error, the end time must be bigger than the beginning time :: end/begin :: ",f13.5," ns / ",f13.5," ns")', t_end_ns, t_begin_ns
        stop
     end if

     number_flux_factor_mm2ns = 1.0d-6 * N_scale_part_m3 * delta_x_m / (t_end_ns - t_begin_ns)

     close (19, status = 'keep')   !############# pick up here
     print '("finished reading file ",A34)', filename_ioncoll

     do n = 1, N_of_segments
        do s = 1, N_spec

!---------------------- particle number/energy fluxes to the wall ---------------------

           filename_segionfluxes = '_NNNN_bo_NN_seg_NN_ion_S_fluxes.asc'
                                 !  ----*----I----*----I----*----I----*
           filename_segionfluxes(2:5)   = convert_int_to_txt_string(snap, 4)
           filename_segionfluxes(10:11) = convert_int_to_txt_string(nbo, 2)
           filename_segionfluxes(17:18) = convert_int_to_txt_string(n, 2)
           filename_segionfluxes(24:24) = convert_int_to_txt_string(s, 1)
        
           open (9, file = filename_segionfluxes)

           write (9, '("# beginning time ",f12.5," ns, end time ",f12.5," ns, number of macroparticles ",i6)') t_begin_ns, t_end_ns, N_part_of_spec(s)
           write (9, '("# column 1 is x-index of the middle of spatial box on segment surface [dim-less]")')
           write (9, '("# column 2 is y-index of the middle of spatial box on segment surface [dim-less]")')
           write (9, '("# column 3 is x-coordinate of the middle [mm]")')
           write (9, '("# column 4 is y-coordinate of the middle [mm]")')
           write (9, '("# column 5 is the particle number flux to the wall [1 / mm^2 ns]")')
           write (9, '("# column 6 is the particle energy flux to the wall [eV / mm^2 ns]")')

           if (segment(n)%istart.eq.segment(n)%iend) then
! vertical segment
              do isb = segment(n)%start_indx, segment(n)%end_indx
                 write (9, '(2x,f8.1,2x,f8.1,2x,f12.7,2x,f12.7,2x,e12.5,2x,e12.5)') &
                      & real(segment(n)%istart), &          ! 1
                      & real(isb) + 0.5, &                  ! 2
                      & segment(n)%spatial_box(isb)%x_m * 1.0e3, &          ! 3
                      & segment(n)%spatial_box(isb)%y_m * 1.0e3, &          ! 4
                      & segment(n)%spatial_box(isb)%particle_number_flux(s) * number_flux_factor_mm2ns, &                     ! 5
                      & segment(n)%spatial_box(isb)%particle_energy_flux(s) * number_flux_factor_mm2ns * energy_factor_eV     ! 6
              end do
           else
! horizontal segment
              do isb = segment(n)%start_indx, segment(n)%end_indx
                 write (9, '(2x,f8.1,2x,f8.1,2x,f12.7,2x,f12.7,2x,e12.5,2x,e12.5)') &
                      & real(isb) + 0.5, &                  ! 1
                      & real(segment(n)%jstart), &          ! 2
                      & segment(n)%spatial_box(isb)%x_m * 1.0e3, &          ! 3
                      & segment(n)%spatial_box(isb)%y_m * 1.0e3, &          ! 4
                      & segment(n)%spatial_box(isb)%particle_number_flux(s) * number_flux_factor_mm2ns, &                     ! 5
                      & segment(n)%spatial_box(isb)%particle_energy_flux(s) * number_flux_factor_mm2ns * energy_factor_eV     ! 6
              end do
           end if

           close (9, status = 'keep')

           print '("created file ",A35)', filename_segionfluxes

!---------------------- energy spectrum of particles collided with the wall ---------------------

           filename_segionspectr = '_NNNN_bo_NN_seg_NN_ion_S_energy_spectrum.asc'
                                 !  ----*----I----*----I----*----I----*----I----
           filename_segionspectr(2:5)   = convert_int_to_txt_string(snap, 4)
           filename_segionspectr(10:11) = convert_int_to_txt_string(nbo, 2)
           filename_segionspectr(17:18) = convert_int_to_txt_string(n, 2)
           filename_segionspectr(24:24) = convert_int_to_txt_string(s, 1)
        
           open (9, file = filename_segionspectr)

           write (9, '("# beginning time ",f12.5," ns, end time ",f12.5," ns, number of macroparticles ",i6)') t_begin_ns, t_end_ns, N_part_of_spec(s)
           write (9, '("# column 1 is x-index of the middle of spatial box on segment surface [dim-less, may be half-integer]")')
           write (9, '("# column 2 is y-index of the middle of spatial box on segment surface [dim-less, may be half-integer]")')
           write (9, '("# column 3 is x-coordinate of the middle [mm]")')
           write (9, '("# column 4 is y-coordinate of the middle [mm]")')
           write (9, '("# column 5 is the index of the energy bin [dim-less]")')
           write (9, '("# column 6 is the energy value of the middle of energy bin [eV]")')
           write (9, '("# column 7 is the number of macroparticles in the bin")')

           if (segment(n)%istart.eq.segment(n)%iend) then
! vertical segment
              do isb = segment(n)%start_indx, segment(n)%end_indx
                 do iwb = 1, N_of_energy_boxes
                    write (9, '(2x,f8.1,2x,f8.1,2x,f12.7,2x,f12.7,2x,i6,2x,f10.4,2x,i6)') &
                         & real(segment(n)%istart), &          ! 1
                         & real(isb) + 0.5, &                  ! 2
                         & segment(n)%spatial_box(isb)%x_m * 1.0e3, &          ! 3
                         & segment(n)%spatial_box(isb)%y_m * 1.0e3, &          ! 4
                         & iwb, &                                                     ! 5
                         & w_min_eV + (real(iwb) - 0.5) * dw_eV, &                           ! 6
                         & segment(n)%spatial_box(isb)%N_of_parts_in_energy_box(iwb,s)       ! 7
                 end do
                 write (9, '(" ")')
              end do
           else
! horizontal segment
              do isb = segment(n)%start_indx, segment(n)%end_indx
                 do iwb = 1, N_of_energy_boxes
                    write (9, '(2x,f8.1,2x,f8.1,2x,f12.7,2x,f12.7,2x,i6,2x,f10.4,2x,i6)') &
                         & real(isb) + 0.5, &                  ! 1
                         & real(segment(n)%jstart), &          ! 2
                         & segment(n)%spatial_box(isb)%x_m * 1.0e3, &          ! 3
                         & segment(n)%spatial_box(isb)%y_m * 1.0e3, &          ! 4
                         & iwb, &                                                     ! 5
                         & w_min_eV + (real(iwb) - 0.5) * dw_eV, &                           ! 6
                         & segment(n)%spatial_box(isb)%N_of_parts_in_energy_box(iwb,s)       ! 7
                 end do
                 write (9, '(" ")')
              end do
           end if

           close (9, status = 'keep')

           print '("created file ",A44)', filename_segionspectr

        end do   !###   do s = 1, N_spec
     end do   !###   do n = 1, N_of_segments

! integrate number fluxes over the surface of the object in order to get average electric current component
! which can be compared with history_bo_NN_avg.dat

     if (save_avg_currents_vs_time_flag.gt.0) then
        avgJ_A_of_ion_spec = 0.0
        do s = 1, N_spec
           do n = 1, N_of_segments
              do isb = segment(n)%start_indx, segment(n)%end_indx
                 avgJ_A_of_ion_spec(s) = avgJ_A_of_ion_spec(s) + segment(n)%spatial_box(isb)%particle_number_flux(s)
              end do
           end do
        end do
        avgJ_A_of_ion_spec = avgJ_A_of_ion_spec * (1.0d9 * e_Cl * N_scale_part_m3 * delta_x_m * delta_x_m / (t_end_ns - t_begin_ns))   ! gives current in Amperes for single charged ions
        write (10, '(2x,f14.6,2x,f14.6,10(2x,e12.5))') t_begin_ns, t_end_ns, avgJ_A_of_ion_spec
     end if

! cleanup
     if (allocated(M_i_amu)) deallocate(M_i_amu, stat = alloc_err)
     if (allocated(N_part_of_spec)) deallocate(N_part_of_spec, stat = alloc_err)
     if (allocated(avgJ_A_of_ion_spec)) deallocate(avgJ_A_of_ion_spec, stat = alloc_err)
     do n = 1, N_of_segments
        do isb = segment(n)%start_indx, segment(n)%end_indx
           if (allocated(segment(n)%spatial_box(isb)%particle_number_flux))     deallocate(segment(n)%spatial_box(isb)%particle_number_flux, stat = alloc_err)
           if (allocated(segment(n)%spatial_box(isb)%particle_energy_flux))     deallocate(segment(n)%spatial_box(isb)%particle_energy_flux, stat = alloc_err)
           if (allocated(segment(n)%spatial_box(isb)%N_of_parts_in_energy_box)) deallocate(segment(n)%spatial_box(isb)%N_of_parts_in_energy_box, stat = alloc_err)
        end do
        if (allocated(segment(n)%spatial_box)) deallocate(segment(n)%spatial_box, stat = alloc_err)
     end do
     deallocate(segment, stat = alloc_err)

  end do   !###   do snap = first_snap, last_snap

  print '("done...")'

  if (save_avg_currents_vs_time_flag.gt.0) then
     close (10, status = 'keep')
     print '("file ",A47," is ready")', Javg_filename
  end if

end program get_energy_spectrum_of_ions_collided_with_walls



!-----------------------------------
! creates a string of length "length_of_string" out of an integer number "int_number"
!
function convert_int_to_txt_string(int_number, length_of_string)

  implicit none

  integer int_number

  integer length_of_string
  character*(length_of_string) convert_int_to_txt_string

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
