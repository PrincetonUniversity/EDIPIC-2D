!----------------------------------------------------
! This program reads files with electron particles data (X/Y/VX/VY/VZ/tag) in binary format
! and saves the particle data in ASCII format. 
!
! During a simulation run, the code creates (if requested) one binary file with electron particles, 
! the particles saved are inside one of several rectangular regions specified before the simulation.
!
! This program saves particles from each of these regions into a separate file.
!
! For example, if the simulation had 3 regions for phase planes in snapshot 5, the binary file would be _0005_epp.bin
! {here "epp" stands for "electron phase plane"}
! and this program will create 3 ASCII files: 
! _0005_epp_01.asc , _0005_epp_02.asc , _0005_epp_03.asc .
!   region 1 ^         region 2 ^         region 3 ^
!
! The binary file has all the information necessary to run this program, 
! no additional input files from the original simulation is required.
!
! The program will ask the user to specify the first and the last snapshot in the sequence to be processed.
!
! The program will also ask whether it is necessary to save the aforementioned rectangular regions as separate data files.
! These files may be quite useful, for example one may plot a 2d colormap of the electric field,
! and on top of it plot the region of the phase plane to show where the selected phase plane (set of particles) is situated.
! The most probable situation is that a sequence of snapshots uses same set of regions, so it is not necessary
! to save the boxes for each snapshot. One can process only the first snapshot in the sequence with saving the boxes,
! and then process the rest of the sequence without saving the boxes.
!
program convert_epp_bin_into_ascii

  implicit none

  integer first_snap, last_snap, snap
  integer save_pp_boxes_flag

!                                   ! ----x----I---
  character(13) filename_epp        ! _NNNN_epp.bin
  logical exists

  integer N_pp_boxes

  TYPE index_limits
     INTEGER imin
     INTEGER jmin
     INTEGER imax
     INTEGER jmax
  END TYPE index_limits

  type(index_limits), ALLOCATABLE :: pp_box(:)
  integer alloc_err

  integer i, j, k, n
  integer N_of_processes
  integer total, itemp

  real delta_x_m

  type particle
     real x
     real y
     real vxms
     real vyms
     real vzms
     integer tag
  end type particle

  type(particle), allocatable :: electron(:)
  real rtemp

!                                         ! ----x----I----x-
  character(16) filename_epp_ascii        ! _NNNN_epp_NN.asc
!                                         ! ----x----I----x----I
  character(20) filename_epp_box          ! _NNNN_epp_box_NN.asc

  integer count

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE
 
  print *, "enter the first snapshot number"
  read *, first_snap

  print *, "enter the last snapshot number"
  read *, last_snap

  print *, "do you need to save spatial boxes of the phase planes? (1/0 = Yes/No)"
  read *, save_pp_boxes_flag

  do snap = first_snap, last_snap

     filename_epp = '_NNNN_epp.bin'
     filename_epp(2:5) = convert_int_to_txt_string(snap, 4)
     inquire (file = filename_epp, exist = exists)

     if (.not.exists) then
        print '("file ",A13," not found, skip")', filename_epp
        cycle
     end if
  
     open (19, file = filename_epp, access = 'stream', form = 'unformatted')
     print '("processing file ",A13)', filename_epp

     read (19) N_pp_boxes
     print '("N_pp_boxes = ",i3)', N_pp_boxes
     if (allocated(pp_box)) deallocate(pp_box, stat = alloc_err)
     allocate(pp_box(1:N_pp_boxes), stat = alloc_err)

     do n = 1, N_pp_boxes
        read (19) pp_box(n)%imin
        read (19) pp_box(n)%jmin
        read (19) pp_box(n)%imax
        read (19) pp_box(n)%jmax
        print '("box ",i2," imin/jmin/imax/jmax :: ",4(2x,i4))', n, pp_box(n)%imin, pp_box(n)%jmin, pp_box(n)%imax, pp_box(n)%jmax
     end do

     read (19) N_of_processes
     print '("N_of_processes = ",i3)', N_of_processes

     total = 0
     do n = 0, N_of_processes-1
        read (19) itemp
        print '("proc ",i4," provided ",i7," particles")', n, itemp
        total = total + itemp
     end do
     print '("total number of particles = ",i8)', total

     read (19) delta_x_m
     print '("delta_x_m = ",f12.10," m")', delta_x_m

     allocate(electron(1:total), stat = alloc_err)

     do k = 1, total
        read (19) electron(k)%x
        read (19) electron(k)%y
        read (19) electron(k)%vxms
        read (19) electron(k)%vyms
        read (19) electron(k)%vzms
        read (19) rtemp
        electron(k)%tag = int(rtemp)
     end do

     close (19, status = 'keep')
     print '("finished reading from file ",A13)', filename_epp

     do n = 1, N_pp_boxes

        filename_epp_ascii = '_NNNN_epp_NN.asc'
        filename_epp_ascii(2:5)   = convert_int_to_txt_string(snap, 4)
        filename_epp_ascii(11:12) = convert_int_to_txt_string(n, 2)

        open (9, file = filename_epp_ascii)
        count = 0
        do k = 1, total
           i = INT(electron(k)%x)
           j = INT(electron(k)%y)
           IF ((i.GE.pp_box(n)%imin).AND.(i.LT.pp_box(n)%imax).AND.(j.GE.pp_box(n)%jmin).AND.(j.LT.pp_box(n)%jmax)) THEN
              count = count + 1
              write (9, '(5(2x,e14.7),2x,i3)') electron(k)%x * delta_x_m, electron(k)%y * delta_x_m, electron(k)%vxms, electron(k)%vyms, electron(k)%vzms, electron(k)%tag
           end IF
        end do
        close (9, status = 'keep')
        print '("created file ",A16," with ",i8," particles")', filename_epp_ascii, count

        if (save_pp_boxes_flag.ne.0) then
           filename_epp_box = '_NNNN_epp_box_NN.asc'
           filename_epp_box(2:5)   = convert_int_to_txt_string(snap, 4)
           filename_epp_box(15:16) = convert_int_to_txt_string(n, 2)
           open (9, file = filename_epp_box)
           write (9, '("# spatial box of phase plane ",i2)') n
           write (9, '("# column 1 is x-index of the grid node in the box vertice [dim-less]")')
           write (9, '("# column 2 is y-index of the grid node in the box vertice [dim-less]")')
           write (9, '("# column 3 is x-coordinate of the box vertice [m]")')
           write (9, '("# column 4 is y-coordinate of the box vertice [m]")')
! left-right along X, up-down along Y
           write (9, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') pp_box(n)%imin, pp_box(n)%jmin, pp_box(n)%imin * delta_x_m, pp_box(n)%jmin * delta_x_m  ! left bottom
           write (9, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') pp_box(n)%imin, pp_box(n)%jmax, pp_box(n)%imin * delta_x_m, pp_box(n)%jmax * delta_x_m  ! left top
           write (9, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') pp_box(n)%imax, pp_box(n)%jmax, pp_box(n)%imax * delta_x_m, pp_box(n)%jmax * delta_x_m  ! right top
           write (9, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') pp_box(n)%imax, pp_box(n)%jmin, pp_box(n)%imax * delta_x_m, pp_box(n)%jmin * delta_x_m  ! right top
           write (9, '(2x,i4,2x,i4,2x,f12.9,2x,f12.9)') pp_box(n)%imin, pp_box(n)%jmin, pp_box(n)%imin * delta_x_m, pp_box(n)%jmin * delta_x_m  ! left bottom
           close (9, status = 'keep')
           print '("created phase plane box file ",A20)', filename_epp_box
        end if

     end do

! cleanup
     if (allocated(pp_box)) deallocate(pp_box, stat = alloc_err)
     if (allocated(electron)) deallocate(electron, stat = alloc_err)

  end do

  print '("done...")'

end program convert_epp_bin_into_ascii



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
