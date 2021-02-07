!------------------------------------------
! This program reads binary snapshot files 
! _NNNN_vdf2d.bin, where NNNN is the snapshot number like 0001, 0002, etc.,
! containing multiple 2d electron velocity distribution functions (EVDF) fe(vx,vy) 
! created using particles from different locations (boxes).
! For each location, the program saves one ASCII file the EVDF as function of vx and vy.
! EVDFS are saved into files with names _NNNN_LLLL_evxvydf.asc,
! boxes are saved into files with names _NNNN_vdf_boxes.asc,
! here NNNN is snapshot number and LLLL is the box number.
! The ASCII files can be plotted with gnuplot, see the header in these files for more details.
!
program convert_vdf2d_bin_to_ascii

  implicit none

  real, parameter :: e_Cl     = 1.602189e-19
  real, parameter :: m_e_kg   = 9.109534e-31
  real, parameter :: amu_kg   = 1.660565e-27

  integer first_snap, last_snap, snap
  character(15) filename                 ! _NNNN_vdf2d.bin
  logical exists

  integer N_boxes_x
  integer N_boxes_y
  integer indx_v_min_e
  integer indx_v_max_e
  integer N_vbins_e

  integer total_N_of_boxes

  real T_e_eV
  integer alloc_err

  real, allocatable :: rbufer(:)

  integer, allocatable :: evxvydf(:,:,:)

  integer i, j, nbox, iv, jv

  type vdfbox
     integer i
     integer j
     real xmin_m
     real xmax_m
     real ymin_m
     real ymax_m
  end type vdfbox

  type(vdfbox), allocatable :: box(:)

                                          ! ----x----I----x----
  character(19) filename_vdfbox           ! _NNNN_vdf_boxes.asc
  character(22) filename_evdf             ! _NNNN_LLLL_evxvydf.asc
                                          ! ----x----I----x----I--

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  print '("enter the number of the first snapshot")'
  read *, first_snap

  print '("enter the number of the last snapshot")'
  read *, last_snap

  do snap = first_snap, last_snap

! read data
     filename = '_NNNN_vdf2d.bin'
     filename(2:5) = convert_int_to_txt_string(snap, 4)

     inquire(file = filename, exist = exists)

     if (.not.exists) then
        print '("file ",A15," not found, skip")', filename
        cycle
     end if

     open (9, file = filename, access = 'stream', form = 'unformatted') 

     read (9) N_boxes_x
     read (9) N_boxes_y
     read (9) indx_v_min_e
     read (9) indx_v_max_e
     read (9) N_vbins_e

     total_N_of_boxes = N_boxes_x * N_boxes_y

     allocate(evxvydf(indx_v_min_e:indx_v_max_e, indx_v_min_e:indx_v_max_e, 1:total_N_of_boxes), stat = alloc_err)
     
     allocate(box(1:total_N_of_boxes), stat = alloc_err)

     read (9) T_e_eV

     nbox=0

     do j = 1, N_boxes_y
        do i = 1, N_boxes_x

           nbox = nbox+1

           box(nbox)%i = i
           box(nbox)%j = j

           read (9) box(nbox)%xmin_m
           read (9) box(nbox)%xmax_m
           read (9) box(nbox)%ymin_m          
           read (9) box(nbox)%ymax_m

           do jv = indx_v_min_e, indx_v_max_e
              do iv = indx_v_min_e, indx_v_max_e
                 read (9) evxvydf(iv, jv, nbox)
              end do
           end do

        end do   ! end of cycle over boxes in the horizontal (x) direction
     end do      ! end of cycle over boxes in the vertical (y) direction

     close (9, status = 'keep')

     print '("reading binary file ",A15," complete")', filename

! save data

! spatial boxes =============================================================

     filename_vdfbox = '_NNNN_vdf_boxes.asc'
     filename_vdfbox(2:5) = convert_int_to_txt_string(snap, 4)

     open (99, file = filename_vdfbox)
! save column description
     do nbox = 1, total_N_of_boxes
        write (99, '("# columns ",i4,2x,i4," are box ",i4," x-index ",i4," y-index ",i4," | xmin ",f9.6," [m] xmax ",f9.6," [m] | ymin ",f9.6," ymax ",f9.6," [m]")') &
             & 2*nbox-1, 2*nbox, nbox, box(nbox)%i, box(nbox)%j, box(nbox)%xmin_m, box(nbox)%xmax_m, box(nbox)%ymin_m, box(nbox)%ymax_m
     end do

     allocate(rbufer(1:2*total_N_of_boxes), stat = alloc_err)

! save left top corner
     do nbox = 1, total_N_of_boxes
        rbufer(2*nbox-1) = box(nbox)%xmin_m
        rbufer(2*nbox)   = box(nbox)%ymax_m
     end do
     write (99, '(1024(3x,f9.6,2x,f9.6))') rbufer

! save right top corner
     do nbox = 1, total_N_of_boxes
        rbufer(2*nbox-1) = box(nbox)%xmax_m
        rbufer(2*nbox)   = box(nbox)%ymax_m
     end do
     write (99, '(1024(3x,f9.6,2x,f9.6))') rbufer

! save right bottom corner
     do nbox = 1, total_N_of_boxes
        rbufer(2*nbox-1) = box(nbox)%xmax_m
        rbufer(2*nbox)   = box(nbox)%ymin_m
     end do
     write (99, '(1024(3x,f9.6,2x,f9.6))') rbufer

! save left bottom corner
     do nbox = 1, total_N_of_boxes
        rbufer(2*nbox-1) = box(nbox)%xmin_m
        rbufer(2*nbox)   = box(nbox)%ymin_m
     end do
     write (99, '(1024(3x,f9.6,2x,f9.6))') rbufer

! save left top corner
     do nbox = 1, total_N_of_boxes
        rbufer(2*nbox-1) = box(nbox)%xmin_m
        rbufer(2*nbox)   = box(nbox)%ymax_m
     end do
     write (99, '(1024(3x,f9.6,2x,f9.6))') rbufer

     close (99, status = 'keep')
     print '("created file ",A19)', filename_vdfbox

! EVDF(VX,VY)=============================================================

     do nbox = 1, total_N_of_boxes
        
        filename_evdf = '_NNNN_LLLL_evxvydf.asc'
        filename_evdf(2:5) = convert_int_to_txt_string(snap, 4)
        filename_evdf(7:10) = convert_int_to_txt_string(nbox, 4)

        open (99, file = filename_evdf)
! save column description
        write (99, '("# col 1 is the X-velocity [middle of the bin, units of v_th_e for Te = ",f8.3," eV which is ",f10.1," m/s]")') T_e_eV, SQRT(2.0 * e_Cl * T_e_eV / m_e_kg)
        write (99, '("# col 2 is the Y-velocity [as above]")')
        write (99, '("# col 3 contains [dim-less] number of macroparticles in the velocity bin calculated in box ",i4)') nbox
        write (99, '("# with x-index ",i4," y-index ",i4," | xmin ",f9.6," [m] xmax ",f9.6," [m] | ymin ",f9.6," ymax ",f9.6," [m]")') &
             & box(nbox)%i, box(nbox)%j, box(nbox)%xmin_m, box(nbox)%xmax_m, box(nbox)%ymin_m, box(nbox)%ymax_m
     
! save the evdf and the velocity
        do jv = indx_v_min_e, indx_v_max_e
           do iv = indx_v_min_e, indx_v_max_e
              write (99, '(2x,f8.4,2x,f8.4,2x,i9)') (real(iv)+0.5) / real(N_vbins_e), (real(jv)+0.5) / real(N_vbins_e), evxvydf(iv,jv,nbox)
           end do
           write (99, '(" ")')
        end do
        close (99, status = 'keep')
        print '("created file ",A22)', filename_evdf

     end do

! cleanup
     if (allocated(evxvydf)) deallocate(evxvydf, stat = alloc_err)
     if (allocated(box)) deallocate(box, stat = alloc_err)
     if (allocated(rbufer)) deallocate(rbufer, stat = alloc_err)

  end do ! end of cycle over snapshots

end program convert_vdf2d_bin_to_ascii

!-----------------------------------
! creates a string of length "length_of_string" out of an integer number "int_number"
!
function convert_int_to_txt_string(int_number, length_of_string)

  implicit none

  integer length_of_string
  character*(length_of_string) convert_int_to_txt_string
  integer int_number
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
