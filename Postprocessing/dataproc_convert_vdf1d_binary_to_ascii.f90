!------------------------------------------
! This program reads binary snapshot files 
! _NNNN_vdf1d.bin, where NNNN is the snapshot number like 0001, 0002, etc.,
! containing multiple 1d electron velocity distributions functions (EVDF) fe(vx), fe(vy), fe(vz)
! and ion velocity distribution functions (IVDF) fis(vx), fis(vy), fis(vz) for each ion species s
! created using particles from different locations (boxes).
! The program saves one ASCII file which contains EVDF fe(vx) from all boxes, _NNNN_evxdf.asc,
!                   one ASCII file which contains EVDF fe(vy) from all boxes, _NNNN_evydf.asc,
!                   one ASCII file which contains EVDF fe(vz) from all boxes, _NNNN_evzdf.asc,
!
! and then for each ion species s, one ASCII file which contains IVDF fis(vx) from all boxes, _NNNN_i_S_vxdf.asc,
!                                  one ASCII file which contains IVDF fis(vy) from all boxes, _NNNN_i_S_vydf.asc,
!                                  one ASCII file which contains IVDF fis(vz) from all boxes, _NNNN_i_S_vzdf.asc,
! here NNNN is snapshot number and S is the ion species number.
! Also, the boxes are saved into files with names _NNNN_vdf_boxes.asc.
! The ASCII files can be plotted with gnuplot, see the header in these files for more details.
!
program convert_vdf1d_bin_to_ascii

  implicit none

  real, parameter :: e_Cl     = 1.602189e-19
  real, parameter :: m_e_kg   = 9.109534e-31
  real, parameter :: amu_kg   = 1.660565e-27

  integer first_snap, last_snap, snap
  character(15) filename                 ! _NNNN_vdf1d.bin
  logical exists

  integer N_boxes_x
  integer N_boxes_y
  integer indx_v_min_e
  integer indx_v_max_e
  integer indx_v_min_i
  integer indx_v_max_i
  integer N_vbins_e
  integer N_vbins_i
  integer N_spec

  integer total_N_of_boxes

  real T_e_eV
  real, allocatable :: M_i_amu(:)
  integer alloc_err

  real, allocatable :: rbufer(:)

  integer, allocatable :: evxdf(:,:)
  integer, allocatable :: evydf(:,:)
  integer, allocatable :: evzdf(:,:)

  integer, allocatable :: isvxdf(:,:,:)
  integer, allocatable :: isvydf(:,:,:)
  integer, allocatable :: isvzdf(:,:,:)

  integer s, i, j, nbox, k

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
  character(15) filename_evdf             ! _NNNN_ev*df.asc      , where *=xyz
  character(18) filename_ivdf             ! _NNNN_i_S_v*df.asc

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
     filename = '_NNNN_vdf1d.bin'
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
     read (9) indx_v_min_i
     read (9) indx_v_max_i
     read (9) N_vbins_e
     read (9) N_vbins_i
     read (9) N_spec

     total_N_of_boxes = N_boxes_x * N_boxes_y

     allocate(M_i_amu(1:N_spec), stat = alloc_err)

     allocate(evxdf(indx_v_min_e:indx_v_max_e, 1:total_N_of_boxes), stat = alloc_err)
     allocate(evydf(indx_v_min_e:indx_v_max_e, 1:total_N_of_boxes), stat = alloc_err)
     allocate(evzdf(indx_v_min_e:indx_v_max_e, 1:total_N_of_boxes), stat = alloc_err)

     allocate(isvxdf(indx_v_min_i:indx_v_max_i, 1:total_N_of_boxes, 1:N_spec), stat = alloc_err)
     allocate(isvydf(indx_v_min_i:indx_v_max_i, 1:total_N_of_boxes, 1:N_spec), stat = alloc_err)
     allocate(isvzdf(indx_v_min_i:indx_v_max_i, 1:total_N_of_boxes, 1:N_spec), stat = alloc_err)
     
     allocate(box(1:total_N_of_boxes), stat = alloc_err)

     read (9) T_e_eV

     do s = 1, N_spec
        read (9) M_i_amu(s)
     end do

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

           do k = indx_v_min_e, indx_v_max_e
              read (9) evxdf(k, nbox)
           end do

           do k = indx_v_min_e, indx_v_max_e
              read (9) evydf(k, nbox)
           end do

           do k = indx_v_min_e, indx_v_max_e
              read (9) evzdf(k, nbox)
           end do

           do s = 1, N_spec
              do k = indx_v_min_i, indx_v_max_i
                 read (9) isvxdf(k, nbox, s)
              end do

              do k = indx_v_min_i, indx_v_max_i
                 read (9) isvydf(k, nbox, s)
              end do

              do k = indx_v_min_i, indx_v_max_i
                 read (9) isvzdf(k, nbox, s)
              end do
           end do  ! end of cycle over ion species

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

! EVDF(VX)=============================================================

     filename_evdf = '_NNNN_evxdf.asc'
     filename_evdf(2:5) = convert_int_to_txt_string(snap, 4)
     open (99, file = filename_evdf)
! save column description
     write (99, '("# col    1 is the X-velocity [middle of the bin, units of v_th_e for Te = ",f8.3," eV which is ",f10.1," m/s]")') T_e_eV, SQRT(2.0 * e_Cl * T_e_eV / m_e_kg)
     write (99, '("# the columns listed below contain [dim-less] number of macroparticles in the velocity bin calculated in the following box :")')
     do nbox = 1, total_N_of_boxes
        write (99, '("# col ",i4," is for box ",i4," x-index ",i4," y-index ",i4," | xmin ",f9.6," [m] xmax ",f9.6," [m] | ymin ",f9.6," ymax ",f9.6," [m]")') &
             & nbox+1, nbox, box(nbox)%i, box(nbox)%j, box(nbox)%xmin_m, box(nbox)%xmax_m, box(nbox)%ymin_m, box(nbox)%ymax_m
     end do
! save the evdf and the velocity
     do k = indx_v_min_e, indx_v_max_e
        write (99, '(2x,f8.4,1024(2x,i9))') (real(k)+0.5) / real(N_vbins_e), evxdf(k, 1:total_N_of_boxes)
     end do
     close (99, status = 'keep')
     print '("created file ",A15)', filename_evdf

! EVDF(VY)=============================================================

     filename_evdf = '_NNNN_evydf.asc'
     filename_evdf(2:5) = convert_int_to_txt_string(snap, 4)
     open (99, file = filename_evdf)
! save column description
     write (99, '("# col    1 is the Y-velocity [middle of the bin, units of v_th_e for Te = ",f8.3," eV which is ",f10.1," m/s]")') T_e_eV, SQRT(2.0 * e_Cl * T_e_eV / m_e_kg)
     write (99, '("# the columns listed below contain [dim-less] number of macroparticles in the velocity bin calculated in the following box :")')
     do nbox = 1, total_N_of_boxes
        write (99, '("# col ",i4," is for box ",i4," x-index ",i4," y-index ",i4," | xmin ",f9.6," [m] xmax ",f9.6," [m] | ymin ",f9.6," ymax ",f9.6," [m]")') &
             & nbox+1, nbox, box(nbox)%i, box(nbox)%j, box(nbox)%xmin_m, box(nbox)%xmax_m, box(nbox)%ymin_m, box(nbox)%ymax_m
     end do
! save the evdf and the velocity
     do k = indx_v_min_e, indx_v_max_e
        write (99, '(2x,f8.4,1024(2x,i9))') (real(k)+0.5) / real(N_vbins_e), evydf(k, 1:total_N_of_boxes)
     end do
     close (99, status = 'keep')
     print '("created file ",A15)', filename_evdf

! EVDF(VZ)=============================================================

     filename_evdf = '_NNNN_evzdf.asc'
     filename_evdf(2:5) = convert_int_to_txt_string(snap, 4)
     open (99, file = filename_evdf)
! save column description
     write (99, '("# col    1 is the Z-velocity [middle of the bin, units of v_th_e for Te = ",f8.3," eV which is ",f10.1," m/s]")') T_e_eV, SQRT(2.0 * e_Cl * T_e_eV / m_e_kg)
     write (99, '("# the columns listed below contain [dim-less] number of macroparticles in the velocity bin calculated in the following box :")')
     do nbox = 1, total_N_of_boxes
        write (99, '("# col ",i4," is for box ",i4," x-index ",i4," y-index ",i4," | xmin ",f9.6," [m] xmax ",f9.6," [m] | ymin ",f9.6," ymax ",f9.6," [m]")') &
             & nbox+1, nbox, box(nbox)%i, box(nbox)%j, box(nbox)%xmin_m, box(nbox)%xmax_m, box(nbox)%ymin_m, box(nbox)%ymax_m
     end do
! save the evdf and the velocity
     do k = indx_v_min_e, indx_v_max_e
        write (99, '(2x,f8.4,1024(2x,i9))') (real(k)+0.5) / real(N_vbins_e), evzdf(k, 1:total_N_of_boxes)
     end do
     close (99, status = 'keep')
     print '("created file ",A15)', filename_evdf

     do s = 1, N_spec

! ion VDF(VX)=============================================================

        filename_ivdf = '_NNNN_i_S_vxdf.asc'
        filename_ivdf(2:5) = convert_int_to_txt_string(snap, 4)
        filename_ivdf(9:9) = convert_int_to_txt_string(s, 1)
        open (99, file = filename_ivdf)
! save column description
        write (99, '("# col    1 is the X-velocity [middle of the bin, units of v_th_i for Ti = ",f8.3," eV and Ms = ",f6.2," (amu) which is ",f10.1," m/s]")') &
             & T_e_eV, M_i_amu(s), SQRT(2.0 * e_Cl * T_e_eV / (amu_kg * M_i_amu(s)))
        write (99, '("# the columns listed below contain [dim-less] number of macroparticles in the velocity bin calculated in the following box :")')
        do nbox = 1, total_N_of_boxes
           write (99, '("# col ",i4," is for box ",i4," x-index ",i4," y-index ",i4," | xmin ",f9.6," [m] xmax ",f9.6," [m] | ymin ",f9.6," ymax ",f9.6," [m]")') &
                & nbox+1, nbox, box(nbox)%i, box(nbox)%j, box(nbox)%xmin_m, box(nbox)%xmax_m, box(nbox)%ymin_m, box(nbox)%ymax_m
        end do
! save the evdf and the velocity
        do k = indx_v_min_i, indx_v_max_i
           write (99, '(2x,f8.4,1024(2x,i9))') (real(k)+0.5) / real(N_vbins_i), isvxdf(k, 1:total_N_of_boxes, s)
        end do
        close (99, status = 'keep')
        print '("created file ",A18)', filename_ivdf

! ion VDF(VY)=============================================================

        filename_ivdf = '_NNNN_i_S_vydf.asc'
        filename_ivdf(2:5) = convert_int_to_txt_string(snap, 4)
        filename_ivdf(9:9) = convert_int_to_txt_string(s, 1)
        open (99, file = filename_ivdf)
! save column description
        write (99, '("# col    1 is the Y-velocity [middle of the bin, units of v_th_i for Ti = ",f8.3," eV and Ms = ",f6.2," (amu) which is ",f10.1," m/s]")') &
             & T_e_eV, M_i_amu(s), SQRT(2.0 * e_Cl * T_e_eV / (amu_kg * M_i_amu(s)))
        write (99, '("# the columns listed below contain [dim-less] number of macroparticles in the velocity bin calculated in the following box :")')
        do nbox = 1, total_N_of_boxes
           write (99, '("# col ",i4," is for box ",i4," x-index ",i4," y-index ",i4," | xmin ",f9.6," [m] xmax ",f9.6," [m] | ymin ",f9.6," ymax ",f9.6," [m]")') &
                & nbox+1, nbox, box(nbox)%i, box(nbox)%j, box(nbox)%xmin_m, box(nbox)%xmax_m, box(nbox)%ymin_m, box(nbox)%ymax_m
        end do
! save the evdf and the velocity
        do k = indx_v_min_i, indx_v_max_i
           write (99, '(2x,f8.4,1024(2x,i9))') (real(k)+0.5) / real(N_vbins_i), isvydf(k, 1:total_N_of_boxes, s)
        end do
        close (99, status = 'keep')
        print '("created file ",A18)', filename_ivdf

! ion VDF(VZ)=============================================================

        filename_ivdf = '_NNNN_i_S_vzdf.asc'
        filename_ivdf(2:5) = convert_int_to_txt_string(snap, 4)
        filename_ivdf(9:9) = convert_int_to_txt_string(s, 1)
        open (99, file = filename_ivdf)
! save column description
        write (99, '("# col    1 is the Z-velocity [middle of the bin, units of v_th_i for Ti = ",f8.3," eV and Ms = ",f6.2," (amu) which is ",f10.1," m/s]")') &
             & T_e_eV, M_i_amu(s), SQRT(2.0 * e_Cl * T_e_eV / (amu_kg * M_i_amu(s)))
        write (99, '("# the columns listed below contain [dim-less] number of macroparticles in the velocity bin calculated in the following box :")')
        do nbox = 1, total_N_of_boxes
           write (99, '("# col ",i4," is for box ",i4," x-index ",i4," y-index ",i4," | xmin ",f9.6," [m] xmax ",f9.6," [m] | ymin ",f9.6," ymax ",f9.6," [m]")') &
                & nbox+1, nbox, box(nbox)%i, box(nbox)%j, box(nbox)%xmin_m, box(nbox)%xmax_m, box(nbox)%ymin_m, box(nbox)%ymax_m
        end do
! save the evdf and the velocity
        do k = indx_v_min_i, indx_v_max_i
           write (99, '(2x,f8.4,1024(2x,i9))') (real(k)+0.5) / real(N_vbins_i), isvzdf(k, 1:total_N_of_boxes, s)
        end do
        close (99, status = 'keep')
        print '("created file ",A18)', filename_ivdf

     end do ! end of cycle  over ion species

! cleanup
     if (allocated(M_i_amu)) deallocate(M_i_amu, stat = alloc_err)

     if (allocated(evxdf)) deallocate(evxdf, stat = alloc_err)
     if (allocated(evydf)) deallocate(evydf, stat = alloc_err)
     if (allocated(evzdf)) deallocate(evzdf, stat = alloc_err)

     if (allocated(isvxdf)) deallocate(isvxdf, stat = alloc_err)
     if (allocated(isvydf)) deallocate(isvydf, stat = alloc_err)
     if (allocated(isvzdf)) deallocate(isvzdf, stat = alloc_err)

     if (allocated(box)) deallocate(box, stat = alloc_err)

     if (allocated(rbufer)) deallocate(rbufer, stat = alloc_err)

  end do ! end of cycle over snapshots

end program convert_vdf1d_bin_to_ascii

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
