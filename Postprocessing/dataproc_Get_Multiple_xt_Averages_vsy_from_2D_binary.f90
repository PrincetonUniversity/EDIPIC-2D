!----------------------------------------------------
! If the simulated domain is a rectangle and the system is periodic along the x direction
! but is not periodic along the y direction, it may be useful to have values of typical parameters 
! (potential, fields, densities, currents, temperatures, velocities, energies) averaged along the x direction
! and represented as functions of y.
! 
! This subroutine reads snapshot binary files with these parameters and averages them along the x direction. 
! It also combines the spatial averaging with the averaging over multiple snapshots, if requested.
!
! Note that in periodic systems, two last nodes (in the periodic direction) are copies of two first nodes (due to periodicity),
! which is why the two last nodes are omitted during the averaging.
!
! The result is saved into a single file with name avgxt_values_vsy_NNNN_NNNN.txt where NNNN_NNNN is replaced
! by numbers of the first and the last snapshots in the range of snapshot numbers used for averaging.
! For example, if you average over snapshots 100 to 200, the file name will be avgxt_values_vsy_0100_0200.txt .
! 
!
program create_multi_xt_avgs_vsy_from_2d_bin

  implicit none

  integer first_snap, last_snap

  integer probestatus

!                                 ! ----x----I----x-
  character(16) filename_F        ! _NNNN_F_V_2D.bin
!                                 ! ----x----I----x---
  character(18) filename_E        ! _NNNN_EX_Vm_2D.bin      ! similar for EY
!                                 ! ----x----I----x---
  character(18) filename_Ne       ! _NNNN_Ne_m3_2D.bin
!                                 ! ----x----I----x----I
  character(20) filename_Je       ! _NNNN_JXe_Am2_2D.bin    ! similar for JY
!                                 ! ----x----I----x----
  character(19) filename_Te       ! _NNNN_TXe_eV_2D.bin     ! similar for TYe 
!                                 ! ----x----I----x----
  character(19) filename_Ve       ! _NNNN_VXe_ms_2D.bin     ! similar for VYe 
!                                 ! ----x----I----x----
  character(19) filename_We       ! _NNNN_WXe_eV_2D.bin     ! similar for WYe 
!                                 ! ----x----I----x----I
  character(20) filename_Ni       ! _NNNN_Ni_1_m3_2D.bin
!                                 ! ----x----I----x----I--
  character(22) filename_Ji       ! _NNNN_JXi_1_Am2_2D.bin  ! similar for JYi
!                                 ! ----x----I----x----I-
  character(21) filename_Ti       ! _NNNN_TXi_1_eV_2D.bin   ! similar for TYi
!                                 ! ----x----I----x----I-
  character(21) filename_Vi       ! _NNNN_VXi_1_ms_2D.bin   ! similar for VYi
!                                 ! ----x----I----x----I-
  character(21) filename_Wi       ! _NNNN_WXi_1_eV_2D.bin   ! similar for WYi

  integer N_points_x
  integer N_points_y

  integer alloc_err
  real, allocatable :: y_m(:)

  real, allocatable :: avgF(:)
  real, allocatable :: avgEX(:)
  real, allocatable :: avgEY(:)

  real, allocatable :: avgNe(:)
  real, allocatable :: avgJeX(:)
  real, allocatable :: avgJeY(:)
  real, allocatable :: avgTeX(:)
  real, allocatable :: avgTeY(:)
  real, allocatable :: avgVeX(:)
  real, allocatable :: avgVeY(:)
  real, allocatable :: avgWeX(:)
  real, allocatable :: avgWeY(:)

  real, allocatable :: avgNi(:)
  real, allocatable :: avgJiX(:)
  real, allocatable :: avgJiY(:)
  real, allocatable :: avgTiX(:)
  real, allocatable :: avgTiY(:)
  real, allocatable :: avgViX(:)
  real, allocatable :: avgViY(:)
  real, allocatable :: avgWiX(:)
  real, allocatable :: avgWiY(:)

  integer proc_status
  integer snap

!                                 ! ----x----I----x----I----x----I
  character(30) filename_out      ! avgxt_values_vsy_NNNN_NNNN.txt

  integer j

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


! identify the size of the domain
! try to access first snapshots with all possible names until success

  print '("identifying the size of the domain...")'
  probestatus  = -1

!---------------------------------------------------- ELECTROSTATIC POTENTIAL -------------------------------------------------!

  if (probestatus.ne.0) then
     filename_F = '_NNNN_F_V_2D.bin'
     filename_F(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_F, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRIC FIELD EX -------------------------------------------------------!

  if (probestatus.ne.0) then
     filename_E = '_NNNN_EX_Vm_2D.bin'
     filename_E(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_E, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRIC FIELD EY -------------------------------------------------------!

  if (probestatus.ne.0) then
     filename_E = '_NNNN_EY_Vm_2D.bin'
     filename_E(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_E, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRON DENSITY Ne -----------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Ne = '_NNNN_Ne_m3_2D.bin'
     filename_Ne(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Ne, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRON CURRENT JXe ----------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Je = '_NNNN_JXe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Je, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRON CURRENT JYe ----------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Je = '_NNNN_JYe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Je, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRON TEMPERATURE TXe ------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Te = '_NNNN_TXe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Te, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRON TEMPERATURE TYe ------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Te = '_NNNN_TYe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Te, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRON VELOCITY VXe ---------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Ve = '_NNNN_VXe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Ve, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRON VELOCITY VYe ---------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Ve = '_NNNN_VYe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Ve, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRON ENERGY WXe -----------------------------------------------------!

  if (probestatus.ne.0) then
     filename_We = '_NNNN_WXe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_We, N_points_x, N_points_y, probestatus )
  end if

!---------------------------------------------------- ELECTRON ENERGY WYe -----------------------------------------------------!

  if (probestatus.ne.0) then
     filename_We = '_NNNN_WYe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_We, N_points_x, N_points_y, probestatus )
  end if

!------------------------------------------------------ ION DENSITY Ni --------------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Ni = '_NNNN_Ni_1_m3_2D.bin'
     filename_Ni(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Ni, N_points_x, N_points_y, probestatus )
  end if

!------------------------------------------------------ ION CURRENT JXi -------------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Ji = '_NNNN_JXi_1_Am2_2D.bin'
     filename_Ji(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Ji, N_points_x, N_points_y, probestatus )
  end if

!------------------------------------------------------ ION CURRENT JYi -------------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Ji = '_NNNN_JYi_1_Am2_2D.bin'
     filename_Ji(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Ji, N_points_x, N_points_y, probestatus )
  end if

!------------------------------------------------------ ION TEMPERATURE TXi ---------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Ti = '_NNNN_TXi_1_eV_2D.bin'
     filename_Ti(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Ti, N_points_x, N_points_y, probestatus )
  end if

!------------------------------------------------------ ION TEMPERATURE TYi ---------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Ti = '_NNNN_TYi_1_eV_2D.bin'
     filename_Ti(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Ti, N_points_x, N_points_y, probestatus )
  end if

!------------------------------------------------------ ION VELOCITY VXi ------------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Vi = '_NNNN_VXi_1_ms_2D.bin'
     filename_Vi(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Vi, N_points_x, N_points_y, probestatus )
  end if

!------------------------------------------------------ ION VELOCITY VYi ------------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Vi = '_NNNN_VYi_1_ms_2D.bin'
     filename_Vi(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Vi, N_points_x, N_points_y, probestatus )
  end if

!------------------------------------------------------ ION ENERGY WXi --------------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Wi = '_NNNN_WXi_1_eV_2D.bin'
     filename_Wi(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Wi, N_points_x, N_points_y, probestatus )
  end if

!------------------------------------------------------ ION ENERGY WYi --------------------------------------------------------!

  if (probestatus.ne.0) then
     filename_Wi = '_NNNN_WYi_1_eV_2D.bin'
     filename_Wi(2:5) = convert_int_to_txt_string(first_snap, 4)
     call probe_binary_data_2( filename_Wi, N_points_x, N_points_y, probestatus )
  end if

! quit if none of the probing attempts above was successful

  if (probestatus.ne.0) then
     print '("could not define the size of the domain, probably no data files found, program terminated")'
     stop
  end if

! since we are here, the size of the domain has been identified eventually

  print '("expected domain size (XxY) :: ",i6," x ",i6)', N_points_x, N_points_y

  allocate(y_m(0:N_points_y-1), stat = alloc_err)

  allocate(avgF(0:N_points_y-1), stat = alloc_err)
  allocate(avgEX(0:N_points_y-1), stat = alloc_err)
  allocate(avgEY(0:N_points_y-1), stat = alloc_err)

  allocate(avgNe(0:N_points_y-1), stat = alloc_err)
  allocate(avgJeX(0:N_points_y-1), stat = alloc_err)
  allocate(avgJeY(0:N_points_y-1), stat = alloc_err)
  allocate(avgTeX(0:N_points_y-1), stat = alloc_err)
  allocate(avgTeY(0:N_points_y-1), stat = alloc_err)
  allocate(avgVeX(0:N_points_y-1), stat = alloc_err)
  allocate(avgVeY(0:N_points_y-1), stat = alloc_err)
  allocate(avgWeX(0:N_points_y-1), stat = alloc_err)
  allocate(avgWeY(0:N_points_y-1), stat = alloc_err)

  allocate(avgNi(0:N_points_y-1), stat = alloc_err)
  allocate(avgJiX(0:N_points_y-1), stat = alloc_err)
  allocate(avgJiY(0:N_points_y-1), stat = alloc_err)
  allocate(avgTiX(0:N_points_y-1), stat = alloc_err)
  allocate(avgTiY(0:N_points_y-1), stat = alloc_err)
  allocate(avgViX(0:N_points_y-1), stat = alloc_err)
  allocate(avgViY(0:N_points_y-1), stat = alloc_err)
  allocate(avgWiX(0:N_points_y-1), stat = alloc_err)
  allocate(avgWiY(0:N_points_y-1), stat = alloc_err)
  
! process the data

!---------------------------------------------------- ELECTROSTATIC POTENTIAL -------------------------------------------------!

  avgF = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_F = '_NNNN_F_V_2D.bin'
     filename_F(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_F, N_points_x, N_points_y, y_m, avgF, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgF = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgF = avgF / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRIC FIELD EX -------------------------------------------------------!

  avgEX = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_E = '_NNNN_EX_Vm_2D.bin'
     filename_E(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_E, N_points_x, N_points_y, y_m, avgEX, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgEX = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgEX = avgEX / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRIC FIELD EY -------------------------------------------------------!

  avgEY = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_E = '_NNNN_EY_Vm_2D.bin'
     filename_E(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_E, N_points_x, N_points_y, y_m, avgEY, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgEY = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgEY = avgEY / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRON DENSITY Ne -----------------------------------------------------!

  avgNe = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Ne = '_NNNN_Ne_m3_2D.bin'
     filename_Ne(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Ne, N_points_x, N_points_y, y_m, avgNe, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgNe = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgNe = avgNe / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRON CURRENT JXe ----------------------------------------------------!

  avgJeX = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Je = '_NNNN_JXe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Je, N_points_x, N_points_y, y_m, avgJeX, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgJeX = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgJeX = avgJeX / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRON CURRENT JYe ----------------------------------------------------!

  avgJeY = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Je = '_NNNN_JYe_Am2_2D.bin'
     filename_Je(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Je, N_points_x, N_points_y, y_m, avgJeY, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgJeY = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgJeY = avgJeY / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRON TEMPERATURE TXe ------------------------------------------------!

  avgTeX = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Te = '_NNNN_TXe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Te, N_points_x, N_points_y, y_m, avgTeX, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgTeX = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgTeX = avgTeX / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRON TEMPERATURE TYe ------------------------------------------------!

  avgTeY = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Te = '_NNNN_TYe_eV_2D.bin'
     filename_Te(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Te, N_points_x, N_points_y, y_m, avgTeY, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgTeY = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgTeY = avgTeY / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRON VELOCITY VXe ---------------------------------------------------!

  avgVeX = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Ve = '_NNNN_VXe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Ve, N_points_x, N_points_y, y_m, avgVeX, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgVeX = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgVeX = avgVeX / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRON VELOCITY VYe ---------------------------------------------------!

  avgVeY = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Ve = '_NNNN_VYe_ms_2D.bin'
     filename_Ve(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Ve, N_points_x, N_points_y, y_m, avgVeY, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgVeY = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgVeY = avgVeY / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRON ENERGY WXe -----------------------------------------------------!

  avgWeX = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_We = '_NNNN_WXe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_We, N_points_x, N_points_y, y_m, avgWeX, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgWeX = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgWeX = avgWeX / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!---------------------------------------------------- ELECTRON ENERGY WYe -----------------------------------------------------!

  avgWeY = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_We = '_NNNN_WYe_eV_2D.bin'
     filename_We(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_We, N_points_x, N_points_y, y_m, avgWeY, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgWeY = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgWeY = avgWeY / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!------------------------------------------------------ ION DENSITY Ni --------------------------------------------------------!

  avgNi = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Ni = '_NNNN_Ni_1_m3_2D.bin'
     filename_Ni(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Ni, N_points_x, N_points_y, y_m, avgNi, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgNi = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgNi = avgNi / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!------------------------------------------------------ ION CURRENT JXi -------------------------------------------------------!

  avgJiX = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Ji = '_NNNN_JXi_1_Am2_2D.bin'
     filename_Ji(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Ji, N_points_x, N_points_y, y_m, avgJiX, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgJiX = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgJiX = avgJiX / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!------------------------------------------------------ ION CURRENT JYi -------------------------------------------------------!

  avgJiY = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Ji = '_NNNN_JYi_1_Am2_2D.bin'
     filename_Ji(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Ji, N_points_x, N_points_y, y_m, avgJiY, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgJiY = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgJiY = avgJiY / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!------------------------------------------------------ ION TEMPERATURE TXi ---------------------------------------------------!

  avgTiX = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Ti = '_NNNN_TXi_1_eV_2D.bin'
     filename_Ti(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Ti, N_points_x, N_points_y, y_m, avgTiX, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgTiX = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgTiX = avgTiX / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!------------------------------------------------------ ION TEMPERATURE TYi ---------------------------------------------------!

  avgTiY = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Ti = '_NNNN_TYi_1_eV_2D.bin'
     filename_Ti(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Ti, N_points_x, N_points_y, y_m, avgTiY, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgTiY = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgTiY = avgTiY / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!------------------------------------------------------ ION VELOCITY VXi ------------------------------------------------------!

  avgViX = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Vi = '_NNNN_VXi_1_ms_2D.bin'
     filename_Vi(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Vi, N_points_x, N_points_y, y_m, avgViX, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgViX = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgViX = avgViX / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!------------------------------------------------------ ION VELOCITY VYi ------------------------------------------------------!

  avgViY = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Vi = '_NNNN_VYi_1_ms_2D.bin'
     filename_Vi(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Vi, N_points_x, N_points_y, y_m, avgViY, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgViY = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgViY = avgViY / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!------------------------------------------------------ ION ENERGY WXi --------------------------------------------------------!

  avgWiX = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Wi = '_NNNN_WXi_1_eV_2D.bin'
     filename_Wi(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Wi, N_points_x, N_points_y, y_m, avgWiX, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgWiX = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgWiX = avgWiX / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!------------------------------------------------------ ION ENERGY WYi --------------------------------------------------------!

  avgWiY = 0.0
  proc_status = -1
  do snap = first_snap, last_snap
     filename_Wi = '_NNNN_WYi_1_eV_2D.bin'
     filename_Wi(2:5) = convert_int_to_txt_string(snap, 4)
     call process_binary_data_4( filename_Wi, N_points_x, N_points_y, y_m, avgWiY, proc_status )
     if (proc_status.ne.0) exit
  end do
  if (proc_status.ne.0) then
     avgWiY = 0.0
     print '("%%%   aborted, sorry   %%%")'
  else
     avgWiY = avgWiY / real((last_snap - first_snap + 1) * (N_points_x-2))
     print '("######   success   ######")'
  end if

!------------------------------------------------------
! save the data

  filename_out = 'avgxt_values_vsy_NNNN_NNNN.txt'
!                 ----x----I----x----I----x----I
  filename_out(18:21) = convert_int_to_txt_string(first_snap, 4)
  filename_out(23:26) = convert_int_to_txt_string(last_snap, 4)

  open (9, file = filename_out)

  write (9, '("# column  1 is Y-coordinate [m]")')

  write (9, '("# values in columns listed below are averaged over X-coordinate and time")')

  write (9, '("# column  2 is electrostatic potential [V]")')
  write (9, '("# column  3 is electric field EX [V/m]")')
  write (9, '("# column  4 is electric field EY [V/m]")')

  write (9, '("# column  5 is electron density [m^-3]")')
  write (9, '("# column  6 is electron current JX [A/m^2]")')
  write (9, '("# column  7 is electron current JY [A/m^2]")')
  write (9, '("# column  8 is electron temperature TX [eV]")')
  write (9, '("# column  9 is electron temperature TY [eV]")')
  write (9, '("# column 10 is electron average velocity VX [m/s]")')
  write (9, '("# column 11 is electron average velocity VY [m/s]")')
  write (9, '("# column 12 is electron average energy WX [eV]")')
  write (9, '("# column 13 is electron average energy WY [eV]")')

  write (9, '("# column 14 is ion density [m^-3]")')
  write (9, '("# column 15 is ion current JX [A/m^2]")')
  write (9, '("# column 16 is ion current JY [A/m^2]")')
  write (9, '("# column 17 is ion temperature TX [eV]")')
  write (9, '("# column 18 is ion temperature TY [eV]")')
  write (9, '("# column 19 is ion average velocity VX [m/s]")')
  write (9, '("# column 20 is ion average velocity VY [m/s]")')
  write (9, '("# column 21 is ion average energy WX [eV]")')
  write (9, '("# column 22 is ion average energy WY [eV]")')

  do j = 0, N_points_y-1
     write (9, '(2x,f12.9,21(2x,e14.7))') &
             & y_m(j), &        ! 1
             & avgF(j), &       ! 2
             & avgEX(j), &      ! 3
             & avgEY(j), &           ! 4
             & avgNe(j), &           ! 5
             & avgJeX(j), &          ! 6
             & avgJeY(j), &          ! 7
             & avgTeX(j), &          ! 8
             & avgTeY(j), &          ! 9
             & avgVeX(j), &          ! 10
             & avgVeY(j), &          ! 11
             & avgWeX(j), &          ! 12
             & avgWeY(j), &          ! 13
             & avgNi(j) , &     ! 14
             & avgJiX(j), &     ! 15
             & avgJiY(j), &     ! 16
             & avgTiX(j), &     ! 17
             & avgTiY(j), &     ! 18
             & avgViX(j), &     ! 19
             & avgViY(j), &     ! 20
             & avgWiX(j), &     ! 21
             & avgWiY(j)        ! 22
  end do

  close (9, status = 'keep')
  print '("file ",A30," is ready, all done")', filename_out

! cleanup
  deallocate(y_m, stat = alloc_err)

  deallocate(avgF, stat = alloc_err)
  deallocate(avgEX, stat = alloc_err)
  deallocate(avgEY, stat = alloc_err)

  deallocate(avgNe, stat = alloc_err)
  deallocate(avgJeX, stat = alloc_err)
  deallocate(avgJeY, stat = alloc_err)
  deallocate(avgTeX, stat = alloc_err)
  deallocate(avgTeY, stat = alloc_err)
  deallocate(avgVeX, stat = alloc_err)
  deallocate(avgVeY, stat = alloc_err)
  deallocate(avgWeX, stat = alloc_err)
  deallocate(avgWeY, stat = alloc_err)

  deallocate(avgNi, stat = alloc_err)
  deallocate(avgJiX, stat = alloc_err)
  deallocate(avgJiY, stat = alloc_err)
  deallocate(avgTiX, stat = alloc_err)
  deallocate(avgTiY, stat = alloc_err)
  deallocate(avgViX, stat = alloc_err)
  deallocate(avgViY, stat = alloc_err)
  deallocate(avgWiX, stat = alloc_err)
  deallocate(avgWiY, stat = alloc_err)

end program create_multi_xt_avgs_vsy_from_2d_bin

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!
subroutine probe_binary_data_2( filename_bin, &
                              & N_points_x, &
                              & N_points_y, &
                              & mystatus )

  implicit none

  character*(*) filename_bin

  integer, intent(out) :: N_points_x
  integer, intent(out) :: N_points_y
  integer, intent(out) :: mystatus

  logical exists

  integer filesize_bytes

  real r_N_points_x

  mystatus = 0
 
! read data
  inquire(file = filename_bin, exist = exists)

  if (.not.exists) then
     print '("probe_binary_data_2 :: ERROR, file ",A30," not found")', filename_bin
     mystatus = -1
     return
  end if

  inquire(file = filename_bin, size = filesize_bytes)
     
  open (19, file = filename_bin, access = 'stream', form = 'unformatted')

  read (19) r_N_points_x

  N_points_x = int(r_N_points_x)

  N_points_y = (filesize_bytes/4) / (N_points_x+1) - 1

  close (19, status = 'keep')

  print '("probe_binary_data_2 :: accessed file ",A30)', filename_bin

end subroutine probe_binary_data_2

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!
subroutine process_binary_data_4( filename_bin, &
                                & exp_N_points_x, &
                                & exp_N_points_y, &
                                & y_m, &
                                & integr_val, &
                                & mystatus )

  implicit none

  character*(*) filename_bin
  integer, intent(in)  :: exp_N_points_x  ! expected number of points along the x direction
  integer, intent(in)  :: exp_N_points_y  ! expected number of points along the y direction
  real,    intent(inout) :: y_m(0:(exp_N_points_y-1))
  real,    intent(inout) :: integr_val(0:(exp_N_points_y-1))  ! must be set to all zeros before the first call
  integer, intent(out) :: mystatus

  logical exists

  integer filesize_bytes

  real r_N_points_x

  integer N_points_x, N_points_y

  real rdummy
! here we don't need x_m array
  real, allocatable :: val(:,:)
  integer alloc_err

  integer i, j

  mystatus = 0

! read data
  inquire(file = filename_bin, exist = exists)

  if (.not.exists) then
     print '("process_binary_data_4 :: ERROR, file ",A30," not found")', filename_bin
     mystatus = -1
     return
  end if

  inquire(file = filename_bin, size = filesize_bytes)
  
  open (19, file = filename_bin, access = 'stream', form = 'unformatted') !, recl=1) 

  read (19) r_N_points_x

  N_points_x = int(r_N_points_x)

  N_points_y = (filesize_bytes/4) / (N_points_x+1) - 1

! foolproof
  if (N_points_x.ne.exp_N_points_x) then
     close (19, status = 'keep')
     print '("process_binary_data_4 :: ERROR, wrong number of points along x in file ",A30," :: ",i6," vs ",i6)', &
          & filename_bin, N_points_x, exp_N_points_x
     mystatus = -2
     return
!     stop
  end if

  if (N_points_y.ne.exp_N_points_y) then
     close (19, status = 'keep')
     print '("process_binary_data_4 :: ERROR, wrong number of points along y in file ",A30," :: ",i6," vs ",i6)', &
          & filename_bin, N_points_y, exp_N_points_y
     mystatus = -3
     return
!     stop
  end if

! note that in the three error situations above, array y_m is not involved
! therefore if this array has been previously defined
! an unsuccessful call of process_binary_data_4 does not change it

  allocate (val(0:N_points_x-1,0:N_points_y-1), stat = alloc_err)

  do i = 0, N_points_x-1
     read (19) rdummy !x_m(i)
  end do

  do j = 0, N_points_y-1
     read (19) y_m(j)
     do i = 0, N_points_x-1
        read (19) val(i,j)
     end do
  end do

  close (19, status = 'keep')

  print '("reading binary file ",A30," complete")', filename_bin

! calculate value integrated over the x-coordinate
! note that in a periodic system, two last saved x-nodes (N_points_x-2, N_points_x-1) correspond to x-nodes 0 and 1
! and therefore must be excluded when integration is performed over the x-period
  do j = 0, N_points_y-1
     do i = 0, N_points_x-3
        integr_val(j) = integr_val(j) + val(i,j)  ! initially, integr_val=0, which is set in the external calling procedure 
     end do
  end do

  if (allocated(val)) deallocate(val, stat = alloc_err)

end subroutine process_binary_data_4



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
