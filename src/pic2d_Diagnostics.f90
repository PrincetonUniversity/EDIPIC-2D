!----------------------------------------
!
subroutine save_electrons

  USE ParallelOperationValues
  USE ElectronParticles

  implicit none

  character(18) snapproc_filename     ! snap_proc_NNNN.dat
                                      ! ----x----I----x---
  integer n
  integer devid

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface

  snapproc_filename = 'snap_proc_NNNN.dat'
  snapproc_filename(11:14) = convert_int_to_txt_string(Rank_of_process, 4)

  devid = 9+Rank_of_process

  open (devid, file=snapproc_filename)
  do n = 1, N_electrons
     write (devid, '(5(2x,e14.7),2x,i3)') &
          & electron(n)%X, &
          & electron(n)%Y, &
          & electron(n)%VX, &
          & electron(n)%VY, &
          & electron(n)%VZ, &
          & electron(n)%tag
  end do
  close (devid, status = 'keep')

  print '("file ",A18," is ready")', snapproc_filename
  
end subroutine save_electrons

!------------------------------------
!
subroutine save_electrons_suf(suf)

  USE ParallelOperationValues
  USE ElectronParticles

  implicit none

  character(3) suf
  character(22) snapproc_filename     ! snap_proc_NNNN_SUF.dat
                                      ! ----x----I----x----I--
  integer n
  integer devid

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface

  snapproc_filename = 'snap_proc_NNNN_SUF.dat'
  snapproc_filename(11:14) = convert_int_to_txt_string(Rank_of_process, 4)
  snapproc_filename(16:18) = suf

  devid = 9+Rank_of_process

  open (devid, file=snapproc_filename)

print '("save_electrons_suf :: proc ",i4," sees N_electrons = ",i8," max_N_electrons = ",i8)', Rank_of_process, N_electrons, max_N_electrons

  do n = 1, N_electrons
     write (devid, '(5(2x,e14.7),2x,i3)') &
          & electron(n)%X, &
          & electron(n)%Y, &
          & electron(n)%VX, &
          & electron(n)%VY, &
          & electron(n)%VZ, &
          & electron(n)%tag
  end do
  close (devid, status = 'keep')
  
  print '("file ",A22," is ready")', snapproc_filename

end subroutine save_electrons_suf

!------------------------------------
!
subroutine save_rho_ei(suf)

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries

  implicit none

  character(3) suf
  character(22) rhoeproc_filename     ! rhoe_proc_NNNN_SUF.dat
                                      ! ----x----I----x----I--
  integer i, j
  integer devid

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface

  rhoeproc_filename = 'rhoe_proc_NNNN_SUF.dat'
  rhoeproc_filename(11:14) = convert_int_to_txt_string(Rank_of_process, 4)
  rhoeproc_filename(16:18) = suf

  devid = 9+Rank_of_process

  open (devid, file=rhoeproc_filename)

  do j = indx_y_min+1, indx_y_max-1
     do i = indx_x_min+1, indx_x_max-1
        write (devid, '(2x,i5,2x,i5,2(2x,e14.7))') &
          & i, &
          & j, &
          & rho_e(i,j), &
          & rho_i(i,j)
     end do
     write (devid, '(" ")')
  end do
  close (devid, status = 'keep')
  
  print '("file ",A22," is ready")', rhoeproc_filename

end subroutine save_rho_ei

!------------------------------------
!
subroutine save_phi(suf)

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries

  implicit none

  character(3) suf
  character(21) phiproc_filename     ! phi_proc_NNNN_SUF.dat
                                     ! ----x----I----x----I-
  integer i, j
  integer devid

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface

!real(8) phi_vs_ij         !############# removeme ################

  phiproc_filename = 'phi_proc_NNNN_SUF.dat'
  phiproc_filename(10:13) = convert_int_to_txt_string(Rank_of_process, 4)
  phiproc_filename(15:17) = suf

  devid = 9+Rank_of_process

  open (devid, file=phiproc_filename)

  do j = indx_y_min, indx_y_max
     do i = indx_x_min, indx_x_max
!        write (devid, '(2x,i5,2x,i5,2x,e14.7,2x,e14.7)') &
        write (devid, '(2x,i5,2x,i5,2x,e14.7)') &
          & i, &
          & j, &
          & phi(i,j) !, &
!          & phi_vs_ij(i,j)        !############# removeme ################
     end do
     write (devid, '(" ")')
  end do
  close (devid, status = 'keep')
  
  print '("file ",A21," is ready")', phiproc_filename

end subroutine save_phi

!------------------------------------
!
subroutine save_E(suf)

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  implicit none

  character(3) suf
  character(21) EXYproc_filename     ! EXY_proc_NNNN_SUF.dat
                                     ! ----x----I----x----I-
  integer i, j
  integer devid

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface

  if (cluster_rank_key.NE.0) return

  EXYproc_filename = 'EXY_proc_NNNN_SUF.dat'
  EXYproc_filename(10:13) = convert_int_to_txt_string(Rank_of_process, 4)
  EXYproc_filename(15:17) = suf

  devid = 9+Rank_of_process

  open (devid, file=EXYproc_filename)

  do j = c_indx_y_min, c_indx_y_max
     do i = c_indx_x_min, c_indx_x_max
        write (devid, '(2x,i5,2x,i5,2(2x,e14.7))') &
          & i, &
          & j, &
          & EX(i,j), &
          & EY(i,j)  
     end do
     write (devid, '(" ")')
  end do
  close (devid, status = 'keep')
  
  print '("file ",A21," is ready")', EXYproc_filename

end subroutine save_E


!------------------------------
!
subroutine report_total_number_of_particles

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles, ONLY : N_electrons, electron
  USE IonParticles, ONLY : N_spec, N_ions, ion

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER ibufer(0:N_spec)
  REAL(8), allocatable :: rbufer(:), pwbufer(:), totpwbufer(:)
  INTEGER ALLOC_ERR
  INTEGER k, s, k1, k2, k3, k4
  INTEGER pos1, pos2

  INTEGER N_particles_cluster(0:N_spec), N_particles_total(0:N_spec)

  INTEGER n
  real(8) surf_char

  ALLOCATE(    rbufer(1:4*(N_spec+1)), STAT = ALLOC_ERR)
  ALLOCATE(   pwbufer(1:4*(N_spec+1)), STAT = ALLOC_ERR)
  ALLOCATE(totpwbufer(1:4*(N_spec+1)), STAT = ALLOC_ERR)

!print '("report_total_number_of_particles :: Rank_of_process ",i4)', Rank_of_process

  rbufer = 0.0_8
  pwbufer = 0.0_8
  totpwbufer = 0.0_8

  DO k = 1, N_electrons
     rbufer(1) = rbufer(1) + electron(k)%VX
     rbufer(2) = rbufer(2) + electron(k)%VY
     rbufer(3) = rbufer(3) + electron(k)%VZ
     rbufer(4) = rbufer(4) + electron(k)%VX**2 + electron(k)%VY**2 + electron(k)%VZ**2
  END DO

  DO s = 1, N_spec
     k1=4*s+1
     k2=4*s+2
     k3=4*s+3
     k4=4*s+4
     DO k = 1, N_ions(s)
        rbufer(k1) = rbufer(k1) + ion(s)%part(k)%VX
        rbufer(k2) = rbufer(k2) + ion(s)%part(k)%VY
        rbufer(k3) = rbufer(k3) + ion(s)%part(k)%VZ
        rbufer(k4) = rbufer(k4) + ion(s)%part(k)%VX**2 + ion(s)%part(k)%VY**2 + ion(s)%part(k)%VZ**2
     END DO
  END DO

  ibufer(0) = N_electrons
  ibufer(1:N_spec) = N_ions(1:N_spec)

! collect number of particles of all species from all processes in a cluster
  CALL MPI_REDUCE(ibufer, N_particles_cluster(0:N_spec), N_spec+1, MPI_INTEGER, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     
! collect number of particles from all processes in a cluster
  CALL MPI_REDUCE(rbufer, pwbufer, 4*(N_spec+1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
     
  IF (cluster_rank_key.EQ.0) THEN

     CALL MPI_REDUCE(N_particles_cluster, N_particles_total, N_spec+1, MPI_INTEGER, MPI_SUM, 0, COMM_HORIZONTAL, ierr)

     CALL MPI_REDUCE(pwbufer, totpwbufer, 4*(N_spec+1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_HORIZONTAL, ierr)

     IF (Rank_horizontal.EQ.0) THEN 
        PRINT '("Total : number of electron particles = ",i10," momentum X/Y/Z = ",3(2x,e16.9)," energy = ",e16.9)', N_particles_total(0), totpwbufer(1:4)
        pos1=5
        pos2=8
        DO s = 1, N_Spec
           PRINT '("Total : number of ion  ",i2,"  particles = ",i10," momentum X/Y/Z = ",3(2x,e16.9)," energy = ",e16.9)', s, N_particles_total(s), totpwbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+4
        END DO
        
        surf_char = 0.0_8
        DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (whole_object(n)%object_type.NE.DIELECTRIC) CYCLE
           DO k = 1, whole_object(n)%N_boundary_nodes
              surf_char = surf_char + whole_object(n)%surface_charge(k)
           END DO
        END DO
        PRINT '("Surface charge on inner dielectric objects ",f12.3," total charge ",f12.3)', surf_char, dble(-N_particles_total(0)+N_particles_total(1))+surf_char

     END IF
  END IF

!print '("report_total_number_of_particles :: Rank_of_process ",i4," cluster_rank_key ",i4," done")', Rank_of_process, cluster_rank_key

end subroutine report_total_number_of_particles

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
