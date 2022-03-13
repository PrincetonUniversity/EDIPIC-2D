!--------------------------------------
!
SUBROUTINE PREPARE_EXTERNAL_FIELDS

  USE ParallelOperationValues
  USE ExternalFields
  USE CurrentProblemValues, ONLY : E_scale_Vm, B_scale_T, delta_x_m, global_maximal_j, pi, mu_0_Hm
  USE IonParticles, ONLY : ions_sense_magnetic_field, ions_sense_EZ

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  LOGICAL exists

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  CHARACTER(1) buf
  INTEGER ions_sense_magnetic_field_flag, ions_sense_EZ_flag
  INTEGER j

  INTEGER ALLOC_ERR
  INTEGER n

! functions
  REAL(8) Bx, By, Bz, Ez

  INQUIRE (FILE = 'init_extfields.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (exists) THEN

     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i5," : init_extfields.dat is found. Reading the data file...")', Rank_of_process
     END IF

     OPEN (9, FILE = 'init_extfields.dat')

     READ (9, '(A1)') buf !"---ddddd.ddd--- X-magnetic field [Gauss]")')
     READ (9, '(3x,f9.3)') Bx_ext
     READ (9, '(A1)') buf !"---ddddd.ddd--- Y-magnetic field [Gauss]")')
     READ (9, '(3x,f9.3)') By_ext

!###     READ (9, '(A1)') buf !"---ddddd.ddd--- Z-magnetic field [Gauss]")')
!###     READ (9, '(3x,f9.3)') Bz_ext

     READ (9, '(A1)') buf !"--------------- parameters of Z-magnetic field as used by Boeuf and Garrigues")')
     READ (9, '(A1)') buf !"---ddddd.ddd--- Y-coordinate of the BZ maximum, y_Bmax [cm]")')
     READ (9, '(3x,f9.3)') y_Bmax
     READ (9, '(A1)') buf !"---ddddd.ddd--- BZ at y=0 [Gauss]")')
     READ (9, '(3x,f9.3)') Bz_0
     READ (9, '(A1)') buf !"---ddddd.ddd--- BZ at y_Bmax [Gauss]")')
     READ (9, '(3x,f9.3)') Bz_max
     READ (9, '(A1)') buf !"---ddddd.ddd--- BZ at y=Lsys [Gauss]")')
     READ (9, '(3x,f9.3)') Bz_Lsys
     READ (9, '(A1)') buf !"---ddddd.ddd--- characteristic length of decay for y<y_Bmax [cm]")')
     READ (9, '(3x,f9.3)') half_over_sigma2_1
     READ (9, '(A1)') buf !"---ddddd.ddd--- characteristic length of decay for y>y_Bmax [cm]")')
     READ (9, '(3x,f9.3)') half_over_sigma2_2

     READ (9, '(A1)') buf !"---ddddd.ddd--- Z-electric field [V/cm]")')
     READ (9, '(3x,f9.3)') Ez_ext
     READ (9, '(A1)') buf !"-------d------- ions sense magnetic field [1=Yes, 0=No]")')
     READ (9, '(7x,i1)') ions_sense_magnetic_field_flag
     READ (9, '(A1)') buf !"-------d------- ions sense Z-electric field [1=Yes, 0=No]")')
     READ (9, '(7x,i1)') ions_sense_EZ_flag

     CLOSE (9, STATUS = 'KEEP')

!B_ext = 1.0d-4 * 100.0_8 / B_scale_T !################  100 Gauss = 0.01 Tesla
     Bx_ext = Bx_ext * 1.0d-4 / B_scale_T  
     By_ext = By_ext * 1.0d-4 / B_scale_T  
!###  Bz_ext = Bz_ext * 1.0d-4 / B_scale_T

     Bz_0    = Bz_0    * 1.0d-4 / B_scale_T

     Bz_max  = Bz_max  * 1.0d-4 / B_scale_T
!### warning ###, nonuniform Bz profile is disabled, to set thte constant external Bz use Bz_max
!### the way how the external fields are specified needs to be revised to make it more universal 
     Bz_ext = Bz_max   !###
     IF (Rank_of_process.EQ.0) PRINT '("### Uniform Bz = ",e12.5," (Gauss) is assumed, set via BZ at y_Bmax in init_extfields.dat ###")', Bz_ext * B_scale_T * 1.0d4

     Bz_Lsys = Bz_Lsys * 1.0d-4 / B_scale_T
     y_Bmax = y_Bmax * 0.01_8 / delta_x_m

     half_over_sigma2_1 = 0.5_8 * (delta_x_m * 100.0_8 / half_over_sigma2_1)**2
     half_over_sigma2_2 = 0.5_8 * (delta_x_m * 100.0_8 / half_over_sigma2_2)**2

     a1 = (Bz_max - Bz_0) / (1.0_8 - EXP(-half_over_sigma2_1 * y_Bmax**2))
     a2 = (Bz_max - Bz_Lsys) / (1.0_8 - EXP(-half_over_sigma2_2 * (DBLE(global_maximal_j)-y_Bmax)**2))

!  b1 = Bz_0 - a1 * EXP(-half_over_sigma2_1 * y_Bmax**2)
!  b2 = Bz_Lsys - a2 * EXP(-half_over_sigma2_2 * (DBLE(global_maximal_j)-y_Bmax)**2)

     b1 = Bz_max - a1
     b2 = Bz_max - a2


     IF (Rank_of_process.EQ.0) THEN
        OPEN (10, FILE = 'external_Bz_vs_y.dat')
        WRITE (10, '("# column 1 is the y-node number [dim-less]")')
        WRITE (10, '("# column 2 is the y-node coordinate [cm]")')
        WRITE (10, '("# column 3 is the BZ [Gauss]")')
        DO j = 0, global_maximal_j
           WRITE (10, '(2x,i6,2x,f12.9,2x,f10.4)') j, j*delta_x_m*100.0_8, Bz(0.0_8, DBLE(j)) * B_scale_T * 1.0d4     ! save magnetic field in Gauss
        END DO
        CLOSE (10, STATUS = 'KEEP')
        PRINT '("Process 0 created file external_Bz_vs_y.dat")'
     END IF

     Ez_ext = Ez_ext * 100.0_8 / E_scale_Vm

     IF (ions_sense_magnetic_field_flag.EQ.0) THEN
        ions_sense_magnetic_field = .FALSE.
     ELSE
        ions_sense_magnetic_field = .TRUE.
     END IF

     IF (ions_sense_EZ_flag.EQ.0) THEN
        ions_sense_EZ = .FALSE.
     ELSE
        ions_sense_EZ = .TRUE.
     END IF

  ELSE
     
     IF (Rank_of_process.EQ.0) PRINT '("### init_extfields.dat not found, use default values :: external fields Bx,By,Bz, and Ez are zero, ions do not sense magnetic and Z-electric field")'

     Bx_ext = 0.0_8
     By_ext = 0.0_8

     y_Bmax = 0.0_8
     a1 = 0.0_8
     b1 = 0.0_8
     a2 = 0.0_8
     b2 = 0.0_8
     half_over_sigma2_1 = 0.0_8
     half_over_sigma2_2 = 0.0_8

     Ez_ext = 0.0_8

     ions_sense_magnetic_field = .FALSE.
     ions_sense_EZ = .FALSE.

  END IF

  N_JZ_wires = 0

  INQUIRE (FILE = 'init_extmagfieldsBxBy.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (exists) THEN
     
     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i5," : init_extmagfieldsBxBy.dat is found. Reading the data file...")', Rank_of_process
     END IF

     OPEN (9, FILE = 'init_extmagfieldsBxBy.dat')

     READ (9, '(A1)') buf !"---dd--- number of wires with the electric current along the Z-direction")')
     READ (9, '(3x,i2)') N_JZ_wires

     IF (N_JZ_wires.GT.0) THEN
        ALLOCATE(JZwire_X(1:N_JZ_wires), STAT=ALLOC_ERR)   ! x-coordinate of the wire in the laboratory CS
        ALLOCATE(JZwire_Y(1:N_JZ_wires), STAT=ALLOC_ERR)   ! y-coordinate
        ALLOCATE(JZwire_JZ(1:N_JZ_wires), STAT=ALLOC_ERR)   ! electric current (JZ) 

        READ (9, '(A1)') buf !"--- provide below for each wire its X coordinate [mm] Y coordinate [mm] and electric current JZ [A]")')
        READ (9, '(A1)') buf !"---ddddd.ddd---ddddd.ddd---ddddd.ddd---")')

        DO n = 1, N_JZ_wires
           READ (9, '(3(3x,f9.3))') JZwire_X(n), JZwire_Y(n), JZwire_JZ(n) 
! make values dimensionless
           JZwire_X(n) = JZwire_X(n) * 0.001_8 / delta_x_m
           JZwire_Y(n) = JZwire_Y(n) * 0.001_8 / delta_x_m
           JZwire_JZ(n) = (mu_0_Hm * JZwire_JZ(n) / (2.0_8 * pi)) / (delta_x_m * B_scale_T)
        END DO
     ELSE
        IF (Rank_of_process.EQ.0) PRINT '("### init_extmagfieldsBxBy.dat found but NO external currents requested, corresponding Bx, By are zero")'
     END IF
     CLOSE (9, STATUS = 'KEEP')

  ELSE

     IF (Rank_of_process.EQ.0) PRINT '("### init_extmagfieldsBxBy.dat not found, magnetic fields Bx, By due to external currents are zero")'

  END IF

END SUBROUTINE PREPARE_EXTERNAL_FIELDS


!----------------------
!
REAL(8) FUNCTION Bx(x, y)

  USE ExternalFields

  IMPLICIT NONE

  REAL(8) x, y
  INTEGER n

  Bx = Bx_ext

  DO n = 1, N_JZ_wires
     Bx = Bx - JZwire_JZ(n) * (y - JZwire_Y(n)) / ((x - JZwire_X(n))**2 + (y - JZwire_Y(n))**2)
  END DO

END FUNCTION Bx

!----------------------
!
REAL(8) FUNCTION By(x, y)

  USE ExternalFields

  IMPLICIT NONE

  REAL(8) x, y
  INTEGER n

  By = By_ext

  DO n = 1, N_JZ_wires
     By = By + JZwire_JZ(n) * (x - JZwire_X(n)) / ((x - JZwire_X(n))**2 + (y - JZwire_Y(n))**2)
  END DO

END FUNCTION By

!----------------------
!
REAL(8) FUNCTION Bz(x, y)

  USE ExternalFields

  IMPLICIT NONE

  REAL(8) x, y
  REAL(8) h

  Bz = Bz_ext

!### new - the nonuniform magnetic field profile is disabled for now, constant Bz assumed
!### 
!###  h = y - y_Bmax
!###  IF (h.LT.0.0_8) THEN
!###     Bz = a1 * EXP(-half_over_sigma2_1 * h * h) + b1
!###  ELSE
!###     Bz = a2 * EXP(-half_over_sigma2_2 * h * h) + b2
!###  END IF

END FUNCTION Bz

!----------------------
!
REAL(8) FUNCTION Ez(x, y)

  USE ExternalFields

  IMPLICIT NONE

  REAL(8) x, y

  Ez = Ez_ext

END FUNCTION Ez
