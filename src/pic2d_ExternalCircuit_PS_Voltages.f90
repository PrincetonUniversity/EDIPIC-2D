
!--------------------------------------------
!
SUBROUTINE PREPARE_ECPS_WAVEFORMS

  USE ParallelOperationValues, ONLY : Rank_of_process
  USE ExternalCircuit
  USE CurrentProblemValues, ONLY : delta_t_s, F_scale_V

  IMPLICIT NONE

  INTEGER n

  CHARACTER(25) initpswf_filename   ! init_ecps_NN_waveform.dat   ! "ecps" stands for "external circuit power supply"
                                    ! ----x----I----x----I----x
  LOGICAL exists
  CHARACTER(1) buf
  INTEGER iostatus
  REAL rdummy

  REAL(8) wf_phi_V, wf_nu_Hz
  INTEGER wf_period

  INTEGER ALLOC_ERR
  INTEGER i
  
  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  DO n = 1, N_of_power_supplies

! default values
     EC_power_supply(n)%N_wf_points = 0
     EC_power_supply(n)%use_waveform = .FALSE.

     initpswf_filename = 'init_ecps_NN_waveform.dat'
     initpswf_filename(11:12) = convert_int_to_txt_string(n, 2)

     INQUIRE (FILE = initpswf_filename, EXIST = exists)
     IF (.NOT.exists) CYCLE

     IF (Rank_of_process.EQ.0) PRINT '("### PREPARE_ECPS_WAVEFORMS :: found file ",A25," for external circuit power supply ",i2," analyzing...")', initpswf_filename, n

     OPEN (11, FILE = initpswf_filename)
     READ (11, '(A1)') buf   ! potential amplitude [V], skip for now
     READ (11, '(A1)') buf   ! frequency [Hz], skip for now
     READ (11, '(A1)') buf   ! comment line, skip
     READ (11, '(A1)') buf   ! comment line, skip
     
     DO 
        READ (11, *, iostat = iostatus) rdummy, rdummy
        IF (iostatus.NE.0) EXIT
        EC_power_supply(n)%N_wf_points = EC_power_supply(n)%N_wf_points + 1
     END DO
     CLOSE (11, STATUS = 'KEEP')

     IF (EC_power_supply(n)%N_wf_points.LE.1) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### PREPARE_ECPS_WAVEFORMS :: WARNING-1 :: not enough (",i4,") valid data points in file ",A25," , waveform for external circuit power supply ",i2," is off ###")', &
             & EC_power_supply(n)%N_wf_points, initpswf_filename,  n
        CYCLE
     END IF

     OPEN (11, FILE = initpswf_filename)
     READ (11, *) wf_phi_V                ! potential amplitude [V]
     READ (11, *) wf_nu_Hz                ! frequency [Hz]

     IF (wf_phi_V.EQ.0.0) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### PREPARE_ECPS_WAVEFORMS :: WARNING-2 :: zero waveform amplitude in file ",A23," , waveform for external circuit power supply ",i2," is off ###")', initpswf_filename,  n
        CYCLE
     END IF

     wf_period = 1.0_8 / (wf_nu_Hz * delta_t_s) ! in units of time steps
     IF (wf_period.LT.(EC_power_supply(n)%N_wf_points-1)) THEN
        IF (Rank_of_process.EQ.0) &
             & PRINT '("### PREPARE_ECPS_WAVEFORMS :: WARNING-3 :: in file ",A25," too many points (",i6,") for given waveform period of ",i6," time steps, waveform for external circuit power supply ",i2," is off ###")', &
             & initpswf_filename, EC_power_supply(n)%N_wf_points, wf_period, n
        CYCLE
     END IF

     ALLOCATE (EC_power_supply(n)%wf_T_cntr(1:EC_power_supply(n)%N_wf_points), STAT = ALLOC_ERR)
     ALLOCATE (EC_power_supply(n)%wf_phi(   1:EC_power_supply(n)%N_wf_points), STAT = ALLOC_ERR)

     READ (11, '(A1)') buf   ! comment line, skip
     READ (11, '(A1)') buf   ! comment line, skip
     DO i = 1, EC_power_supply(n)%N_wf_points
        READ (11, *) rdummy, EC_power_supply(n)%wf_phi(i)
        EC_power_supply(n)%wf_T_cntr(i) = INT(rdummy * wf_period)
        EC_power_supply(n)%wf_phi(i) = (wf_phi_V / F_scale_V) * EC_power_supply(n)%wf_phi(i)
     END DO
     CLOSE (11, STATUS = 'KEEP')

! enforce the ends
     EC_power_supply(n)%wf_T_cntr(1) = 0
     EC_power_supply(n)%wf_T_cntr(EC_power_supply(n)%N_wf_points) = wf_period

! enforce increasing times
     DO i = 2, EC_power_supply(n)%N_wf_points-1
         EC_power_supply(n)%wf_T_cntr(i) = MAX(EC_power_supply(n)%wf_T_cntr(i), EC_power_supply(n)%wf_T_cntr(i-1)+1)
     END DO

! final check
     i = EC_power_supply(n)%N_wf_points
     IF (EC_power_supply(n)%wf_T_cntr(i).GT.EC_power_supply(n)%wf_T_cntr(i-1)) THEN
! passed
        IF (Rank_of_process.EQ.0) &
             & PRINT '("### PREPARE_ECPS_WAVEFORMS :: voltage of external circuit power supply ",i2," will be calculated with the waveform  ###")', n
        EC_power_supply(n)%use_waveform = .TRUE.
     ELSE
! did not pass
        IF (Rank_of_process.EQ.0) &
             & PRINT '("### PREPARE_ECPS_WAVEFORMS :: WARNING-4 :: inconsistent data in file ",A25," , waveform for external circuit power supply ",i2," is off ###")', &
             & initpswf_filename, n
! cleanup
        IF (ALLOCATED(EC_power_supply(n)%wf_T_cntr)) DEALLOCATE(EC_power_supply(n)%wf_T_cntr, STAT = ALLOC_ERR)
        IF (ALLOCATED(EC_power_supply(n)%wf_phi)) DEALLOCATE(EC_power_supply(n)%wf_phi, STAT = ALLOC_ERR)
     END IF

  END DO   !###   DO n = 1, N_of_power_supplies

END SUBROUTINE PREPARE_ECPS_WAVEFORMS

!--------------------------------------------
!
SUBROUTINE PREPARE_ECPS_OSCILLATIONS_AMPLITUDE_PROFILE

  USE ParallelOperationValues, ONLY : Rank_of_process
  USE ExternalCircuit
  USE CurrentProblemValues, ONLY : delta_t_s

  IMPLICIT NONE

  INTEGER n

  CHARACTER(34) initpsap_filename   ! init_ecps_NN_amplitude_profile.dat
                                    ! ----x----I----x----I----x----I----
  LOGICAL exists
  CHARACTER(1) buf
  INTEGER iostatus
  REAL rdummy

  INTEGER wf_period
  INTEGER delta_T_cntr

  INTEGER ALLOC_ERR
  INTEGER i
  
  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  DO n = 1, N_of_power_supplies

! default values
     EC_power_supply(n)%N_ap_points = 0
     EC_power_supply(n)%use_amplitude_profile = .FALSE.

     initpsap_filename = 'init_ecps_NN_amplitude_profile.dat'
     initpsap_filename(11:12) = convert_int_to_txt_string(n, 2)

     INQUIRE (FILE = initpsap_filename, EXIST = exists)
     IF (.NOT.exists) CYCLE

     IF (Rank_of_process.EQ.0) PRINT '("### PREPARE_ECPS_OSCILLATIONS_AMPLITUDE_PROFILE :: found file ",A34," for external circuit power supply ",i2," analyzing...")', initpsap_filename, n

     OPEN (11, FILE = initpsap_filename)
     READ (11, '(A1)') buf   ! column 1 is time (ns), column 2 is amplitude factor (dimensionless)
     DO 
        READ (11, *, iostat = iostatus) rdummy, rdummy
        IF (iostatus.NE.0) EXIT
        EC_power_supply(n)%N_ap_points = EC_power_supply(n)%N_ap_points + 1
     END DO
     CLOSE (11, STATUS = 'KEEP')

     IF (EC_power_supply(n)%N_ap_points.LE.1) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### PREPARE_ECPS_OSCILLATIONS_AMPLITUDE_PROFILE :: WARNING-1 :: not enough (",i4,") valid data points in file ",A34," , oscillations amplitude profile for external circuit power supply ",i2," is off ###")', &
             & EC_power_supply(n)%N_ap_points, initpsap_filename,  n
        CYCLE
     END IF

     OPEN (11, FILE = initpsap_filename)

     ALLOCATE (EC_power_supply(n)%ap_T_cntr(1:EC_power_supply(n)%N_ap_points), STAT = ALLOC_ERR)
     ALLOCATE (EC_power_supply(n)%ap_factor(1:EC_power_supply(n)%N_ap_points), STAT = ALLOC_ERR)

     READ (11, '(A1)') buf   ! column 1 is time (ns), column 2 is amplitude factor (dimensionless)
     DO i = 1, EC_power_supply(n)%N_ap_points
        READ (11, *) rdummy, EC_power_supply(n)%ap_factor(i)
        EC_power_supply(n)%ap_T_cntr(i) = INT(rdummy * 1.0d-9 / delta_t_s)
     END DO
     CLOSE (11, STATUS = 'KEEP')

! enforce the very first point
     EC_power_supply(n)%ap_T_cntr(1) = 0

! enforce increasing times
     DO i = 2, EC_power_supply(n)%N_ap_points-1
         EC_power_supply(n)%ap_T_cntr(i) = MAX(EC_power_supply(n)%ap_T_cntr(i), EC_power_supply(n)%ap_T_cntr(i-1)+1)
     END DO

! if the object uses waveforms, adjust ends of non-zero-amplitude-factor intervals to an integer number of waveform periods
     IF (EC_power_supply(n)%use_waveform) THEN

        wf_period = EC_power_supply(n)%wf_T_cntr(EC_power_supply(n)%N_wf_points)

        DO i = 1, EC_power_supply(n)%N_ap_points-1
           IF ((EC_power_supply(n)%ap_factor(i).EQ.0.0_8).AND.(EC_power_supply(n)%ap_factor(i+1).NE.0.0_8)) THEN
! the non-zero interval begins
              delta_T_cntr = EC_power_supply(n)%ap_T_cntr(i+1) - EC_power_supply(n)%ap_T_cntr(i)
              EC_power_supply(n)%ap_T_cntr(i) = wf_period * INT(EC_power_supply(n)%ap_T_cntr(i) / wf_period)
              EC_power_supply(n)%ap_T_cntr(i+1) = EC_power_supply(n)%ap_T_cntr(i) + delta_T_cntr
           END IF
        END DO

        DO i = 2, EC_power_supply(n)%N_ap_points
           IF ((EC_power_supply(n)%ap_factor(i-1).NE.0.0_8).AND.(EC_power_supply(n)%ap_factor(i).EQ.0.0_8)) THEN
! the non-zero interval ends
              delta_T_cntr = EC_power_supply(n)%ap_T_cntr(i) - EC_power_supply(n)%ap_T_cntr(i-1)
              EC_power_supply(n)%ap_T_cntr(i) = wf_period * INT(EC_power_supply(n)%ap_T_cntr(i) / wf_period)
              EC_power_supply(n)%ap_T_cntr(i-1) = EC_power_supply(n)%ap_T_cntr(i) - delta_T_cntr
           END IF
        END DO

        i = EC_power_supply(n)%N_ap_points
        EC_power_supply(n)%ap_T_cntr(i) = wf_period * INT(EC_power_supply(n)%ap_T_cntr(i) / wf_period)

     END IF

! final check
     EC_power_supply(n)%use_amplitude_profile = .TRUE.
     DO i = 1, EC_power_supply(n)%N_ap_points-1
        IF (EC_power_supply(n)%ap_T_cntr(i+1).GT.EC_power_supply(n)%ap_T_cntr(i)) CYCLE
        EC_power_supply(n)%use_amplitude_profile = .FALSE.
        EXIT
     END DO

     IF (EC_power_supply(n)%use_amplitude_profile) THEN
! passed
        IF (Rank_of_process.EQ.0) THEN
           PRINT '("### PREPARE_ECPS_OSCILLATIONS_AMPLITUDE_PROFILE :: oscillatory voltage of external circuit power supply ",i2," will be calculated with the variable amplitude ###")', n
           DO i = 1, EC_power_supply(n)%N_ap_points
              PRINT '(2x,i3,2x,i4,2x,i10,2x,f8.3)', n, i, EC_power_supply(n)%ap_T_cntr(i), EC_power_supply(n)%ap_factor(i)
           END DO
        END IF
     ELSE
! did not pass
        IF (Rank_of_process.EQ.0) THEN
           PRINT '("### PREPARE_ECPS_OSCILLATIONS_AMPLITUDE_PROFILE :: WARNING-4 :: inconsistent data in file ",A34," , oscillation amplitude profile for external circuit power supply ",i2," is off ###")', &
             & initpsap_filename, n
           DO i = 1, EC_power_supply(n)%N_ap_points
              PRINT '(12x,i3,2x,i4,2x,i10,2x,f8.3)', n, i, EC_power_supply(n)%ap_T_cntr(i), EC_power_supply(n)%ap_factor(i)
           END DO
        END IF
! cleanup
        IF (ALLOCATED(EC_power_supply(n)%ap_T_cntr)) DEALLOCATE(EC_power_supply(n)%ap_T_cntr, STAT = ALLOC_ERR)
        IF (ALLOCATED(EC_power_supply(n)%ap_factor)) DEALLOCATE(EC_power_supply(n)%ap_factor, STAT = ALLOC_ERR)
     END IF

  END DO

END SUBROUTINE PREPARE_ECPS_OSCILLATIONS_AMPLITUDE_PROFILE

!-------------------------------
!
REAL(8) FUNCTION ECPS_Voltage(n_ps, T_cntr)

  USE ExternalCircuit

  IMPLICIT NONE

  INTEGER n_ps        ! number of the external circuit power supply
  INTEGER T_cntr

! functions
  REAL(8) ECPS_Harmonic_Oscillation
  REAL(8) ECPS_Custom_Waveform
  REAL(8) ECPS_Amplitude_profile

  ECPS_Voltage = 0.0_8

  IF (n_ps.LT.1) RETURN
  IF (n_ps.GT.N_of_power_supplies) RETURN 

  ECPS_Voltage = EC_power_supply(n_ps)%phi_const + ( ECPS_Harmonic_Oscillation(n_ps, T_cntr) + ECPS_Custom_Waveform(n_ps, T_cntr) ) * ECPS_Amplitude_Profile(n_ps, T_cntr)

END FUNCTION ECPS_Voltage


!-------------------------------------
!
REAL(8) FUNCTION ECPS_Harmonic_Oscillation(n_ps, T_cntr)

  USE ExternalCircuit
  USE CurrentProblemValues, ONLY : pi

  IMPLICIT NONE

  INTEGER n_ps        ! number of the external circuit power supply
  INTEGER T_cntr

  INTEGER ap_period
  INTEGER n_periods

  ECPS_Harmonic_Oscillation = 0.0_8

  IF (EC_power_supply(n_ps)%use_amplitude_profile) THEN

     ap_period = EC_power_supply(n_ps)%ap_T_cntr(EC_power_supply(n_ps)%N_ap_points)
     IF (T_cntr.EQ.(INT(T_cntr / ap_period) * ap_period)) THEN

!EC_power_supply(n)%omega * T_cntr + EC_power_supply(n)%phase_adjusted = integer * 2 * pi + EC_power_supply(n)%phase

        n_periods = INT((EC_power_supply(n_ps)%omega * T_cntr) / (2.0_8 * pi))
        EC_power_supply(n_ps)%phase_adjusted = n_periods * 2.0_8 * pi + EC_power_supply(n_ps)%phase - EC_power_supply(n_ps)%omega * T_cntr
     END IF

     ECPS_Harmonic_Oscillation = EC_power_supply(n_ps)%phi_var * SIN(EC_power_supply(n_ps)%omega * T_cntr + EC_power_supply(n_ps)%phase_adjusted)

  ELSE

     ECPS_Harmonic_Oscillation = EC_power_supply(n_ps)%phi_var * SIN(EC_power_supply(n_ps)%omega * T_cntr + EC_power_supply(n_ps)%phase)

  END IF

END FUNCTION ECPS_Harmonic_Oscillation

!-------------------------------------
!
REAL(8) FUNCTION ECPS_Custom_Waveform(n_ps, T_cntr)

  USE ExternalCircuit

  IMPLICIT NONE

  INTEGER n_ps        ! number of the external circuit power supply
  INTEGER T_cntr

  INTEGER wf_period
  INTEGER my_wf_T_cntr
  INTEGER i
  REAL ai

  ECPS_Custom_Waveform = 0.0_8

  IF (.NOT.EC_power_supply(n_ps)%use_waveform) RETURN

  wf_period = EC_power_supply(n_ps)%wf_T_cntr(EC_power_supply(n_ps)%N_wf_points)

  my_wf_T_cntr = T_cntr - INT(T_cntr / wf_period) * wf_period

  DO i = 1, EC_power_supply(n_ps)%N_wf_points-1
     IF ( (my_wf_T_cntr.GE.EC_power_supply(n_ps)%wf_T_cntr(i)).AND.(my_wf_T_cntr.LE.EC_power_supply(n_ps)%wf_T_cntr(i+1)) ) THEN

        ai = DBLE(EC_power_supply(n_ps)%wf_T_cntr(i+1) - my_wf_T_cntr) / DBLE(EC_power_supply(n_ps)%wf_T_cntr(i+1) - EC_power_supply(n_ps)%wf_T_cntr(i))

        ECPS_Custom_Waveform = DBLE(EC_power_supply(n_ps)%wf_phi(i) * ai + EC_power_supply(n_ps)%wf_phi(i+1) * (1.0_8 - ai))

        RETURN
     END IF
  END DO

END FUNCTION ECPS_Custom_Waveform

!-------------------------------------
!
REAL(8) FUNCTION ECPS_Amplitude_Profile(n_ps, T_cntr)

  USE ExternalCircuit

  IMPLICIT NONE

  INTEGER n_ps        ! number of the external circuit power supply
  INTEGER T_cntr

  INTEGER ap_period
  INTEGER my_ap_T_cntr
  INTEGER i
  REAL ai

  ECPS_Amplitude_Profile = 1.0_8

  IF (.NOT.EC_power_supply(n_ps)%use_amplitude_profile) RETURN

  ap_period = EC_power_supply(n_ps)%ap_T_cntr(EC_power_supply(n_ps)%N_ap_points)

  my_ap_T_cntr = T_cntr - INT(T_cntr / ap_period) * ap_period

  DO i = 1, EC_power_supply(n_ps)%N_ap_points-1
     IF ( (my_ap_T_cntr.GE.EC_power_supply(n_ps)%ap_T_cntr(i)).AND.(my_ap_T_cntr.LE.EC_power_supply(n_ps)%ap_T_cntr(i+1)) ) THEN

        ai = DBLE(EC_power_supply(n_ps)%ap_T_cntr(i+1) - my_ap_T_cntr) / DBLE(EC_power_supply(n_ps)%ap_T_cntr(i+1) - EC_power_supply(n_ps)%ap_T_cntr(i))

        ECPS_Amplitude_Profile = DBLE(EC_power_supply(n_ps)%ap_factor(i) * ai + EC_power_supply(n_ps)%ap_factor(i+1) * (1.0_8 - ai))

        RETURN
     END IF
  END DO

END FUNCTION ECPS_Amplitude_Profile

!-------------------------------------------
! This subroutine adjusts phase of harmonic oscillations of voltage of external circuit power supplies.
! It is similar to subroutine ADJUST_HARMONIC_OSCILLATIONS_PHASE which is applied to voltages of boundary/inner objects.
! The phase adjustment is necessary when the harmonically oscillating voltage is applied in pulses (that is Amplitude_Profile is not always 1).
! This subroutine does what is normally executed in the beginning of each pulse by the ECPS_Harmonic_Oscillation function.
! Note that when the simulation is restarted from a checkpoint, the beginning of the pulse may be missed,
! and the proper value of the phase will not be specified.
!
SUBROUTINE ADJUST_ECPS_HARMONIC_OSCILLATIONS_PHASE

  USE ExternalCircuit
  USE CurrentProblemValues, ONLY : pi, Start_T_cntr

  IMPLICIT NONE

  INTEGER n_ps        ! number of the whole object
  INTEGER my_T_cntr

  INTEGER ap_period
  INTEGER n_periods

  IF (Start_T_cntr.EQ.0) RETURN

! we are here if the checkpoint was requested and Start_T_cntr is not zero

  DO n_ps = 1, N_of_power_supplies
     IF (EC_power_supply(n_ps)%use_amplitude_profile) THEN

        ap_period = EC_power_supply(n_ps)%ap_T_cntr(EC_power_supply(n_ps)%N_ap_points)   ! pulse period

        my_T_cntr = INT(Start_T_cntr / ap_period) * ap_period
        n_periods = INT((EC_power_supply(n_ps)%omega * my_T_cntr) / (2.0_8 * pi))
 
        EC_power_supply(n_ps)%phase_adjusted = n_periods * 2.0_8 * pi + EC_power_supply(n_ps)%phase - EC_power_supply(n_ps)%omega * my_T_cntr
 
     END IF
  END DO

END SUBROUTINE ADJUST_ECPS_HARMONIC_OSCILLATIONS_PHASE
  
