SUBROUTINE UPDATE_WALL_POTENTIALS(T_cntr)

  USE CurrentProblemValues, ONLY : whole_object, N_of_boundary_objects, N_of_boundary_and_inner_objects, METAL_WALL, VACUUM_GAP

  IMPLICIT NONE

  INTEGER T_cntr
  INTEGER n, i, j

! function
  REAL(8) Harmonic_Oscillation
  REAL(8) Custom_Waveform
  REAL(8) Amplitude_profile

! set potentials of metal electrodes
  DO n = 1, N_of_boundary_and_inner_objects
     IF (whole_object(n)%object_type.NE.METAL_WALL) CYCLE   
     IF (whole_object(n)%potential_must_be_solved) THEN
        whole_object(n)%phi = 0.0_8
     ELSE
        whole_object(n)%phi = whole_object(n)%phi_const + &
!                            & ( whole_object(n)%phi_var * SIN(whole_object(n)%omega * T_cntr + whole_object(n)%phase) + &
                            & ( Harmonic_Oscillation(n, T_cntr) + &
                            &   Custom_Waveform(n, T_cntr) ) * &
                            & Amplitude_Profile(n, T_cntr)
     END IF
  END DO

! set potentials across vacuum gaps
  DO n = 1, N_of_boundary_objects
     IF (whole_object(n)%object_type.NE.VACUUM_GAP) CYCLE
     
     IF (whole_object(n)%segment(1)%jstart.NE.whole_object(n)%segment(1)%jend) THEN
! the gap is vertical (stretched along y)
        DO j = whole_object(n)%segment(1)%jstart+1, whole_object(n)%segment(1)%jend-1
           whole_object(n)%phi_profile(j) = whole_object(whole_object(n)%n_connected_to_start)%phi +  &
                                          & (whole_object(whole_object(n)%n_connected_to_end)%phi - whole_object(whole_object(n)%n_connected_to_start)%phi) * &
                                          & DBLE(j - whole_object(n)%segment(1)%jstart) / DBLE(whole_object(n)%segment(1)%jend - whole_object(n)%segment(1)%jstart)
        END DO
        whole_object(n)%phi_profile(whole_object(n)%segment(1)%jstart) = whole_object(whole_object(n)%n_connected_to_start)%phi
        whole_object(n)%phi_profile(whole_object(n)%segment(1)%jend)   = whole_object(whole_object(n)%n_connected_to_end)%phi
     ELSE IF (whole_object(n)%segment(1)%istart.NE.whole_object(n)%segment(1)%iend) THEN
! the gap is horizontal (stretched along x)
        DO i = whole_object(n)%segment(1)%istart+1, whole_object(n)%segment(1)%iend-1
           whole_object(n)%phi_profile(i) = whole_object(whole_object(n)%n_connected_to_start)%phi +  &
                                          & (whole_object(whole_object(n)%n_connected_to_end)%phi - whole_object(whole_object(n)%n_connected_to_start)%phi) * &
                                          & DBLE(i - whole_object(n)%segment(1)%istart) / DBLE(whole_object(n)%segment(1)%iend - whole_object(n)%segment(1)%istart)
        END DO
        whole_object(n)%phi_profile(whole_object(n)%segment(1)%istart) = whole_object(whole_object(n)%n_connected_to_start)%phi
        whole_object(n)%phi_profile(whole_object(n)%segment(1)%iend)   = whole_object(whole_object(n)%n_connected_to_end)%phi
     END IF
  END DO

END SUBROUTINE UPDATE_WALL_POTENTIALS

!-------------------------------------
!
REAL(8) FUNCTION Harmonic_Oscillation(nwo, T_cntr)

  USE CurrentProblemValues, ONLY : whole_object, pi

  IMPLICIT NONE

  INTEGER nwo        ! number of the whole object
  INTEGER T_cntr

  INTEGER ap_period
  INTEGER n_periods

  Harmonic_Oscillation = 0.0_8

  IF (whole_object(nwo)%use_amplitude_profile) THEN

     ap_period = whole_object(nwo)%ap_T_cntr(whole_object(nwo)%N_ap_points)
     IF (T_cntr.EQ.(INT(T_cntr / ap_period) * ap_period)) THEN

!whole_object(n)%omega * T_cntr + whole_object(n)%phase_adjusted = integer * 2 * pi + whole_object(n)%phase

        n_periods = INT((whole_object(nwo)%omega * T_cntr) / (2.0_8 * pi))
        whole_object(nwo)%phase_adjusted = n_periods * 2.0_8 * pi + whole_object(nwo)%phase - whole_object(nwo)%omega * T_cntr
     END IF

     Harmonic_Oscillation = whole_object(nwo)%phi_var * SIN(whole_object(nwo)%omega * T_cntr + whole_object(nwo)%phase_adjusted)

  ELSE

     Harmonic_Oscillation = whole_object(nwo)%phi_var * SIN(whole_object(nwo)%omega * T_cntr + whole_object(nwo)%phase)

  END IF

END FUNCTION Harmonic_Oscillation

!-------------------------------------
!
REAL(8) FUNCTION Custom_Waveform(nwo, T_cntr)

  USE CurrentProblemValues, ONLY : whole_object

  IMPLICIT NONE

  INTEGER nwo        ! number of the whole object
  INTEGER T_cntr

  INTEGER wf_period
  INTEGER my_wf_T_cntr
  INTEGER i
  REAL ai

  Custom_Waveform = 0.0_8

  IF (.NOT.whole_object(nwo)%use_waveform) RETURN

  wf_period = whole_object(nwo)%wf_T_cntr(whole_object(nwo)%N_wf_points)

  my_wf_T_cntr = T_cntr - INT(T_cntr / wf_period) * wf_period

  DO i = 1, whole_object(nwo)%N_wf_points-1
     IF ( (my_wf_T_cntr.GE.whole_object(nwo)%wf_T_cntr(i)).AND.(my_wf_T_cntr.LE.whole_object(nwo)%wf_T_cntr(i+1)) ) THEN

        ai = DBLE(whole_object(nwo)%wf_T_cntr(i+1) - my_wf_T_cntr) / DBLE(whole_object(nwo)%wf_T_cntr(i+1) - whole_object(nwo)%wf_T_cntr(i))

        Custom_Waveform = DBLE(whole_object(nwo)%wf_phi(i) * ai + whole_object(nwo)%wf_phi(i+1) * (1.0_8 - ai))

        RETURN
     END IF
  END DO

END FUNCTION Custom_Waveform
  
!-------------------------------------
!
REAL(8) FUNCTION Amplitude_Profile(nwo, T_cntr)

  USE CurrentProblemValues, ONLY : whole_object

  IMPLICIT NONE

  INTEGER nwo        ! number of the whole object
  INTEGER T_cntr

  INTEGER ap_period
  INTEGER my_ap_T_cntr
  INTEGER i
  REAL ai

  Amplitude_Profile = 1.0_8

  IF (.NOT.whole_object(nwo)%use_amplitude_profile) RETURN

  ap_period = whole_object(nwo)%ap_T_cntr(whole_object(nwo)%N_ap_points)

  my_ap_T_cntr = T_cntr - INT(T_cntr / ap_period) * ap_period

  DO i = 1, whole_object(nwo)%N_ap_points-1
     IF ( (my_ap_T_cntr.GE.whole_object(nwo)%ap_T_cntr(i)).AND.(my_ap_T_cntr.LE.whole_object(nwo)%ap_T_cntr(i+1)) ) THEN

        ai = DBLE(whole_object(nwo)%ap_T_cntr(i+1) - my_ap_T_cntr) / DBLE(whole_object(nwo)%ap_T_cntr(i+1) - whole_object(nwo)%ap_T_cntr(i))

        Amplitude_Profile = DBLE(whole_object(nwo)%ap_factor(i) * ai + whole_object(nwo)%ap_factor(i+1) * (1.0_8 - ai))

        RETURN
     END IF
  END DO

END FUNCTION Amplitude_Profile

!-------------------------------------------
! This subroutine adjusts phase of harmonic oscillations of voltage applied to metal boundary objects. 
! The phase adjustment is necessary when the harmonically oscillating voltage is applied in pulses (that is Amplitude_Profile is not always 1).
! This subroutine does what is normally executed in the beginning of each pulse by the Harmonic_Oscillation function.
! Note that when the simulation is restarted from a checkpoint, the beginning of the pulse may be missed,
! and the proper value of the phase will not be specified.
!
SUBROUTINE ADJUST_HARMONIC_OSCILLATIONS_PHASE

  USE CurrentProblemValues, ONLY : whole_object, pi, Start_T_cntr, METAL_WALL, N_of_boundary_and_inner_objects

  IMPLICIT NONE

  INTEGER n        ! number of the whole object
  INTEGER my_T_cntr

  INTEGER ap_period
  INTEGER n_periods

  IF (Start_T_cntr.EQ.0) RETURN

! we are here if the checkpoint was requested and Start_T_cntr is not zero

  DO n = 1, N_of_boundary_and_inner_objects
     IF (whole_object(n)%object_type.NE.METAL_WALL) CYCLE   
     IF (whole_object(n)%potential_must_be_solved) CYCLE
     IF (whole_object(n)%use_amplitude_profile) THEN

        ap_period = whole_object(n)%ap_T_cntr(whole_object(n)%N_ap_points)   ! pulse period

        my_T_cntr = INT(Start_T_cntr / ap_period) * ap_period
        n_periods = INT((whole_object(n)%omega * my_T_cntr) / (2.0_8 * pi))
 
        whole_object(n)%phase_adjusted = n_periods * 2.0_8 * pi + whole_object(n)%phase - whole_object(n)%omega * my_T_cntr
 
     END IF
  END DO

END SUBROUTINE ADJUST_HARMONIC_OSCILLATIONS_PHASE
  
