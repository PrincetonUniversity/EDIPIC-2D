
!------------------------------------------------
!
REAL(8) FUNCTION frequency_of_en_collision(energy_eV, indx_neutral, colproc_id)

  USE MCCollisions
  USE CurrentProblemValues, ONLY : V_scale_ms, m_e_kg, e_Cl, N_subcycles, delta_t_s 

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: energy_eV
  INTEGER, INTENT(IN) :: indx_neutral
  INTEGER, INTENT(IN) :: colproc_id

  INTEGER N_crsect_points
  REAL(8) f_temp
  INTEGER j
  REAL(8) energy_j_eV, energy_jp1_eV, f_j, f_jp1

  N_crsect_points = neutral(indx_neutral)%en_colproc(colproc_id)%N_crsect_points

  IF (energy_eV.GE.neutral(indx_neutral)%en_colproc(colproc_id)%energy_eV(N_crsect_points)) THEN

     f_temp = neutral(indx_neutral)%en_colproc(colproc_id)%crsect_m2(N_crsect_points) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)

  ELSE IF (energy_eV.LT.neutral(indx_neutral)%en_colproc(colproc_id)%energy_eV(1)) THEN

     IF (neutral(indx_neutral)%en_colproc(colproc_id)%type.LT.20) THEN
! for elastic collisions only
        f_temp = neutral(indx_neutral)%en_colproc(colproc_id)%crsect_m2(1) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
     ELSE
! no inelastic and ionization collisions if energy below threshold
        f_temp = 0.0_8
     END IF

  ELSE
     
     DO j = 1, N_crsect_points-1
        IF ( (energy_eV.GE.neutral(indx_neutral)%en_colproc(colproc_id)%energy_eV(j)).AND. &
           & (energy_eV.LT.neutral(indx_neutral)%en_colproc(colproc_id)%energy_eV(j+1)) ) EXIT
     END DO

     j = MIN(j, N_crsect_points-1)

     energy_j_eV   = neutral(indx_neutral)%en_colproc(colproc_id)%energy_eV(j)
     energy_jp1_eV = neutral(indx_neutral)%en_colproc(colproc_id)%energy_eV(j+1)

     f_j   = neutral(indx_neutral)%en_colproc(colproc_id)%crsect_m2(j)   * SQRT(2.0_8 * energy_j_eV   * e_Cl / m_e_kg) 
     f_jp1 = neutral(indx_neutral)%en_colproc(colproc_id)%crsect_m2(j+1) * SQRT(2.0_8 * energy_jp1_eV * e_Cl / m_e_kg)

     f_temp = f_j + (f_jp1 - f_j) * (energy_eV - energy_j_eV) / (energy_jp1_eV - energy_j_eV)

  END IF

  frequency_of_en_collision = f_temp * neutral(indx_neutral)%N_m3 * N_subcycles * delta_t_s 

END FUNCTION frequency_of_en_collision

!-------------------------------------------------------------------------
!
SUBROUTINE en_Collision_Elastic_10(indx_neutral, indx_particle, energy_eV, counter)

  USE MCCollisions
  USE ElectronParticles
!???E CurrentProblemValues, ONLY : 
!???  USE ParallelOperationValues
!???  USE Diagnostics, ONLY : Rate_energy_coll
!???  USE ElectronInjection, ONLY : UseSmartTagsFlag

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER indx_neutral          ! ordering numer of neutral species
  INTEGER indx_particle         ! ordering number of particle 
  REAL(8) energy_eV             ! particle energy [eV]
  INTEGER counter               ! counter of collisions

  REAL(8) R                ! random number             
  REAL(8) Ksi              ! scattering angle (relative to the initial direction)
  REAL(8) CosKsi, SinKsi     
  REAL(8) Fi               ! azimuthal scattering angle 
  REAL(8) CosFi, SinFi     
  REAL(8) Vx, Vy, Vz       ! velocity components, before scattering
  REAL(8) Vx_s, Vy_s, Vz_s ! velocity components, after scattering
  REAL(8) V, V_xy, a, b            
  REAL(8) delta_energy     ! electron energy drop due to collision
  REAL(8) alpha            ! coefficient, accounting the electron energy drop

!  REAL(8) energy_change    ! change of energy of colliding electron

!  IF (energy_eV.LT.0.0_8) THEN 
!     PRINT '(/2x,"Process ",i3," : Potential ERROR in CollideElectron_1 (elastic e-n collisions):")', Rank_of_process
!     PRINT  '(2x,"Particle energy ",f10.3,"(eV) is EXACT ZERO")', energy_eV
!     PRINT  '(2x,"******* ******* THIS COLLISION EVENT WILL BE SKIPPED ******* *******")'
!     e_n_1_count = e_n_1_count - 1
!     RETURN
!  END IF

!?????? instead of a mild safety rule above, let us not scatter particles with energy below, say, 0.000001 eV
!  IF (energy_eV.LT.1.0d-6) THEN
!     e_n_1_count = e_n_1_count - 1
!     RETURN
!  END IF

! Calculate the scattering angle relative to the initial direction of electron
  R = well_random_number()

! #####  CosKsi = (2.0_8 + energy_eV - 2.0_8 * (1.0_8 + energy_eV)**R) / energy_eV #####
! the formula above was in the older code and it was based on Surendra's differential cross section
! below is the corrected expression from Okhrimovsky et al., Phys.Rev.E, 65, 037402 (2002).
  CosKsi = 1.0_8 - 2.0_8 * R / (1.0_8 + 8.0_8 * (energy_eV / 27.21_8) * (1.0_8 - R))

  CosKsi = MAX(MIN(0.999999999999_8, CosKsi), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi|>1
  Ksi = ACOS(CosKsi)
  SinKsi = SIN(Ksi)
! Calculate the azimuthal scattering angle
  R = well_random_number()
  Fi = R * 6.28318530718_8
  CosFi = COS(Fi)
  SinFi = SIN(Fi)
! Take the velocity 
  Vx = electron(indx_particle)%VX !VX_of_spec(1)%part(num)
  Vy = electron(indx_particle)%VY !VY_of_spec(1)%part(num)
  Vz = electron(indx_particle)%VZ !VZ_of_spec(1)%part(num)
! Turn the velocity  
  V    = SQRT(Vx*Vx + Vy*Vy + Vz*Vz) 
  V_xy = SQRT(Vx*Vx + Vy*Vy)
!  IF (V_xy.GT.0.0_8) THEN                                  ! was like this, changed for a new version (below)
!     a    = SinKsi * SinFi * V / V_xy                      ! to avoid possible error when Vx, Vy are very small but non-zero
!     b    = SinKsi * CosFi * Vz / V_xy                     ! (the error may ??? appear in terms with V_x,y/V_xy)
!     Vx_s = Vx * CosKsi + Vy * a + Vx * b                  !
!     Vy_s = Vy * CosKsi - Vx * a + Vy * b                  !
!     Vz_s = Vz * CosKsi          - V_xy * SinKsi * CosFi   !
  IF (V_xy.GT.1.0d-20) THEN
     a = Vx / V_xy
     b = Vy / V_xy
     Vx_s = Vx * CosKsi + (SinFi * V * b + CosFi * Vz * a) * SinKsi
     Vy_s = Vy * CosKsi - (SinFi * V * a - CosFi * Vz * b) * SinKsi
     Vz_s = Vz * CosKsi - V_xy * CosFi * SinKsi
  ELSE
     Vx_s = ABS(Vz) * SinKsi * CosFi
     Vy_s = ABS(Vz) * SinKsi * SinFi
     Vz_s = Vz * CosKsi
  END IF

!print *, Vx_s, Vy_s, Vz_s

! Calculate the energy drop [eV]
  delta_energy = (0.001097161_8 / neutral(indx_neutral)%M_amu) * (1.0_8 - CosKsi) ! * energy_eV 
!! 0.001097161 = 2 * m_e_kg / 1_amu_kg = 2 * 9.109534e-31 / 1.660565e-27 
!! ####### NOTE, VAHEDI'S EQ.(12) EITHER HAS A MISTAKE IN THE ABOVE EXPRESSION - MISSED " * ENERGY_EV"
!! ####### OR IS GIVEN ALREADY FOR THE RELATIVE ENERGY LOSS ...
!  alpha    = SQRT(1.0_8 - delta_eV / energy_eV)
  alpha = 1.0_8-delta_energy
! IF (alpha.LT.0.0_8) THEN
!    PRINT '(/2x,"Process ",i3," : ERROR in CollideElectron_1 (elastic e-n collisions):")', Rank_of_process
!    PRINT  '(2x,"scattering angle cosine: CosKsi = ",e16.9)', CosKsi
!    PRINT  '(2x,"factor for energy drop begore sqrt (should be positive) alpha = ",e16.9)', alpha
!    PRINT  '(2x,"terminating the program")'
!    STOP
! END IF
  alpha = SQRT(alpha)
! Renormalize the velocity in order to account the energy drop
  Vx_s = Vx_s * alpha
  Vy_s = Vy_s * alpha
  Vz_s = Vz_s * alpha
  electron(indx_particle)%VX = Vx_s !  VX_of_spec(1)%part(num) = Vx_s
  electron(indx_particle)%VY = Vy_s !  VY_of_spec(1)%part(num) = Vy_s
  electron(indx_particle)%VZ = Vz_s !  VZ_of_spec(1)%part(num) = Vz_s

  counter = counter + 1

!??? Mark the electron as the one collided with a neutral atom ONLY if it changed the direction of its velocity
!  IF (((Vx*Vx_s).LE.0.0_8).AND.(UseSmartTagsFlag.EQ.1)) THEN
!     species(1)%part(num)%Tag = 0 ! Tag_of_spec(1)%part(num) = 0 !eTag_Coll_Neutral
!  END IF

!??? calculate and save the change of energy
!  energy_change = (Vx_s*Vx_s + Vy_s*Vy_s + Vz_s*Vz_s) - (Vx*Vx + Vy*Vy + Vz*Vz)
!  Rate_energy_coll(1) = Rate_energy_coll(1) + energy_change

END SUBROUTINE en_Collision_Elastic_10

!-------------------------------------------------------------------------
!
SUBROUTINE en_Collision_Inelastic_20(indx_neutral, indx_particle, energy_eV, threshold_energy_eV, counter)

!???  USE MCCollisions
  USE ElectronParticles
!???  USE CurrentProblemValues
!???  USE ParallelOperationValues
!???  USE Diagnostics, ONLY : Rate_energy_coll
!???  USE ElectronInjection, ONLY : UseSmartTagsFlag

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER indx_neutral          ! ordering numer of neutral species
  INTEGER indx_particle         ! ordering number of particle 
  REAL(8) energy_eV             ! particle energy [eV]
  REAL(8) threshold_energy_eV   ! threshold energy for the collision to occur [eV]
  INTEGER counter               ! counter of collisions

  REAL(8) energy_sc_eV     ! particle energy after scattering [eV]
  REAL(8) R                ! random number             
  REAL(8) Ksi              ! scattering angle (relative to the initial direction)
  REAL(8) CosKsi, SinKsi     
  REAL(8) Fi               ! azimuthal scattering angle 
  REAL(8) CosFi, SinFi     
  REAL(8) Vx, Vy, Vz       ! velocity components, before scattering
  REAL(8) Vx_s, Vy_s, Vz_s ! velocity components, after scattering
  REAL(8) V, V_xy, a, b            
  REAL(8) alpha            ! coefficient, accounting the electron energy drop

!  REAL(8) energy_change    ! change of energy of colliding electron

!  REAL(8) threshold_eV                        ! the excitation threshold, Argon, approximate
!  threshold_eV = 11.55_8                      ! the excitation threshold, Argon, approximate

!print *, 'excitation, enter'

!  IF (energy_eV.LT.(Thresh_en_excit_eV)) THEN 
!     PRINT '(/2x,"Process ",i3," : Potential ERROR in CollideElectron_2 (excitation e-n collisions):")', Rank_of_process
!     PRINT  '(2x,"Particle energy ",f10.3,"(eV) is not more than the excitation threshold ",f6.2,"(eV)")', energy_eV, Thresh_en_excit_eV
!     PRINT  '(2x,"Such low energy particle cannot take part in this kind of collisions")'
!     PRINT  '(2x,"******* ******* THIS COLLISION EVENT WILL BE SKIPPED ******* *******")'
!     e_n_2_count = e_n_2_count - 1
!     RETURN
!!     PRINT  '(2x,"The program will be terminated now :(")'
!!     STOP
!  END IF

!??? instead of a mild safety rule above, let us not scatter particles with energy less than 0.000001 eV above the threshold
  IF (energy_eV.LE.threshold_energy_eV) THEN 
!     e_n_2_count = e_n_2_count - 1
     RETURN
  END IF

! Calculate the energy of the scattered electron 
  energy_sc_eV = MAX(0.0_8, energy_eV - threshold_energy_eV)

! Calculate the scattering angle relative to the initial direction of electron ! [Vahedi]: Use the modified energy "energy_sc_eV" here 
  R = well_random_number()
  
! ##### CosKsi = (2.0_8 + energy_sc_eV - 2.0_8 * (1.0_8 + energy_sc_eV)**R) / energy_sc_eV ##### 
  CosKsi = 1.0_8 - 2.0_8 * R / (1.0_8 + 8.0_8 * (energy_sc_eV / 27.21_8) * (1.0_8 - R))

  CosKsi = MAX(MIN(0.999999999999_8, CosKsi), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi|>1
  Ksi = ACOS(CosKsi)
  SinKsi = SIN(Ksi)
! Calculate the azimuthal scattering angle
  R = well_random_number()
  Fi = R * 6.28318530718_8
  CosFi = COS(Fi)
  SinFi = SIN(Fi)
! Take the velocity 
  Vx = electron(indx_particle)%VX !VX_of_spec(1)%part(num)
  Vy = electron(indx_particle)%VY !VY_of_spec(1)%part(num)
  Vz = electron(indx_particle)%VZ !VZ_of_spec(1)%part(num)
!print *, Vx, Vy, Vz
! Turn the velocity  
  V    = SQRT(Vx*Vx + Vy*Vy + Vz*Vz) 
  V_xy = SQRT(Vx*Vx + Vy*Vy)
!  IF (V_xy.GT.0.0_8) THEN
!     a    = SinKsi * SinFi * V / V_xy
!     b    = SinKsi * CosFi * Vz / V_xy 
!     Vx_s = Vx * CosKsi + Vy * a + Vx * b
!     Vy_s = Vy * CosKsi - Vx * a + Vy * b
!     Vz_s = Vz * CosKsi          - V_xy * SinKsi * CosFi
  IF (V_xy.GT.1.0d-20) THEN
     a = Vx / V_xy
     b = Vy / V_xy
     Vx_s = Vx * CosKsi + (SinFi * V * b + CosFi * Vz * a) * SinKsi
     Vy_s = Vy * CosKsi - (SinFi * V * a - CosFi * Vz * b) * SinKsi
     Vz_s = Vz * CosKsi - V_xy * CosFi * SinKsi
  ELSE
     Vx_s = ABS(Vz) * SinKsi * CosFi
     Vy_s = ABS(Vz) * SinKsi * SinFi
     Vz_s = Vz * CosKsi
  END IF

! Calculate the energy drop 
!  alpha    = SQRT(1.0_8 - Thresh_en_excit_eV / energy_eV)
  alpha    = SQRT(energy_sc_eV / energy_eV)
! Renormalize the velocity in order to account the energy drop
  Vx_s = Vx_s * alpha    
  Vy_s = Vy_s * alpha    
  Vz_s = Vz_s * alpha    
  electron(indx_particle)%VX = Vx_s !  VX_of_spec(1)%part(num) = Vx_s 
  electron(indx_particle)%VY = Vy_s !  VY_of_spec(1)%part(num) = Vy_s 
  electron(indx_particle)%VZ = Vz_s !  VZ_of_spec(1)%part(num) = Vz_s 

  counter = counter + 1

!??? Mark the electron as the one collided with a neutral atom ONLY if it changed the direction of its velocity 
!  IF (((Vx*Vx_s).LE.0.0_8).AND.(UseSmartTagsFlag.EQ.1)) THEN
!     species(1)%part(num)%Tag = 0 ! Tag_of_spec(1)%part(num) = 0 !eTag_Coll_Neutral
!  END IF

!??? calculate and save the change of energy
!  energy_change = (Vx_s*Vx_s + Vy_s*Vy_s + Vz_s*Vz_s) - (Vx*Vx + Vy*Vy + Vz*Vz)
!  Rate_energy_coll(1) = Rate_energy_coll(1) + energy_change

!print *, 'excitation, exit'

END SUBROUTINE en_Collision_Inelastic_20

!-------------------------------------------------------------------------
!
SUBROUTINE en_Collision_Ionization_30(indx_neutral, indx_particle, energy_inc_eV, threshold_energy_eV, ion_species_produced, ion_velocity_factor, counter)

!  USE MCCollisions
  USE ElectronParticles
  USE SetupValues, ONLY : factor_convert_vion_i
!???  USE CurrentProblemValues
!  USE ParallelOperationValues
!???  USE Diagnostics, ONLY : Rate_energy_coll
!???  USE ElectronInjection, ONLY : UseSmartTagsFlag

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER indx_neutral          ! ordering numer of neutral species
  INTEGER indx_particle         ! ordering number of particle 
  REAL(8) energy_inc_eV         ! particle energy [eV]
  REAL(8) threshold_energy_eV   ! threshold energy for the collision to occur [eV]
  INTEGER ion_species_produced  ! index of ion species produced (1<N_spec)
  REAL(8) ion_velocity_factor   ! factor to get ion velocity, dependes on specific neutral temperature
  INTEGER counter               ! counter of collisions

  REAL(8) B_of_Einc_eV                       ! some known function of incident electron energy

  REAL(8) energy_sc_eV     ! energy of incident electron after scattering [eV]
  REAL(8) Ksi_sc
  REAL(8) CosKsi_sc, SinKsi_sc     
  REAL(8) Fi_sc
  REAL(8) CosFi_sc, SinFi_sc     

  REAL(8) energy_ej_eV     ! energy of ejected electron [eV] 
  REAL(8) Ksi_ej
  REAL(8) CosKsi_ej, SinKsi_ej     
  REAL(8) Fi_ej
  REAL(8) CosFi_ej, SinFi_ej     

  REAL(8) V, V_xy, a, b            

  REAL(8) Vx, Vy, Vz          ! velocity components of incident electron, before scattering
  REAL(8) Vx_sc, Vy_sc, Vz_sc ! velocity components of incident electron, after scattering
  REAL(8) Vx_ej, Vy_ej, Vz_ej ! velocity components of ejected electron
  REAL(8) Vx_i, Vy_i, Vz_i    ! velocity components of produced ion

  REAL(8) alpha            

  REAL(8)   R                ! random number             

  REAL(8) energy_change_e    ! change of kinetic energy of electrons, incident and ejected
  REAL(8) energy_change_i    ! change of kinetic energy due to produced ion

  INTEGER left_node, right_node

!  REAL threshold_eV                       ! the ionization threshold, Argon, approximate
  B_of_Einc_eV = 10.0_8                   ! is constant for Einc < 70 eV, Argon, approximate
!############################################################################################      ##    ##
!###################### NOTE: MUST BE: "B_of_Einc_eV = 8.7" ###################################      ##    ##
!############################################################################################      ##    ##
!  threshold_eV = 15.0_8                   ! the ionization threshold, Argon, approximate

!print *, 'ionization'

!  IF (energy_inc_eV.LE.Thresh_en_ioniz_eV) THEN 
!     PRINT '(/2x,"Process ",i3," : Potential ERROR in CollideElectron_3 (ionization e-n collision):")', Rank_of_process
!     PRINT  '(2x,"Particle energy ",f10.3,"(eV) is not more than the ionization threshold ",f6.2,"(eV)")',&
!                                                                    & energy_inc_eV, Thresh_en_ioniz_eV
!     PRINT  '(2x,"Such low energy particle cannot take part in this kind of collisions")'
!     PRINT  '(2x,"******* ******* THIS COLLISION EVENT WILL BE SKIPPED ******* *******")'
!     e_n_3_count = e_n_3_count - 1
!     RETURN
!!     PRINT  '(2x,"The program will be terminated now :(")'
!!     STOP
!  END IF

!  IF (energy_inc_eV.LT.(Thresh_en_ioniz_eV+1.0d-6)) THEN 
  IF (energy_inc_eV.LE.threshold_energy_eV) THEN 
!     e_n_3_count = e_n_3_count - 1
     RETURN
  END IF

! Calculate the energy of the ejected electron
  R = well_random_number()
  energy_ej_eV = B_of_Einc_eV * TAN(R * ATAN( 0.5_8 * (energy_inc_eV - threshold_energy_eV) / B_of_Einc_eV ))  
! note, theoretically, for R<=1, the value of energy_ej_eV should never exceed (energy_inc_eV - Thresh_en_ioniz_eV)/2
! but below we shall enforce this manually

! I added the line blow to be sure that the energy of ejected electron is always
! (a) non-zero and (b) below half of energy difference [incident - threshold]
! note, 5d-7 = 1d-6/2, where 1d-6 is the minimal positive energy difference [incident - threshold]
  energy_ej_eV = MAX(0.0_8, MIN(energy_ej_eV, 0.5_8 * (energy_inc_eV - threshold_energy_eV)))

!  if (energy_ej_eV.gt.(0.5_8 * (energy_inc_eV - Thresh_en_ioniz_eV))) then
!     PRINT '(/2x,"Process ",i3," : Error in CollideElectron_3:")', Rank_of_process
!     PRINT  '(2x,"energy of  ",f10.3,"(eV) is greater than ",f6.2,"(eV)")', energy_ej_eV, &
!                                              & 0.5_8 * (energy_inc_eV - Thresh_en_ioniz_eV)
!     PRINT  '(2x,"The energy conservation law can be violated")'
!     PRINT  '(2x,"The program will be terminated now :(")'
!     STOP
!  END IF

! Calculate the energy of the scattered (incident) electron
  energy_sc_eV = MAX(0.0_8, energy_inc_eV - threshold_energy_eV - energy_ej_eV)

!##### it is expected that here neither energy_ej_eV nor energy_sc_eV are zeros (both should be nor less than 5d-7 eV) #####

!print *, energy_inc_eV, energy_sc_eV, energy_ej_eV

! Calculate the scattering angle Ksi for the incident electron ! [Vahedi]: use the modified energy "energy_sc_eV" here
  R = well_random_number()
!####  CosKsi_sc = (2.0_8 + energy_sc_eV - 2.0_8 * (1.0_8 + energy_sc_eV)**R) / energy_sc_eV ####
  CosKsi_sc = 1.0_8 - 2.0_8 * R / (1.0_8 + 8.0_8 * (energy_sc_eV / 27.21_8) * (1.0_8 - R))
  CosKsi_sc = MAX(MIN(0.999999999999_8, CosKsi_sc), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi_sc|>1
  Ksi_sc    = ACOS(CosKsi_sc)
  SinKsi_sc = SIN(Ksi_sc)

! Calculate the azimuthal scattering angle for the incident electron
  R = well_random_number()
  Fi_sc    = R * 6.28318530718_8
  CosFi_sc = COS(Fi_sc)
  SinFi_sc = SIN(Fi_sc)

! Calculate the scattering angle Ksi for the ejected electron
  R = well_random_number()
!####  CosKsi_ej = (2.0_8 + energy_ej_eV - 2.0_8 * (1.0_8 + energy_ej_eV)**R) / energy_ej_eV ####
  CosKsi_ej = 1.0_8 - 2.0_8 * R / (1.0_8 + 8.0_8 * (energy_ej_eV / 27.21_8) * (1.0_8 - R))
  CosKsi_ej = MAX(MIN(0.999999999999_8, CosKsi_ej), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi_ej|>1
  Ksi_ej    = ACOS(CosKsi_ej)
  SinKsi_ej = SIN(Ksi_ej)

! Calculate the azimuthal scattering angle for the ejected electron
  R = well_random_number()
  Fi_ej    = R * 6.28318530718_8
  CosFi_ej = COS(Fi_ej)
  SinFi_ej = SIN(Fi_ej)

! Take the velocity of the incident electron before the scattering
  Vx = electron(indx_particle)%VX !VX_of_spec(1)%part(num)
  Vy = electron(indx_particle)%VY !VY_of_spec(1)%part(num)
  Vz = electron(indx_particle)%VZ !VZ_of_spec(1)%part(num)

! Precalculate   
  V    = SQRT(Vx*Vx + Vy*Vy + Vz*Vz) 
  V_xy = SQRT(Vx*Vx + Vy*Vy)

! Calculate the velocity of the INCIDENT ELECTRON after scattering (turn the velocity)
!  IF (V_xy.GT.0.0_8) THEN
!     a     = SinKsi_sc * SinFi_sc * V / V_xy
!     b     = SinKsi_sc * CosFi_sc * Vz / V_xy 
!     Vx_sc = Vx * CosKsi_sc + Vy * a + Vx * b
!     Vy_sc = Vy * CosKsi_sc - Vx * a + Vy * b
!     Vz_sc = Vz * CosKsi_sc - V_xy * SinKsi_sc * CosFi_sc
  IF (V_xy.GT.1.0d-20) THEN
     a = Vx / V_xy
     b = Vy / V_xy
     Vx_sc = Vx * CosKsi_sc + (SinFi_sc * V * b + CosFi_sc * Vz * a) * SinKsi_sc
     Vy_sc = Vy * CosKsi_sc - (SinFi_sc * V * a - CosFi_sc * Vz * b) * SinKsi_sc
     Vz_sc = Vz * CosKsi_sc - V_xy * CosFi_sc * SinKsi_sc
  ELSE
     Vx_sc = ABS(Vz) * SinKsi_sc * CosFi_sc
     Vy_sc = ABS(Vz) * SinKsi_sc * SinFi_sc
     Vz_sc = Vz * CosKsi_sc
  END IF

! Calculate the relative energy drop 
  alpha = SQRT(energy_sc_eV / energy_inc_eV)
! Renormalize the velocity of incident electron in order to account the energy drop
  Vx_sc = Vx_sc * alpha
  Vy_sc = Vy_sc * alpha
  Vz_sc = Vz_sc * alpha
  electron(indx_particle)%VX = Vx_sc !  VX_of_spec(1)%part(num) = Vx_sc
  electron(indx_particle)%VY = Vy_sc !  VY_of_spec(1)%part(num) = Vy_sc
  electron(indx_particle)%VZ = Vz_sc !  VZ_of_spec(1)%part(num) = Vz_sc

!??? Mark the incident electron as the one collided with a neutral atom ONLY if it changed the direction of its velocity
!  IF (((Vx*Vx_sc).LE.0.0_8).AND.(UseSmartTagsFlag.EQ.1)) THEN
!     species(1)%part(num)%Tag = 0 !Tag_of_spec(1)%part(num) = 0 !eTag_Coll_Neutral
!  END IF

! Calculate the velocity of the EJECTED ELECTRON after scattering (turn the velocity)
!  IF (V_xy.GT.0.0_8) THEN
!     a     = SinKsi_ej * SinFi_ej * V / V_xy
!     b     = SinKsi_ej * CosFi_ej * Vz / V_xy 
!     Vx_ej = Vx * CosKsi_ej + Vy * a + Vx * b
!     Vy_ej = Vy * CosKsi_ej - Vx * a + Vy * b
!     Vz_ej = Vz * CosKsi_ej          - V_xy * SinKsi_ej * CosFi_ej
  IF (V_xy.GT.1.0d-20) THEN
     a = Vx / V_xy
     b = Vy / V_xy
     Vx_ej = Vx * CosKsi_ej + (SinFi_ej * V * b + CosFi_ej * Vz * a) * SinKsi_ej
     Vy_ej = Vy * CosKsi_ej - (SinFi_ej * V * a - CosFi_ej * Vz * b) * SinKsi_ej
     Vz_ej = Vz * CosKsi_ej - V_xy * CosFi_ej * SinKsi_ej
  ELSE
     Vx_ej = ABS(Vz) * SinKsi_ej * CosFi_ej
     Vy_ej = ABS(Vz) * SinKsi_ej * SinFi_ej
     Vz_ej = Vz * CosKsi_ej
  END IF

!print *, Vx**2 + Vy**2 + Vz**2, Vx_ej**2 + Vy_ej**2 + Vz_ej**2   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate the relative energy drop 
  alpha = SQRT(energy_ej_eV / energy_inc_eV)
! Renormalize the velocity of EJECTED ELECTRON in order to account the energy drop
  Vx_ej = Vx_ej * alpha
  Vy_ej = Vy_ej * alpha
  Vz_ej = Vz_ej * alpha

!print *, energy_inc_eV, ((Vx_sc**2+Vy_sc**2+Vz_sc**2)*T_e_eV/N_box_vel**2)/energy_sc_eV, ((Vx_ej**2+Vy_ej**2+Vz_ej**2)*T_e_eV/N_box_vel**2)/energy_ej_eV
!print *, "en_Collision_Ionization_30 :: before ADD_ELECTRON_TO_ADD_LIST"
  CALL ADD_ELECTRON_TO_ADD_LIST(electron(indx_particle)%X, electron(indx_particle)%Y, Vx_ej, Vy_ej, Vz_ej, 0)   ! tag=0
!print *, "en_Collision_Ionization_30 :: after ADD_ELECTRON_TO_ADD_LIST"
!??? diagnostics
!     N_inject(1) = N_inject(1) + 1


!??? calculate and save the change of energy for electrons
!  energy_change_e = (Vx_ej**2 + Vy_ej**2 + Vz_ej**2) + &    ! energy of ejected electron PLUS
!                  & (Vx_sc**2 + Vy_sc**2 + Vz_sc**2) - &    ! energy of incident electron AFTER collision (scattered) MINUS
!                  & (Vx**2    + Vy**2    + Vz**2)           ! energy of incident electron BEFORE collision
!  Rate_energy_coll(1) = Rate_energy_coll(1) + energy_change_e

! Take random velocity from normalized maxwell distribution  
  CALL GetMaxwellVelocity(Vx_i) 
  CALL GetMaxwellVelocity(Vy_i) 
  CALL GetMaxwellVelocity(Vz_i) 

! Use the factor above to obtain the dim-less velocity (V * N_box_vel / V_te) of the produced ion
  Vx_i = Vx_i * ion_velocity_factor !factor_convert_vion_i(ion_species_produced)  !alpha_Vscl
  Vy_i = Vy_i * ion_velocity_factor !factor_convert_vion_i(ion_species_produced)  !alpha_Vscl
  Vz_i = Vz_i * ion_velocity_factor !factor_convert_vion_i(ion_species_produced)  !alpha_Vscl

!print *, "en_Collision_Ionization_30 :: before ADD_ION_TO_ADD_LIST"

  CALL ADD_ION_TO_ADD_LIST(ion_species_produced, electron(indx_particle)%X, electron(indx_particle)%Y, Vx_i, Vy_i, Vz_i, 0)   ! tag=0
!??? diagnostics
!print *, "en_Collision_Ionization_30 :: after ADD_ION_TO_ADD_LIST"

  counter = counter + 1

!??? calculate and save the change of energy for ion
!  energy_change_i = Vx_i*Vx_i + Vy_i*Vy_i + Vz_i*Vz_i        
!  Rate_energy_coll(2) = Rate_energy_coll(2) + energy_change_i

END SUBROUTINE en_Collision_Ionization_30
!=====================================================================================================
