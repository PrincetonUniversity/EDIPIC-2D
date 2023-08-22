
!-------------------------------------------------------------------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION(x, y, vx, vy, vz, tag, myobject, m, dirflag)

  USE ParallelOperationValues
  USE CurrentProblemValues

  USE rng_wrapper

  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr

  REAL(8) x, y         ! primary electron coordinates
  REAL(8) vx, vy, vz   ! primary electron velocities
  INTEGER tag          ! primary electron tag 
!  INTEGER nwo          ! number of the whole object with which the primary electron collided
  TYPE(boundary_object) myobject
  INTEGER m            ! number of the local [belonging to this particular cluster] part of the whole object with which the primary electron collided
  
  INTEGER dirflag      ! direction flag, here left/right/below/above are -x/+x/-y/+y directions
                       ! 1 = collision with wall on the left
                       ! 2 = collision with wall above
                       ! 3 = collision with wall on the right
                       ! 4 = collision with wall below

  REAL(8) energy_inc      ! energy of primary electron
  REAL(8) v               ! speed of primary electron
  REAL(8) theta_inc       ! angle of incidence with respect to wall normal 

  REAL(8) coef_elastic    !
  REAL(8) coef_inelastic  ! corresponding emission coefficients at current energy and angle of incidence
  REAL(8) coef_true       !
  REAL(8) coef_total      !

  REAL(8) R               ! random number 

! functions
  REAL(8) Coeff_SEE_Elastic
  REAL(8) Coeff_SEE_Inelastic
  REAL(8) Coeff_SEE_True

! calculate the dim-less energy of incident electron
  energy_inc = vx**2 + vy**2 + vz**2
  IF (energy_inc.LE.myobject%lowest_energy_for_see) RETURN

! calculate the angle of incidence
  v = SQRT(energy_inc)

  IF (v.EQ.0.0_8) THEN
! we should not be here, but just in case
     theta_inc = 0.0_8
  ELSE
     IF (dirflag.EQ.1) THEN
        theta_inc = ACOS(MIN(1.0_8, ABS(vx)/v)) 
     ELSE IF (dirflag.EQ.2) THEN
        theta_inc = ACOS(MIN(1.0_8, ABS(vy)/v))
     ELSE IF (dirflag.EQ.3) THEN
        theta_inc = ACOS(MIN(1.0_8, ABS(vx)/v))
     ELSE IF (dirflag.EQ.4) THEN
        theta_inc = ACOS(MIN(1.0_8, ABS(vy)/v))
     END IF
  END IF

! calculate the coefficients of elastic/inelastic scattering and true secondary emission
  coef_elastic   = Coeff_SEE_Elastic(energy_inc, theta_inc, myobject)
  coef_inelastic = Coeff_SEE_Inelastic(energy_inc, theta_inc, myobject)
  IF ((coef_elastic + coef_inelastic).GT.1.0_8) THEN
     PRINT '("Process ",i3," : Error in secondary emission model!")', Rank_of_process
     PRINT '("The total coefficient of elastic/inelastic backscattering is greater than 1 !!!")'
     PRINT '("elastic: ",f5.2," inelastic: ",f5.2)', coef_elastic, coef_inelastic
     PRINT '("Program will be terminated now :(")'
     errcode=370
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF  
  coef_true = Coeff_SEE_True(energy_inc, theta_inc, myobject)

! then we combine the rest (fractional part) of the true secondary emission coefficient 
! with the coefficients of elastic/inelastic backscattering to complete the total emission
  coef_total = coef_elastic + coef_inelastic + coef_true

! quit subroutine if SEE is already exhausted / turned off
  IF (coef_total.EQ.0.0_8) RETURN

! We can always treat coefficients (relative yields) of elastic/inelastic backscattering as probabilities.
! If the total coefficient is less or equal than 1, we simply treat ALL coefficients as probabilities.
! If at this stage the total coefficient is greater than 1, this means that, on one hand, 
! not all electrons will produce true secondary electrons, 
! on the other hand, those incident electrons, which produce the true secondary electrons, 
! will (in average) produce more than 1 secondary electron per 1 incident electron. 
! We cannot treat the true secondary coefficient as the probability yet,
! (the difference from the other coefficients is in that only true secondary emission 
! can produce more than 1 electron per incident electron).
! Also, we cannot start here with the injection of the true secondaries, because afterwards it is 
! possible that the elastically or inelastically reflected particle will be injected.
! Therefore we must start from the attempts to inject a reflected (elastically or inelastically) particle.

! take a random number
  R = well_random_number()

! try  to inject the elastically reflected electron
  IF (R.LE.coef_elastic) THEN
     CALL INJECT_ELASTIC_REFLECTED_ELECTRON(x, y, vx, vy, vz, v, tag, myobject, m, dirflag)       
     RETURN
  END IF

! try to inject the inelasticaly backscattered electron
  IF ((R.GT.coef_elastic).AND.(R.LE.(coef_elastic + coef_inelastic))) THEN
     CALL INJECT_INELASTIC_BACKSCATTERED_ELECTRON(x, y, v, tag, myobject, m, dirflag)
     RETURN
  END IF

  IF (coef_true.LE.0.0_8) RETURN   ! this should take care of myobject%Emitted_model(3)=0

! we want to have some safety mechanism to prevent emission of more particles than the energy of the primary electron can afford

! reduce the primary electron energy by the threshold for the true secondary electron emission
  IF (myobject%Emitted_model(3).EQ.1) THEN
     energy_inc = MAX(0.0_8, energy_inc - myobject%minE_see_true)
  ELSE IF (myobject%Emitted_model(3).EQ.2) THEN
     energy_inc = MAX(0.0_8, energy_inc - myobject%E_see_0)
  END IF  

  IF (coef_total.LE.1.0_8) THEN   ! in this case we can treat the true secondary yield as the probability
! try to inject the true secondary electron
     IF ((R.GT.(coef_elastic + coef_inelastic)).AND.(R.LE.coef_total)) THEN
        CALL INJECT_TRUE_SECONDARY_ELECTRON(x, y, energy_inc, tag, myobject, m, dirflag)   ! this call reduces energy_inc
     END IF
! if nothing was injected (R.GT.coef_total), the particle is attached to the wall and we leave the subroutine
     RETURN
  END IF

! if we are here the total SEE yield is greater than 1 (coef_total.GT.1) and at least one true secondary electron must be injected
! we calculate the ratio of the emitted true secondary flux 
! to the portion of the incident flux, which does not produce elastically/inelastically reflected electrons:
  coef_true = coef_true / (1.0_8 - (coef_elastic + coef_inelastic))   ! note, we cannot be here if coef_elastic + coef_inelastic = 1

! then we inject true secondary electrons until we exhaust the integer part of the ratio
  DO WHILE (coef_true.GE.1.0_8)
     CALL INJECT_TRUE_SECONDARY_ELECTRON(x, y, energy_inc, tag, myobject, m, dirflag)   ! this call reduces energy_inc
     coef_true = coef_true - 1.0_8
  END DO
! the remaining part gives us the probability for statistical injection of fractional part of true secondary electron

! take a random number
  R = well_random_number()
  IF (R.LT.coef_true) CALL INJECT_TRUE_SECONDARY_ELECTRON(x, y, energy_inc, tag, myobject, m, dirflag)   ! this call reduces energy_inc

END SUBROUTINE PROCESS_ELECTRON_INDUCED_ELECTRON_EMISSION

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE INJECT_ELASTIC_REFLECTED_ELECTRON(x, y, vx, vy, vz, v, tag, myobject, m, dirflag)

!  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues

  USE rng_wrapper

  use mpi

  IMPLICIT NONE

  REAL(8) x, y         ! primary electron coordinates
  REAL(8) vx, vy, vz   ! primary electron velocity components
  REAL(8) v            ! absolute value of the primary electron velocity
  INTEGER tag          ! primary electron tag 
  
!  INTEGER nwo          ! number of the whole object with which the primary electron collided
  TYPE(boundary_object) myobject
  INTEGER m            ! number of the local [belonging to this particular cluster] part of the whole object with which the primary electron collided
  
  INTEGER dirflag      ! direction flag, here left/right/below/above are -x/+x/-y/+y directions
                       ! 1 = collision with wall on the left
                       ! 2 = collision with wall above
                       ! 3 = collision with wall on the right
                       ! 4 = collision with wall below

  REAL(8) x_new, y_new             ! scattered electron coordinates
  REAL(8) vx_new, vy_new, vz_new   ! scattered electron velocity components
  INTEGER tag_new                  ! scattered electron tag 

  INTEGER jbelow, jabove
  REAL(8) dqbelow, dqabove

  INTEGER ileft, iright
  REAL(8) dqleft, dqright

  REAL(8) theta, fi               ! scattering angles

! for inner object
  INTEGER i_left_top, i_right_top, i_right_bottom, i_left_bottom_bis, ip1
  INTEGER i
  REAL(8) dqi, dqip1

  IF (m.GT.0) THEN
! domain boundary 

     IF (myobject%Elast_refl_type.EQ.0) THEN    ! if elastic reflection occurs specularly

        IF (dirflag.EQ.1) THEN
! collision with the wall on the left
           x_new = DBLE(c_indx_x_min) + 1.0d-6
           y_new = y
           vx_new = -vx
           vy_new =  vy
           vz_new =  vz
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = y - jbelow
              dqbelow = 1.0_8 - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) + dqbelow   ! note the plus. remember ::
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) + dqabove   ! we subtract electron charge, it's negative
           END IF

        ELSE IF (dirflag.EQ.2) THEN
!  collision with the wall above
           x_new = x
           y_new = DBLE(c_indx_y_max) - 1.0d-6 
           vx_new =  vx
           vy_new = -vy
           vz_new =  vz
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = x - ileft
              dqleft = 1.0_8 - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft)   + dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) + dqright
           END IF

        ELSE IF (dirflag.EQ.3) THEN
! collision with the wall on the right
           x_new = DBLE(c_indx_x_max) - 1.0d-6
           y_new = y
           vx_new = -vx
           vy_new =  vy
           vz_new =  vz
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = y - jbelow
              dqbelow = 1.0_8 - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) + dqbelow   ! note the plus. remember ::
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) + dqabove   ! we subtract electron charge, it's negative
           END IF

        ELSE IF (dirflag.EQ.4) THEN
! collision with the wall below
           x_new = x
           y_new = DBLE(c_indx_y_min) + 1.0d-6 
           vx_new =  vx
           vy_new = -vy
           vz_new =  vz
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = x - ileft
              dqleft = 1.0_8 - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft)   + dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) + dqright
           END IF

        END IF   ! IF (dirflag.EQ.1) THEN

     ELSE                             ! if elastic reflection occurs at a random angle

! get the angles of reflection, you have to choose between the following:

! the distribution over angle theta is uniform 
!     theta = 1.5707963_8 * well_random_number()                 ! note that 1.5707963 is sligtly less than the exact pi/2=1.5707963267949

 ! the distribution over angle theta angle is f(theta) = COS(theta)   
        theta = ASIN(MIN(well_random_number(),1.0_8))              
        fi   = 6.283185307179_8 * well_random_number()

        IF (dirflag.EQ.1) THEN
! collision with the wall on the left
           x_new = DBLE(c_indx_x_min) + 1.0d-6
           y_new = y
           vx_new = v * COS(theta)
           vy_new = v * SIN(theta) * SIN(fi)
           vz_new = v * SIN(theta) * COS(fi)
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = y - jbelow
              dqbelow = 1.0_8 - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) + dqbelow   ! note the plus. remember ::
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) + dqabove   ! we subtract electron charge, it's negative
           END IF
        
        ELSE IF (dirflag.EQ.2) THEN
!  collision with the wall above
           x_new = x
           y_new = DBLE(c_indx_y_max) - 1.0d-6 
           vx_new =  v * SIN(theta) * SIN(fi) 
           vy_new = -v * COS(theta)
           vz_new =  v * SIN(theta) * COS(fi)
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = x - ileft
              dqleft = 1.0_8 - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft)   + dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) + dqright
           END IF

        ELSE IF (dirflag.EQ.3) THEN
! collision with the wall on the right
           x_new = DBLE(c_indx_x_max) - 1.0d-6
           y_new = y
           vx_new = -v * COS(theta)
           vy_new =  v * SIN(theta) * SIN(fi)
           vz_new =  v * SIN(theta) * COS(fi)
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
              jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
              dqabove = y - jbelow
              dqbelow = 1.0_8 - dqabove
              c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) + dqbelow   ! note the plus. remember ::
              c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) + dqabove   ! we subtract electron charge, it's negative
           END IF

        ELSE IF (dirflag.EQ.4) THEN
! collision with the wall below
           x_new = x
           y_new = DBLE(c_indx_y_min) + 1.0d-6 
           vx_new = v * SIN(theta) * SIN(fi) 
           vy_new = v * COS(theta)
           vz_new = v * SIN(theta) * COS(fi)
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              ileft = MAX(INT(x), c_local_object_part(m)%istart)
              iright = MIN(ileft + 1, c_local_object_part(m)%iend)
              dqright = x - ileft
              dqleft = 1.0_8 - dqright
              c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft)   + dqleft
              c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) + dqright
           END IF

        END IF      ! IF (dirflag.EQ.1) THEN
      
     END IF    ! IF (myobject%Elast_refl_type.EQ.0) THEN

  ELSE
! inner object

! ilt   -  ---- irt
!  |             |
!  |             |
!  1 ilbb ---- irb
!
     i_left_top        =                  myobject%jtop   - myobject%jbottom + 1
     i_right_top       = i_left_top     + myobject%iright - myobject%ileft
     i_right_bottom    = i_right_top    + myobject%jtop   - myobject%jbottom
     i_left_bottom_bis = i_right_bottom + myobject%iright - myobject%ileft - 1

     IF (myobject%Elast_refl_type.EQ.0) THEN    ! if elastic reflection occurs specularly

        IF (dirflag.EQ.1) THEN
! collision with the wall on the left
           x_new = myobject%xmax + 1.0d-6
           y_new = y
           vx_new = -vx
           vy_new =  vy
           vz_new =  vz
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              i = MIN(INT(myobject%ymax - y) + i_right_top, i_right_bottom - 1)
              dqi = y - INT(y)
              dqip1 = 1.0_8 - dqi
              myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi     ! note the plus. remember ::
              myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1   ! we subtract electron charge, it's negative
           END IF

        ELSE IF (dirflag.EQ.2) THEN
!  collision with the wall above
           x_new = x
           y_new = myobject%ymin - 1.0d-6 
           vx_new =  vx
           vy_new = -vy
           vz_new =  vz
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              i = MIN(INT(myobject%xmax - x) + i_right_bottom, i_left_bottom_bis)
              dqi = x - INT(x)
              dqip1 = 1.0_8 - dqi
              ip1 = i+1
              IF (i.EQ.i_left_bottom_bis) ip1 = 1
              myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
              myobject%surface_charge_variation(ip1) = myobject%surface_charge_variation(ip1) + dqip1
           END IF

        ELSE IF (dirflag.EQ.3) THEN
! collision with the wall on the right
           x_new = myobject%xmin - 1.0d-6
           y_new = y
           vx_new = -vx
           vy_new =  vy
           vz_new =  vz
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              i = MIN(INT(y - myobject%ymin) + 1, i_left_top - 1)
              dqip1 = y - INT(y)
              dqi = 1.0_8 - dqip1
              myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
              myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1
           END IF

        ELSE IF (dirflag.EQ.4) THEN
! collision with the wall below
           x_new = x
           y_new = myobject%ymax + 1.0d-6 
           vx_new =  vx
           vy_new = -vy
           vz_new =  vz
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              i = MIN(INT(x - myobject%xmin) + i_left_top, i_right_top - 1)
              dqip1 = x - INT(x)
              dqi = 1.0_8 - dqip1
              myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
              myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1
           END IF

        END IF   ! IF (dirflag.EQ.1) THEN

     ELSE                             ! if elastic reflection occurs at a random angle

! get the angles of reflection, you have to choose between the following:

! the distribution over angle theta is uniform 
!     theta = 1.5707963_8 * well_random_number()                 ! note that 1.5707963 is sligtly less than the exact pi/2=1.5707963267949

 ! the distribution over angle theta angle is f(theta) = COS(theta)   
        theta = ASIN(MIN(well_random_number(),1.0_8))              
        fi   = 6.283185307179_8 * well_random_number()

        IF (dirflag.EQ.1) THEN
! collision with the wall on the left
           x_new = myobject%xmax + 1.0d-6
           y_new = y
           vx_new = v * COS(theta)
           vy_new = v * SIN(theta) * SIN(fi)
           vz_new = v * SIN(theta) * COS(fi)
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              i = MIN(INT(myobject%ymax - y) + i_right_top, i_right_bottom - 1)
              dqi = y - INT(y)
              dqip1 = 1.0_8 - dqi
              myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi     ! note the plus. remember ::
              myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1   ! we subtract electron charge, it's negative
           END IF
        
        ELSE IF (dirflag.EQ.2) THEN
!  collision with the wall above
           x_new = x
           y_new = myobject%ymin - 1.0d-6 
           vx_new =  v * SIN(theta) * SIN(fi) 
           vy_new = -v * COS(theta)
           vz_new =  v * SIN(theta) * COS(fi)
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              i = MIN(INT(myobject%xmax - x) + i_right_bottom, i_left_bottom_bis)
              dqi = x - INT(x)
              dqip1 = 1.0_8 - dqi
              ip1 = i+1
              IF (i.EQ.i_left_bottom_bis) ip1 = 1
              myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
              myobject%surface_charge_variation(ip1) = myobject%surface_charge_variation(ip1) + dqip1
           END IF

        ELSE IF (dirflag.EQ.3) THEN
! collision with the wall on the right
           x_new = myobject%xmin - 1.0d-6
           y_new = y
           vx_new = -v * COS(theta)
           vy_new =  v * SIN(theta) * SIN(fi)
           vz_new =  v * SIN(theta) * COS(fi)
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              i = MIN(INT(y - myobject%ymin) + 1, i_left_top - 1)
              dqip1 = y - INT(y)
              dqi = 1.0_8 - dqip1
              myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
              myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1
            END IF

        ELSE IF (dirflag.EQ.4) THEN
! collision with the wall below
           x_new = x
           y_new = myobject%ymax + 1.0d-6 
           vx_new = v * SIN(theta) * SIN(fi) 
           vy_new = v * COS(theta)
           vz_new = v * SIN(theta) * COS(fi)
           tag_new = myobject%object_id_number

           IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
              i = MIN(INT(x - myobject%xmin) + i_left_top, i_right_top - 1)
              dqip1 = x - INT(x)
              dqi = 1.0_8 - dqip1
              myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
              myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1
           END IF

        END IF      ! IF (dirflag.EQ.1) THEN
      
     END IF    ! IF (myobject%Elast_refl_type.EQ.0) THEN

  END IF   ! IF (m.GT.0) THEN

  CALL ADD_ELECTRON_TO_ADD_LIST(x_new, y_new, vx_new, vy_new, vz_new, tag_new)    ! the particle gets tag equal to nwo - the number of the whole object that emitted it

  myobject%electron_emit_count = myobject%electron_emit_count + 1

END SUBROUTINE INJECT_ELASTIC_REFLECTED_ELECTRON

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE INJECT_INELASTIC_BACKSCATTERED_ELECTRON(x, y, v, tag, myobject, m, dirflag)

!  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues

  USE rng_wrapper

  use mpi

  IMPLICIT NONE

  REAL(8) x, y         ! primary electron coordinates
  REAL(8) v            ! absolute value of the primary electron velocity
  INTEGER tag          ! primary electron tag 
!  INTEGER nwo          ! number of the whole object with which the primary electron collided
  TYPE(boundary_object) myobject
  INTEGER m            ! number of the local [belonging to this particular cluster] part of the whole object with which the primary electron collided
  
  INTEGER dirflag      ! direction flag, here left/right/below/above are -x/+x/-y/+y directions
                       ! 1 = collision with wall on the left
                       ! 2 = collision with wall above
                       ! 3 = collision with wall on the right
                       ! 4 = collision with wall below

  REAL(8) theta, fi               ! scattering angles

  REAL(8) v_new                    ! scattered electron absolute velocity value
  REAL(8) x_new, y_new             ! scattered electron coordinates
  REAL(8) vx_new, vy_new, vz_new   ! scattered electron velocity components
  INTEGER tag_new                  ! scattered electron tag 

  INTEGER jbelow, jabove
  REAL(8) dqbelow, dqabove

  INTEGER ileft, iright
  REAL(8) dqleft, dqright

! for inner object
  INTEGER i_left_top, i_right_top, i_right_bottom, i_left_bottom_bis, ip1
  INTEGER i
  REAL(8) dqi, dqip1

! get the angles of reflection, you have to choose between the following:

! the distribution over angle theta is uniform 
!     theta = 1.5707963_8 * well_random_number()                 ! note that 1.5707963 is sligtly less than the exact pi/2=1.5707963267949

 ! the distribution over angle theta angle is f(theta) = COS(theta)   
  theta = ASIN(MIN(well_random_number(),1.0_8))              
  fi   = 6.283185307179_8 * well_random_number()

! get the absolute velocity of backscattered electron
  v_new = v * SQRT(well_random_number())           ! assume that the energy (NOT the velocity) is uniformly distributed between 0 and initial_energy

  IF (m.GT.0) THEN
! domain boundary 

     IF (dirflag.EQ.1) THEN
! collision with the wall on the left
        x_new = DBLE(c_indx_x_min) + 1.0d-6
        y_new = y
        vx_new = v_new * COS(theta)
        vy_new = v_new * SIN(theta) * SIN(fi)
        vz_new = v_new * SIN(theta) * COS(fi)
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
           jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
           dqabove = y - jbelow
           dqbelow = 1.0_8 - dqabove
           c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) + dqbelow   ! note the plus. remember ::
           c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) + dqabove   ! we subtract electron charge, it's negative
        END IF

     ELSE IF (dirflag.EQ.2) THEN
!  collision with the wall above
        x_new = x
        y_new = DBLE(c_indx_y_max) - 1.0d-6 
        vx_new =  v_new * SIN(theta) * SIN(fi) 
        vy_new = -v_new * COS(theta)
        vz_new =  v_new * SIN(theta) * COS(fi)
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           ileft = MAX(INT(x), c_local_object_part(m)%istart)
           iright = MIN(ileft + 1, c_local_object_part(m)%iend)
           dqright = x - ileft
           dqleft = 1.0_8 - dqright
           c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft)   + dqleft
           c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) + dqright
        END IF

     ELSE IF (dirflag.EQ.3) THEN
! collision with the wall on the right
        x_new = DBLE(c_indx_x_max) - 1.0d-6
        y_new = y
        vx_new = -v_new * COS(theta)
        vy_new =  v_new * SIN(theta) * SIN(fi)
        vz_new =  v_new * SIN(theta) * COS(fi)
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
           jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
           dqabove = y - jbelow
           dqbelow = 1.0_8 - dqabove
           c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) + dqbelow   ! note the plus. remember ::
           c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) + dqabove   ! we subtract electron charge, it's negative
        END IF

     ELSE IF (dirflag.EQ.4) THEN
! collision with the wall below
        x_new = x
        y_new = DBLE(c_indx_y_min) + 1.0d-6 
        vx_new = v_new * SIN(theta) * SIN(fi) 
        vy_new = v_new * COS(theta)
        vz_new = v_new * SIN(theta) * COS(fi)
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           ileft = MAX(INT(x), c_local_object_part(m)%istart)
           iright = MIN(ileft + 1, c_local_object_part(m)%iend)
           dqright = x - ileft
           dqleft = 1.0_8 - dqright
           c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft)   + dqleft
           c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) + dqright
        END IF

     END IF  ! IF (dirflag.EQ.1) THEN

  ELSE
! inner object

! ilt   -  ---- irt
!  |             |
!  |             |
!  1 ilbb ---- irb
!
     i_left_top        =                  myobject%jtop   - myobject%jbottom + 1
     i_right_top       = i_left_top     + myobject%iright - myobject%ileft
     i_right_bottom    = i_right_top    + myobject%jtop   - myobject%jbottom
     i_left_bottom_bis = i_right_bottom + myobject%iright - myobject%ileft - 1

     IF (dirflag.EQ.1) THEN
! collision with the wall on the left
        x_new = myobject%xmax + 1.0d-6
        y_new = y
        vx_new = v_new * COS(theta)
        vy_new = v_new * SIN(theta) * SIN(fi)
        vz_new = v_new * SIN(theta) * COS(fi)
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           i = MIN(INT(myobject%ymax - y) + i_right_top, i_right_bottom - 1)
           dqi = y - INT(y)
           dqip1 = 1.0_8 - dqi
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi     ! note the plus. remember ::
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1   ! we subtract electron charge, it's negative
        END IF

     ELSE IF (dirflag.EQ.2) THEN
!  collision with the wall above
        x_new = x
        y_new = myobject%ymin - 1.0d-6 
        vx_new =  v_new * SIN(theta) * SIN(fi) 
        vy_new = -v_new * COS(theta)
        vz_new =  v_new * SIN(theta) * COS(fi)
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           i = MIN(INT(myobject%xmax - x) + i_right_bottom, i_left_bottom_bis)
           dqi = x - INT(x)
           dqip1 = 1.0_8 - dqi
           ip1 = i+1
           IF (i.EQ.i_left_bottom_bis) ip1 = 1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
           myobject%surface_charge_variation(ip1) = myobject%surface_charge_variation(ip1) + dqip1
        END IF

     ELSE IF (dirflag.EQ.3) THEN
! collision with the wall on the right
        x_new = myobject%xmin - 1.0d-6
        y_new = y
        vx_new = -v_new * COS(theta)
        vy_new =  v_new * SIN(theta) * SIN(fi)
        vz_new =  v_new * SIN(theta) * COS(fi)
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           i = MIN(INT(y - myobject%ymin) + 1, i_left_top - 1)
           dqip1 = y - INT(y)
           dqi = 1.0_8 - dqip1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1
        END IF

     ELSE IF (dirflag.EQ.4) THEN
! collision with the wall below
        x_new = x
        y_new = myobject%ymax + 1.0d-6 
        vx_new = v_new * SIN(theta) * SIN(fi) 
        vy_new = v_new * COS(theta)
        vz_new = v_new * SIN(theta) * COS(fi)
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           i = MIN(INT(x - myobject%xmin) + i_left_top, i_right_top - 1)
           dqip1 = x - INT(x)
           dqi = 1.0_8 - dqip1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1
        END IF

     END IF  ! IF (dirflag.EQ.1) THEN

  END IF   ! IF (m.GT.0) THEN

  CALL ADD_ELECTRON_TO_ADD_LIST(x_new, y_new, vx_new, vy_new, vz_new, tag_new)    ! the particle gets tag equal to nwo - the number of the whole object that emitted it

  myobject%electron_emit_count = myobject%electron_emit_count + 1

END SUBROUTINE INJECT_INELASTIC_BACKSCATTERED_ELECTRON

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE INJECT_TRUE_SECONDARY_ELECTRON(x, y, energy_inc, tag, myobject, m, dirflag)

!  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues

  USE rng_wrapper

  use mpi

  IMPLICIT NONE

  REAL(8) x, y         ! primary electron coordinates
  REAL(8) energy_inc   ! energy leftover of the primary electron velocity
  INTEGER tag          ! primary electron tag 
!  INTEGER nwo          ! number of the whole object with which the primary electron collided
  TYPE(boundary_object) myobject
  INTEGER m            ! number of the local [belonging to this particular cluster] part of the whole object with which the primary electron collided
  
  INTEGER dirflag      ! direction flag, here left/right/below/above are -x/+x/-y/+y directions
                       ! 1 = collision with wall on the left
                       ! 2 = collision with wall above
                       ! 3 = collision with wall on the right
                       ! 4 = collision with wall below

  REAL(8) x_new, y_new             ! emitted electron coordinates
  REAL(8) vx_new, vy_new, vz_new   ! emitted electron velocity components
  INTEGER tag_new                  ! emitted electron tag 

  INTEGER jbelow, jabove
  REAL(8) dqbelow, dqabove

  INTEGER ileft, iright
  REAL(8) dqleft, dqright

  REAL(8) energy_new               ! scattered electron energy
  REAL(8) alpha                    ! factor used to correct velocity of emitted electron

! for inner object
  INTEGER i_left_top, i_right_top, i_right_bottom, i_left_bottom_bis, ip1
  INTEGER i
  REAL(8) dqi, dqip1

! check the energy leftover of the primary electron to make sure that the energy conservation is not violated
  IF (energy_inc.LE.0.0_8) RETURN

  IF (m.GT.0) THEN
! domain boundary

     IF (dirflag.EQ.1) THEN
! collision with the wall on the left
        x_new = DBLE(c_indx_x_min) + 1.0d-6
        y_new = y
        CALL GetInjMaxwellVelocity(vx_new)
        CALL GetMaxwellVelocity(vy_new)
        CALL GetMaxwellVelocity(vz_new)
        vx_new = vx_new * myobject%factor_convert_seetrue_vinj
        vy_new = vy_new * myobject%factor_convert_seetrue_vinj
        vz_new = vz_new * myobject%factor_convert_seetrue_vinj
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
           jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
           dqabove = y - jbelow
           dqbelow = 1.0_8 - dqabove
           c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) + dqbelow   ! note the plus. remember ::
           c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) + dqabove   ! we subtract electron charge, it's negative
        END IF

     ELSE IF (dirflag.EQ.2) THEN
!  collision with the wall above
        x_new = x
        y_new = DBLE(c_indx_y_max) - 1.0d-6 
        CALL GetMaxwellVelocity(vx_new)
        CALL GetInjMaxwellVelocity(vy_new)
        CALL GetMaxwellVelocity(vz_new)
        vx_new =  vx_new * myobject%factor_convert_seetrue_vinj
        vy_new = -vy_new * myobject%factor_convert_seetrue_vinj
        vz_new =  vz_new * myobject%factor_convert_seetrue_vinj
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           ileft = MAX(INT(x), c_local_object_part(m)%istart)
           iright = MIN(ileft + 1, c_local_object_part(m)%iend)
           dqright = x - ileft
           dqleft = 1.0_8 - dqright
           c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft)   + dqleft
           c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) + dqright
        END IF

     ELSE IF (dirflag.EQ.3) THEN
! collision with the wall on the right
        x_new = DBLE(c_indx_x_max) - 1.0d-6
        y_new = y
        CALL GetInjMaxwellVelocity(vx_new)
        CALL GetMaxwellVelocity(vy_new)
        CALL GetMaxwellVelocity(vz_new)
        vx_new = -vx_new * myobject%factor_convert_seetrue_vinj
        vy_new =  vy_new * myobject%factor_convert_seetrue_vinj
        vz_new =  vz_new * myobject%factor_convert_seetrue_vinj
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           jbelow = MAX(INT(y), c_local_object_part(m)%jstart)
           jabove = MIN(jbelow + 1, c_local_object_part(m)%jend)
           dqabove = y - jbelow
           dqbelow = 1.0_8 - dqabove
           c_local_object_part(m)%surface_charge(jbelow) = c_local_object_part(m)%surface_charge(jbelow) + dqbelow   ! note the plus. remember ::
           c_local_object_part(m)%surface_charge(jabove) = c_local_object_part(m)%surface_charge(jabove) + dqabove   ! we subtract electron charge, it's negative
        END IF

     ELSE IF (dirflag.EQ.4) THEN
! collision with the wall below
        x_new = x
        y_new = DBLE(c_indx_y_min) + 1.0d-6 
        CALL GetMaxwellVelocity(vx_new)
        CALL GetInjMaxwellVelocity(vy_new)
        CALL GetMaxwellVelocity(vz_new)
        vx_new = vx_new * myobject%factor_convert_seetrue_vinj
        vy_new = vy_new * myobject%factor_convert_seetrue_vinj
        vz_new = vz_new * myobject%factor_convert_seetrue_vinj
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           ileft = MAX(INT(x), c_local_object_part(m)%istart)
           iright = MIN(ileft + 1, c_local_object_part(m)%iend)
           dqright = x - ileft
           dqleft = 1.0_8 - dqright
           c_local_object_part(m)%surface_charge(ileft) = c_local_object_part(m)%surface_charge(ileft)   + dqleft
           c_local_object_part(m)%surface_charge(iright) = c_local_object_part(m)%surface_charge(iright) + dqright
        END IF

     END IF  ! IF (dirflag.EQ.1) THEN

  ELSE
! inner object

! ilt   -  ---- irt
!  |             |
!  |             |
!  1 ilbb ---- irb
!
     i_left_top        =                  myobject%jtop   - myobject%jbottom + 1
     i_right_top       = i_left_top     + myobject%iright - myobject%ileft
     i_right_bottom    = i_right_top    + myobject%jtop   - myobject%jbottom
     i_left_bottom_bis = i_right_bottom + myobject%iright - myobject%ileft - 1

     IF (dirflag.EQ.1) THEN
! collision with the wall on the left
        x_new = myobject%xmax + 1.0d-6
        y_new = y
        CALL GetInjMaxwellVelocity(vx_new)
        CALL GetMaxwellVelocity(vy_new)
        CALL GetMaxwellVelocity(vz_new)
        vx_new = vx_new * myobject%factor_convert_seetrue_vinj
        vy_new = vy_new * myobject%factor_convert_seetrue_vinj
        vz_new = vz_new * myobject%factor_convert_seetrue_vinj
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           i = MIN(INT(myobject%ymax - y) + i_right_top, i_right_bottom - 1)
           dqi = y - INT(y)
           dqip1 = 1.0_8 - dqi
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi     ! note the plus. remember ::
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1   ! we subtract electron charge, it's negative
        END IF

     ELSE IF (dirflag.EQ.2) THEN
!  collision with the wall above
        x_new = x
        y_new = myobject%ymin - 1.0d-6 
        CALL GetMaxwellVelocity(vx_new)
        CALL GetInjMaxwellVelocity(vy_new)
        CALL GetMaxwellVelocity(vz_new)
        vx_new =  vx_new * myobject%factor_convert_seetrue_vinj
        vy_new = -vy_new * myobject%factor_convert_seetrue_vinj
        vz_new =  vz_new * myobject%factor_convert_seetrue_vinj
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           i = MIN(INT(myobject%xmax - x) + i_right_bottom, i_left_bottom_bis)
           dqi = x - INT(x)
           dqip1 = 1.0_8 - dqi
           ip1 = i+1
           IF (i.EQ.i_left_bottom_bis) ip1 = 1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
           myobject%surface_charge_variation(ip1) = myobject%surface_charge_variation(ip1) + dqip1
        END IF

     ELSE IF (dirflag.EQ.3) THEN
! collision with the wall on the right
        x_new = myobject%xmin - 1.0d-6
        y_new = y
        CALL GetInjMaxwellVelocity(vx_new)
        CALL GetMaxwellVelocity(vy_new)
        CALL GetMaxwellVelocity(vz_new)
        vx_new = -vx_new * myobject%factor_convert_seetrue_vinj
        vy_new =  vy_new * myobject%factor_convert_seetrue_vinj
        vz_new =  vz_new * myobject%factor_convert_seetrue_vinj
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           i = MIN(INT(y - myobject%ymin) + 1, i_left_top - 1)
           dqip1 = y - INT(y)
           dqi = 1.0_8 - dqip1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1
        END IF

     ELSE IF (dirflag.EQ.4) THEN
! collision with the wall below
        x_new = x
        y_new = myobject%ymax + 1.0d-6 
        CALL GetMaxwellVelocity(vx_new)
        CALL GetInjMaxwellVelocity(vy_new)
        CALL GetMaxwellVelocity(vz_new)
        vx_new = vx_new * myobject%factor_convert_seetrue_vinj
        vy_new = vy_new * myobject%factor_convert_seetrue_vinj
        vz_new = vz_new * myobject%factor_convert_seetrue_vinj
        tag_new = myobject%object_id_number

        IF (myobject%object_type.EQ.DIELECTRIC) THEN
! update the surface charge
           i = MIN(INT(x - myobject%xmin) + i_left_top, i_right_top - 1)
           dqip1 = x - INT(x)
           dqi = 1.0_8 - dqip1
           myobject%surface_charge_variation(i)   = myobject%surface_charge_variation(i)   + dqi
           myobject%surface_charge_variation(i+1) = myobject%surface_charge_variation(i+1) + dqip1
        END IF

     END IF  ! IF (dirflag.EQ.1) THEN

  END IF  ! IF (m.GT.0) THEN

! make sure that the energy of the emitted electron does not exceed the present leftover of the primary electron energy
  energy_new = vx_new * vx_new + vy_new * vy_new + vz_new * vz_new
  IF (energy_new.GT.energy_inc) THEN
     alpha = SQRT(energy_inc/energy_new)
     vx_new = vx_new * alpha
     vy_new = vy_new * alpha
     vz_new = vz_new * alpha
     energy_inc = 0.0_8
  ELSE
     energy_inc = energy_inc - energy_new
  END IF

  CALL ADD_ELECTRON_TO_ADD_LIST(x_new, y_new, vx_new, vy_new, vz_new, tag_new)    ! the particle gets tag equal to nwo - the number of the whole object that emitted it

  myobject%electron_emit_count = myobject%electron_emit_count + 1

END SUBROUTINE INJECT_TRUE_SECONDARY_ELECTRON

!-------------------------------------------------------------------------------------------------
!
REAL(8) FUNCTION Coeff_SEE_Elastic(energy, theta, myobject)

  USE CurrentProblemValues !, ONLY : whole_object
  use mpi

  IMPLICIT NONE

  REAL(8) energy       ! dim-less energy of incident electron
  REAL(8) theta        ! angle of incidence 
  TYPE(boundary_object) myobject
!  INTEGER nwo          ! number of the whole object with which the primary electron collided

  REAL(8) coeff_max_theta
  REAL(8) Coeff_SEE_Classic

  SELECT CASE (myobject%Emitted_model(1))

     CASE (0)
        Coeff_SEE_Elastic = 0.0_8

     CASE (1)
        IF ((energy.GT.myobject%minE_see_elastic).AND.(energy.LE.myobject%maxE_see_elastic)) THEN
           Coeff_SEE_Elastic = myobject%setD_see_elastic
        END IF

     CASE (2) 
        IF (energy.GE.myobject%E_elast_max) THEN
           coeff_max_theta    = myobject%maxD_elast !!!* (1.0_8 + 0.159154943_8 * k_smooth * theta**2) ! 0.159154943 = 1 / 2*pi
           Coeff_SEE_Elastic = coeff_max_theta * &
                             & EXP( - (energy - myobject%E_elast_max) / myobject%dE_elast) * &
                             & (energy - myobject%E_elast_max + myobject%dE_elast) / myobject%dE_elast + &
                             & myobject%Frac_elast_highenergy * Coeff_SEE_Classic(energy, theta, myobject)

        ELSE IF (energy.GT.myobject%E_elast_0) THEN
           coeff_max_theta    = myobject%maxD_elast !!!* (1.0_8 + 0.159154943_8 * k_smooth * theta**2) ! 0.159154943 = 1 / 2*pi         
           Coeff_SEE_Elastic = coeff_max_theta * & 
                             & EXP(-(energy - myobject%E_elast_max) / (myobject%E_elast_max - myobject%E_elast_0) ) * &
                             & (energy - myobject%E_elast_0) / (myobject%E_elast_max - myobject%E_elast_0) + & 
                             & myobject%Frac_elast_highenergy * Coeff_SEE_Classic(energy, theta, myobject)

        ELSE
           Coeff_SEE_Elastic = 0.0_8

        END IF
        
  END SELECT

END FUNCTION Coeff_SEE_Elastic

!-------------------------------------------------------------------------------------------------
!
REAL(8) FUNCTION Coeff_SEE_Inelastic(energy, theta, myobject)

  USE CurrentProblemValues !, ONLY : whole_object
  use mpi

  IMPLICIT NONE

  REAL(8) energy       ! dim-less energy of incident electron
  REAL(8) theta        ! angle of incidence 
!  INTEGER nwo          ! number of the whole object with which the primary electron collided
  TYPE(boundary_object) myobject

  REAL(8) Coeff_SEE_Classic

!print '("I am here with ",e12.5,2x,e12.5,2x,i3)', energy, theta, nwo

  SELECT CASE (myobject%Emitted_model(2))

     CASE (0)
        Coeff_SEE_Inelastic = 0.0_8

     CASE (1)
        IF ((energy.GT.myobject%minE_see_inelastic).AND.(energy.LE.myobject%maxE_see_inelastic)) THEN
!print '("I should not be here")'
           Coeff_SEE_Inelastic = myobject%setD_see_inelastic
        END IF

     CASE (2) 
        Coeff_SEE_Inelastic = myobject%Frac_inelastic * Coeff_SEE_Classic(energy, theta, myobject)
!print '(2(2x,e10.3))', energy, Coeff_SEE_Inelastic
        
  END SELECT

END FUNCTION Coeff_SEE_Inelastic


!-------------------------------------------------------------------------------------------------
!
REAL(8) FUNCTION Coeff_SEE_True(energy, theta, myobject)

  USE CurrentProblemValues !, ONLY : whole_object
  use mpi

  IMPLICIT NONE

  REAL(8) energy       ! dim-less energy of incident electron
  REAL(8) theta        ! angle of incidence, with respect to the surface normal  (0 = normal incidence)
!  INTEGER nwo          ! number of the whole object with which the primary electron collided
  TYPE(boundary_object) myobject

  REAL(8) Coeff_SEE_Classic

  SELECT CASE (myobject%Emitted_model(3))

     CASE (0)
        Coeff_SEE_True = 0.0_8

     CASE (1)
        IF ((energy.GT.myobject%minE_see_true).AND.(energy.LE.myobject%maxE_see_true)) THEN
           Coeff_SEE_True = myobject%setD_see_true
        END IF

     CASE (2) 
        Coeff_SEE_True = Coeff_SEE_Classic(energy, theta, myobject) * (1.0_8 - myobject%Frac_elast_highenergy - myobject%Frac_inelastic)
        IF (Coeff_SEE_True.LT.0.0_8) Coeff_SEE_True = 0.0_8        
        
  END SELECT

END FUNCTION Coeff_SEE_True

!-------------------------------------------------------------------------------------------------
!
REAL(8) FUNCTION Coeff_SEE_Classic(energy, theta, myobject)

  USE CurrentProblemValues !, ONLY : whole_object
  use mpi

  IMPLICIT NONE

  REAL(8) energy       ! dim-less energy of incident electron
  REAL(8) theta        ! angle of incidence, with respect to the surface normal  (0 = normal incidence)
!  INTEGER nwo          ! number of the whole object with which the primary electron collided
  TYPE(boundary_object) myobject

  REAL(8) energy_max_theta
  REAL(8) v
  REAL(8) coeff_max_theta
  REAL k

  IF (energy.GT.myobject%E_see_0) THEN

     energy_max_theta = myobject%E_see_max * (1.0_8 + 0.318309886_8 * myobject%k_smooth * theta**2)     ! 0.318309886 = 1 / pi

     v = (energy - myobject%E_see_0) / (energy_max_theta - myobject%E_see_0) 

     coeff_max_theta  = myobject%maxD_see_classic * (1.0_8 + 0.159154943_8 * myobject%k_smooth * theta**2) ! 0.159154943 = 1 / 2*pi
     
     IF (v.LT.1.0_8) THEN 
        k = 0.62
     ELSE
        k = 0.25
     END IF
     
     Coeff_SEE_Classic = coeff_max_theta * (v * EXP(1.0_8 - v)) ** k
     
  ELSE
     
     Coeff_SEE_Classic = 0.0_8

  END IF

END FUNCTION Coeff_SEE_Classic

