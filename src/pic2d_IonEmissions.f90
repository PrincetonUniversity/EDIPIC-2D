
!-------------------------------------------------------------------------------------------------
!
SUBROUTINE PROCESS_ION_INDUCED_ELECTRON_EMISSION(s, x, y, vx, vy, vz, tag, myobject, m, dirflag)

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles, ONLY : Ms

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER s            ! ion species index
  REAL(8) x, y         ! primary ion coordinates
  REAL(8) vx, vy, vz   ! primary ion velocities
  INTEGER tag          ! primary ion tag 
!  INTEGER nwo          ! number of the whole object with which the primary electron collided
  TYPE(boundary_object) myobject
  INTEGER m            ! number of the local [belonging to this particular cluster] part of the whole object with which the primary electron collided
  
  INTEGER dirflag      ! direction flag, here left/right/below/above are -x/+x/-y/+y directions
                       ! 1 = collision with wall on the left
                       ! 2 = collision with wall above
                       ! 3 = collision with wall on the right
                       ! 4 = collision with wall below

  REAL(8) energy_inc      ! energy of primary ion
  REAL(8) v               ! speed of primary ion
  REAL(8) theta_inc       ! angle of incidence with respect to wall normal 

! presently there is only one type of electron emission caused by ion impact, similar to the true secondary electron emission

  REAL(8) coef_total      ! corresponding emission coefficient at current energy and angle of incidence

  REAL(8) R               ! random number 

! functions
  REAL(8) Coeff_ii_EE_True

! calculate the dim-less energy of incident ion
  v = vx**2 + vy**2 + vz**2
  energy_inc = Ms(s) * v
  IF (energy_inc.LE.myobject%minE_ii_ee_true(s)) RETURN
  IF (energy_inc.GE.myobject%maxE_ii_ee_true(s)) RETURN

! calculate the angle of incidence
  v = SQRT(v)

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

  coef_total = Coeff_ii_EE_True(s, energy_inc, theta_inc, myobject)

! quit subroutine if ion-induced EE is turned off
  IF (coef_total.EQ.0.0_8) RETURN

! we want to have some safety mechanism to prevent emission of more particles than the energy of the primary electron can afford
! reduce the primary electron energy by the threshold for the true secondary electron emission
  energy_inc = MAX(0.0_8, energy_inc - myobject%minE_ii_ee_true(s))

! then we inject true secondary electrons until we exhaust the integer part of the ratio
  DO WHILE (coef_total.GE.1.0_8)
     CALL INJECT_ION_IMPACT_EXTRACTED_ELECTRON(s, x, y, energy_inc, tag, myobject, m, dirflag)   ! this call reduces energy_inc
     coef_total = coef_total - 1.0_8
  END DO
! the remaining part gives us the probability for statistical injection of fractional part of true secondary electron

  IF (coef_total.EQ.0.0_8) RETURN

! take a random number
  R = well_random_number()
  IF (R.LT.coef_total) CALL INJECT_ION_IMPACT_EXTRACTED_ELECTRON(s, x, y, energy_inc, tag, myobject, m, dirflag)   ! this call reduces energy_inc

END SUBROUTINE PROCESS_ION_INDUCED_ELECTRON_EMISSION

!-------------------------------------------------------------------------------------------------
!
SUBROUTINE INJECT_ION_IMPACT_EXTRACTED_ELECTRON(s, x, y, energy_inc, tag, myobject, m, dirflag)

!  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER s            ! ion species index
  REAL(8) x, y         ! primary ion coordinates
  REAL(8) energy_inc   ! energy leftover of the primary ion
  INTEGER tag          ! primary ion tag 
  TYPE(boundary_object) myobject
  INTEGER m            ! number of the local [belonging to this particular cluster] part of the whole object with which the primary electron collided
                       ! must be -1 for innner objects
  
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
        vx_new = vx_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vy_new = vy_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vz_new = vz_new * myobject%factor_convert_ii_ee_true_vinj(s)
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
        vx_new =  vx_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vy_new = -vy_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vz_new =  vz_new * myobject%factor_convert_ii_ee_true_vinj(s)
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
        vx_new = -vx_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vy_new =  vy_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vz_new =  vz_new * myobject%factor_convert_ii_ee_true_vinj(s)
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
        vx_new = vx_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vy_new = vy_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vz_new = vz_new * myobject%factor_convert_ii_ee_true_vinj(s)
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
        vx_new = vx_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vy_new = vy_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vz_new = vz_new * myobject%factor_convert_ii_ee_true_vinj(s)
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
        vx_new =  vx_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vy_new = -vy_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vz_new =  vz_new * myobject%factor_convert_ii_ee_true_vinj(s)
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
        vx_new = -vx_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vy_new =  vy_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vz_new =  vz_new * myobject%factor_convert_ii_ee_true_vinj(s)
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
        vx_new = vx_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vy_new = vy_new * myobject%factor_convert_ii_ee_true_vinj(s)
        vz_new = vz_new * myobject%factor_convert_ii_ee_true_vinj(s)
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

END SUBROUTINE INJECT_ION_IMPACT_EXTRACTED_ELECTRON

!-------------------------------------------------------------------------------------------------
!
REAL(8) FUNCTION Coeff_ii_EE_True(s, energy, theta, myobject)

  USE CurrentProblemValues !, ONLY : whole_object
  IMPLICIT NONE

  INTEGER s            ! ion species index
  REAL(8) energy       ! dim-less energy of incident ion
  REAL(8) theta        ! angle of incidence, with respect to the surface normal  (0 = normal incidence)
  TYPE(boundary_object) myobject

  Coeff_ii_EE_True = 0.0_8

  IF ((energy.GT.myobject%minE_ii_ee_true(s)).AND.(energy.LE.myobject%maxE_ii_ee_true(s))) THEN
     Coeff_ii_EE_True = myobject%setD_ii_ee_true(s)
  END IF

END FUNCTION Coeff_ii_EE_True

!-------------------------------------------------------------------------------------------------
! when 100% ion reflection occurs, charge density on dielectric objects is not changed
! correspondingly, this procedure does not update the charge density 
!
SUBROUTINE INJECT_REFLECTED_ION(s, x, y, vx, vy, vz, tag, myobject, m, dirflag)

!  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE CurrentProblemValues

  USE rng_wrapper

  IMPLICIT NONE

  INTEGER s            ! index of ion species
  REAL(8) x, y         ! primary ion coordinates
  REAL(8) vx, vy, vz   ! primary ion velocity components
  INTEGER tag          ! primary electron tag 
  
  TYPE(boundary_object) myobject
  INTEGER m            ! number of the local [belonging to this particular cluster] part of the whole object with which the primary ion collided
  
  INTEGER dirflag      ! direction flag, here left/right/below/above are -x/+x/-y/+y directions
                       ! 1 = collision with wall on the left
                       ! 2 = collision with wall above
                       ! 3 = collision with wall on the right
                       ! 4 = collision with wall below

  REAL(8) x_new, y_new             ! reflected ion coordinates
  REAL(8) vx_new, vy_new, vz_new   ! reflected ion velocity components
  INTEGER tag_new                  ! reflected ion tag 

  IF (m.GT.0) THEN
! domain boundary 

! only specular reflection for now

     IF (dirflag.EQ.1) THEN
! collision with the wall on the left
        x_new = DBLE(c_indx_x_min) + 1.0d-6
        y_new = y
        vx_new = -vx
        vy_new =  vy
        vz_new =  vz
        tag_new = myobject%object_id_number

     ELSE IF (dirflag.EQ.2) THEN
!  collision with the wall above
        x_new = x
        y_new = DBLE(c_indx_y_max) - 1.0d-6 
        vx_new =  vx
        vy_new = -vy
        vz_new =  vz
        tag_new = myobject%object_id_number

     ELSE IF (dirflag.EQ.3) THEN
! collision with the wall on the right
        x_new = DBLE(c_indx_x_max) - 1.0d-6
        y_new = y
        vx_new = -vx
        vy_new =  vy
        vz_new =  vz
        tag_new = myobject%object_id_number

     ELSE IF (dirflag.EQ.4) THEN
! collision with the wall below
        x_new = x
        y_new = DBLE(c_indx_y_min) + 1.0d-6 
        vx_new =  vx
        vy_new = -vy
        vz_new =  vz
        tag_new = myobject%object_id_number

     END IF   ! IF (dirflag.EQ.1) THEN

  ELSE
! inner object

     IF (dirflag.EQ.1) THEN
! collision with the wall on the left
        x_new = myobject%xmax + 1.0d-6
        y_new = y
        vx_new = -vx
        vy_new =  vy
        vz_new =  vz
        tag_new = myobject%object_id_number

     ELSE IF (dirflag.EQ.2) THEN
!  collision with the wall above
        x_new = x
        y_new = myobject%ymin - 1.0d-6 
        vx_new =  vx
        vy_new = -vy
        vz_new =  vz
        tag_new = myobject%object_id_number

     ELSE IF (dirflag.EQ.3) THEN
! collision with the wall on the right
        x_new = myobject%xmin - 1.0d-6
        y_new = y
        vx_new = -vx
        vy_new =  vy
        vz_new =  vz
        tag_new = myobject%object_id_number

     ELSE IF (dirflag.EQ.4) THEN
! collision with the wall below
        x_new = x
        y_new = myobject%ymax + 1.0d-6 
        vx_new =  vx
        vy_new = -vy
        vz_new =  vz
        tag_new = myobject%object_id_number

     END IF   ! IF (dirflag.EQ.1) THEN

  END IF   ! IF (m.GT.0) THEN

  CALL ADD_ION_TO_ADD_LIST(s, x_new, y_new, vx_new, vy_new, vz_new, tag_new)    ! the particle gets tag equal to the number of the whole object that reflected it

!???  myobject%ion_emit_count(s) = myobject%ion_emit_count(s) + 1

END SUBROUTINE INJECT_REFLECTED_ION
