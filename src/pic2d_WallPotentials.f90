SUBROUTINE UPDATE_WALL_POTENTIALS(T_cntr)

  USE CurrentProblemValues, ONLY : whole_object, N_of_boundary_objects, METAL_WALL, VACUUM_GAP

  IMPLICIT NONE

  INTEGER T_cntr
  INTEGER n, i, j

! set potentials of metal electrodes
  DO n = 1, N_of_boundary_objects
     IF (whole_object(n)%object_type.NE.METAL_WALL) CYCLE   
     whole_object(n)%phi = whole_object(n)%phi_const + whole_object(n)%phi_var * SIN(whole_object(n)%omega * T_cntr + whole_object(n)%phase)
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
