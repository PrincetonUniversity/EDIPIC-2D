!-----------------------------------------
!
SUBROUTINE PREPARE_WALL_MATERIALS

  USE ParallelOperationValues
  USE CurrentProblemValues

  IMPLICIT NONE

  INTEGER n

  CHARACTER(24) initmaterial_filename     ! init_material_AAAAAA.dat
                                          ! ----x----I----x----I----
  LOGICAl exists
  
  CHARACTER(1) buf 

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  DO n = 1, N_of_boundary_objects
     whole_object(n)%eps_diel = 1.0_8
     whole_object(n)%SEE_enabled = .FALSE.
     whole_object(n)%ion_induced_EE_enabled = .FALSE.

     IF (whole_object(n)%object_type.EQ.VACUUM_GAP) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### PREPARE_WALL_MATERIALS :: no material properties data file for VACUUM_GAP boundary object ",i2)', n
        CYCLE
     END IF

     IF (whole_object(n)%object_type.EQ.PERIODIC_PIPELINE_X) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### PREPARE_WALL_MATERIALS :: no material properties data file for PERIODIC_PIPELINE_X boundary object ",i2)', n
        CYCLE
     END IF

     IF (whole_object(n)%object_type.EQ.PERIODIC_PIPELINE_Y) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### PREPARE_WALL_MATERIALS :: no material properties data file for PERIODIC_PIPELINE_Y boundary object ",i2)', n
        CYCLE
     END IF

     initmaterial_filename = 'init_material_AAAAAA.dat'
     initmaterial_filename(15:20) = whole_object(n)%material

     INQUIRE (FILE = initmaterial_filename, EXIST = exists)
     IF (.NOT.exists) THEN
        IF (Rank_of_process.EQ.0) PRINT '("### PREPARE_WALL_MATERIALS :: file ",A24," not found, use default material properties for boundary object ",i2)', initmaterial_filename, n
        CYCLE
     END IF

  END DO

END SUBROUTINE PREPARE_WALL_MATERIALS
