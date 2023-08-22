!-----------------------------------
!
! This module initially was written by Janhunen Jani Salomon
! A number of changes listed below was made by DS (me)
!
module PETSc_Solver
  use ParallelOperationValues
  use CurrentProblemValues
  use BlockAndItsBoundaries
  use petsc
  implicit none
  private

#include <petsc/finclude/petsc.h>

  KSP ksp     ! solver object
  PC  pc      ! pre-conditioner object
  
  Mat Amat                ! Matrix type
  Vec bvec                ! Right hand side
  Vec xvec                ! Solution
  PetscInt ntoteq         ! total number of equations
  PetscInt nloceq         ! local number of equations
  
  integer :: parComm
  
  public SolverInitialization, Solve, InsertData, FetchData, SolverDestroy
  
contains
  
  subroutine Solve
    PetscErrorCode :: ierr

! sets KSP options from the options database
    call KSPSetFromOptions(ksp,ierr)
! solves linear system
    call KSPSolve(ksp, bvec, xvec, ierr)

    return
  end subroutine Solve
  
  subroutine SolverInitialization(comm)

    implicit none
    
    integer, intent(IN) ::  comm  !??? not used anywhere???

    PetscInt :: one, four, five

    character(20) :: petsc_config='petsc.rc'   ! PETSc database file
    PetscErrorCode :: ierr

    integer jbegin, jend, ibegin, iend
    PetscInt :: irow_global
    PetscInt :: jcolumn_global(1:5)
    PetscScalar :: value_at_jcol(1:5)
    integer i, j

    INTEGER nio, position_flag

    REAL(8), ALLOCATABLE :: eps_ishifted_j(:,:)
    REAL(8), ALLOCATABLE :: eps_i_jshifted(:,:)
    INTEGER ALLOC_ERR

!    integer            :: m, n, nx, ny, i, j, k, ix, jy
!    PetscInt :: nrows, ncols, one=1, five=5, temp
!    PetscInt :: irow(5), jcol(5), tmpcol(5)
!    PetscScalar,parameter  :: wv(5)=(/ 0.25, 0.25, -1.0, 0.25, 0.25 /)
!    PetscScalar :: v(5)
!    PetscViewer :: viewer
!    PetscScalar :: tol
!    INTEGER local_x_max,local_x_min,local_y_min,local_y_max
!    INTEGER i_local,j_local

    one=1
    four=4
    five=5
    
! Initializes the petsc database and mpi
    call PetscInitialize(petsc_config, ierr)

! Remember the communicator for future reference
    parComm=PETSC_COMM_WORLD

! Creates a matrix where the type is determined from either a call to MatSetType() 
! or from the options database with a call to MatSetFromOptions()
! The default matrix type is AIJ
    call MatCreate(parComm, Amat, ierr)

    ntoteq=N_to_solve_total
    nloceq=block_N_of_nodes_to_solve   !######## ??????????????????????????? NEW ##########

! Set MPI distribution for A matrix [by S.J.]

! Sets the local and global sizes and checks to determine compatibility
! below we use PETSC_DECIDE instead of number of local rows/columns
!####??????    call MatSetSizes(Amat, PETSC_DECIDE, PETSC_DECIDE, ntoteq, ntoteq, ierr)
    call MatSetSizes(Amat, nloceq, nloceq, ntoteq, ntoteq, ierr)    !######## ??????????????????????????? NEW ##########

! Builds matrix object for a particular matrix type
! MATAIJ = "aij" - a matrix type to be used for sparce matrices
    call MatSetType(Amat, MATAIJ, ierr)

! Creates a matrix where the type is determined from the options database
    call MatSetFromOptions(Amat, ierr)
    five=5    !?????  why???

! For good matrix assembly performance the user should preallocate the matrix storage
! the second argument is the number of nonzeros per row (same for all rows)
    call MatSeqAIJSetPreallocation(Amat, five, PETSC_NULL_INTEGER, ierr)
    five=5    !????? why???

! Preallocates memory for a sparce parallel matrix in AIJ format (the default parallel petsc format)
! the second argument is the number of nonzeros per row in DIAGONAL portion of local submatrix
! the fourth argument is the number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix
    call MatMPIAIJSetPreallocation(Amat, five, PETSC_NULL_INTEGER, five, PETSC_NULL_INTEGER, ierr)
    one=1    !????  why???

! ny - the number of columns
! jcol - global indices of columns
! v - a logically two-dimensional array of values

!----------------------------------------
! Initialization of matrix coefficients written by DS
!

    jbegin = indx_y_min+1
    jend   = indx_y_max-1
    ibegin = indx_x_min+1
    iend   = indx_x_max-1

    IF (Rank_of_process_left.LT.0)  ibegin = indx_x_min
    IF (Rank_of_process_right.LT.0) iend   = indx_x_max
    IF (Rank_of_process_below.LT.0) jbegin = indx_y_min
    IF (Rank_of_process_above.LT.0) jend   = indx_y_max

    ALLOCATE(eps_ishifted_j(ibegin:iend+1, jbegin:jend), STAT=ALLOC_ERR)   ! eps_ishifted_j(i,j) is between nodes {i-1,j} and {i,j}
    ALLOCATE(eps_i_jshifted(ibegin:iend, jbegin:jend+1), STAT=ALLOC_ERR)   ! eps_i_jshifted(i,j) is between nodes {i,j-1} and {i,j}

    eps_ishifted_j = 1.0_8
    eps_i_jshifted = 1.0_8

    DO j = jbegin, jend
       DO i = ibegin, iend+1
          CALL SET_EPS_ISHIFTED(i, j, eps_ishifted_j(i,j))
       END DO
    END DO

    DO j = jbegin, jend+1
       DO i = ibegin, iend
          CALL SET_EPS_JSHIFTED(i, j, eps_i_jshifted(i,j))
       END DO
    END DO

    irow_global = global_offset

!    j = indx_y_min !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
    IF (jbegin.EQ.indx_y_min) THEN
! boundary object along bottom border
       DO i = ibegin, iend
          irow_global = irow_global + 1
!          number_of_columns = 1
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       END DO
    END IF

!    j = indx_y_min+1 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    j = indx_y_min+1

    i = indx_x_min

    IF (ibegin.EQ.indx_x_min) THEN
! boundary object along left border
       irow_global = irow_global + 1
       IF (.NOT.block_has_symmetry_plane_X_left) THEN
! Dirichlet (given potential) boundary
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       ELSE
! the left border is a symmetry plane
          IF (jbegin.EQ.indx_y_min) THEN                       ! BELOW
! boundary object along the bottom border
             jcolumn_global(1) = irow_global - (iend-ibegin+1)                          ! use the own node
          ELSE
! use a node from neighbor below
             jcolumn_global(1) = process_below_left_top_inner_node-1                    ! use a node from the neighbor below
          END IF
          jcolumn_global(2) = irow_global                      ! CENTER
          jcolumn_global(3) = irow_global+1                    ! RIGHT
          jcolumn_global(4) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
          CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min, indx_y_min+1, nio, position_flag)

! note that we are at the left edge of the domain and the symmetry is applied here,
! so the boundary object must be symmetric relative to x=0 as well
! therefore we are either at the bottom surface, or inside, or at the top surface of the inner object

          SELECT CASE (position_flag)
             CASE (9)
! metal
                jcolumn_global(1) = irow_global
                value_at_jcol(1) = 1.0_8
                call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!             CASE (1,2)
!! dielectric surface above
!                value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
!             CASE (6,7)
!! dielectric surface below
!                value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr)              
             CASE DEFAULT 
! inside dielectric or plasma
                value_at_jcol(1) =   eps_i_jshifted(i,j)
                value_at_jcol(2) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j))
                value_at_jcol(3) =   eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j)
                value_at_jcol(4) =   eps_i_jshifted(i,j+1)
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.25_8
                call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
          END SELECT

       END IF   !### IF (.NOT.block_has_symmetry_plane_X_left) THEN
    END IF      !### IF (ibegin.EQ.indx_x_min) THEN

!    i = indx_x_min+1

    i = indx_x_min+1
    irow_global = irow_global + 1

    IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
    ELSE

       IF (jbegin.EQ.indx_y_min) THEN                       ! BELOW
! boundary object along the bottom border
          jcolumn_global(1) = irow_global - (iend-ibegin+1)                          ! use the own node
       ELSE
! use a node from neighbor below
          jcolumn_global(1) = process_below_left_top_inner_node                      ! use a node from the neighbor below
       END IF
       IF (ibegin.EQ.indx_x_min) THEN                       ! LEFT
! boundary object along the left border
          jcolumn_global(2) = irow_global-1                                          ! use the own node
       ELSE
          jcolumn_global(2) = process_left_bottom_right_inner_node                   ! use a node from the left neighbor
       END IF
       jcolumn_global(3) = irow_global                      ! CENTER
       jcolumn_global(4) = irow_global+1                    ! RIGHT
       jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
       CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min+1, indx_y_min+1, nio, position_flag)

       SELECT CASE (position_flag)
          CASE (9)
! metal
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!          CASE (8)
!! dielectric surface on the right
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (4)
!! dielectric surface on the left
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (2)
!! dielectric surface above
!             value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (6)
!! dielectric surface below
!             value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
             value_at_jcol(1) =   eps_i_jshifted(i,j)
             value_at_jcol(2) =   eps_ishifted_j(i,j)
             value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
             value_at_jcol(4) =   eps_ishifted_j(i+1,j)
             value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.25_8
             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END SELECT
    END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

    IF (jbegin.EQ.indx_y_min) THEN
! boundary object along the bottom border
       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

             jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
             jcolumn_global(2) = irow_global-1                    ! LEFT
             jcolumn_global(3) = irow_global                      ! CENTER
             jcolumn_global(4) = irow_global+1                    ! RIGHT
             jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, indx_y_min+1, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (8)
!! dielectric surface on the right
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (4)
!! dielectric surface on the left
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (6)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) =   eps_ishifted_j(i,j)
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             END SELECT
          END IF    !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       END DO    !### DO i = indx_x_min+2, indx_x_max-2

    ELSE   !### IF (jbegin.EQ.indx_y_min) THEN

! use a node from neighbor below
      DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

             jcolumn_global(1) = process_below_left_top_inner_node + (i-indx_x_min-1)  ! BELOW
             jcolumn_global(2) = irow_global-1                                         ! LEFT
             jcolumn_global(3) = irow_global                                           ! CENTER
             jcolumn_global(4) = irow_global+1                                         ! RIGHT
             jcolumn_global(5) = irow_global + (iend-ibegin+1)                         ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, indx_y_min+1, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (8)
!! dielectric surface on the right
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (4)
!! dielectric surface on the left
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (6)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) =   eps_ishifted_j(i,j)
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             END SELECT
          END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       END DO   !### DO i = indx_x_min+2, indx_x_max-2
    END IF   !### IF (jbegin.EQ.indx_y_min) THEN
 
    i = indx_x_max-1
    irow_global = irow_global + 1

    IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
    ELSE

       IF (jbegin.EQ.indx_y_min) THEN                       ! BELOW
! boundary object along the bottom border
          jcolumn_global(1) = irow_global - (iend-ibegin+1)                         ! use the own node
       ELSE
! use a node from neighbor below
          jcolumn_global(1) = process_below_left_top_inner_node + (i-indx_x_min-1)  ! use a node from the neighbor below
       END IF
       jcolumn_global(2) = irow_global-1                    ! LEFT
       jcolumn_global(3) = irow_global                      ! CENTER
       IF (iend.EQ.indx_x_max) THEN                         ! RIGHT
! boundary object along the right border
          jcolumn_global(4) = irow_global+1                                         ! use the own node
       ELSE
          jcolumn_global(4) = process_right_bottom_left_inner_node                  ! use a node from the right neighbor
       END IF
       jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
       CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_max-1, indx_y_min+1, nio, position_flag)

       SELECT CASE (position_flag)
          CASE (9)
! metal
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!          CASE (8)
!! dielectric surface on the right
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (4)
!! dielectric surface on the left
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (2)
!! dielectric surface above
!             value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (6)
!! dielectric surface below
!             value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
             value_at_jcol(1) =   eps_i_jshifted(i,j)
             value_at_jcol(2) =   eps_ishifted_j(i,j)
             value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
             value_at_jcol(4) =   eps_ishifted_j(i+1,j)
             value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.25_8
             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END SELECT
    END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

    i = indx_x_max

    IF (iend.EQ.indx_x_max) THEN
! boundary object along right border
       irow_global = irow_global + 1
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
    END IF

    DO j = indx_y_min+2, indx_y_max-2 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

       i = indx_x_min

       IF (ibegin.EQ.indx_x_min) THEN
! boundary object along left border
          irow_global = irow_global + 1
          IF (.NOT.block_has_symmetry_plane_X_left) THEN
! Dirichlet (given potential) boundary
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
          ELSE
! the left border is a symmetry plane
! use own nodes
             jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
             jcolumn_global(2) = irow_global                      ! CENTER
             jcolumn_global(3) = irow_global+1                    ! RIGHT
             jcolumn_global(4) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min, j, nio, position_flag)

! note that we are at the left edge of the domain and the symmetry is applied here,
! so the boundary object must be symmetric relative to x=0 as well
! therefore we are either at the bottom surface, or inside, or at the top surface of the inner object

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (1,2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = -1.0_8
!                   value_at_jcol(3) = 0.5_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
!                CASE (6,7)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = -1.0_8
!                   value_at_jcol(3) = 0.5_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr)              
                CASE DEFAULT 
! inside dielectric or plasma
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j))
                   value_at_jcol(3) =   eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j)
                   value_at_jcol(4) =   eps_i_jshifted(i,j+1)              
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = -1.0_8
!                   value_at_jcol(3) = 0.5_8
!                   value_at_jcol(4) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
             END SELECT

          END IF   !### IF (.NOT.block_has_symmetry_plane_X_left) THEN
       END IF      !### IF (ibegin.EQ.indx_x_min) THEN

!       i = indx_x_min+1

       i = indx_x_min+1
       irow_global = irow_global + 1

       IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
       ELSE

          jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
          IF (ibegin.EQ.indx_x_min) THEN                       ! LEFT
! boundary object along the left border
             jcolumn_global(2) = irow_global-1                                                                                     ! use the own node
          ELSE
             jcolumn_global(2) = process_left_bottom_right_inner_node + (j-indx_y_min-1) * process_left_solved_nodes_row_length     ! use a node from the left neighbor
          END IF
          jcolumn_global(3) = irow_global                      ! CENTER
          jcolumn_global(4) = irow_global+1                    ! RIGHT
          jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
          CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min+1, j, nio, position_flag)

          SELECT CASE (position_flag)
             CASE (9)
! metal
                jcolumn_global(1) = irow_global
                value_at_jcol(1) = 1.0_8
                call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!             CASE (8)
!! dielectric surface on the right
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(5) = 0.25_8
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (4)
!! dielectric surface on the left
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(5) = 0.25_8
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (2)
!! dielectric surface above
!                value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (6)
!! dielectric surface below
!                value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                value_at_jcol(1) =   eps_i_jshifted(i,j)
                value_at_jcol(2) =   eps_ishifted_j(i,j)
                value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
                value_at_jcol(4) =   eps_ishifted_j(i+1,j)
                value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.25_8
                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          END SELECT
       END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

             jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
             jcolumn_global(2) = irow_global-1                    ! LEFT
             jcolumn_global(3) = irow_global                      ! CENTER
             jcolumn_global(4) = irow_global+1                    ! RIGHT
             jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, j, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (8)
!! dielectric surface on the right
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (4)
!! dielectric surface on the left
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (6)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) =   eps_ishifted_j(i,j)
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             END SELECT
          END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       END DO   !### DO i = indx_x_min+2, indx_x_max-2

       i = indx_x_max-1
       irow_global = irow_global + 1

       IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
       ELSE

          jcolumn_global(1) = irow_global - (iend-ibegin+1)       ! BELOW
          jcolumn_global(2) = irow_global-1                       ! LEFT
          jcolumn_global(3) = irow_global                         ! CENTER
          IF (iend.EQ.indx_x_max) THEN                            ! RIGHT
! boundary object along the right border
             jcolumn_global(4) = irow_global+1                                                                                   ! use the own node
          ELSE
             jcolumn_global(4) = process_right_bottom_left_inner_node + (j-indx_y_min-1) * process_right_solved_nodes_row_length  ! use a node from the right neighbor
          END IF
          jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
          CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_max-1, j, nio, position_flag)

          SELECT CASE (position_flag)
             CASE (9)
! metal
                jcolumn_global(1) = irow_global
                value_at_jcol(1) = 1.0_8
                call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!             CASE (8)
!! dielectric surface on the right
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(5) = 0.25_8
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (4)
!! dielectric surface on the left
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(5) = 0.25_8
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (2)
!! dielectric surface above
!                value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!             CASE (6)
!! dielectric surface below
!                value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                value_at_jcol(1) =   eps_i_jshifted(i,j)
                value_at_jcol(2) =   eps_ishifted_j(i,j)
                value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
                value_at_jcol(4) =   eps_ishifted_j(i+1,j)
                value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = 0.25_8
!                value_at_jcol(3) = -1.0_8
!                value_at_jcol(4) = 0.25_8
!                value_at_jcol(5) = 0.25_8
                call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          END SELECT
       END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

       i = indx_x_max

       IF (iend.EQ.indx_x_max) THEN
! boundary object along right border
          irow_global = irow_global + 1
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       END IF

    END DO !### DO j = indx_y_min+2, indx_y_max-2

    j = indx_y_max-1 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    i = indx_x_min

    IF (ibegin.EQ.indx_x_min) THEN
! boundary object along left border
       irow_global = irow_global + 1
       IF (.NOT.block_has_symmetry_plane_X_left) THEN
! Dirichlet (given potential) boundary
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       ELSE
! the left border is a symmetry plane
          jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
          jcolumn_global(2) = irow_global                      ! CENTER
          jcolumn_global(3) = irow_global+1                    ! RIGHT
          IF (jend.EQ.indx_y_max) THEN                         ! ABOVE
! boundary object along the top border
             jcolumn_global(4) = irow_global + (iend-ibegin+1)                         ! use the own node
          ELSE
             jcolumn_global(4) = process_above_left_bottom_inner_node-1                ! use a node from the neighbor above
          END IF

! check whether the point is inside or at the surface of any inner object
          CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min, indx_y_max-1, nio, position_flag)

! note that we are at the left edge of the domain and the symmetry is applied here,
! so the boundary object must be symmetric relative to x=0 as well
! therefore we are either at the bottom surface, or inside, or at the top surface of the inner object

          SELECT CASE (position_flag)
             CASE (9)
! metal
                jcolumn_global(1) = irow_global
                value_at_jcol(1) = 1.0_8
                call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!             CASE (1,2)
!! dielectric surface above
!                value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
!             CASE (6,7)
!! dielectric surface below
!                value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr)              
             CASE DEFAULT 
! inside dielectric or plasma
                value_at_jcol(1) =   eps_i_jshifted(i,j)
                value_at_jcol(2) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j))
                value_at_jcol(3) =   eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j)
                value_at_jcol(4) =   eps_i_jshifted(i,j+1)
!                value_at_jcol(1) = 0.25_8
!                value_at_jcol(2) = -1.0_8
!                value_at_jcol(3) = 0.5_8
!                value_at_jcol(4) = 0.25_8
                call MatSetValues(Amat, one, irow_global, four, jcolumn_global(1:4), value_at_jcol(1:4), INSERT_VALUES, ierr) 
          END SELECT

       END IF
    END IF

!    i = indx_x_min+1

    i = indx_x_min+1
    irow_global = irow_global + 1

    IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
    ELSE

       jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
       IF (ibegin.EQ.indx_x_min) THEN                       ! LEFT
! boundary object along the left border
          jcolumn_global(2) = irow_global-1                                                                                      ! use the own node
       ELSE
          jcolumn_global(2) = process_left_bottom_right_inner_node + (j-indx_y_min-1) * process_left_solved_nodes_row_length     ! use a node from the left neighbor
       END IF
       jcolumn_global(3) = irow_global                      ! CENTER
       jcolumn_global(4) = irow_global+1                    ! RIGHT
       IF (jend.EQ.indx_y_max) THEN                         ! ABOVE
! boundary object along the top border
          jcolumn_global(5) = irow_global + (iend-ibegin+1)                         ! use the own node
       ELSE
          jcolumn_global(5) = process_above_left_bottom_inner_node                  ! use a node from the neighbor above
       END IF

! check whether the point is inside or at the surface of any inner object
       CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_min+1, indx_y_max-1, nio, position_flag)

       SELECT CASE (position_flag)
          CASE (9)
! metal
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!          CASE (8)
!! dielectric surface on the right
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (4)
!! dielectric surface on the left
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (2)
!! dielectric surface above
!             value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (6)
!! dielectric surface below
!             value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
             value_at_jcol(1) =   eps_i_jshifted(i,j)
             value_at_jcol(2) =   eps_ishifted_j(i,j)
             value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
             value_at_jcol(4) =   eps_ishifted_j(i+1,j)
             value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.25_8
             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END SELECT
    END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

    IF (jend.EQ.indx_y_max) THEN
! boundary object along the top border
       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

             jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
             jcolumn_global(2) = irow_global-1                    ! LEFT
             jcolumn_global(3) = irow_global                      ! CENTER
             jcolumn_global(4) = irow_global+1                    ! RIGHT
             jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, indx_y_max-1, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (8)
!! dielectric surface on the right
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (4)
!! dielectric surface on the left
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (6)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) =   eps_ishifted_j(i,j)
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             END SELECT
          END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       END DO   !### DO i = indx_x_min+2, indx_x_max-2
    ELSE   !### IF (jend.EQ.indx_y_max) THEN
! use a node from the neighbor above
       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

             jcolumn_global(1) = irow_global - (iend-ibegin+1)                            ! BELOW
             jcolumn_global(2) = irow_global-1                                            ! LEFT
             jcolumn_global(3) = irow_global                                              ! CENTER
             jcolumn_global(4) = irow_global+1                                            ! RIGHT
             jcolumn_global(5) = process_above_left_bottom_inner_node + (i-indx_x_min-1)  ! ABOVE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, indx_y_max-1, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = irow_global
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!                CASE (8)
!! dielectric surface on the right
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (4)
!! dielectric surface on the left
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(5) = 0.25_8
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (2)
!! dielectric surface above
!                   value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!                CASE (6)
!! dielectric surface below
!                   value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
                   value_at_jcol(1) =   eps_i_jshifted(i,j)
                   value_at_jcol(2) =   eps_ishifted_j(i,j)
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!                   value_at_jcol(1) = 0.25_8
!                   value_at_jcol(2) = 0.25_8
!                   value_at_jcol(3) = -1.0_8
!                   value_at_jcol(4) = 0.25_8
!                   value_at_jcol(5) = 0.25_8
                   call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
             END SELECT
          END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       END DO   !### DO i = indx_x_min+2, indx_x_max-2
    END IF   !### IF (jend.EQ.indx_y_max) THEN
 
    i = indx_x_max-1
    irow_global = irow_global + 1

    IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
    ELSE

       jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
       jcolumn_global(2) = irow_global-1                    ! LEFT
       jcolumn_global(3) = irow_global                      ! CENTER
       IF (iend.EQ.indx_x_max) THEN                         ! RIGHT
! boundary object along the right border
          jcolumn_global(4) = irow_global+1                                                                                    ! use the own node
       ELSE
          jcolumn_global(4) = process_right_bottom_left_inner_node + (j-indx_y_min-1) * process_right_solved_nodes_row_length  ! use a node from the right neighbor
       END IF
       IF (jend.EQ.indx_y_max) THEN                         ! ABOVE
! boundary object along the top border
          jcolumn_global(5) = irow_global + (iend-ibegin+1)                            ! use the own node
       ELSE
          jcolumn_global(5) = process_above_left_bottom_inner_node + (i-indx_x_min-1)  ! use a node from the neighbor above
       END IF

! check whether the point is inside or at the surface of any inner object
       CALL FIND_INNER_OBJECT_CONTAINING_POINT(indx_x_max-1, indx_y_max-1, nio, position_flag)

       SELECT CASE (position_flag)
          CASE (9)
! metal
             jcolumn_global(1) = irow_global
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
!          CASE (8)
!! dielectric surface on the right
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (4)
!! dielectric surface on the left
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(5) = 0.25_8
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (2)
!! dielectric surface above
!             value_at_jcol(1) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
!          CASE (6)
!! dielectric surface below
!             value_at_jcol(1) = 0.5_8 * whole_object(nio)%eps_diel / (1.0_8 + whole_object(nio)%eps_diel)
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.5_8 / (1.0_8 + whole_object(nio)%eps_diel)
!             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
          CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object
             value_at_jcol(1) =   eps_i_jshifted(i,j)
             value_at_jcol(2) =   eps_ishifted_j(i,j)
             value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))
             value_at_jcol(4) =   eps_ishifted_j(i+1,j)
             value_at_jcol(5) =   eps_i_jshifted(i,j+1)
!             value_at_jcol(1) = 0.25_8
!             value_at_jcol(2) = 0.25_8
!             value_at_jcol(3) = -1.0_8
!             value_at_jcol(4) = 0.25_8
!             value_at_jcol(5) = 0.25_8
             call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END SELECT
    END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

    i = indx_x_max

    IF (iend.EQ.indx_x_max) THEN
! boundary object along right border
       irow_global = irow_global + 1
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
    END IF

!    j = indx_y_max !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.

    j = indx_y_max

    IF (jend.EQ.indx_y_max) THEN
! boundary object along top border
       DO i = ibegin, iend
          irow_global = irow_global + 1
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       END DO
    END IF

! end of initialization of matrix coefficients written by DS ---------------

    DEALLOCATE(eps_ishifted_j, STAT=ALLOC_ERR)   ! eps_ishifted_j(i,j) is between nodes {i-1,j} and {i,j}
    DEALLOCATE(eps_i_jshifted, STAT=ALLOC_ERR)   ! eps_i_jshifted(i,j) is between nodes {i,j-1} and {i,j}

! call MatSetValues(Amat, one, ix, ny, jcol, v, INSERT_VALUES,ierr) 
! inserts or adds a block of values into a matrix
! MatAssemblyBegin() and MatAssemblyEnd() must be called after all calls to matSetValues have been completed.
! Amat - the matrix
! one - the number of rows
! ix - global indices of rows
! ny - the number of columns
! jcol - global indices of columns
! v - a logically two-dimensional array of values
! INSERT_VALUES replaces existing entries with new values

! Commit changes to the matrix [S.J.]

! Begins assembling the matrix
    call MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY, ierr)
! Completes assembling the matrix
    call MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY, ierr)
    
! Save matrix for testing ....
!    if(.true.) then
!       call  PetscViewerBinaryOpen(parComm, 'Amat.bin', FILE_MODE_WRITE, viewer, ierr)
!       call  MatView(Amat, viewer, ierr)
!       call  PetscViewerDestroy(viewer,ierr)        
!    end if
    
! Create solver [S.J.]

! KSP is abstract petsc object that manages all Krylov methods. 
! It manages the linear solves (even when krylov accelerators are not used like in direct solvers

! Creates the default KSP context
! here ksp is location where to put KSP context
    call KSPCreate(parComm, ksp, ierr)

! Sets the matrix associated with the linear system and a (possibly) different one associated with the preconditioner
! argument #3 is the matrix to be used in constructing the preconditioner, 
! usually the same as argument #2 (the matrix that defines the linear system)
    call KSPSetOperators(ksp, Amat, Amat, ierr)

! . . Set method (this set-up uses the LU-decomposition as the solver) [S.J.]

! returns a pointer ot the preconditioner context set with KSPSetPc
    call KSPGetPc(ksp, pc, ierr)

! sets the preconditioner to be used to calculate the application of the preconditioner on a vector
    call KSPSetUp(ksp, ierr)

! Create B and X vectors [S.J.]

! note that ### bvec ### is the vector of the right hand side and ### xvec ### is the solution vector

! get vector(s) compatible with the matrix, i.e. with the same parallel layout
    call MatCreateVecs(Amat, PETSC_NULL_VEC, bvec, ierr)

! configures the vector from the options database
    call VecSetFromOptions(bvec, ierr)

! creates a new vector of the same type as an existing vector
    call VecDuplicate(bvec, xvec, ierr)
    
    ! Done with setting up vectors and matrices!
    return
  end subroutine SolverInitialization
  
!-------------------------------------
!
! this subroutine inserts dvec which in the actual call is the real(8) array rhs
! into vector bvec
  subroutine InsertData(nnodes,dvec)

    PetscInt, intent(IN) :: nnodes
    PetscScalar, intent(IN) :: dvec(:)

    PetscInt, save :: n, val
    PetscInt, allocatable, save :: jvec(:)

    integer ibegin, jbegin, iend, jend

    integer :: i, j !, k
!    PetscViewer :: viewer
!    character(20) :: filename='bvector.dat'
    PetscErrorCode :: ierr
    
    if(.not.allocated(jvec)) then

       allocate(jvec(nnodes))

! calculation of jvec modified by DS

       jbegin = indx_y_min+1
       jend   = indx_y_max-1
       ibegin = indx_x_min+1
       iend   = indx_x_max-1

       IF (Rank_of_process_left.LT.0)  ibegin = indx_x_min
       IF (Rank_of_process_right.LT.0) iend   = indx_x_max
       IF (Rank_of_process_below.LT.0) jbegin = indx_y_min
       IF (Rank_of_process_above.LT.0) jend   = indx_y_max

       n=0
       do j = jbegin, jend  !indx_y_min+1, indx_y_max-1
          do i = ibegin, iend  !indx_x_min+1, indx_x_max-1
!          if(i<1.or.i>=global_maximal_i.or.j<1.or.j>=global_maximal_j) cycle
             n = n+1
             jvec(n)=global_offset + n  !(N_grid_block_x*N_grid_block_y)*Rank_of_process+N_grid_block_x*(j-indx_y_min-1)+i-indx_x_min-1 
          end do
       end do

    end if

    val = n 
    if(.not.(nnodes.eq.n)) print *,'Not conformant!'
    if(size(dvec)<n) print *,'dvec too small!', size(dvec), n

! inserts or adds values into certain locations of a vector
! bvec - vector to insert in
! val - number of elements to add
! jvec - indices where to add
! dvec - array of values
! INSERT_VALUES - replaces existing entries with new values
    call VecSetValues(bvec, val, jvec, dvec, INSERT_VALUES, ierr)
    
! begins assembling the vector, should be called after completing all calls to VecSetValues
    call VecAssemblyBegin(bvec, ierr)

! completes assembling the vector
    call VecAssemblyEnd(bvec, ierr)

!    call PetscViewerCreate(parComm,viewer,ierr)
!    call PetscViewerASCIIOpen(parComm,filename,viewer,ierr)
!!    call PetscViewerSetType(viewer,PETSCVIEWERMATLAB,ierr)
!    call VecView(bvec,viewer,ierr)
!    call PetscViewerDestroy(parComm,viewer,ierr)

    return
  end subroutine InsertData

  subroutine FetchData(nnodes, dvec)

    PetscInt, intent(IN) :: nnodes
    PetscScalar, intent(OUT) :: dvec(:)

    integer ibegin, jbegin, iend, jend

    integer :: n, i, j   !, k, val
    PetscInt, allocatable, save :: jvec(:)
!    PetscViewer :: viewer
!    character(20) :: filename='soln.mat'
    PetscErrorCode :: ierr
    
    if(.not.allocated(jvec)) then
       allocate(jvec(nnodes))

! calculation of jvec modified by DS

       jbegin = indx_y_min+1
       jend   = indx_y_max-1
       ibegin = indx_x_min+1
       iend   = indx_x_max-1

       IF (Rank_of_process_left.LT.0)  ibegin = indx_x_min
       IF (Rank_of_process_right.LT.0) iend   = indx_x_max
       IF (Rank_of_process_below.LT.0) jbegin = indx_y_min
       IF (Rank_of_process_above.LT.0) jend   = indx_y_max

       n=0
       do j = jbegin, jend  !indx_y_min+1, indx_y_max-1
          do i = ibegin, iend  !indx_x_min+1, indx_x_max-1
!            if(i<1.or.i>=global_maximal_i.or.j<1.or.j>=global_maximal_j) cycle
             n = n+1
             jvec(n)=global_offset + n  !(N_grid_block_x*N_grid_block_y)*Rank_of_process+N_grid_block_x*(j-indx_y_min-1)+i-indx_x_min-1
          end do
       end do
    end if

! gets values from certain locations of a vector
    call VecGetValues(xvec, nnodes, jvec, dvec, ierr)
    
! begins assembling the vector
    call VecAssemblyBegin(xvec, ierr)
! completes assembling the vector
    call VecAssemblyEnd(xvec, ierr)

!  Get the local part of the solution
!    call PetscViewerCreate(parComm,viewer,ierr)
!    call PetscViewerASCIIOpen(parComm,filename,viewer,ierr)
!!    call PetscViewerSetType(viewer,PETSCVIEWERMATLAB,ierr)
!    call VecView(bvec,viewer,ierr)
!    call PetscViewerDestroy(parComm,viewer,ierr)
    return
  end subroutine FetchData
  
  subroutine SolverDestroy
    PetscErrorCode :: ierr
    call KSPDestroy(ksp,ierr)
    call MatDestroy(Amat,ierr)
    return
  end subroutine SolverDestroy
  
end module PETSc_Solver

!-------------------------------
! 7---6---5
! |       |
! 8   0   4
! |       |
! 1---2---3
!
SUBROUTINE FIND_INNER_OBJECT_CONTAINING_POINT(i, j, nio, position_flag)

  USE CurrentProblemValues

  use mpi

  IMPLICIT NONE

  INTEGER, INTENT(IN) ::  i, j  ! x,y indices of the point of interest
  INTEGER, INTENT(OUT) :: nio   ! number of inner object containing the point
                                ! it is meaningful only if position_flag>=0
  INTEGER, INTENT(OUT) :: position_flag   ! 9 for a metal object
                                          ! 0 inside a dielectric object
                                          ! 1 to 8 at the surface of the dielectric object
                                          ! -1 if there is no object containing the given point
  INTEGER n

  nio = -1
  position_flag = -1

  IF (N_of_inner_objects.EQ.0) RETURN

  DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
       IF (i.LT.whole_object(n)%ileft) CYCLE
       IF (i.GT.whole_object(n)%iright) CYCLE
       IF (j.LT.whole_object(n)%jbottom) CYCLE
       IF (j.GT.whole_object(n)%jtop) CYCLE

! since we are here, the point is either at the surface or inside inner material object n

       nio = n  ! save the number of the inner object

       IF (whole_object(n)%object_type.EQ.METAL_WALL) THEN
! the material object is metal, stop scanning inner objects, dielectric properties are omitted (relevant if metal is inside dielectric)
          position_flag = 9
          RETURN
       END IF

       position_flag = 0
       IF (i.EQ.whole_object(n)%ileft) THEN
! point on the left side of a dielectric object
          IF (j.EQ.whole_object(n)%jbottom) THEN
! left bottom corner
             position_flag = 1
          ELSE IF (j.EQ.whole_object(n)%jtop) THEN
! left top corner
             position_flag = 7
          ELSE
! left side, not a corner
             position_flag = 8
          END IF
       ELSE IF (i.EQ.whole_object(n)%iright) THEN
! point on the right side of a dielectric object
          IF (j.EQ.whole_object(n)%jbottom) THEN
! right bottom corner
             position_flag = 3
          ELSE IF (j.EQ.whole_object(n)%jtop) THEN
! right top corner
             position_flag = 5
          ELSE
! right side, not a corner
             position_flag = 4
          END IF
       ELSE IF (j.EQ.whole_object(n)%jbottom) THEN
! point on the bottom side of a dielectric object
          position_flag = 2
       ELSE IF (j.EQ.whole_object(n)%jtop) THEN
! point on the top side of a dielectric object
          position_flag = 6
       END IF
! note that we do not stop scanning inner material objects here since there may be metal objects inside the dielectric
    END DO

    RETURN

END SUBROUTINE FIND_INNER_OBJECT_CONTAINING_POINT

!-------------------------------
! 7---6---5
! |       |
! 8   0   4
! |       |
! 1---2---3
!
SUBROUTINE CHECK_IF_INNER_OBJECT_CONTAINS_POINT(myobject, i, j, position_flag)

  USE CurrentProblemValues

  use mpi

  IMPLICIT NONE

  TYPE(boundary_object), INTENT(IN) :: myobject
  INTEGER, INTENT(IN) :: i, j  ! x,y indices of the point of interest
  INTEGER, INTENT(OUT) :: position_flag   ! 9 for a metal object
                                          ! 0 inside a dielectric object
                                          ! 1 to 8 at the surface of the dielectric object
                                          ! -1 if there is no object containing the given point
  INTEGER n

  position_flag = -1

!  IF (N_of_inner_objects.EQ.0) RETURN

  IF (i.LT.myobject%ileft) RETURN
  IF (i.GT.myobject%iright) RETURN
  IF (j.LT.myobject%jbottom) RETURN
  IF (j.GT.myobject%jtop) RETURN

! since we are here, the point is either at the surface or inside inner material object n

  IF (myobject%object_type.EQ.METAL_WALL) THEN
! the material object is metal, stop scanning inner objects, dielectric properties are omitted (relevant if metal is inside dielectric)
     position_flag = 9
     RETURN
  END IF

  position_flag = 0
  IF (i.EQ.myobject%ileft) THEN
! point on the left side of a dielectric object
     IF (j.EQ.myobject%jbottom) THEN
! left bottom corner
        position_flag = 1
     ELSE IF (j.EQ.myobject%jtop) THEN
! left top corner
        position_flag = 7
     ELSE
! left side, not a corner
        position_flag = 8
     END IF
  ELSE IF (i.EQ.myobject%iright) THEN
! point on the right side of a dielectric object
     IF (j.EQ.myobject%jbottom) THEN
! right bottom corner
        position_flag = 3
     ELSE IF (j.EQ.myobject%jtop) THEN
! right top corner
        position_flag = 5
     ELSE
! right side, not a corner
        position_flag = 4
     END IF
  ELSE IF (j.EQ.myobject%jbottom) THEN
! point on the bottom side of a dielectric object
     position_flag = 2
  ELSE IF (j.EQ.myobject%jtop) THEN
! point on the top side of a dielectric object
     position_flag = 6
  END IF
! note that we do not stop scanning inner material objects here since there may be metal objects inside the dielectric

  RETURN

END SUBROUTINE CHECK_IF_INNER_OBJECT_CONTAINS_POINT

!--------------------------------------------------
!
REAL(8) FUNCTION Get_Surface_Charge_Inner_Object(i,j,position_flag, myobject)

  USE CurrentProblemValues
  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr
 
  INTEGER i, j            ! x,y indices of the point of interest
  INTEGER position_flag   ! 9 for a metal object
                          ! 0 inside a dielectric object
                          ! 1 to 8 at the surface of the dielectric object
                          ! -1 if there is no object containing the given point
  TYPE(boundary_object) myobject

  INTEGER pos

  INTEGER i_left_top, i_right_top, i_right_bottom

  Get_Surface_Charge_Inner_Object = 0.0_8 

  IF (N_of_inner_objects.EQ.0) RETURN
!  IF (nio.LT.N_of_boundary_objects+1) RETURN
!  IF (nio.GT.N_of_boundary_and_inner_objects) RETURN
  IF (position_flag.LT.1) RETURN
  IF (position_flag.GT.8) RETURN

  pos=-1

! ilt   -  ---- irt
!  |             |
!  |             |
!  1 ilbb ---- irb
!
  i_left_top        =                  myobject%jtop   - myobject%jbottom + 1
  i_right_top       = i_left_top     + myobject%iright - myobject%ileft
  i_right_bottom    = i_right_top    + myobject%jtop   - myobject%jbottom
!  i_left_bottom_bis = i_right_bottom + myobject%iright - myobject%ileft - 1

  SELECT CASE (position_flag)
     CASE (1)
! left bottom corner
        pos=1
     CASE (8)
! left side
        pos = j - myobject%jbottom  + 1
     CASE (7)
! left top corner
        pos = i_left_top
     CASE (6)
! top side
        pos = i_left_top + i - myobject%ileft
     CASE (5) 
! right top corner
        pos = i_right_top
     CASE (4)
! right side
        pos = i_right_top + myobject%jtop - j
     CASE (3)
! right bottom corner
        pos = i_right_bottom
     CASE (2)
  ! bottom side
        pos = i_right_bottom + myobject%iright - i
  END SELECT

  IF ((pos.GE.1).AND.(pos.LE.myobject%N_boundary_nodes)) THEN
     Get_Surface_Charge_Inner_Object = myobject%surface_charge(pos)
     RETURN
  ELSE
     print '("Error in Get_Surface_Charge_Inner_Object")'
     errcode=350
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

END FUNCTION Get_Surface_Charge_Inner_Object

!-----------------------------------------------
!
SUBROUTINE SET_EPS_ISHIFTED(i, j, eps)  ! here point {i,j} is between nodes {i-1,j} and {i,j}

  USE CurrentProblemValues

  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr
 
  INTEGER, INTENT(IN) :: i, j
  REAL(8), INTENT(OUT) :: eps

  INTEGER count  ! counts dielectric objects
  LOGICAL segment_is_inside_dielectric_object
  LOGICAL segment_is_on_surface_of_metal_object
  INTEGER n1

! find all inner objects owning segment {i-1,j}-{i,j}
! assume that only two dielectric objects may own a common segment
! allow a metal object to be added on top of that

  eps = 0.0_8
  count = 0
  segment_is_inside_dielectric_object = .FALSE.

! first, look through the dielectric objects only
  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.EQ.METAL_WALL) CYCLE

! find an inner object containing point {i-1,j}

     IF (i-1.LT.whole_object(n1)%ileft) CYCLE
     IF (i-1.GT.whole_object(n1)%iright) CYCLE
     IF (j.LT.whole_object(n1)%jbottom) CYCLE
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i-1,j} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (i.GT.whole_object(n1)%iright) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i-1,j}-{i,j} is owned by the object
     
! find whether segment {i-1,j}-{i,j} is at the surface or inside
! it only can be on the surface if this is the top or the bottom edge of the object

     IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is on the surface

        count = count + 1
        eps = eps + whole_object(n1)%eps_diel

     ELSE   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is inside
        segment_is_inside_dielectric_object = .TRUE.
        eps = whole_object(n1)%eps_diel
        count = 1
        EXIT

     END IF   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! at this point
! if count==0 then the segment is either in/at metal object or in vacuum
! if count==1 then the segment either inside a dielectric object or on the surface of ONE (not two) dielectric object
! if count==2 then the segment is at the interface between two dielectrics
! if count==3 this is an error

! second, look through the metal objects only
! if segment belongs to a metal object 
! return zero if it is inside
! if it is on the surface
! return (#1) epsilon of adjacent dielectric object or (#2) 1 if the metal object is in vacuum or (#3) 0 if the metal object is attached to another metal object

  segment_is_on_surface_of_metal_object = .FALSE.

  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.NE.METAL_WALL) CYCLE

! find an inner object containing point {i-1,j}

     IF (i-1.LT.whole_object(n1)%ileft) CYCLE
     IF (i-1.GT.whole_object(n1)%iright) CYCLE
     IF (j.LT.whole_object(n1)%jbottom) CYCLE
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i-1,j} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (i.GT.whole_object(n1)%iright) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i-1,j}-{i,j} is owned by the object
 
! find whether segment {i-1,j}-{i,j} is at the surface or inside
! it only can be on the surface if this is the top or the bottom edge of the object
     
     IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is on the surface

        IF (segment_is_on_surface_of_metal_object) THEN
! #3
! another metal object already has this segment on its border, so the segment is between two metal objects
! return zero epsilon
           eps = 0.0_8
           RETURN
        END IF
! remember that a metal object was found
        segment_is_on_surface_of_metal_object = .TRUE.

     ELSE   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is inside

        eps = 0.0_8
        RETURN

     END IF   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! note that cases when segment is inside a metal object or in between two attached metal objects are already processed above, followed by RETURN calls

  IF (segment_is_on_surface_of_metal_object) THEN
     IF (count.EQ.0) THEN
! segment on metal surface facing vacuum
        eps = 1.0_8
     ELSE
! segment on metal surface facing dielectric
        eps = eps / count
     END IF
     RETURN
  END IF

  IF (count.EQ.0) THEN
! segment is in vacuum
     eps = 1.0_8
     RETURN
  ELSE IF (count.EQ.1) THEN
! segment is inside a dielectric object
     IF (segment_is_inside_dielectric_object) RETURN
! segment is on the surface of a dielectric object, facing vacuum
     eps = (1.0_8 + eps) / 2.0_8
     RETURN
  ELSE IF (count.EQ.2) THEN
     eps = eps / DBLE(count)
     RETURN
  ELSE
     PRINT '("Error-3 in SET_EPS_ISHIFTED for i/j ",2x,i4,2x,i4)', i, j
     errcode=351
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

END SUBROUTINE SET_EPS_ISHIFTED


!-----------------------------------------------
!
SUBROUTINE SET_EPS_JSHIFTED(i, j, eps)  ! here point {i,j} is between nodes {i,j-1} and {i,j}

  USE CurrentProblemValues

  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr
 
  INTEGER, INTENT(IN) :: i, j
  REAL(8), INTENT(OUT) :: eps

  INTEGER count  ! counts dielectric objects
  LOGICAL segment_is_inside_dielectric_object
  LOGICAL segment_is_on_surface_of_metal_object
  INTEGER n1

! find all inner objects owning segment {i,j-1}-{i,j}
! assume that only two dielectric objects may own a common segment
! allow a metal object to be added on top of that

  eps = 0.0_8
  count = 0
  segment_is_inside_dielectric_object = .FALSE.

! first, look through the dielectric objects only
  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.EQ.METAL_WALL) CYCLE

! find an inner object containing point {i,j-1}

     IF (i.LT.whole_object(n1)%ileft) CYCLE
     IF (i.GT.whole_object(n1)%iright) CYCLE
     IF (j-1.LT.whole_object(n1)%jbottom) CYCLE
     IF (j-1.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i,j-1} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i,j-1}-{i,j} is owned by the object
     
! find whether segment {i,j-1}-{i,j} is at the surface or inside
! it only can be on the surface if this is the left or the right edge of the object

     IF ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN
! segment is on the surface

        count = count + 1
        eps = eps + whole_object(n1)%eps_diel

     ELSE   !### IF ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN
! segment is inside
        segment_is_inside_dielectric_object = .TRUE.
        eps = whole_object(n1)%eps_diel
        count = 1
        EXIT

     END IF   !### ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! at this point
! if count==0 then the segment is either in/on metal object or in vacuum
! if count==1 then the segment is either inside a dielectric object or on the surface of ONE (not two) dielectric object
! if count==2 then the segment is at the interface between two dielectrics
! if count>=3 this is an error

! second, look through the metal objects only
! if segment belongs to a metal object 
! return zero if it is inside
! if it is on the surface
! return (#1) epsilon of adjacent dielectric object or (#2) 1 if the metal object is in vacuum or (#3) 0 if the metal object is attached to another metal object

  segment_is_on_surface_of_metal_object = .FALSE.

  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.NE.METAL_WALL) CYCLE

! find an inner object containing point {i,j-1}

     IF (i.LT.whole_object(n1)%ileft) CYCLE
     IF (i.GT.whole_object(n1)%iright) CYCLE
     IF (j-1.LT.whole_object(n1)%jbottom) CYCLE
     IF (j-1.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i-1,j} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i,j-1}-{i,j} is owned by the object
 
! find whether segment {i-1,j}-{i,j} is at the surface or inside
! it only can be on the surface if this is the top or the bottom edge of the object
     
     IF ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN
! segment is on the surface

        IF (segment_is_on_surface_of_metal_object) THEN
! #3
! another metal object already has this segment on its border, so the segment is between two metal objects
! return zero epsilon
           eps = 0.0_8
           RETURN
        END IF
! remember that a metal object was found
        segment_is_on_surface_of_metal_object = .TRUE.

     ELSE   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is inside

        eps = 0.0_8
        RETURN

     END IF   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! note that cases when segment is inside a metal object or in between two attached metal objects are already processed above, followed by RETURN calls

  IF (segment_is_on_surface_of_metal_object) THEN
     IF (count.EQ.0) THEN
! segment on metal surface facing vacuum
        eps = 1.0_8
     ELSE
! segment on metal surface facing dielectric
        eps = eps / count
     END IF
     RETURN
  END IF

  IF (count.EQ.0) THEN
! segment is in vacuum
     eps = 1.0_8
     RETURN
  ELSE IF (count.EQ.1) THEN
! segment is inside a dielectric object
     IF (segment_is_inside_dielectric_object) RETURN
! segment is on the surface of a dielectric object, facing vacuum
     eps = (1.0_8 + eps) / 2.0_8
     RETURN
  ELSE IF (count.EQ.2) THEN
     eps = eps / DBLE(count)
     RETURN
  ELSE
     PRINT '("Error-3 in SET_EPS_JSHIFTED for i/j ",2x,i4,2x,i4)', i, j
     errcode=352
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

END SUBROUTINE SET_EPS_JSHIFTED

!-----------------------------------------------
! this is a simplified verions
! it is expected that point x is INSIDE some inner object or OUTSIDE, but not at the border of the object
!
SUBROUTINE GET_EPS_IN_POINT(x, y, eps)

  USE CurrentProblemValues

  use mpi

  IMPLICIT NONE


  INTEGER errcode,ierr
 
  REAL(8), INTENT(IN) :: x, y
  REAL(8), INTENT(OUT) :: eps

  INTEGER count  ! counts dielectric objects
  LOGICAL point_is_inside_dielectric_object
  LOGICAL point_is_on_surface_of_metal_object
  INTEGER n1

! find all inner objects owning segment {i,j-1}-{i,j}
! assume that only two dielectric objects may own a common segment
! allow a metal object to be added on top of that

  eps = 0.0_8
  count = 0
  point_is_inside_dielectric_object = .FALSE.

! first, look through the dielectric objects only
  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.EQ.METAL_WALL) CYCLE

! find an inner object containing point x, y

     IF (x.LT.whole_object(n1)%Xmin) CYCLE
     IF (x.GT.whole_object(n1)%Xmax) CYCLE
     IF (y.LT.whole_object(n1)%Ymin) CYCLE
     IF (y.GT.whole_object(n1)%Ymax) CYCLE

! since we are here, point x,y is either at the surface or inside inner object n1

     count = count + 1
     eps = eps + whole_object(n1)%eps_diel

     IF ( (x.GT.whole_object(n1)%Xmin).AND. &
        & (x.LT.whole_object(n1)%Xmax).AND. &
        & (y.GT.whole_object(n1)%Ymin).AND. &
        & (y.LT.whole_object(n1)%Ymax) ) THEN
! point inside
        point_is_inside_dielectric_object = .TRUE.
        EXIT
     END IF

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! at this point
! if count==0 then the point is either in/on metal object or in vacuum
! if count==1 then the point is either inside a dielectric object or on the surface of ONE (not two) dielectric object
! if count==2,3,4 then the point is at the interface between 2,3,4 dielectrics
! if count>=5 this is an error

! second, look through the metal objects only
! if point belongs to a metal object 
! return zero if it is inside
! if it is on the surface
! return (#1) epsilon of adjacent dielectric object or (#2) 1 if the metal object is in vacuum or (#3) 0 if the metal object is attached to another metal object

  point_is_on_surface_of_metal_object = .FALSE.

  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.NE.METAL_WALL) CYCLE

! find an inner object containing point x, y

     IF (x.LT.whole_object(n1)%Xmin) CYCLE
     IF (x.GT.whole_object(n1)%Xmax) CYCLE
     IF (y.LT.whole_object(n1)%Ymin) CYCLE
     IF (y.GT.whole_object(n1)%Ymax) CYCLE

! since we are here, point x,y is either at the surface or inside inner object n1

     IF ( (x.EQ.whole_object(n1)%Xmin).OR. &
        & (x.EQ.whole_object(n1)%Xmax).OR. &
        & (y.EQ.whole_object(n1)%Ymin).OR. &
        & (y.EQ.whole_object(n1)%Ymax) ) THEN
! point is on the surface
        IF (point_is_on_surface_of_metal_object) THEN
! #3
! another metal object already has this segment on its border, so the segment is between two metal objects
! return zero epsilon
           eps = 0.0_8
           RETURN
        END IF
! remember that a metal object was found
        point_is_on_surface_of_metal_object = .TRUE.
     ELSE
! point is inside
        eps = 0.0_8
        RETURN
     END IF

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! note that cases when segment is inside a metal object or in between two attached metal objects are already processed above, followed by RETURN calls

  IF (point_is_on_surface_of_metal_object) THEN
     IF (count.EQ.0) THEN
! point on metal surface facing vacuum
        eps = 1.0_8
     ELSE
! point on metal surface facing dielectric
        eps = eps / count
     END IF
     RETURN
  END IF

  IF (count.EQ.0) THEN
! point is in vacuum
     eps = 1.0_8
     RETURN
  ELSE IF (count.EQ.1) THEN
! point is inside a dielectric object
     IF (point_is_inside_dielectric_object) RETURN
! segment is on the surface of a dielectric object, facing vacuum
     eps = (1.0_8 + eps) / 2.0_8
     RETURN
  ELSE IF (count.LT.4) THEN
     eps = eps / DBLE(count)
     RETURN
  ELSE
     PRINT '("Error-3 in GET_EPS_IN_POINT for x/y ",2x,f11.3,2x,f11.3)', x, y
     errcode=353
     CALL MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
  END IF

END SUBROUTINE GET_EPS_IN_POINT

