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

    PetscInt :: one, five

    character(20) :: petsc_config='petsc.rc'   ! PETSc database file
    PetscErrorCode :: ierr

    integer jbegin, jend, ibegin, iend
    PetscInt :: irow_global
    PetscInt :: jcolumn_global(1:5)
    PetscScalar :: value_at_jcol(1:5)
    integer i, j

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

!    i = indx_x_min

    IF (ibegin.EQ.indx_x_min) THEN
! boundary object along left border
       irow_global = irow_global + 1
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
    END IF

!    i = indx_x_min+1

    irow_global = irow_global + 1
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
    value_at_jcol(1) = 0.25_8
    value_at_jcol(2) = 0.25_8
    value_at_jcol(3) = -1.0_8
    value_at_jcol(4) = 0.25_8
    value_at_jcol(5) = 0.25_8
    call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 

    IF (jbegin.EQ.indx_y_min) THEN
! boundary object along the bottom border
       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1
          jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
          jcolumn_global(2) = irow_global-1                    ! LEFT
          jcolumn_global(3) = irow_global                      ! CENTER
          jcolumn_global(4) = irow_global+1                    ! RIGHT
          jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE
          value_at_jcol(1) = 0.25_8
          value_at_jcol(2) = 0.25_8
          value_at_jcol(3) = -1.0_8
          value_at_jcol(4) = 0.25_8
          value_at_jcol(5) = 0.25_8
          call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END DO
    ELSE
! use a node from neighbor below
      DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1
          jcolumn_global(1) = process_below_left_top_inner_node + (i-indx_x_min-1)  ! BELOW
          jcolumn_global(2) = irow_global-1                                         ! LEFT
          jcolumn_global(3) = irow_global                                           ! CENTER
          jcolumn_global(4) = irow_global+1                                         ! RIGHT
          jcolumn_global(5) = irow_global + (iend-ibegin+1)                         ! ABOVE
          value_at_jcol(1) = 0.25_8
          value_at_jcol(2) = 0.25_8
          value_at_jcol(3) = -1.0_8
          value_at_jcol(4) = 0.25_8
          value_at_jcol(5) = 0.25_8
          call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END DO
    END IF
 
    i = indx_x_max-1
    irow_global = irow_global + 1
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
    value_at_jcol(1) = 0.25_8
    value_at_jcol(2) = 0.25_8
    value_at_jcol(3) = -1.0_8
    value_at_jcol(4) = 0.25_8
    value_at_jcol(5) = 0.25_8
    call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 

!    i = indx_x_max

    IF (iend.EQ.indx_x_max) THEN
! boundary object along right border
       irow_global = irow_global + 1
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
    END IF

    DO j = indx_y_min+2, indx_y_max-2 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!       i = indx_x_min

       IF (ibegin.EQ.indx_x_min) THEN
! boundary object along left border
          irow_global = irow_global + 1
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       END IF

!       i = indx_x_min+1

       irow_global = irow_global + 1
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
       value_at_jcol(1) = 0.25_8
       value_at_jcol(2) = 0.25_8
       value_at_jcol(3) = -1.0_8
       value_at_jcol(4) = 0.25_8
       value_at_jcol(5) = 0.25_8
       call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 

       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1
          jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
          jcolumn_global(2) = irow_global-1                    ! LEFT
          jcolumn_global(3) = irow_global                      ! CENTER
          jcolumn_global(4) = irow_global+1                    ! RIGHT
          jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE
          value_at_jcol(1) = 0.25_8
          value_at_jcol(2) = 0.25_8
          value_at_jcol(3) = -1.0_8
          value_at_jcol(4) = 0.25_8
          value_at_jcol(5) = 0.25_8
          call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END DO

       i = indx_x_max-1
       irow_global = irow_global + 1
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
       value_at_jcol(1) = 0.25_8
       value_at_jcol(2) = 0.25_8
       value_at_jcol(3) = -1.0_8
       value_at_jcol(4) = 0.25_8
       value_at_jcol(5) = 0.25_8
       call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 

!       i = indx_x_max

       IF (iend.EQ.indx_x_max) THEN
! boundary object along right border
          irow_global = irow_global + 1
          jcolumn_global(1) = irow_global
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       END IF

    END DO !### DO j = indx_y_min+2, indx_y_max-2

    j = indx_y_max-1 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!    i = indx_x_min

    IF (ibegin.EQ.indx_x_min) THEN
! boundary object along left border
       irow_global = irow_global + 1
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
    END IF

!    i = indx_x_min+1

    irow_global = irow_global + 1
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
    value_at_jcol(1) = 0.25_8
    value_at_jcol(2) = 0.25_8
    value_at_jcol(3) = -1.0_8
    value_at_jcol(4) = 0.25_8
    value_at_jcol(5) = 0.25_8
    call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 

    IF (jend.EQ.indx_y_max) THEN
! boundary object along the top border
       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1
          jcolumn_global(1) = irow_global - (iend-ibegin+1)    ! BELOW
          jcolumn_global(2) = irow_global-1                    ! LEFT
          jcolumn_global(3) = irow_global                      ! CENTER
          jcolumn_global(4) = irow_global+1                    ! RIGHT
          jcolumn_global(5) = irow_global + (iend-ibegin+1)    ! ABOVE
          value_at_jcol(1) = 0.25_8
          value_at_jcol(2) = 0.25_8
          value_at_jcol(3) = -1.0_8
          value_at_jcol(4) = 0.25_8
          value_at_jcol(5) = 0.25_8
          call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END DO
    ELSE
! use a node from the neighbor above
       DO i = indx_x_min+2, indx_x_max-2
          irow_global = irow_global + 1
          jcolumn_global(1) = irow_global - (iend-ibegin+1)                            ! BELOW
          jcolumn_global(2) = irow_global-1                                            ! LEFT
          jcolumn_global(3) = irow_global                                              ! CENTER
          jcolumn_global(4) = irow_global+1                                            ! RIGHT
          jcolumn_global(5) = process_above_left_bottom_inner_node + (i-indx_x_min-1)  ! ABOVE
          value_at_jcol(1) = 0.25_8
          value_at_jcol(2) = 0.25_8
          value_at_jcol(3) = -1.0_8
          value_at_jcol(4) = 0.25_8
          value_at_jcol(5) = 0.25_8
          call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 
       END DO
    END IF
 
    i = indx_x_max-1
    irow_global = irow_global + 1
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
    value_at_jcol(1) = 0.25_8
    value_at_jcol(2) = 0.25_8
    value_at_jcol(3) = -1.0_8
    value_at_jcol(4) = 0.25_8
    value_at_jcol(5) = 0.25_8
    call MatSetValues(Amat, one, irow_global, five, jcolumn_global, value_at_jcol, INSERT_VALUES, ierr) 

!    i = indx_x_max

    IF (iend.EQ.indx_x_max) THEN
! boundary object along right border
       irow_global = irow_global + 1
       jcolumn_global(1) = irow_global
       value_at_jcol(1) = 1.0_8
       call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
    END IF

!    j = indx_y_max !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
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
