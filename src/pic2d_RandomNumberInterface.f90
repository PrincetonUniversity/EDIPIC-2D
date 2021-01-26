!-------------------------------------------------------
! this interface was written by Janhunen Jani Salomon
!
module rng_wrapper

  use, intrinsic :: iso_c_binding

  interface
     subroutine set_rng_state(init_state,state_ind,func) bind(C,name="set_state")
       use iso_c_binding, only : c_int
       integer(KIND=C_INT), intent(INOUT), dimension(624) :: init_state
       integer(KIND=C_INT), intent(INOUT)                 :: state_ind, func
     end subroutine set_rng_state

     subroutine get_rng_state(init_state,state_ind,func) bind(C,name="get_state")
       use iso_c_binding, only : c_int
       integer(KIND=C_INT), intent(INOUT), dimension(624) :: init_state
       integer(KIND=C_INT), intent(INOUT)                 :: state_ind, func
     end subroutine get_rng_state
  end interface

  abstract interface
     real(KIND=C_DOUBLE) function WELLRNG_c()
       use, intrinsic :: ISO_C_BINDING
     end function WELLRNG_c
  end interface

  type(c_funptr), bind(C,name="WELLRNG19937") :: WELLRNG19937

contains

  real(8) function well_random_number()
    procedure(WELLRNG_c), POINTER :: WELLRNG_f
    CALL C_F_PROCPOINTER(WELLRNG19937,WELLRNG_f)
    well_random_number=WELLRNG_f()
    return
  end function well_random_number

  subroutine well_random_seed(seed)
    use iso_c_binding, only : c_int
    integer(KIND=C_INT) :: init_state(624)
    integer(KIND=C_INT) :: seed
    integer(KIND=C_INT) :: j, k
    integer*8           :: tmp

!     init_state(1)=iand(seed,4294967295)   ! AND 0xFFFFFFFF
    init_state(1)=iand(seed,X'FFFFFFFF')   ! AND 0xFFFFFFFF
    do j=2,624
       tmp=(1812433253_8*(ieor(init_state(j-1),ishft(init_state(j-1),30)))+j-1)
!       tmp=(X'1812433253'*(ieor(init_state(j-1),ishft(init_state(j-1),30)))+j-1)  
! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
       init_state(j)=int(iand(tmp,X'FFFFFFFF'),kind(i4)) ! AND 0xFFFFFFFF
    end do
    j=0; k=1; 
    call set_rng_state(init_state,j,k)
    return
  end subroutine well_random_seed

end module rng_wrapper  

