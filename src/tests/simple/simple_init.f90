!< FriVolous, Finite Volume block-structured Fortran test.
!<
!< A simple *use* and *init* test: it initializes a simple, Cartesian block.

program simple_init
!< FriVolous, Finite Volume block-structured Fortran test.
!<
!< A simple *use* and *init* test: it initializes a simple, Cartesian block.
use frivolous, only : block_object

implicit none
type(block_object) :: ablock

call ablock%initialize(cartesian=.true., Ni=3, Nj=4, Nk=5)
print "(A, L1)", "Is the block Cartesian? ", ablock%is_cartesian()
endprogram simple_init
