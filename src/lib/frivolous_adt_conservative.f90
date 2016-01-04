!< Define the abstract type *conservative* for building FriVolous finite-volume conservative residuals.
module frivolous_adt_conservative
!-----------------------------------------------------------------------------------------------------------------------------------
!< Define the abstract type *conservative* for building FriVolous finite-volume conservative residuals.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use frivolous_kinds, only : R_P
use vecfor, only : vector
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: conservative
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: conservative
  private
  !< Abstract type for building FriVolous finite-volume conservative residuals.
#ifdef CAF
  class(*), allocatable :: dummy_to_allow_extensions[:] !< Dummy member to allow concrete extensions with coarray members.
#endif
  contains
    private
    ! public deferred procedures that concrete integrand-field must implement
    ! operators
    procedure(symmetric_operator),      pass(lhs), deferred, public :: conservative_add_conservative !< Conservative+conservative.
    procedure(symmetric_operator),      pass(lhs), deferred, public :: conservative_sub_conservative !< Conservative-conservative.
    procedure(conservative_op_real),    pass(lhs), deferred, public :: conservative_multiply_real    !< Conservative*real.
    procedure(conservative_op_real),    pass(lhs), deferred, public :: conservative_divide_real      !< Conservative/real.
    procedure(conservative_op_vector),  pass(lhs), deferred, public :: conservative_dot_vector       !< Conservative.dot.vector.
    procedure(assignment_conservative), pass(lhs), deferred, public :: assign_conservative           !< Conservative=Conservative.
    ! operators overloading
    generic, public :: operator(+) => conservative_add_conservative !< Overloading + operator.
    generic, public :: operator(-) => conservative_sub_conservative !< Overloading - operator.
    generic, public :: operator(*) => conservative_multiply_real    !< Overloading * operator.
    generic, public :: operator(/) => conservative_divide_real      !< Overloading / operator.
    generic, public :: operator(.dot.) => conservative_dot_vector   !< Overloading .dot. operator.
    generic, public :: assignment(=) => assign_conservative         !< Overloading = assignament.
endtype conservative

abstract interface
  !< Abstract type bound procedures necessary for implementing a concrete extension of the class(conservative).
  function symmetric_operator(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Symmetric type operator integrand.op.integrand.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: conservative
  class(conservative), intent(IN)  :: lhs             !< Left hand side.
  class(conservative), intent(IN)  :: rhs             !< Right hand side.
  class(conservative), allocatable :: operator_result !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction symmetric_operator

  function conservative_op_real(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Asymmetric type operator conservative.op.real.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: conservative, R_P
  class(conservative), intent(IN)  :: lhs              !< Left hand side.
  real(R_P),           intent(IN)  :: rhs              !< Right hand side.
  class(conservative), allocatable :: operator_result  !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction conservative_op_real

  function conservative_op_vector(lhs, rhs) result(operator_result)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Asymmetric type operator conservative.op.vector.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: conservative, vector
  class(conservative), intent(IN)  :: lhs              !< Left hand side.
  type(vector),        intent(IN)  :: rhs              !< Right hand side.
  class(conservative), allocatable :: operator_result  !< Operator result.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction conservative_op_vector

  pure subroutine assignment_conservative(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Symmetric assignment conservative = conservative.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: conservative
  class(conservative), intent(INOUT) :: lhs !< Left hand side.
  class(conservative), intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assignment_conservative
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule frivolous_adt_conservative
